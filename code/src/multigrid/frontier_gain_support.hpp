#pragma once

#include "multigrid/residual_budget_support.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace tgi {

struct FrontierGainSupportOptions {
    int base_patch_layers = 2;
    int maximum_rounds = 10;
    int maximum_extra_nodes_per_column = 128;
    int maximum_total_nodes_per_round = 256;
    int maximum_nodes_per_column_per_round = 16;
    int minimum_nodes_per_marked_column = 2;
    int maximum_stagnant_rounds = 3;
    double column_marking_fraction = 0.80;
    double candidate_marking_fraction = 0.90;
    double target_residual_ratio = 0.22;
    double absolute_residual_tolerance = 1.0e-12;
    double minimum_round_relative_gain = 0.002;
    double refinement_tolerance = 1.0e-6;
    StrengthScaling strength_scaling = StrengthScaling::SymmetricDiagonal;
    double strong_edge_fraction = 0.10;
    int thread_count = 1;
};

struct FrontierGainRoundReport {
    int round = 0;
    int active_columns = 0;
    int marked_columns = 0;
    int expanded_columns = 0;
    int frozen_columns = 0;
    int added_nodes = 0;
    double mean_scaled_residual = 0.0;
    double maximum_scaled_residual = 0.0;
    double mean_observed_gain = 0.0;
    double selection_ms = 0.0;
    double refinement_ms = 0.0;
};

struct FrontierGainSupportReport {
    int rounds = 0;
    int expanded_columns = 0;
    int frozen_columns = 0;
    int total_extra_nodes = 0;
    int maximum_extra_nodes = 0;
    double initial_mean_scaled_residual = 0.0;
    double initial_maximum_scaled_residual = 0.0;
    double final_mean_scaled_residual = 0.0;
    double final_maximum_scaled_residual = 0.0;
    double mean_final_to_initial_ratio = 0.0;
    double graph_ms = 0.0;
    double residual_ms = 0.0;
    double selection_ms = 0.0;
    double refinement_ms = 0.0;
    double total_ms = 0.0;
    std::vector<FrontierGainRoundReport> history;
};

struct FrontierGainSupportResult {
    SparseMatrix prolongation;
    std::vector<std::vector<int>> supports;
    FrontierGainSupportReport report;
};

namespace frontier_gain_support_detail {

using Clock = std::chrono::steady_clock;

inline double milliseconds(Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

struct Candidate {
    double score = 0.0;
    int coarse = -1;
    int node = -1;
};

} // namespace frontier_gain_support_detail

inline FrontierGainSupportResult build_frontier_gain_interpolation(
    const StructuredGrid& grid, const SparseMatrix& a,
    const SparseMatrix& initial_prolongation,
    const InterpolationOptions& interpolation_options,
    const FrontierGainSupportOptions& options = {}) {
    if (a.rows() != grid.fine_size() || a.cols() != grid.fine_size() ||
        initial_prolongation.rows() != grid.fine_size() ||
        initial_prolongation.cols() != grid.coarse_size()) {
        throw std::invalid_argument(
            "build_frontier_gain_interpolation: incompatible dimensions");
    }
    if (options.base_patch_layers <= 0 || options.maximum_rounds < 0 ||
        options.maximum_extra_nodes_per_column < 0 ||
        options.maximum_total_nodes_per_round <= 0 ||
        options.maximum_nodes_per_column_per_round <= 0 ||
        options.minimum_nodes_per_marked_column < 0 ||
        options.minimum_nodes_per_marked_column >
            options.maximum_nodes_per_column_per_round ||
        options.maximum_stagnant_rounds <= 0 ||
        !(options.column_marking_fraction > 0.0 &&
          options.column_marking_fraction <= 1.0) ||
        !(options.candidate_marking_fraction > 0.0 &&
          options.candidate_marking_fraction <= 1.0) ||
        !(options.target_residual_ratio >= 0.0) ||
        !(options.minimum_round_relative_gain >= 0.0) ||
        !(options.refinement_tolerance > 0.0) ||
        !(options.strong_edge_fraction >= 0.0) ||
        options.thread_count <= 0) {
        throw std::invalid_argument(
            "build_frontier_gain_interpolation: invalid options");
    }

    using namespace frontier_gain_support_detail;
    using residual_budget_support_detail::ResidualSnapshot;
    using residual_budget_support_detail::build_strength_graph;
    using residual_budget_support_detail::mean_and_maximum;
    using residual_budget_support_detail::scaled_f_residual;

    const auto total_begin = Clock::now();
    FrontierGainSupportResult result;
    result.prolongation = initial_prolongation;
    result.supports.resize(static_cast<std::size_t>(grid.coarse_size()));
    for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
        grid.patch_f_nodes(
            coarse, options.base_patch_layers,
            result.supports[static_cast<std::size_t>(coarse)]);
    }

    const auto graph_begin = Clock::now();
    const SparseMatrix strong_graph = build_strength_graph(
        a, options.strength_scaling, options.strong_edge_fraction);
    result.report.graph_ms = milliseconds(graph_begin, Clock::now());

    const auto initial_residual_begin = Clock::now();
    ResidualSnapshot snapshot = scaled_f_residual(
        grid, a, result.prolongation, result.supports,
        options.thread_count);
    result.report.residual_ms +=
        milliseconds(initial_residual_begin, Clock::now());
    const std::vector<double> initial_norm = snapshot.norm;
    const auto [initial_mean, initial_maximum] =
        mean_and_maximum(initial_norm);
    result.report.initial_mean_scaled_residual = initial_mean;
    result.report.initial_maximum_scaled_residual = initial_maximum;

    std::vector<int> extras(
        static_cast<std::size_t>(grid.coarse_size()), 0);
    std::vector<int> stagnant(
        static_cast<std::size_t>(grid.coarse_size()), 0);
    std::vector<unsigned char> ever_expanded(
        static_cast<std::size_t>(grid.coarse_size()), 0U);
    std::vector<unsigned char> frozen(
        static_cast<std::size_t>(grid.coarse_size()), 0U);
    std::vector<unsigned char> refine_column(
        static_cast<std::size_t>(grid.coarse_size()), 0U);
    std::vector<unsigned char> in_support(
        static_cast<std::size_t>(grid.fine_size()), 0U);
    Vector residual_value(
        static_cast<std::size_t>(grid.fine_size()), 0.0);
    Vector best_score(
        static_cast<std::size_t>(grid.fine_size()), 0.0);
    const Vector diagonal = a.diagonal();
    std::vector<int> touched;
    std::vector<int> added_this_round(
        static_cast<std::size_t>(grid.coarse_size()), 0);

    for (int round = 1; round <= options.maximum_rounds; ++round) {
        FrontierGainRoundReport round_report;
        round_report.round = round;
        const auto [round_mean, round_maximum] =
            mean_and_maximum(snapshot.norm);
        round_report.mean_scaled_residual = round_mean;
        round_report.maximum_scaled_residual = round_maximum;
        std::fill(refine_column.begin(), refine_column.end(), 0U);
        std::fill(added_this_round.begin(), added_this_round.end(), 0);

        const auto selection_begin = Clock::now();
        std::vector<std::pair<double, int>> active;
        double active_squared = 0.0;
        for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
            const std::size_t offset = static_cast<std::size_t>(coarse);
            const double target = std::max(
                options.absolute_residual_tolerance,
                options.target_residual_ratio * initial_norm[offset]);
            if (frozen[offset] != 0U || snapshot.norm[offset] <= target ||
                extras[offset] >= options.maximum_extra_nodes_per_column) {
                continue;
            }
            const double squared = snapshot.norm[offset] * snapshot.norm[offset];
            active.push_back({squared, coarse});
            active_squared += squared;
        }
        round_report.active_columns = static_cast<int>(active.size());
        std::sort(active.begin(), active.end(),
            [](const auto& lhs, const auto& rhs) {
                if (lhs.first != rhs.first) return lhs.first > rhs.first;
                return lhs.second < rhs.second;
            });

        std::vector<int> marked;
        double marked_squared = 0.0;
        for (const auto& [squared, coarse] : active) {
            marked.push_back(coarse);
            marked_squared += squared;
            if (marked_squared >=
                options.column_marking_fraction * active_squared) break;
        }
        round_report.marked_columns = static_cast<int>(marked.size());

        std::vector<Candidate> candidates;
        for (int coarse : marked) {
            const std::size_t coarse_offset = static_cast<std::size_t>(coarse);
            auto& support = result.supports[coarse_offset];
            for (int node : support) {
                in_support[static_cast<std::size_t>(node)] = 1U;
            }
            for (int position = snapshot.transpose.row_ptr()[coarse_offset];
                 position < snapshot.transpose.row_ptr()[coarse_offset + 1U];
                 ++position) {
                const int node = snapshot.transpose.col_idx()[
                    static_cast<std::size_t>(position)];
                residual_value[static_cast<std::size_t>(node)] =
                    std::abs(snapshot.transpose.values()[
                        static_cast<std::size_t>(position)]) /
                    std::sqrt(diagonal[static_cast<std::size_t>(node)]);
            }
            for (int source : support) {
                const std::size_t source_offset = static_cast<std::size_t>(source);
                for (int position = strong_graph.row_ptr()[source_offset];
                     position < strong_graph.row_ptr()[source_offset + 1U];
                     ++position) {
                    const int node = strong_graph.col_idx()[
                        static_cast<std::size_t>(position)];
                    const std::size_t node_offset = static_cast<std::size_t>(node);
                    if (grid.is_coarse_node(node) ||
                        in_support[node_offset] != 0U) continue;
                    if (best_score[node_offset] == 0.0) touched.push_back(node);
                    const double strength = strong_graph.values()[
                        static_cast<std::size_t>(position)];
                    // Residual payoff divided by one added degree of freedom.
                    // A small source-leakage term stabilizes ties when the
                    // frontier residual is initially tiny.
                    const double score = strength *
                        (residual_value[node_offset] +
                         0.20 * residual_value[source_offset] + 1.0e-30);
                    best_score[node_offset] =
                        std::max(best_score[node_offset], score);
                }
            }
            for (int node : touched) {
                candidates.push_back(
                    {best_score[static_cast<std::size_t>(node)], coarse, node});
                best_score[static_cast<std::size_t>(node)] = 0.0;
            }
            touched.clear();
            for (int node : support) {
                in_support[static_cast<std::size_t>(node)] = 0U;
            }
            for (int position = snapshot.transpose.row_ptr()[coarse_offset];
                 position < snapshot.transpose.row_ptr()[coarse_offset + 1U];
                 ++position) {
                residual_value[static_cast<std::size_t>(
                    snapshot.transpose.col_idx()[static_cast<std::size_t>(position)])]
                    = 0.0;
            }
        }
        std::sort(candidates.begin(), candidates.end(),
            [](const Candidate& lhs, const Candidate& rhs) {
                if (lhs.score != rhs.score) return lhs.score > rhs.score;
                if (lhs.coarse != rhs.coarse) return lhs.coarse < rhs.coarse;
                return lhs.node < rhs.node;
            });
        double total_candidate_score = 0.0;
        for (const Candidate& candidate : candidates) {
            total_candidate_score += candidate.score;
        }
        double selected_score = 0.0;
        std::vector<unsigned char> selected(candidates.size(), 0U);
        // A small fairness floor prevents a few very large columns from
        // consuming the complete round budget and starving other marked
        // channel segments.
        for (std::size_t index = 0; index < candidates.size(); ++index) {
            const Candidate& candidate = candidates[index];
            if (round_report.added_nodes >=
                options.maximum_total_nodes_per_round) break;
            const std::size_t coarse_offset =
                static_cast<std::size_t>(candidate.coarse);
            if (added_this_round[coarse_offset] >=
                    options.minimum_nodes_per_marked_column ||
                extras[coarse_offset] + added_this_round[coarse_offset] >=
                    options.maximum_extra_nodes_per_column) {
                continue;
            }
            result.supports[coarse_offset].push_back(candidate.node);
            selected[index] = 1U;
            ++added_this_round[coarse_offset];
            ++round_report.added_nodes;
            selected_score += candidate.score;
        }
        for (std::size_t index = 0; index < candidates.size(); ++index) {
            if (selected[index] != 0U) continue;
            const Candidate& candidate = candidates[index];
            if (round_report.added_nodes >=
                options.maximum_total_nodes_per_round) break;
            const std::size_t coarse_offset =
                static_cast<std::size_t>(candidate.coarse);
            if (added_this_round[coarse_offset] >=
                    options.maximum_nodes_per_column_per_round ||
                extras[coarse_offset] + added_this_round[coarse_offset] >=
                    options.maximum_extra_nodes_per_column) {
                continue;
            }
            result.supports[coarse_offset].push_back(candidate.node);
            ++added_this_round[coarse_offset];
            ++round_report.added_nodes;
            selected_score += candidate.score;
            if (selected_score >= options.candidate_marking_fraction *
                    total_candidate_score) break;
        }
        for (int coarse : marked) {
            const std::size_t offset = static_cast<std::size_t>(coarse);
            if (added_this_round[offset] == 0) continue;
            auto& support = result.supports[offset];
            std::sort(support.begin(), support.end());
            support.erase(std::unique(support.begin(), support.end()), support.end());
            extras[offset] += added_this_round[offset];
            refine_column[offset] = 1U;
            ever_expanded[offset] = 1U;
            ++round_report.expanded_columns;
        }
        round_report.selection_ms = milliseconds(selection_begin, Clock::now());
        result.report.selection_ms += round_report.selection_ms;
        if (round_report.expanded_columns == 0) {
            result.report.history.push_back(round_report);
            break;
        }

        const std::vector<double> previous_norm = snapshot.norm;
        const auto refinement_begin = Clock::now();
        InterpolationOptions warm_options = interpolation_options;
        warm_options.thread_count = options.thread_count;
        warm_options.local_tolerance = std::min(
            warm_options.local_tolerance, options.refinement_tolerance);
        const InterpolationResult refined =
            refine_selected_energy_interpolation_on_supports(
                grid, a, result.supports, refine_column,
                result.prolongation, warm_options);
        result.prolongation = refined.prolongation;
        round_report.refinement_ms =
            milliseconds(refinement_begin, Clock::now());
        result.report.refinement_ms += round_report.refinement_ms;

        const auto residual_begin = Clock::now();
        snapshot = scaled_f_residual(
            grid, a, result.prolongation, result.supports,
            options.thread_count);
        result.report.residual_ms += milliseconds(residual_begin, Clock::now());
        double gain_sum = 0.0;
        int gain_count = 0;
        for (int coarse : marked) {
            const std::size_t offset = static_cast<std::size_t>(coarse);
            if (added_this_round[offset] == 0) continue;
            const double gain = previous_norm[offset] > 0.0
                ? std::max(0.0, 1.0 - snapshot.norm[offset] /
                    previous_norm[offset])
                : 0.0;
            gain_sum += gain;
            ++gain_count;
            if (gain < options.minimum_round_relative_gain) {
                ++stagnant[offset];
                if (stagnant[offset] >= options.maximum_stagnant_rounds) {
                    frozen[offset] = 1U;
                    ++round_report.frozen_columns;
                }
            } else {
                stagnant[offset] = 0;
            }
        }
        round_report.mean_observed_gain = gain_count > 0
            ? gain_sum / static_cast<double>(gain_count) : 0.0;
        result.report.history.push_back(round_report);
        result.report.rounds = round;
    }

    const auto [final_mean, final_maximum] = mean_and_maximum(snapshot.norm);
    result.report.final_mean_scaled_residual = final_mean;
    result.report.final_maximum_scaled_residual = final_maximum;
    double ratio_sum = 0.0;
    for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
        const std::size_t offset = static_cast<std::size_t>(coarse);
        if (ever_expanded[offset] != 0U) ++result.report.expanded_columns;
        if (frozen[offset] != 0U) ++result.report.frozen_columns;
        result.report.total_extra_nodes += extras[offset];
        result.report.maximum_extra_nodes = std::max(
            result.report.maximum_extra_nodes, extras[offset]);
        ratio_sum += initial_norm[offset] > 0.0
            ? snapshot.norm[offset] / initial_norm[offset] : 0.0;
    }
    result.report.mean_final_to_initial_ratio = ratio_sum /
        static_cast<double>(grid.coarse_size());
    result.report.total_ms = milliseconds(total_begin, Clock::now());
    return result;
}

} // namespace tgi
