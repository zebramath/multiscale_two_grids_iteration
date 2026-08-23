#pragma once

#include "multigrid/energy_interpolation.hpp"
#include "multigrid/two_grid_solver.hpp"
#include "pde/diffusion_problem.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <exception>
#include <functional>
#include <limits>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

namespace tgi {

struct ResidualStrongSupportOptions {
    int base_patch_layers = 2;
    int maximum_extra_nodes_per_column = 64;
    int maximum_graph_hops = 256;
    double strong_edge_fraction = 0.25;
    int thread_count = 1;
};

struct FixedSupportExpansionReport {
    double selection_ms = 0.0;
    double maximum_off_diagonal = 0.0;
    double strong_edge_threshold = 0.0;
    int expanded_columns = 0;
    int total_extra_nodes = 0;
    int maximum_extra_nodes = 0;
    double mean_scaled_column_residual = 0.0;
    double maximum_scaled_column_residual = 0.0;
};

struct FixedSupportExpansionResult {
    std::vector<std::vector<int>> supports;
    FixedSupportExpansionReport report;
};

inline FixedSupportExpansionResult build_residual_strong_supports(
    const StructuredGrid& grid, const SparseMatrix& a,
    const SparseMatrix& base_prolongation,
    const ResidualStrongSupportOptions& options = {});

namespace support_expansion_detail {

using Clock = std::chrono::steady_clock;

inline double milliseconds(Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

struct QueueEntry {
    double score = 0.0;
    int hops = 0;
    int node = -1;
};

struct QueueLess {
    bool operator()(const QueueEntry& lhs,
                    const QueueEntry& rhs) const {
        if (lhs.score != rhs.score) return lhs.score < rhs.score;
        if (lhs.hops != rhs.hops) return lhs.hops > rhs.hops;
        return lhs.node > rhs.node;
    }
};

inline double maximum_off_diagonal(const SparseMatrix& a) {
    double maximum = 0.0;
    for (int row = 0; row < a.rows(); ++row) {
        for (int position =
                 a.row_ptr()[static_cast<std::size_t>(row)];
             position <
                 a.row_ptr()[static_cast<std::size_t>(row) + 1U];
             ++position) {
            if (a.col_idx()[static_cast<std::size_t>(position)] != row) {
                maximum = std::max(
                    maximum,
                    std::abs(a.values()[static_cast<std::size_t>(position)]));
            }
        }
    }
    return maximum;
}

inline SparseMatrix build_strong_graph(
    const SparseMatrix& a, double threshold) {
    std::vector<int> row_ptr(
        static_cast<std::size_t>(a.rows()) + 1U, 0);
    std::vector<int> col_idx;
    Vector weights;
    col_idx.reserve(a.nnz());
    weights.reserve(a.nnz());
    for (int row = 0; row < a.rows(); ++row) {
        for (int position =
                 a.row_ptr()[static_cast<std::size_t>(row)];
             position <
                 a.row_ptr()[static_cast<std::size_t>(row) + 1U];
             ++position) {
            const int col =
                a.col_idx()[static_cast<std::size_t>(position)];
            const double edge = std::abs(
                a.values()[static_cast<std::size_t>(position)]);
            if (col != row && edge >= threshold) {
                col_idx.push_back(col);
                weights.push_back(edge);
            }
        }
        row_ptr[static_cast<std::size_t>(row) + 1U] =
            static_cast<int>(weights.size());
    }
    return SparseMatrix(
        a.rows(), a.cols(), std::move(row_ptr),
        std::move(col_idx), std::move(weights));
}

}


inline FixedSupportExpansionResult build_residual_strong_supports(
    const StructuredGrid& grid, const SparseMatrix& a,
    const SparseMatrix& base_prolongation,
    const ResidualStrongSupportOptions& options) {
    FixedSupportExpansionResult result;
    result.supports.resize(
        static_cast<std::size_t>(grid.coarse_size()));
    const auto begin = support_expansion_detail::Clock::now();
    result.report.maximum_off_diagonal =
        support_expansion_detail::maximum_off_diagonal(a);
    result.report.strong_edge_threshold =
        options.strong_edge_fraction *
        result.report.maximum_off_diagonal;
    if (!(result.report.maximum_off_diagonal > 0.0)) {
        throw std::runtime_error(
            "build_residual_strong_supports: matrix has no graph edges");
    }
    const SparseMatrix strong_graph =
        support_expansion_detail::build_strong_graph(
            a, result.report.strong_edge_threshold);

    const SparseMatrix ap = sparse_multiply(
        a, base_prolongation, 0.0, options.thread_count);
    const SparseMatrix ap_transpose =
        ap.transpose(options.thread_count);
    const Vector diagonal = a.diagonal();
    std::vector<double> column_residuals(
        static_cast<std::size_t>(grid.coarse_size()), 0.0);
    for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
        double squared = 0.0;
        for (int position =
                 ap_transpose.row_ptr()[static_cast<std::size_t>(coarse)];
             position <
                 ap_transpose.row_ptr()[
                     static_cast<std::size_t>(coarse) + 1U];
             ++position) {
            const int node =
                ap_transpose.col_idx()[static_cast<std::size_t>(position)];
            if (grid.is_coarse_node(node)) continue;
            const double scaled =
                ap_transpose.values()[static_cast<std::size_t>(position)] /
                diagonal[static_cast<std::size_t>(node)];
            squared += scaled * scaled;
        }
        const double value = std::sqrt(squared);
        column_residuals[static_cast<std::size_t>(coarse)] = value;
        result.report.mean_scaled_column_residual += value;
        result.report.maximum_scaled_column_residual = std::max(
            result.report.maximum_scaled_column_residual, value);
    }
    result.report.mean_scaled_column_residual /=
        static_cast<double>(grid.coarse_size());
    std::vector<unsigned char> coarse_mask(
        static_cast<std::size_t>(grid.fine_size()), 0U);
    for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
        coarse_mask[static_cast<std::size_t>(
            grid.coarse_fine_id(coarse))] = 1U;
    }
    std::vector<int> extras_per_column(
        static_cast<std::size_t>(grid.coarse_size()), 0);
    const unsigned int requested = options.thread_count > 0
        ? static_cast<unsigned int>(options.thread_count)
        : std::thread::hardware_concurrency();
    const int worker_count = std::max(
        1, std::min(
               grid.coarse_size(),
               static_cast<int>(requested == 0U ? 1U : requested)));
    std::atomic<int> next_coarse{0};
    std::exception_ptr worker_error;
    std::mutex error_mutex;

    auto worker = [&]() {
        Vector residual(
            static_cast<std::size_t>(grid.fine_size()), 0.0);
        std::vector<int> residual_active;
        std::vector<unsigned char> in_support(
            static_cast<std::size_t>(grid.fine_size()), 0U);
        std::vector<unsigned char> visited(
            static_cast<std::size_t>(grid.fine_size()), 0U);
        std::vector<int> visited_nodes;
        try {
            while (true) {
                const int coarse = next_coarse.fetch_add(
                    1, std::memory_order_relaxed);
                if (coarse >= grid.coarse_size()) break;
                auto& support =
                    result.supports[static_cast<std::size_t>(coarse)];
                grid.patch_f_nodes(
                    coarse, options.base_patch_layers, support);
                for (int node : support) {
                    in_support[static_cast<std::size_t>(node)] = 1U;
                }

                residual_active.clear();
                for (int position =
                         ap_transpose.row_ptr()[
                             static_cast<std::size_t>(coarse)];
                     position <
                         ap_transpose.row_ptr()[
                             static_cast<std::size_t>(coarse) + 1U];
                     ++position) {
                    const int node =
                        ap_transpose.col_idx()[
                            static_cast<std::size_t>(position)];
                    if (coarse_mask[
                            static_cast<std::size_t>(node)] != 0U) {
                        continue;
                    }
                    const double scaled =
                        ap_transpose.values()[
                            static_cast<std::size_t>(position)] /
                        diagonal[static_cast<std::size_t>(node)];
                    residual[static_cast<std::size_t>(node)] = scaled;
                    residual_active.push_back(node);
                }

                int extras = 0;
                if (column_residuals[
                        static_cast<std::size_t>(coarse)] >
                        0.0 &&
                    options.maximum_extra_nodes_per_column > 0) {
                    std::priority_queue<
                        support_expansion_detail::QueueEntry,
                        std::vector<
                            support_expansion_detail::QueueEntry>,
                        support_expansion_detail::QueueLess> queue;

                    for (int node : support) {
                        const double source_score = std::abs(
                            residual[static_cast<std::size_t>(node)]);
                        for (int position =
                                 strong_graph.row_ptr()[
                                     static_cast<std::size_t>(node)];
                             position <
                                 strong_graph.row_ptr()[
                                     static_cast<std::size_t>(node)
                                     + 1U];
                             ++position) {
                            const int neighbor =
                                strong_graph.col_idx()[
                                    static_cast<std::size_t>(
                                        position)];
                            if (in_support[
                                    static_cast<std::size_t>(
                                        neighbor)] != 0U) {
                                continue;
                            }
                            const double edge =
                                strong_graph.values()[
                                    static_cast<std::size_t>(
                                        position)];
                            const double neighbor_score = std::abs(
                                residual[static_cast<std::size_t>(
                                    neighbor)]);
                            queue.push({
                                std::max(
                                    source_score, neighbor_score) +
                                    edge /
                                        result.report
                                            .maximum_off_diagonal,
                                1,
                                neighbor
                            });
                        }
                    }

                    while (!queue.empty() &&
                           extras <
                               options
                                   .maximum_extra_nodes_per_column) {
                        const auto entry = queue.top();
                        queue.pop();
                        if (entry.hops >
                            options.maximum_graph_hops) {
                            continue;
                        }
                        const std::size_t offset =
                            static_cast<std::size_t>(entry.node);
                        if (visited[offset] != 0U ||
                            in_support[offset] != 0U) {
                            continue;
                        }
                        visited[offset] = 1U;
                        visited_nodes.push_back(entry.node);
                        if (coarse_mask[offset] == 0U) {
                            support.push_back(entry.node);
                            in_support[offset] = 1U;
                            ++extras;
                        }
                        for (int position =
                                 strong_graph.row_ptr()[offset];
                             position <
                                 strong_graph.row_ptr()[offset + 1U];
                             ++position) {
                            const int neighbor =
                                strong_graph.col_idx()[
                                    static_cast<std::size_t>(
                                        position)];
                            const std::size_t neighbor_offset =
                                static_cast<std::size_t>(neighbor);
                            if (visited[neighbor_offset] != 0U ||
                                in_support[neighbor_offset] != 0U) {
                                continue;
                            }
                            const double edge =
                                strong_graph.values()[
                                    static_cast<std::size_t>(
                                        position)];
                            queue.push({
                                entry.score * 0.999 +
                                    1.0e-3 * edge /
                                        result.report
                                            .maximum_off_diagonal,
                                entry.hops + 1,
                                neighbor
                            });
                        }
                    }
                }

                std::sort(support.begin(), support.end());
                support.erase(
                    std::unique(support.begin(), support.end()),
                    support.end());
                extras_per_column[
                    static_cast<std::size_t>(coarse)] = extras;

                for (int node : support) {
                    in_support[static_cast<std::size_t>(node)] = 0U;
                }
                for (int node : residual_active) {
                    residual[static_cast<std::size_t>(node)] = 0.0;
                }
                for (int node : visited_nodes) {
                    visited[static_cast<std::size_t>(node)] = 0U;
                }
                visited_nodes.clear();
            }
        } catch (...) {
            std::lock_guard<std::mutex> lock(error_mutex);
            if (!worker_error) {
                worker_error = std::current_exception();
            }
            next_coarse.store(grid.coarse_size());
        }
    };

    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(worker_count));
    for (int index = 0; index < worker_count; ++index) {
        workers.emplace_back(worker);
    }
    for (auto& thread : workers) thread.join();
    if (worker_error) std::rethrow_exception(worker_error);

    for (int extras : extras_per_column) {
        if (extras > 0) ++result.report.expanded_columns;
        result.report.total_extra_nodes += extras;
        result.report.maximum_extra_nodes = std::max(
            result.report.maximum_extra_nodes, extras);
    }
    result.report.selection_ms =
        support_expansion_detail::milliseconds(
            begin, support_expansion_detail::Clock::now());
    return result;
}

}

namespace tgi {

enum class StrengthScaling {
    GlobalMaximum,
    RowMaximum,
    SymmetricDiagonal
};

struct ResidualBudgetSupportOptions {
    int base_patch_layers = 2;
    int maximum_rounds = 8;
    int maximum_extra_nodes_per_column = 128;
    int maximum_nodes_per_round = 16;
    double marking_fraction = 0.70;
    double target_residual_ratio = 0.25;
    double absolute_residual_tolerance = 1.0e-12;
    double refinement_tolerance = 1.0e-6;
    StrengthScaling strength_scaling = StrengthScaling::RowMaximum;
    double strong_edge_fraction = 0.25;
    int thread_count = 1;
};

struct ResidualBudgetRoundReport {
    int round = 0;
    int active_columns = 0;
    int expanded_columns = 0;
    int added_nodes = 0;
    double mean_scaled_residual = 0.0;
    double maximum_scaled_residual = 0.0;
    double selection_ms = 0.0;
    double refinement_ms = 0.0;
};

struct ResidualBudgetSupportReport {
    int rounds = 0;
    int expanded_columns = 0;
    int converged_columns = 0;
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
    std::vector<ResidualBudgetRoundReport> history;
};

struct ResidualBudgetSupportResult {
    SparseMatrix prolongation;
    std::vector<std::vector<int>> supports;
    ResidualBudgetSupportReport report;
};

inline ResidualBudgetSupportResult build_residual_budget_interpolation(
    const StructuredGrid& grid, const SparseMatrix& a,
    const SparseMatrix& initial_prolongation,
    const InterpolationOptions& interpolation_options,
    const ResidualBudgetSupportOptions& support_options = {});

namespace residual_budget_support_detail {

using Clock = std::chrono::steady_clock;

inline double milliseconds(Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

struct QueueEntry {
    double score = 0.0;
    int hops = 0;
    int node = -1;
};

struct QueueLess {
    bool operator()(const QueueEntry& lhs, const QueueEntry& rhs) const {
        if (lhs.score != rhs.score) return lhs.score < rhs.score;
        if (lhs.hops != rhs.hops) return lhs.hops > rhs.hops;
        return lhs.node > rhs.node;
    }
};

inline SparseMatrix build_strength_graph(
    const SparseMatrix& a, StrengthScaling scaling, double fraction) {
    const Vector diagonal = a.diagonal();
    Vector row_max(static_cast<std::size_t>(a.rows()), 0.0);
    double global_max = 0.0;
    for (int row = 0; row < a.rows(); ++row) {
        for (int position = a.row_ptr()[static_cast<std::size_t>(row)];
             position < a.row_ptr()[static_cast<std::size_t>(row) + 1U];
             ++position) {
            const int column =
                a.col_idx()[static_cast<std::size_t>(position)];
            if (column == row) continue;
            const double edge =
                std::abs(a.values()[static_cast<std::size_t>(position)]);
            row_max[static_cast<std::size_t>(row)] = std::max(
                row_max[static_cast<std::size_t>(row)], edge);
            global_max = std::max(global_max, edge);
        }
    }
    if (!(global_max > 0.0)) {
        throw std::runtime_error("strength graph has no off-diagonal edges");
    }

    std::vector<int> row_ptr(
        static_cast<std::size_t>(a.rows()) + 1U, 0);
    std::vector<int> col_idx;
    Vector weights;
    col_idx.reserve(a.nnz());
    weights.reserve(a.nnz());
    for (int row = 0; row < a.rows(); ++row) {
        for (int position = a.row_ptr()[static_cast<std::size_t>(row)];
             position < a.row_ptr()[static_cast<std::size_t>(row) + 1U];
             ++position) {
            const int column =
                a.col_idx()[static_cast<std::size_t>(position)];
            if (column == row) continue;
            const double edge =
                std::abs(a.values()[static_cast<std::size_t>(position)]);
            double normalized = 0.0;
            if (scaling == StrengthScaling::GlobalMaximum) {
                normalized = edge / global_max;
            } else if (scaling == StrengthScaling::RowMaximum) {
                const double scale = row_max[static_cast<std::size_t>(row)];
                normalized = scale > 0.0 ? edge / scale : 0.0;
            } else {
                const double scale = std::sqrt(
                    diagonal[static_cast<std::size_t>(row)] *
                    diagonal[static_cast<std::size_t>(column)]);
                normalized = scale > 0.0 ? edge / scale : 0.0;
            }
            if (normalized >= fraction) {
                col_idx.push_back(column);
                weights.push_back(normalized);
            }
        }
        row_ptr[static_cast<std::size_t>(row) + 1U] =
            static_cast<int>(weights.size());
    }
    return SparseMatrix(
        a.rows(), a.cols(), std::move(row_ptr),
        std::move(col_idx), std::move(weights));
}

struct ResidualSnapshot {
    SparseMatrix transpose;
    std::vector<double> norm;
};

inline ResidualSnapshot scaled_f_residual(
    const StructuredGrid& grid, const SparseMatrix& a,
    const SparseMatrix& prolongation,
    const std::vector<std::vector<int>>& supports,
    int thread_count) {
    const SparseMatrix ap = sparse_multiply(
        a, prolongation, 0.0, thread_count);
    SparseMatrix transpose = ap.transpose(thread_count);
    const Vector diagonal = a.diagonal();
    std::vector<double> norm(
        static_cast<std::size_t>(grid.coarse_size()), 0.0);
    std::vector<unsigned char> in_support(
        static_cast<std::size_t>(grid.fine_size()), 0U);
    for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
        const std::size_t coarse_offset = static_cast<std::size_t>(coarse);
        for (int node : supports[coarse_offset]) {
            in_support[static_cast<std::size_t>(node)] = 1U;
        }
        double squared = 0.0;
        for (int position =
                 transpose.row_ptr()[static_cast<std::size_t>(coarse)];
             position < transpose.row_ptr()[
                 static_cast<std::size_t>(coarse) + 1U];
             ++position) {
            const int node = transpose.col_idx()[
                static_cast<std::size_t>(position)];
            if (grid.is_coarse_node(node) ||
                in_support[static_cast<std::size_t>(node)] != 0U) {
                continue;
            }
            const double scaled = transpose.values()[
                static_cast<std::size_t>(position)] /
                std::sqrt(diagonal[static_cast<std::size_t>(node)]);
            squared += scaled * scaled;
        }
        norm[static_cast<std::size_t>(coarse)] = std::sqrt(squared);
        for (int node : supports[coarse_offset]) {
            in_support[static_cast<std::size_t>(node)] = 0U;
        }
    }
    return {std::move(transpose), std::move(norm)};
}

inline std::pair<double, double> mean_and_maximum(
    const std::vector<double>& values) {
    double mean = 0.0;
    double maximum = 0.0;
    for (double value : values) {
        mean += value;
        maximum = std::max(maximum, value);
    }
    if (!values.empty()) mean /= static_cast<double>(values.size());
    return {mean, maximum};
}

}

inline ResidualBudgetSupportResult build_residual_budget_interpolation(
    const StructuredGrid& grid, const SparseMatrix& a,
    const SparseMatrix& initial_prolongation,
    const InterpolationOptions& interpolation_options,
    const ResidualBudgetSupportOptions& options) {
    using namespace residual_budget_support_detail;
    const auto total_begin = Clock::now();
    ResidualBudgetSupportResult result;
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
    result.report.residual_ms += milliseconds(
        initial_residual_begin, Clock::now());
    const std::vector<double> initial_norm = snapshot.norm;
    const auto [initial_mean, initial_maximum] =
        mean_and_maximum(initial_norm);
    result.report.initial_mean_scaled_residual = initial_mean;
    result.report.initial_maximum_scaled_residual = initial_maximum;

    std::vector<int> extras(
        static_cast<std::size_t>(grid.coarse_size()), 0);
    std::vector<unsigned char> ever_expanded(
        static_cast<std::size_t>(grid.coarse_size()), 0U);
    std::vector<unsigned char> refine_column(
        static_cast<std::size_t>(grid.coarse_size()), 0U);
    std::vector<unsigned char> in_support(
        static_cast<std::size_t>(grid.fine_size()), 0U);
    Vector residual_value(
        static_cast<std::size_t>(grid.fine_size()), 0.0);
    const Vector diagonal = a.diagonal();
    std::vector<int> residual_nodes;
    std::vector<unsigned char> visited(
        static_cast<std::size_t>(grid.fine_size()), 0U);
    std::vector<int> visited_nodes;
    std::vector<int> marked_columns;

    for (int round_index = 0;
         round_index < options.maximum_rounds; ++round_index) {
        ResidualBudgetRoundReport round_report;
        round_report.round = round_index + 1;
        const auto [round_mean, round_maximum] =
            mean_and_maximum(snapshot.norm);
        round_report.mean_scaled_residual = round_mean;
        round_report.maximum_scaled_residual = round_maximum;
        std::fill(refine_column.begin(), refine_column.end(), 0U);

        const auto selection_begin = Clock::now();
        std::vector<std::pair<double, int>> active_columns;
        double active_squared = 0.0;
        for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
            const std::size_t coarse_offset =
                static_cast<std::size_t>(coarse);
            const double target = std::max(
                options.absolute_residual_tolerance,
                options.target_residual_ratio * initial_norm[coarse_offset]);
            if (snapshot.norm[coarse_offset] <= target ||
                extras[coarse_offset] >=
                    options.maximum_extra_nodes_per_column) {
                continue;
            }
            const double squared = snapshot.norm[coarse_offset] *
                snapshot.norm[coarse_offset];
            active_columns.push_back({squared, coarse});
            active_squared += squared;
        }
        round_report.active_columns =
            static_cast<int>(active_columns.size());
        std::sort(
            active_columns.begin(), active_columns.end(),
            [](const auto& lhs, const auto& rhs) {
                if (lhs.first != rhs.first) return lhs.first > rhs.first;
                return lhs.second < rhs.second;
            });
        marked_columns.clear();
        double marked_squared = 0.0;
        const double column_target =
            options.marking_fraction * active_squared;
        for (const auto& [squared, coarse] : active_columns) {
            marked_columns.push_back(coarse);
            marked_squared += squared;
            if (marked_squared >= column_target) break;
        }

        for (int coarse : marked_columns) {
            const std::size_t coarse_offset =
                static_cast<std::size_t>(coarse);
            auto& support = result.supports[coarse_offset];
            for (int node : support) {
                in_support[static_cast<std::size_t>(node)] = 1U;
            }

            residual_nodes.clear();
            for (int position = snapshot.transpose.row_ptr()[coarse_offset];
                 position < snapshot.transpose.row_ptr()[coarse_offset + 1U];
                 ++position) {
                const int node = snapshot.transpose.col_idx()[
                    static_cast<std::size_t>(position)];
                if (grid.is_coarse_node(node) ||
                    in_support[static_cast<std::size_t>(node)] != 0U) {
                    continue;
                }
                residual_value[static_cast<std::size_t>(node)] =
                    std::abs(snapshot.transpose.values()[
                        static_cast<std::size_t>(position)]) /
                    std::sqrt(diagonal[static_cast<std::size_t>(node)]);
                residual_nodes.push_back(node);
            }

            std::priority_queue<
                QueueEntry, std::vector<QueueEntry>, QueueLess> queue;
            for (int node : support) {
                const std::size_t node_offset =
                    static_cast<std::size_t>(node);
                for (int position = strong_graph.row_ptr()[node_offset];
                     position < strong_graph.row_ptr()[node_offset + 1U];
                     ++position) {
                    const int neighbor = strong_graph.col_idx()[
                        static_cast<std::size_t>(position)];
                    const double strength = strong_graph.values()[
                        static_cast<std::size_t>(position)];
                    if (in_support[static_cast<std::size_t>(neighbor)] != 0U) {
                        continue;
                    }
                    const double residual_seed =
                        residual_value[node_offset] +
                        residual_value[static_cast<std::size_t>(neighbor)];
                    queue.push({
                        (residual_seed + 1.0e-12) * strength,
                        1,
                        neighbor
                    });
                }
            }

            int added = 0;
            const int available =
                options.maximum_extra_nodes_per_column -
                extras[coarse_offset];
            const int limit = std::min(
                options.maximum_nodes_per_round, available);
            while (!queue.empty() && added < limit) {
                const QueueEntry entry = queue.top();
                queue.pop();
                const std::size_t offset =
                    static_cast<std::size_t>(entry.node);
                if (visited[offset] != 0U || in_support[offset] != 0U) {
                    continue;
                }
                visited[offset] = 1U;
                visited_nodes.push_back(entry.node);
                if (!grid.is_coarse_node(entry.node)) {
                    support.push_back(entry.node);
                    in_support[offset] = 1U;
                    ++added;
                }
                for (int position = strong_graph.row_ptr()[offset];
                     position < strong_graph.row_ptr()[offset + 1U];
                     ++position) {
                    const int neighbor = strong_graph.col_idx()[
                        static_cast<std::size_t>(position)];
                    const std::size_t neighbor_offset =
                        static_cast<std::size_t>(neighbor);
                    if (visited[neighbor_offset] != 0U ||
                        in_support[neighbor_offset] != 0U) {
                        continue;
                    }
                    const double strength = strong_graph.values()[
                        static_cast<std::size_t>(position)];
                    const double residual_seed =
                        residual_value[neighbor_offset];
                    queue.push({
                        entry.score * (0.90 + 0.10 * strength) +
                            residual_seed,
                        entry.hops + 1,
                        neighbor
                    });
                }
            }
            if (added > 0) {
                std::sort(support.begin(), support.end());
                support.erase(
                    std::unique(support.begin(), support.end()),
                    support.end());
                extras[coarse_offset] += added;
                refine_column[coarse_offset] = 1U;
                ever_expanded[coarse_offset] = 1U;
                ++round_report.expanded_columns;
                round_report.added_nodes += added;
            }

            for (int node : support) {
                in_support[static_cast<std::size_t>(node)] = 0U;
            }
            for (int node : residual_nodes) {
                residual_value[static_cast<std::size_t>(node)] = 0.0;
            }
            for (int node : visited_nodes) {
                visited[static_cast<std::size_t>(node)] = 0U;
            }
            visited_nodes.clear();
        }
        round_report.selection_ms = milliseconds(
            selection_begin, Clock::now());
        result.report.selection_ms += round_report.selection_ms;
        if (round_report.expanded_columns == 0) {
            result.report.history.push_back(round_report);
            break;
        }

        const auto refinement_begin = Clock::now();
        InterpolationOptions warm_options = interpolation_options;
        warm_options.thread_count = options.thread_count;
        warm_options.local_tolerance = std::min(
            warm_options.local_tolerance,
            options.refinement_tolerance);
        const InterpolationResult refined =
            refine_selected_energy_interpolation_on_supports(
                grid, a, result.supports, refine_column,
                result.prolongation, warm_options);
        result.prolongation = refined.prolongation;
        round_report.refinement_ms = milliseconds(
            refinement_begin, Clock::now());
        result.report.refinement_ms += round_report.refinement_ms;
        result.report.history.push_back(round_report);
        result.report.rounds = round_index + 1;

        const auto residual_begin = Clock::now();
        snapshot = scaled_f_residual(
            grid, a, result.prolongation, result.supports,
            options.thread_count);
        result.report.residual_ms += milliseconds(
            residual_begin, Clock::now());
    }

    const auto [final_mean, final_maximum] =
        mean_and_maximum(snapshot.norm);
    result.report.final_mean_scaled_residual = final_mean;
    result.report.final_maximum_scaled_residual = final_maximum;
    double ratio_sum = 0.0;
    for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
        const std::size_t offset = static_cast<std::size_t>(coarse);
        if (ever_expanded[offset] != 0U) {
            ++result.report.expanded_columns;
        }
        result.report.total_extra_nodes += extras[offset];
        result.report.maximum_extra_nodes = std::max(
            result.report.maximum_extra_nodes, extras[offset]);
        const double target = std::max(
            options.absolute_residual_tolerance,
            options.target_residual_ratio * initial_norm[offset]);
        if (snapshot.norm[offset] <= target) {
            ++result.report.converged_columns;
        }
        ratio_sum += initial_norm[offset] > 0.0
            ? snapshot.norm[offset] / initial_norm[offset]
            : 0.0;
    }
    result.report.mean_final_to_initial_ratio = ratio_sum /
        static_cast<double>(grid.coarse_size());
    result.report.total_ms = milliseconds(total_begin, Clock::now());
    return result;
}

}
