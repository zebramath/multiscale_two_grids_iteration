#pragma once

#include "multigrid/global_pcg_path.hpp"
#include "multigrid/two_grid_solver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace tgi {

struct AdaptiveGlobalPcgOptions {
    int minimum_steps = 12;
    int maximum_steps = 64;
    int step_quantum = 2;
    int pilot_iterations = 24;
    int tail_window = 8;
    int maximum_candidate_hierarchies = 4;
    double target_forecast_fraction = 0.12;
    double early_accept_forecast_multiple = 2.0;
    double acceptable_cycle_slack = 0.08;
    double forecast_uncertainty_weight = 1.0;
    int maximum_cycles = 40000;
    int smoothing_steps = 1;
    int thread_count = 1;
    bool include_initial_candidate = true;
    double solve_tolerance = 1.0e-6;
    double drop_tolerance = 0.0;
};

struct AdaptiveGlobalPcgCheckpoint {
    int steps = 0;
    std::string phase;
    int pilot_iterations = 0;
    int predicted_cycles = 0;
    int confirmed_cycles = -1;
    double density_percent = 0.0;
    double rho_rhs_pilot = 0.0;
    double pilot_relative_residual = 0.0;
    double forecast_relative_uncertainty = 0.0;
    double maximum_pcg_residual = 0.0;
    double rms_pcg_residual = 0.0;
    double preconditioned_pcg_residual = 0.0;
    double path_ms = 0.0;
    double coarse_setup_ms = 0.0;
    double pilot_ms = 0.0;
    bool selected = false;
};

struct AdaptiveGlobalPcgReport {
    int selected_steps = 0;
    int selected_cycles = 0;
    int estimated_selected_cycles = 0;
    int screened_candidates = 0;
    bool early_accept = false;
    bool selected_cycles_confirmed = false;
    double selection_wall_ms = 0.0;
    std::vector<AdaptiveGlobalPcgCheckpoint> history;
};

struct AdaptiveGlobalPcgResult {
    SparseMatrix prolongation;
    AdaptiveGlobalPcgReport report;
};

namespace adaptive_global_pcg_detail {

using Clock = std::chrono::steady_clock;

inline double milliseconds(Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

inline double median(std::vector<double> values) {
    const std::size_t middle = values.size() / 2U;
    std::nth_element(values.begin(), values.begin() + middle, values.end());
    double value = values[middle];
    if (values.size() % 2U == 0U) {
        std::nth_element(
            values.begin(), values.begin() + middle - 1U, values.end());
        value = 0.5 * (value + values[middle - 1U]);
    }
    return value;
}

struct TailModel {
    double rho = 1.0;
    double relative_uncertainty = 1.0;
};

inline TailModel fit_tail(
    const std::vector<double>& norms, int tail_window,
    double uncertainty_weight) {
    const int completed = static_cast<int>(norms.size()) - 1;
    const int block = std::max(2, tail_window / 2);
    const int begin = std::max(block, completed - tail_window + 1);
    std::vector<double> rates;
    rates.reserve(static_cast<std::size_t>(completed - begin + 1));
    for (int iteration = begin; iteration <= completed; ++iteration) {
        const double current = std::max(
            norms[static_cast<std::size_t>(iteration)], 1.0e-300);
        const double previous = std::max(norms[
            static_cast<std::size_t>(iteration - block)], 1.0e-300);
        rates.push_back(
            std::log(current / previous) / static_cast<double>(block));
    }
    const double center = median(rates);
    std::vector<double> deviations;
    deviations.reserve(rates.size());
    for (double rate : rates) deviations.push_back(std::abs(rate - center));
    const double scale = 1.4826 * median(std::move(deviations));
    const double uncertainty = scale /
        std::sqrt(static_cast<double>(rates.size()));
    const double upper = *std::max_element(rates.begin(), rates.end());
    const double conservative = std::max(
        upper, center + uncertainty_weight * uncertainty);
    TailModel model;
    model.rho = std::exp(std::min(-1.0e-10, conservative));
    model.relative_uncertainty = std::min(
        2.0, uncertainty / std::max(std::abs(center), 1.0e-10));
    return model;
}

inline int forecast_cycles(
    int completed, double relative_residual, double rho,
    double tolerance, int maximum) {
    if (relative_residual <= tolerance) return completed;
    if (rho >= 1.0) return maximum;
    const double remaining = std::ceil(
        std::log(tolerance / relative_residual) / std::log(rho));
    const long long forecast = static_cast<long long>(completed) +
        std::max(0LL, static_cast<long long>(remaining));
    return static_cast<int>(std::min(
        static_cast<long long>(maximum), forecast));
}

struct Candidate {
    std::shared_ptr<SparseMatrix> prolongation;
    std::unique_ptr<TwoGridCycle> cycle;
    AdaptiveGlobalPcgCheckpoint checkpoint;
};

inline std::unique_ptr<Candidate> pilot_candidate(
    const SparseMatrix& a, const Vector& rhs, SparseMatrix prolongation,
    int steps, const std::string& phase,
    const AdaptiveGlobalPcgOptions& options,
    const GlobalPcgPathReport& path_report) {
    auto candidate = std::make_unique<Candidate>();
    candidate->prolongation = std::make_shared<SparseMatrix>(
        std::move(prolongation));
    auto& checkpoint = candidate->checkpoint;
    checkpoint.steps = steps;
    checkpoint.phase = phase;
    checkpoint.density_percent = 100.0 *
        static_cast<double>(candidate->prolongation->nnz()) /
        (static_cast<double>(candidate->prolongation->rows()) *
         candidate->prolongation->cols());
    checkpoint.path_ms = path_report.total_ms;
    checkpoint.maximum_pcg_residual =
        path_report.maximum_relative_residual;
    checkpoint.rms_pcg_residual = path_report.rms_relative_residual;
    checkpoint.preconditioned_pcg_residual =
        path_report.relative_preconditioned_residual;

    const auto begin = Clock::now();
    candidate->cycle = std::make_unique<TwoGridCycle>(
        a, *candidate->prolongation,
        options.smoothing_steps, options.thread_count);
    checkpoint.coarse_setup_ms = candidate->cycle->setup_report().total_ms;
    Vector solution(rhs.size(), 0.0);
    Vector residual = rhs;
    TwoGridCycle::Workspace workspace;
    const double initial = std::max(norm2(rhs), 1.0e-300);
    std::vector<double> norms;
    norms.reserve(static_cast<std::size_t>(options.pilot_iterations) + 1U);
    norms.push_back(initial);
    for (int iteration = 0; iteration < options.pilot_iterations; ++iteration) {
        const double squared = candidate->cycle->iterate(
            rhs, solution, residual, workspace);
        const double current = std::max(std::sqrt(squared), 1.0e-300);
        norms.push_back(current);
        checkpoint.pilot_iterations = iteration + 1;
        if (current / initial <= options.solve_tolerance) break;
    }
    checkpoint.pilot_ms = milliseconds(begin, Clock::now());
    checkpoint.pilot_relative_residual = norms.back() / initial;
    const TailModel model = fit_tail(
        norms, std::min(options.tail_window, checkpoint.pilot_iterations),
        options.forecast_uncertainty_weight);
    checkpoint.rho_rhs_pilot = model.rho;
    checkpoint.forecast_relative_uncertainty = model.relative_uncertainty;
    checkpoint.predicted_cycles = forecast_cycles(
        checkpoint.pilot_iterations,
        checkpoint.pilot_relative_residual,
        checkpoint.rho_rhs_pilot,
        options.solve_tolerance, options.maximum_cycles);
    return candidate;
}

inline int quantize(int steps, int quantum, int lower, int upper) {
    const int rounded = ((steps + quantum / 2) / quantum) * quantum;
    return std::max(lower, std::min(upper, rounded));
}

inline int projected_checkpoint(
    const Candidate& initial, const Candidate& current,
    const AdaptiveGlobalPcgOptions& options) {
    const double first = static_cast<double>(
        std::max(1, initial.checkpoint.predicted_cycles));
    const double second = static_cast<double>(
        std::max(1, current.checkpoint.predicted_cycles));
    const int span = std::max(1,
        current.checkpoint.steps - initial.checkpoint.steps);
    const double slope = std::log(second / first) /
        static_cast<double>(span);
    if (slope >= -1.0e-10) {
        return quantize(
            (current.checkpoint.steps + options.maximum_steps) / 2,
            options.step_quantum,
            current.checkpoint.steps + options.step_quantum,
            options.maximum_steps);
    }
    const double pilot_floor = options.early_accept_forecast_multiple *
        static_cast<double>(options.pilot_iterations);
    const double target = std::max(
        pilot_floor, options.target_forecast_fraction * second);
    const double raw = static_cast<double>(current.checkpoint.steps) +
        std::log(target / second) / slope;
    const int lower = std::min(
        options.maximum_steps,
        current.checkpoint.steps + options.step_quantum);
    const int midpoint = quantize(
        (current.checkpoint.steps + options.maximum_steps) / 2,
        options.step_quantum, lower, options.maximum_steps);
    return quantize(
        static_cast<int>(std::ceil(raw)), options.step_quantum,
        lower, midpoint);
}

inline std::size_t select_candidate(
    const std::vector<std::unique_ptr<Candidate>>& candidates,
    double slack) {
    int best = std::numeric_limits<int>::max();
    for (const auto& candidate : candidates) {
        best = std::min(best, candidate->checkpoint.predicted_cycles);
    }
    if (best == candidates.front()->checkpoint.predicted_cycles &&
        std::all_of(
            candidates.begin(), candidates.end(),
            [&](const auto& candidate) {
                return candidate->checkpoint.predicted_cycles == best;
            })) {
        return static_cast<std::size_t>(std::min_element(
            candidates.begin(), candidates.end(),
            [](const auto& lhs, const auto& rhs) {
                if (lhs->checkpoint.preconditioned_pcg_residual !=
                    rhs->checkpoint.preconditioned_pcg_residual) {
                    return lhs->checkpoint.preconditioned_pcg_residual <
                        rhs->checkpoint.preconditioned_pcg_residual;
                }
                if (lhs->checkpoint.rho_rhs_pilot !=
                    rhs->checkpoint.rho_rhs_pilot) {
                    return lhs->checkpoint.rho_rhs_pilot <
                        rhs->checkpoint.rho_rhs_pilot;
                }
                return lhs->checkpoint.pilot_relative_residual <
                    rhs->checkpoint.pilot_relative_residual;
            }) - candidates.begin());
    }
    const double limit = (1.0 + slack) * static_cast<double>(best);
    std::size_t selected = 0;
    bool found = false;
    for (std::size_t index = 0; index < candidates.size(); ++index) {
        if (static_cast<double>(
                candidates[index]->checkpoint.predicted_cycles) <= limit &&
            (!found || candidates[index]->checkpoint.steps <
                candidates[selected]->checkpoint.steps)) {
            selected = index;
            found = true;
        }
    }
    return selected;
}

}

inline AdaptiveGlobalPcgResult build_adaptive_global_pcg_interpolation(
    const StructuredGrid& grid, const SparseMatrix& a,
    const SparseMatrix& initial_prolongation,
    const AdaptiveGlobalPcgOptions& options,
    const Vector* representative_rhs) {
    using namespace adaptive_global_pcg_detail;
    const auto begin = Clock::now();
    const Vector& rhs = *representative_rhs;
    std::vector<std::unique_ptr<Candidate>> candidates;
    GlobalPcgPathReport initial_report;
    initial_report.maximum_relative_residual = 1.0;
    initial_report.rms_relative_residual = 1.0;
    initial_report.relative_preconditioned_residual = 1.0;
    if (options.include_initial_candidate) {
        candidates.push_back(pilot_candidate(
            a, rhs, initial_prolongation, 0, "initial",
            options, initial_report));
    }

    AdaptiveGlobalPcgResult result;
    const double easy_limit = options.early_accept_forecast_multiple *
        static_cast<double>(options.pilot_iterations);
    if (!candidates.empty() &&
        candidates.front()->checkpoint.predicted_cycles <= easy_limit) {
        result.report.early_accept = true;
    } else {
        GlobalEnergyPcgPath path(
            grid, a, initial_prolongation, options.thread_count);
        path.advance_to(options.minimum_steps);
        GlobalPcgPathReport path_report = path.report();
        candidates.push_back(pilot_candidate(
            a, rhs, path.prolongation(options.drop_tolerance),
            options.minimum_steps, "anchor", options, path_report));

        const bool anchor_is_easy =
            candidates.back()->checkpoint.predicted_cycles <= easy_limit;
        if (anchor_is_easy) {
            result.report.early_accept = true;
        } else if (static_cast<int>(candidates.size()) <
                       options.maximum_candidate_hierarchies) {
            const int projected = projected_checkpoint(
                *candidates.front(), *candidates.back(), options);
            if (projected > path.steps()) {
                const int anchor_steps = path.steps();
                const int refinement_steps = quantize(
                    (anchor_steps + projected) / 2,
                    options.step_quantum,
                    anchor_steps + options.step_quantum,
                    projected);
                SparseMatrix refinement_prolongation;
                GlobalPcgPathReport refinement_report;
                if (refinement_steps < projected) {
                    path.advance_to(refinement_steps);
                    refinement_report = path.report();
                    refinement_prolongation = path.prolongation(
                        options.drop_tolerance);
                }
                path.advance_to(projected);
                path_report = path.report();
                candidates.push_back(pilot_candidate(
                    a, rhs, path.prolongation(options.drop_tolerance),
                    projected, "projected", options, path_report));
                const auto& anchor = *candidates[candidates.size() - 2U];
                const auto& projection = *candidates.back();
                const double energy_ratio =
                    projection.checkpoint.preconditioned_pcg_residual /
                    std::max(
                        anchor.checkpoint.preconditioned_pcg_residual,
                        1.0e-300);
                if (static_cast<int>(candidates.size()) <
                    options.maximum_candidate_hierarchies) {
                    if (projection.checkpoint.predicted_cycles >
                        4.0 * easy_limit && projected < options.maximum_steps) {
                        const int forward = quantize(
                            (projected + options.maximum_steps) / 2,
                            options.step_quantum,
                            projected + options.step_quantum,
                            options.maximum_steps);
                        path.advance_to(forward);
                        path_report = path.report();
                        candidates.push_back(pilot_candidate(
                            a, rhs,
                            path.prolongation(options.drop_tolerance),
                            forward, "forward", options, path_report));
                    } else if (energy_ratio < 0.02 &&
                               refinement_prolongation.rows() > 0) {
                        candidates.push_back(pilot_candidate(
                            a, rhs, std::move(refinement_prolongation),
                            refinement_steps, "refine", options,
                            refinement_report));
                    }
                }
            }
        }
    }

    const std::size_t selected = select_candidate(
        candidates, options.acceptable_cycle_slack);
    candidates[selected]->checkpoint.selected = true;
    result.prolongation = *candidates[selected]->prolongation;
    result.report.selected_steps = candidates[selected]->checkpoint.steps;
    result.report.selected_cycles =
        candidates[selected]->checkpoint.predicted_cycles;
    result.report.estimated_selected_cycles = result.report.selected_cycles;
    result.report.screened_candidates = static_cast<int>(candidates.size());
    for (const auto& candidate : candidates) {
        result.report.history.push_back(candidate->checkpoint);
    }
    result.report.selection_wall_ms = milliseconds(begin, Clock::now());
    return result;
}

}
