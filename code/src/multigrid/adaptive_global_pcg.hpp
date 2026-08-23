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
    double early_accept_forecast_multiple = 2.0;
    double forward_refinement_forecast_multiple = 4.0;
    double acceptable_cycle_slack = 0.08;
    int maximum_cycles = 40000;
    int smoothing_steps = 1;
    int thread_count = 1;
    double solve_tolerance = 1.0e-6;
    double drop_tolerance = 0.0;
};

struct AdaptiveGlobalPcgCheckpoint {
    int steps = 0;
    std::string phase;
    int pilot_iterations = 0;
    int predicted_cycles = 0;
    double density_percent = 0.0;
    double rho_rhs_pilot = 0.0;
    double pilot_relative_residual = 0.0;
    double forecast_relative_uncertainty = 0.0;
    double preconditioned_pcg_residual = 0.0;
    double path_ms = 0.0;
    double pilot_ms = 0.0;
    bool selected = false;
};

struct AdaptiveGlobalPcgReport {
    int selected_steps = 0;
    int estimated_selected_cycles = 0;
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
    const std::vector<double>& norms, int tail_window) {
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
    const double conservative = std::max(upper, center + uncertainty);
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
    checkpoint.preconditioned_pcg_residual =
        path_report.relative_preconditioned_residual;

    const auto begin = Clock::now();
    candidate->cycle = std::make_unique<TwoGridCycle>(
        a, *candidate->prolongation,
        options.smoothing_steps, options.thread_count);
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
        norms, std::min(options.tail_window, checkpoint.pilot_iterations));
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
                return lhs->checkpoint.preconditioned_pcg_residual <
                    rhs->checkpoint.preconditioned_pcg_residual;
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
    initial_report.relative_preconditioned_residual = 1.0;
    candidates.push_back(pilot_candidate(
        a, rhs, initial_prolongation, 0, "initial",
        options, initial_report));

    AdaptiveGlobalPcgResult result;
    const double easy_limit = options.early_accept_forecast_multiple *
        static_cast<double>(options.pilot_iterations);
    if (candidates.front()->checkpoint.predicted_cycles > easy_limit) {
        GlobalEnergyPcgPath path(
            grid, a, initial_prolongation, options.thread_count);
        path.advance_to(options.minimum_steps);
        GlobalPcgPathReport path_report = path.report();
        candidates.push_back(pilot_candidate(
            a, rhs, path.prolongation(options.drop_tolerance),
            options.minimum_steps, "anchor", options, path_report));

        const bool anchor_is_easy =
            candidates.back()->checkpoint.predicted_cycles <= easy_limit;
        if (!anchor_is_easy) {
            const int midpoint = quantize(
                (options.minimum_steps + options.maximum_steps) / 2,
                options.step_quantum,
                options.minimum_steps + options.step_quantum,
                options.maximum_steps);
            const int backward = quantize(
                (options.minimum_steps + midpoint) / 2,
                options.step_quantum,
                options.minimum_steps + options.step_quantum,
                midpoint);
            SparseMatrix backward_prolongation;
            GlobalPcgPathReport backward_report;
            if (backward < midpoint) {
                path.advance_to(backward);
                backward_report = path.report();
                backward_prolongation = path.prolongation(
                    options.drop_tolerance);
            }
            path.advance_to(midpoint);
            path_report = path.report();
            candidates.push_back(pilot_candidate(
                a, rhs, path.prolongation(options.drop_tolerance),
                midpoint, "midpoint", options, path_report));

            const int midpoint_forecast =
                candidates.back()->checkpoint.predicted_cycles;
            if (midpoint_forecast > easy_limit) {
                if (midpoint_forecast >
                        easy_limit + options.pilot_iterations &&
                    candidates.back()->checkpoint
                            .preconditioned_pcg_residual <=
                        std::sqrt(options.solve_tolerance) &&
                    backward_prolongation.rows() > 0) {
                    candidates.push_back(pilot_candidate(
                        a, rhs, std::move(backward_prolongation),
                        backward, "refine-backward", options,
                        backward_report));
                } else if (midpoint_forecast >
                        options.forward_refinement_forecast_multiple *
                            easy_limit &&
                    midpoint < options.maximum_steps) {
                    const int forward = quantize(
                        (midpoint + options.maximum_steps) / 2,
                        options.step_quantum,
                        midpoint + options.step_quantum,
                        options.maximum_steps);
                    path.advance_to(forward);
                    path_report = path.report();
                    candidates.push_back(pilot_candidate(
                        a, rhs, path.prolongation(options.drop_tolerance),
                        forward, "refine-forward", options, path_report));
                }
            }
        }
    }

    const std::size_t selected = select_candidate(
        candidates, options.acceptable_cycle_slack);
    candidates[selected]->checkpoint.selected = true;
    result.prolongation = *candidates[selected]->prolongation;
    result.report.selected_steps = candidates[selected]->checkpoint.steps;
    result.report.estimated_selected_cycles =
        candidates[selected]->checkpoint.predicted_cycles;
    for (const auto& candidate : candidates) {
        result.report.history.push_back(candidate->checkpoint);
    }
    result.report.selection_wall_ms = milliseconds(begin, Clock::now());
    return result;
}

}
