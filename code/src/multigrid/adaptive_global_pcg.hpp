#pragma once

#include "multigrid/global_pcg_path.hpp"
#include "multigrid/two_grid_solver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace tgi {

struct AdaptiveGlobalPcgOptions {
    int minimum_steps = 16;
    int maximum_steps = 64;
    int step_increment = 4;
    int patience = 3;
    int probe_count = 3;
    int power_iterations = 30;
    int rhs_pilot_iterations = 32;
    int rhs_tail_window = 8;
    int smoothing_steps = 1;
    int thread_count = 1;
    int expected_rhs = 8;
    int maximum_predicted_cycles = 40000;
    bool include_initial_candidate = true;
    double solve_tolerance = 1.0e-6;
    double minimum_relative_improvement = 0.005;
    double drop_tolerance = 0.0;
    std::uint64_t probe_seed = 0x51adbeefULL;
};

struct AdaptiveGlobalPcgCheckpoint {
    int steps = 0;
    int predicted_cycles = 0;
    double density_percent = 0.0;
    double rho_hat = 0.0;
    double rho_power = 0.0;
    double rho_rhs_pilot = 0.0;
    double maximum_pcg_residual = 0.0;
    double path_ms = 0.0;
    double coarse_setup_ms = 0.0;
    double probe_ms = 0.0;
    double estimated_cycle_ms = 0.0;
    double predicted_total_ms = 0.0;
    bool improved = false;
};

struct AdaptiveGlobalPcgReport {
    int selected_steps = 0;
    int evaluated_checkpoints = 0;
    bool stopped_by_patience = false;
    double selected_rho_hat = 0.0;
    double selected_score_ms = 0.0;
    double selection_wall_ms = 0.0;
    std::vector<AdaptiveGlobalPcgCheckpoint> history;
};

struct AdaptiveGlobalPcgResult {
    SparseMatrix prolongation;
    AdaptiveGlobalPcgReport report;
};

inline int predicted_stationary_iterations(
    double rho, double tolerance, int maximum) {
    if (!(rho >= 0.0) || !std::isfinite(rho) || rho >= 1.0) {
        return maximum;
    }
    if (rho == 0.0) return 1;
    const double count = std::ceil(std::log(tolerance) / std::log(rho));
    if (!std::isfinite(count)) return maximum;
    return std::max(1, std::min(maximum, static_cast<int>(count)));
}

inline AdaptiveGlobalPcgResult build_adaptive_global_pcg_interpolation(
    const StructuredGrid& grid, const SparseMatrix& a,
    const SparseMatrix& initial_prolongation,
    const AdaptiveGlobalPcgOptions& options = {},
    const Vector* representative_rhs = nullptr) {
    if (options.minimum_steps < 0 ||
        options.maximum_steps < options.minimum_steps ||
        options.step_increment <= 0 || options.patience <= 0 ||
        options.probe_count <= 0 || options.power_iterations <= 0 ||
        options.rhs_pilot_iterations < 0 || options.rhs_tail_window <= 0 ||
        (options.rhs_pilot_iterations > 0 &&
         options.rhs_tail_window > options.rhs_pilot_iterations) ||
        options.smoothing_steps < 0 || options.thread_count <= 0 ||
        options.expected_rhs <= 0 || options.maximum_predicted_cycles <= 0 ||
        !(options.solve_tolerance > 0.0 && options.solve_tolerance < 1.0) ||
        !(options.minimum_relative_improvement >= 0.0 &&
          options.minimum_relative_improvement < 1.0) ||
        !(options.drop_tolerance >= 0.0)) {
        throw std::invalid_argument(
            "build_adaptive_global_pcg_interpolation: invalid options");
    }
    if (representative_rhs != nullptr &&
        representative_rhs->size() != static_cast<std::size_t>(a.rows())) {
        throw std::invalid_argument(
            "build_adaptive_global_pcg_interpolation: wrong pilot RHS size");
    }

    using Clock = std::chrono::steady_clock;
    const auto selection_begin = Clock::now();
    GlobalEnergyPcgPath path(
        grid, a, initial_prolongation, options.thread_count);
    AdaptiveGlobalPcgResult result;
    double best_score = std::numeric_limits<double>::infinity();
    int stale = 0;

    std::vector<int> checkpoints;
    if (options.include_initial_candidate && options.minimum_steps > 0) {
        checkpoints.push_back(0);
    }
    for (int steps = options.minimum_steps;
         steps <= options.maximum_steps; steps += options.step_increment) {
        checkpoints.push_back(steps);
    }
    for (int steps : checkpoints) {
        path.advance_to(steps);
        SparseMatrix candidate = path.prolongation(options.drop_tolerance);
        const GlobalPcgPathReport path_report = path.report();

        const auto probe_begin = Clock::now();
        const TwoGridCycle cycle(
            a, candidate, options.smoothing_steps, options.thread_count);
        double rho_power = 0.0;
        for (int probe = 0; probe < options.probe_count; ++probe) {
            const std::uint64_t seed = options.probe_seed +
                0x9e3779b97f4a7c15ULL *
                    static_cast<std::uint64_t>(probe + 1);
            rho_power = std::max(
                rho_power,
                cycle.estimate_convergence_factor(
                    options.power_iterations, seed));
        }
        double rho_rhs = 0.0;
        if (representative_rhs != nullptr &&
            options.rhs_pilot_iterations > 0) {
            Vector solution(representative_rhs->size(), 0.0);
            Vector residual = *representative_rhs;
            TwoGridCycle::Workspace workspace;
            std::vector<double> norms;
            norms.reserve(static_cast<std::size_t>(
                options.rhs_pilot_iterations + 1));
            norms.push_back(std::max(norm2(residual), 1.0e-300));
            for (int iteration = 0;
                 iteration < options.rhs_pilot_iterations; ++iteration) {
                const double squared = cycle.iterate(
                    *representative_rhs, solution, residual, workspace);
                norms.push_back(std::max(std::sqrt(squared), 1.0e-300));
            }
            const int end = options.rhs_pilot_iterations;
            const int window = options.rhs_tail_window;
            // Use overlapping tail windows. This is less sensitive than one
            // residual ratio but still exposes slowly emerging channel modes.
            for (int finish = end - window + 1; finish <= end; ++finish) {
                if (finish < window) continue;
                const double ratio = norms[static_cast<std::size_t>(finish)] /
                    norms[static_cast<std::size_t>(finish - window)];
                rho_rhs = std::max(
                    rho_rhs,
                    std::pow(ratio, 1.0 / static_cast<double>(window)));
            }
        }
        const double rho_hat = std::max(rho_power, rho_rhs);
        const double probe_ms =
            std::chrono::duration<double, std::milli>(
                Clock::now() - probe_begin).count();
        const int measured_applications =
            options.probe_count * options.power_iterations +
            (representative_rhs != nullptr
                 ? options.rhs_pilot_iterations : 0);
        const double cycle_ms = probe_ms /
            static_cast<double>(
                std::max(1, measured_applications));
        const int predicted_cycles = predicted_stationary_iterations(
            rho_hat, options.solve_tolerance,
            options.maximum_predicted_cycles);
        // This is the predicted cost of choosing this checkpoint a priori.
        // Probe work is reported as selection overhead, not charged repeatedly
        // to every future right-hand side.
        const double charged_path_ms = steps == 0 ? 0.0 : path_report.total_ms;
        const double score = charged_path_ms +
            cycle.setup_report().total_ms +
            static_cast<double>(options.expected_rhs * predicted_cycles) *
                cycle_ms;
        AdaptiveGlobalPcgCheckpoint checkpoint;
        checkpoint.steps = steps;
        checkpoint.predicted_cycles = predicted_cycles;
        checkpoint.density_percent =
            100.0 * static_cast<double>(candidate.nnz()) /
            (static_cast<double>(candidate.rows()) * candidate.cols());
        checkpoint.rho_hat = rho_hat;
        checkpoint.rho_power = rho_power;
        checkpoint.rho_rhs_pilot = rho_rhs;
        checkpoint.maximum_pcg_residual =
            path_report.maximum_relative_residual;
        checkpoint.path_ms = charged_path_ms;
        checkpoint.coarse_setup_ms = cycle.setup_report().total_ms;
        checkpoint.probe_ms = probe_ms;
        checkpoint.estimated_cycle_ms = cycle_ms;
        checkpoint.predicted_total_ms = score;

        const double improvement_threshold = std::isfinite(best_score)
            ? best_score * (1.0 - options.minimum_relative_improvement)
            : best_score;
        if (!std::isfinite(best_score) || score < improvement_threshold) {
            best_score = score;
            result.prolongation = candidate;
            result.report.selected_steps = steps;
            result.report.selected_rho_hat = rho_hat;
            result.report.selected_score_ms = score;
            checkpoint.improved = true;
            stale = 0;
        } else {
            ++stale;
        }
        result.report.history.push_back(checkpoint);
        if (stale >= options.patience &&
            steps < options.maximum_steps) {
            result.report.stopped_by_patience = true;
            break;
        }
    }
    result.report.evaluated_checkpoints =
        static_cast<int>(result.report.history.size());
    result.report.selection_wall_ms =
        std::chrono::duration<double, std::milli>(
            Clock::now() - selection_begin).count();
    if (result.prolongation.rows() == 0) {
        throw std::runtime_error(
            "adaptive global PCG did not produce a candidate");
    }
    return result;
}

} // namespace tgi
