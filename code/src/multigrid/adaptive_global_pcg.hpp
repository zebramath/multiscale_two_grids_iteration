#pragma once

#include "multigrid/global_pcg_path.hpp"
#include "multigrid/two_grid_solver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace tgi {

struct AdaptiveGlobalPcgOptions {
    // v3.2 default: screen at most two positive PCG checkpoints and accept
    // a candidate whose forecast is within the requested cycle slack.  The
    // older screen-refine-confirm policy remains available for diagnostics.
    bool cost_aware_mode = true;
    int budget_pilot_iterations = 16;
    int budget_tail_window = 6;
    bool scale_budget_pilot_with_grid_ratio = true;
    int budget_middle_offset_steps = 12;
    int budget_far_steps = 40;
    bool skip_middle_checkpoint_on_large_grids = true;
    int budget_easy_forecast_cycles = 48;
    double acceptable_cycle_slack = 0.05;
    int minimum_steps = 16;
    int maximum_steps = 64;
    int maximum_screening_steps = 48;
    int screening_increment = 16;
    int screening_pilot_iterations = 32;
    int screening_tail_window = 8;
    int screening_patience = 2;
    int minimum_screened_positive_candidates = 5;
    int refinement_backtrack_steps = 10;
    int refinement_stop_before_anchor_steps = 4;
    int refinement_increment = 2;
    int refinement_pilot_iterations = 64;
    int refinement_tail_window = 16;
    int confirmation_candidates = 2;
    int easy_accept_cycles = 48;
    int medium_accept_cycles = 64;
    int initial_safety_forecast_cycles = 112;
    int maximum_confirmation_cycles = 40000;
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
    double maximum_pcg_residual = 0.0;
    double path_ms = 0.0;
    double coarse_setup_ms = 0.0;
    double pilot_ms = 0.0;
    double confirmation_ms = 0.0;
    bool improved = false;
    bool selected = false;
};

struct AdaptiveGlobalPcgReport {
    int selected_steps = 0;
    int selected_cycles = 0;
    int screening_anchor_steps = 0;
    int screened_candidates = 0;
    int refined_candidates = 0;
    int confirmed_candidates = 0;
    bool quick_accepted_initial = false;
    bool budget_mode = false;
    bool budget_early_accept = false;
    bool selected_cycles_confirmed = true;
    int estimated_selected_cycles = 0;
    bool screening_stopped_early = false;
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

inline int forecast_total_cycles(
    int completed, double relative_residual, double rho,
    double tolerance, int maximum) {
    if (relative_residual <= tolerance) return completed;
    if (!(rho > 0.0) || !std::isfinite(rho) || rho >= 1.0) return maximum;
    const double remaining = std::ceil(
        std::log(tolerance / relative_residual) / std::log(rho));
    if (!std::isfinite(remaining)) return maximum;
    const int count = completed + std::max(0, static_cast<int>(remaining));
    return std::max(completed, std::min(maximum, count));
}

struct CandidateRecord {
    std::shared_ptr<SparseMatrix> prolongation;
    AdaptiveGlobalPcgCheckpoint checkpoint;
    std::unique_ptr<TwoGridCycle> cycle;
    Vector solution;
    Vector residual;
    double initial_residual_norm = 1.0;
};

struct PcgSnapshot {
    int steps = 0;
    SparseMatrix prolongation;
    double path_ms = 0.0;
    double maximum_pcg_residual = 0.0;
};

inline CandidateRecord run_pilot(
    const SparseMatrix& a, const Vector& rhs, SparseMatrix prolongation,
    int steps, const std::string& phase, int pilot_iterations,
    int tail_window, int smoothing_steps, int thread_count,
    double tolerance, int maximum_cycles, double path_ms,
    double maximum_pcg_residual) {
    CandidateRecord record;
    record.prolongation = std::make_shared<SparseMatrix>(
        std::move(prolongation));
    record.checkpoint.steps = steps;
    record.checkpoint.phase = phase;
    record.checkpoint.density_percent =
        100.0 * static_cast<double>(record.prolongation->nnz()) /
        (static_cast<double>(record.prolongation->rows()) *
         record.prolongation->cols());
    record.checkpoint.path_ms = path_ms;
    record.checkpoint.maximum_pcg_residual = maximum_pcg_residual;
    record.solution.assign(rhs.size(), 0.0);
    record.residual = rhs;
    record.initial_residual_norm = std::max(norm2(rhs), 1.0e-300);
    std::vector<double> norms;
    norms.reserve(static_cast<std::size_t>(pilot_iterations) + 1U);
    norms.push_back(record.initial_residual_norm);

    const auto begin = Clock::now();
    record.cycle = std::make_unique<TwoGridCycle>(
        a, *record.prolongation, smoothing_steps, thread_count);
    record.checkpoint.coarse_setup_ms = record.cycle->setup_report().total_ms;
    TwoGridCycle::Workspace workspace;
    for (int iteration = 0; iteration < pilot_iterations; ++iteration) {
        const double squared = record.cycle->iterate(
            rhs, record.solution, record.residual, workspace);
        const double current = std::max(std::sqrt(squared), 1.0e-300);
        norms.push_back(current);
        record.checkpoint.pilot_iterations = iteration + 1;
        if (current / record.initial_residual_norm <= tolerance) break;
    }
    record.checkpoint.pilot_ms = milliseconds(begin, Clock::now());
    const int completed = record.checkpoint.pilot_iterations;
    const double relative = norms.back() / record.initial_residual_norm;
    record.checkpoint.pilot_relative_residual = relative;
    double rho = 0.0;
    if (completed >= tail_window) {
        const int first_finish = std::max(
            tail_window, completed - tail_window + 1);
        const int finish_count = completed - first_finish + 1;
        for (int offset = 0; offset < finish_count; ++offset) {
            const int finish = first_finish + offset;
            const double ratio = norms[static_cast<std::size_t>(finish)] /
                norms[static_cast<std::size_t>(finish - tail_window)];
            rho = std::max(
                rho, std::pow(ratio, 1.0 / static_cast<double>(tail_window)));
        }
    }
    record.checkpoint.rho_rhs_pilot = rho;
    record.checkpoint.predicted_cycles = forecast_total_cycles(
        completed, relative, rho, tolerance, maximum_cycles);
    return record;
}

inline void confirm_candidate(
    const SparseMatrix& a, const Vector& rhs, CandidateRecord& record,
    int smoothing_steps, int thread_count, double tolerance,
    int maximum_cycles) {
    const auto begin = Clock::now();
    if (!record.cycle) {
        record.cycle = std::make_unique<TwoGridCycle>(
            a, *record.prolongation, smoothing_steps, thread_count);
    }
    TwoGridCycle::Workspace workspace;
    int completed = record.checkpoint.pilot_iterations;
    double relative = record.checkpoint.pilot_relative_residual;
    while (completed < maximum_cycles && relative > tolerance) {
        const double squared = record.cycle->iterate(
            rhs, record.solution, record.residual, workspace);
        ++completed;
        relative = std::sqrt(squared) / record.initial_residual_norm;
    }
    record.checkpoint.confirmed_cycles = completed;
    record.checkpoint.confirmation_ms = milliseconds(begin, Clock::now());
}

inline bool forecast_less(
    const std::unique_ptr<CandidateRecord>& lhs,
    const std::unique_ptr<CandidateRecord>& rhs) {
    if (lhs->checkpoint.predicted_cycles !=
        rhs->checkpoint.predicted_cycles) {
        return lhs->checkpoint.predicted_cycles <
            rhs->checkpoint.predicted_cycles;
    }
    return lhs->checkpoint.steps < rhs->checkpoint.steps;
}

inline AdaptiveGlobalPcgResult build_cost_aware_interpolation(
    const StructuredGrid& grid, const SparseMatrix& a,
    const SparseMatrix& initial_prolongation,
    const AdaptiveGlobalPcgOptions& options,
    const Vector& representative_rhs) {
    const auto selection_begin = Clock::now();
    AdaptiveGlobalPcgResult result;
    result.report.budget_mode = true;
    result.report.selected_cycles_confirmed = false;

    std::vector<std::unique_ptr<CandidateRecord>> candidates;
    const long long requested_pilot_iterations =
        static_cast<long long>(options.budget_pilot_iterations) +
        (options.scale_budget_pilot_with_grid_ratio
             ? 4LL * std::max(0, grid.ratio() - 4)
             : 0LL);
    const int pilot_iterations = static_cast<int>(std::min(
        static_cast<long long>(options.maximum_confirmation_cycles),
        requested_pilot_iterations));
    const int tail_window = std::min(
        pilot_iterations,
        options.budget_tail_window +
            (options.scale_budget_pilot_with_grid_ratio
                 ? std::max(0, grid.ratio() - 4) / 2
                 : 0));
    auto add_candidate = [&](SparseMatrix prolongation, int steps,
                             double path_ms, double maximum_pcg_residual) {
        auto candidate = std::make_unique<CandidateRecord>(run_pilot(
            a, representative_rhs, std::move(prolongation), steps,
            "budget", pilot_iterations,
            tail_window, options.smoothing_steps,
            options.thread_count, options.solve_tolerance,
            options.maximum_confirmation_cycles, path_ms,
            maximum_pcg_residual));
        ++result.report.screened_candidates;
        candidates.push_back(std::move(candidate));
    };

    if (options.include_initial_candidate) {
        add_candidate(initial_prolongation, 0, 0.0, 1.0);
        if (candidates.back()->checkpoint.predicted_cycles <=
            options.budget_easy_forecast_cycles) {
            candidates.back()->checkpoint.selected = true;
            result.prolongation = *candidates.back()->prolongation;
            result.report.selected_steps = 0;
            result.report.selected_cycles =
                candidates.back()->checkpoint.predicted_cycles;
            result.report.estimated_selected_cycles =
                result.report.selected_cycles;
            result.report.quick_accepted_initial = true;
            result.report.budget_early_accept = true;
            result.report.history.push_back(candidates.back()->checkpoint);
            result.report.selection_wall_ms = milliseconds(
                selection_begin, Clock::now());
            return result;
        }
    }

    GlobalEnergyPcgPath path(
        grid, a, initial_prolongation, options.thread_count);
    const int near_steps = options.minimum_steps;
    const int middle_steps = static_cast<int>(std::min(
        static_cast<long long>(options.maximum_steps),
        static_cast<long long>(near_steps) +
            options.budget_middle_offset_steps));
    const int far_steps = std::max(
        middle_steps,
        std::min(options.maximum_steps, options.budget_far_steps));
    const int checkpoints[3] = {near_steps, middle_steps, far_steps};
    int previous_steps = -1;
    GlobalPcgPathReport path_report;
    for (int checkpoint_steps : checkpoints) {
        if (checkpoint_steps == middle_steps &&
            options.skip_middle_checkpoint_on_large_grids &&
            grid.intervals() >= 128) {
            continue;
        }
        if (checkpoint_steps == previous_steps) continue;
        previous_steps = checkpoint_steps;
        path.advance_to(checkpoint_steps);
        path_report = path.report();
        add_candidate(
            path.prolongation(options.drop_tolerance), checkpoint_steps,
            path_report.total_ms, path_report.maximum_relative_residual);
        if (checkpoint_steps != far_steps &&
            candidates.back()->checkpoint.predicted_cycles <=
                options.budget_easy_forecast_cycles) {
            candidates.back()->checkpoint.selected = true;
            result.prolongation = *candidates.back()->prolongation;
            result.report.selected_steps = checkpoint_steps;
            result.report.selected_cycles =
                candidates.back()->checkpoint.predicted_cycles;
            result.report.estimated_selected_cycles =
                result.report.selected_cycles;
            result.report.screening_anchor_steps = checkpoint_steps;
            result.report.budget_early_accept = true;
            for (const auto& candidate : candidates) {
                result.report.history.push_back(candidate->checkpoint);
            }
            result.report.selection_wall_ms = milliseconds(
                selection_begin, Clock::now());
            return result;
        }
    }

    int best_forecast = options.maximum_confirmation_cycles;
    for (const auto& candidate : candidates) {
        best_forecast = std::min(
            best_forecast, candidate->checkpoint.predicted_cycles);
    }
    const double allowed = (1.0 + options.acceptable_cycle_slack) *
        static_cast<double>(best_forecast);
    std::size_t selected = candidates.size() - 1U;
    for (std::size_t index = 0; index < candidates.size(); ++index) {
        if (static_cast<double>(
                candidates[index]->checkpoint.predicted_cycles) <= allowed) {
            selected = index;
            break;
        }
    }
    candidates[selected]->checkpoint.selected = true;
    result.prolongation = *candidates[selected]->prolongation;
    result.report.selected_steps = candidates[selected]->checkpoint.steps;
    result.report.selected_cycles =
        candidates[selected]->checkpoint.predicted_cycles;
    result.report.estimated_selected_cycles = result.report.selected_cycles;
    result.report.screening_anchor_steps = result.report.selected_steps;
    for (const auto& candidate : candidates) {
        result.report.history.push_back(candidate->checkpoint);
    }
    result.report.selection_wall_ms = milliseconds(
        selection_begin, Clock::now());
    return result;
}

} // namespace adaptive_global_pcg_detail

inline AdaptiveGlobalPcgResult build_adaptive_global_pcg_interpolation(
    const StructuredGrid& grid, const SparseMatrix& a,
    const SparseMatrix& initial_prolongation,
    const AdaptiveGlobalPcgOptions& options,
    const Vector* representative_rhs) {
    if (representative_rhs == nullptr ||
        representative_rhs->size() != static_cast<std::size_t>(a.rows())) {
        throw std::invalid_argument(
            "staged adaptive PCG requires a representative RHS");
    }
    if (options.cost_aware_mode) {
        if (options.minimum_steps < 0 ||
            options.maximum_steps < options.minimum_steps ||
            options.maximum_confirmation_cycles <= 0 ||
            options.smoothing_steps < 0 || options.thread_count <= 0 ||
            !(options.solve_tolerance > 0.0 &&
              options.solve_tolerance < 1.0) ||
            !(options.drop_tolerance >= 0.0) ||
            options.budget_pilot_iterations <= 0 ||
            options.budget_tail_window <= 0 ||
            options.budget_tail_window > options.budget_pilot_iterations ||
            options.budget_middle_offset_steps <= 0 ||
            options.budget_far_steps <= options.minimum_steps ||
            options.budget_easy_forecast_cycles <= 0 ||
            !(options.acceptable_cycle_slack >= 0.0)) {
            throw std::invalid_argument(
                "build_adaptive_global_pcg_interpolation: invalid budget options");
        }
        return adaptive_global_pcg_detail::build_cost_aware_interpolation(
            grid, a, initial_prolongation, options, *representative_rhs);
    }
    if (options.minimum_steps < 0 ||
        options.maximum_steps < options.minimum_steps ||
        options.maximum_screening_steps < options.minimum_steps ||
        options.maximum_screening_steps > options.maximum_steps ||
        options.screening_increment <= 0 ||
        options.screening_pilot_iterations <= 0 ||
        options.screening_tail_window <= 0 ||
        options.screening_tail_window > options.screening_pilot_iterations ||
        options.screening_patience <= 0 ||
        options.minimum_screened_positive_candidates <= 0 ||
        options.refinement_backtrack_steps < 0 ||
        options.refinement_stop_before_anchor_steps < 0 ||
        options.refinement_increment <= 0 ||
        options.refinement_pilot_iterations <= 0 ||
        options.refinement_tail_window <= 0 ||
        options.refinement_tail_window > options.refinement_pilot_iterations ||
        options.confirmation_candidates <= 0 ||
        options.easy_accept_cycles <= 0 ||
        options.medium_accept_cycles < options.easy_accept_cycles ||
        options.initial_safety_forecast_cycles < options.easy_accept_cycles ||
        options.maximum_confirmation_cycles <= 0 ||
        options.smoothing_steps < 0 || options.thread_count <= 0 ||
        !(options.solve_tolerance > 0.0 && options.solve_tolerance < 1.0) ||
        !(options.drop_tolerance >= 0.0)) {
        throw std::invalid_argument(
            "build_adaptive_global_pcg_interpolation: invalid options");
    }

    using namespace adaptive_global_pcg_detail;
    const auto selection_begin = Clock::now();
    AdaptiveGlobalPcgResult result;
    int best_screen_cycles = options.maximum_confirmation_cycles;
    int best_screen_steps = options.minimum_steps;
    int stale = 0;
    int positive_count = 0;
    std::vector<PcgSnapshot> best_refinement_snapshots;
    SparseMatrix initial_safety_prolongation;
    int initial_safety_cycles = options.maximum_confirmation_cycles;
    std::size_t initial_history_index = std::numeric_limits<std::size_t>::max();

    if (options.include_initial_candidate) {
        CandidateRecord initial = run_pilot(
            a, *representative_rhs, initial_prolongation, 0, "screen",
            options.screening_pilot_iterations,
            options.screening_tail_window, options.smoothing_steps,
            options.thread_count, options.solve_tolerance,
            options.maximum_confirmation_cycles, 0.0, 1.0);
        initial.checkpoint.improved = true;
        best_screen_cycles = initial.checkpoint.predicted_cycles;
        best_screen_steps = 0;
        ++result.report.screened_candidates;
        if (initial.checkpoint.predicted_cycles <=
            options.easy_accept_cycles) {
            confirm_candidate(
                a, *representative_rhs, initial,
                options.smoothing_steps, options.thread_count,
                options.solve_tolerance,
                options.maximum_confirmation_cycles);
            initial.checkpoint.selected = true;
            result.prolongation = *initial.prolongation;
            result.report.selected_steps = 0;
            result.report.selected_cycles = initial.checkpoint.confirmed_cycles;
            result.report.screening_anchor_steps = 0;
            result.report.confirmed_candidates = 1;
            result.report.quick_accepted_initial = true;
            result.report.history.push_back(initial.checkpoint);
            result.report.selection_wall_ms = milliseconds(
                selection_begin, Clock::now());
            return result;
        }
        if (initial.checkpoint.predicted_cycles <=
            options.initial_safety_forecast_cycles) {
            confirm_candidate(
                a, *representative_rhs, initial,
                options.smoothing_steps, options.thread_count,
                options.solve_tolerance,
                options.maximum_confirmation_cycles);
            ++result.report.confirmed_candidates;
            initial_safety_cycles = initial.checkpoint.confirmed_cycles;
            initial_safety_prolongation = *initial.prolongation;
            if (initial_safety_cycles <= options.easy_accept_cycles) {
                initial.checkpoint.selected = true;
                result.prolongation = *initial.prolongation;
                result.report.selected_steps = 0;
                result.report.selected_cycles = initial_safety_cycles;
                result.report.screening_anchor_steps = 0;
                result.report.quick_accepted_initial = true;
                result.report.history.push_back(initial.checkpoint);
                result.report.selection_wall_ms = milliseconds(
                    selection_begin, Clock::now());
                return result;
            }
        }
        initial_history_index = result.report.history.size();
        result.report.history.push_back(initial.checkpoint);
    }

    GlobalEnergyPcgPath screening_path(
        grid, a, initial_prolongation, options.thread_count);
    for (long long step_value = options.minimum_steps;
         step_value <= options.maximum_screening_steps;
         step_value += options.screening_increment) {
        const int steps = static_cast<int>(step_value);
        int tentative_begin = std::max(
            options.minimum_steps,
            steps - options.refinement_backtrack_steps);
        const int tentative_end = std::max(
            options.minimum_steps,
            steps - options.refinement_stop_before_anchor_steps);
        const int tentative_remainder =
            (tentative_begin - options.minimum_steps) %
            options.refinement_increment;
        if (tentative_remainder != 0) {
            tentative_begin +=
                options.refinement_increment - tentative_remainder;
        }
        std::vector<PcgSnapshot> tentative_snapshots;
        for (long long candidate_value = tentative_begin;
             candidate_value <= tentative_end && candidate_value < steps;
             candidate_value += options.refinement_increment) {
            const int candidate_steps = static_cast<int>(candidate_value);
            if (candidate_steps <= screening_path.steps()) continue;
            screening_path.advance_to(candidate_steps);
            const GlobalPcgPathReport candidate_path_report =
                screening_path.report();
            tentative_snapshots.push_back({
                candidate_steps,
                screening_path.prolongation(options.drop_tolerance),
                candidate_path_report.total_ms,
                candidate_path_report.maximum_relative_residual});
        }
        screening_path.advance_to(steps);
        const GlobalPcgPathReport path_report = screening_path.report();
        CandidateRecord candidate = run_pilot(
            a, *representative_rhs,
            screening_path.prolongation(options.drop_tolerance),
            steps, "screen", options.screening_pilot_iterations,
            options.screening_tail_window, options.smoothing_steps,
            options.thread_count, options.solve_tolerance,
            options.maximum_confirmation_cycles, path_report.total_ms,
            path_report.maximum_relative_residual);
        ++positive_count;
        ++result.report.screened_candidates;
        if (candidate.checkpoint.predicted_cycles <=
            options.medium_accept_cycles) {
            confirm_candidate(
                a, *representative_rhs, candidate,
                options.smoothing_steps, options.thread_count,
                options.solve_tolerance,
                options.maximum_confirmation_cycles);
            ++result.report.confirmed_candidates;
            if (candidate.checkpoint.confirmed_cycles <=
                options.medium_accept_cycles) {
                if (candidate.checkpoint.confirmed_cycles <=
                    options.easy_accept_cycles) {
                    candidate.checkpoint.selected = true;
                    result.prolongation = *candidate.prolongation;
                    result.report.selected_steps = steps;
                    result.report.selected_cycles =
                        candidate.checkpoint.confirmed_cycles;
                    result.report.screening_anchor_steps = steps;
                    result.report.history.push_back(candidate.checkpoint);
                    result.report.selection_wall_ms = milliseconds(
                        selection_begin, Clock::now());
                    return result;
                }

                const int lookahead_steps = static_cast<int>(std::min(
                    static_cast<long long>(options.maximum_steps),
                    static_cast<long long>(steps) + std::max(
                        options.refinement_increment,
                        options.screening_increment / 2)));
                screening_path.advance_to(lookahead_steps);
                const GlobalPcgPathReport lookahead_path_report =
                    screening_path.report();
                CandidateRecord lookahead = run_pilot(
                    a, *representative_rhs,
                    screening_path.prolongation(options.drop_tolerance),
                    lookahead_steps, "lookahead",
                    options.screening_pilot_iterations,
                    options.screening_tail_window,
                    options.smoothing_steps, options.thread_count,
                    options.solve_tolerance,
                    options.maximum_confirmation_cycles,
                    lookahead_path_report.total_ms,
                    lookahead_path_report.maximum_relative_residual);
                ++result.report.screened_candidates;
                confirm_candidate(
                    a, *representative_rhs, lookahead,
                    options.smoothing_steps, options.thread_count,
                    options.solve_tolerance,
                    options.maximum_confirmation_cycles);
                ++result.report.confirmed_candidates;
                CandidateRecord* selected = &candidate;
                if (lookahead.checkpoint.confirmed_cycles <
                    candidate.checkpoint.confirmed_cycles) {
                    selected = &lookahead;
                }
                selected->checkpoint.selected = true;
                result.prolongation = *selected->prolongation;
                result.report.selected_steps = selected->checkpoint.steps;
                result.report.selected_cycles =
                    selected->checkpoint.confirmed_cycles;
                result.report.screening_anchor_steps = steps;
                result.report.history.push_back(candidate.checkpoint);
                result.report.history.push_back(lookahead.checkpoint);
                result.report.selection_wall_ms = milliseconds(
                    selection_begin, Clock::now());
                return result;
            }
        }
        if (candidate.checkpoint.predicted_cycles < best_screen_cycles) {
            best_screen_cycles = candidate.checkpoint.predicted_cycles;
            best_screen_steps = steps;
            candidate.checkpoint.improved = true;
            stale = 0;
            if (!tentative_snapshots.empty()) {
                best_refinement_snapshots = std::move(tentative_snapshots);
            }
        } else {
            ++stale;
        }
        result.report.history.push_back(candidate.checkpoint);
        if (positive_count >= options.minimum_screened_positive_candidates &&
            stale >= options.screening_patience &&
            steps < options.maximum_screening_steps) {
            result.report.screening_stopped_early = true;
            break;
        }
    }
    if (best_screen_steps == 0) best_screen_steps = options.minimum_steps;
    result.report.screening_anchor_steps = best_screen_steps;

    if (best_refinement_snapshots.empty()) {
        GlobalEnergyPcgPath fallback_path(
            grid, a, initial_prolongation, options.thread_count);
        fallback_path.advance_to(options.minimum_steps);
        const GlobalPcgPathReport fallback_report = fallback_path.report();
        best_refinement_snapshots.push_back({
            options.minimum_steps,
            fallback_path.prolongation(options.drop_tolerance),
            fallback_report.total_ms,
            fallback_report.maximum_relative_residual});
    }
    std::vector<std::unique_ptr<CandidateRecord>> refined;
    for (PcgSnapshot& snapshot : best_refinement_snapshots) {
        auto candidate = std::make_unique<CandidateRecord>(run_pilot(
            a, *representative_rhs,
            std::move(snapshot.prolongation),
            snapshot.steps, "refine", options.refinement_pilot_iterations,
            options.refinement_tail_window, options.smoothing_steps,
            options.thread_count, options.solve_tolerance,
            options.maximum_confirmation_cycles, snapshot.path_ms,
            snapshot.maximum_pcg_residual));
        refined.push_back(std::move(candidate));
    }
    result.report.refined_candidates = static_cast<int>(refined.size());
    if (refined.empty()) {
        throw std::runtime_error("adaptive PCG refinement set is empty");
    }

    std::vector<std::size_t> confirmation_indices{0U};
    std::vector<std::size_t> order(refined.size());
    for (std::size_t index = 0; index < order.size(); ++index) {
        order[index] = index;
    }
    std::sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
        return forecast_less(refined[lhs], refined[rhs]);
    });
    for (std::size_t index : order) {
        if (static_cast<int>(confirmation_indices.size()) >=
            options.confirmation_candidates) break;
        if (std::find(
                confirmation_indices.begin(), confirmation_indices.end(),
                index) == confirmation_indices.end()) {
            confirmation_indices.push_back(index);
        }
    }

    std::size_t selected_index = confirmation_indices.front();
    int selected_cycles = options.maximum_confirmation_cycles;
    for (std::size_t index : confirmation_indices) {
        confirm_candidate(
            a, *representative_rhs, *refined[index],
            options.smoothing_steps, options.thread_count,
            options.solve_tolerance,
            options.maximum_confirmation_cycles);
        ++result.report.confirmed_candidates;
        const int cycles = refined[index]->checkpoint.confirmed_cycles;
        if (cycles < selected_cycles ||
            (cycles == selected_cycles &&
             refined[index]->checkpoint.steps <
                 refined[selected_index]->checkpoint.steps)) {
            selected_cycles = cycles;
            selected_index = index;
        }
    }
    refined[selected_index]->checkpoint.selected = true;
    result.prolongation = *refined[selected_index]->prolongation;
    result.report.selected_steps = refined[selected_index]->checkpoint.steps;
    result.report.selected_cycles = selected_cycles;
    if (initial_safety_cycles < selected_cycles) {
        refined[selected_index]->checkpoint.selected = false;
        result.prolongation = initial_safety_prolongation;
        result.report.selected_steps = 0;
        result.report.selected_cycles = initial_safety_cycles;
        if (initial_history_index < result.report.history.size()) {
            result.report.history[initial_history_index].selected = true;
        }
    }
    for (const auto& candidate : refined) {
        result.report.history.push_back(candidate->checkpoint);
    }
    result.report.selection_wall_ms = milliseconds(
        selection_begin, Clock::now());
    return result;
}

} // namespace tgi
