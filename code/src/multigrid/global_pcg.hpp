#pragma once

#include "multigrid/energy_interpolation.hpp"
#include "multigrid/two_grid_solver.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace tgi {

struct GlobalPcgPathReport {
    int steps = 0;
    double relative_preconditioned_residual = 0.0;
    double total_ms = 0.0;
};

class GlobalEnergyPcgPath {
public:
    GlobalEnergyPcgPath(
        const StructuredGrid& grid, const SparseMatrix& a,
        const SparseMatrix& initial_prolongation, int thread_count = 1);

    void advance_to(int target_steps);
    SparseMatrix prolongation(double drop_tolerance = 0.0);
    GlobalPcgPathReport report() const;

private:
    using Clock = std::chrono::steady_clock;

    struct ColumnState {
        Vector solution;
        Vector residual;
        Vector direction;
        double initial_rz = 1.0;
        double rz = 0.0;
        int iterations = 0;
        bool active = true;
    };

    const StructuredGrid& grid_;
    energy_interpolation_detail::GlobalFSystem system_;
    std::vector<ColumnState> columns_;
    int thread_count_ = 1;
    int steps_ = 0;
    double solve_ms_ = 0.0;
    Clock::time_point begin_;
};

inline GlobalEnergyPcgPath::GlobalEnergyPcgPath(
    const StructuredGrid& grid, const SparseMatrix& a,
    const SparseMatrix& initial_prolongation, int thread_count)
    : grid_(grid), begin_(Clock::now()) {
    system_ = energy_interpolation_detail::assemble_global_f_system(grid, a);
    thread_count_ = std::max(1, std::min(thread_count, grid.coarse_size()));
    columns_.resize(static_cast<std::size_t>(grid.coarse_size()));
    const SparseMatrix initial_transpose =
        initial_prolongation.transpose(thread_count_);

    for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
        ColumnState& state = columns_[static_cast<std::size_t>(coarse)];
        const std::size_t n = system_.f_nodes.size();
        state.solution.assign(n, 0.0);
        for (int position =
                 initial_transpose.row_ptr()[static_cast<std::size_t>(coarse)];
             position < initial_transpose.row_ptr()[
                 static_cast<std::size_t>(coarse) + 1U]; ++position) {
            const int fine = initial_transpose.col_idx()[
                static_cast<std::size_t>(position)];
            if (grid.is_coarse_node(fine)) continue;
            const int local =
                system_.local_index[static_cast<std::size_t>(fine)];
            if (local >= 0) {
                state.solution[static_cast<std::size_t>(local)] =
                    initial_transpose.values()[static_cast<std::size_t>(position)];
            }
        }

        Vector rhs(n, 0.0);
        for (const auto& [row, value] :
             system_.rhs_entries[static_cast<std::size_t>(coarse)]) {
            rhs[static_cast<std::size_t>(row)] += value;
        }
        state.residual = rhs;
        Vector product;
        system_.matrix.multiply(state.solution, product);
        axpy(-1.0, product, state.residual);
        const double residual_scale = std::max(norm2(rhs), 1.0e-30);
        state.direction.resize(n);
        for (std::size_t index = 0; index < n; ++index) {
            state.direction[index] = system_.inverse_diagonal[index] *
                state.residual[index];
        }
        state.rz = dot(state.residual, state.direction);
        state.initial_rz = std::max(state.rz, 1.0e-300);
        const double threshold = std::numeric_limits<double>::epsilon() *
            residual_scale;
        state.active = std::isfinite(state.rz) && state.rz > 0.0 &&
            norm2(state.residual) > threshold;
    }
    system_.assembly_ms = std::chrono::duration<double, std::milli>(
        Clock::now() - begin_).count();
}

inline void GlobalEnergyPcgPath::advance_to(int target_steps) {
    if (target_steps < steps_) {
        throw std::invalid_argument(
            "GlobalEnergyPcgPath checkpoints must be nondecreasing");
    }
    if (target_steps < 0) {
        throw std::invalid_argument(
            "GlobalEnergyPcgPath target steps must be nonnegative");
    }
    if (target_steps == steps_) return;
    const auto solve_begin = Clock::now();
    std::atomic<int> next_column{0};
    std::exception_ptr worker_error;
    std::mutex error_mutex;
    auto worker = [&]() {
        Vector product;
        Vector z(system_.f_nodes.size(), 0.0);
        try {
            while (true) {
                const int coarse =
                    next_column.fetch_add(1, std::memory_order_relaxed);
                if (coarse >= grid_.coarse_size()) break;
                ColumnState& state = columns_[static_cast<std::size_t>(coarse)];
                while (state.active && state.iterations < target_steps) {
                    system_.matrix.multiply(state.direction, product);
                    const double denominator = dot(state.direction, product);
                    if (!(denominator > 0.0) || !std::isfinite(denominator) ||
                        !(state.rz > 0.0) || !std::isfinite(state.rz)) {
                        state.active = false;
                        break;
                    }
                    const double alpha = state.rz / denominator;
                    axpy(alpha, state.direction, state.solution);
                    axpy(-alpha, product, state.residual);
                    ++state.iterations;
                    for (std::size_t index = 0; index < z.size(); ++index) {
                        z[index] = system_.inverse_diagonal[index] *
                            state.residual[index];
                    }
                    const double rz_new = dot(state.residual, z);
                    if (!(rz_new > 0.0) || !std::isfinite(rz_new)) {
                        state.rz = rz_new;
                        state.active = false;
                        break;
                    }
                    const double beta = rz_new / state.rz;
                    for (std::size_t index = 0;
                         index < state.direction.size(); ++index) {
                        state.direction[index] =
                            z[index] + beta * state.direction[index];
                    }
                    state.rz = rz_new;
                }
            }
        } catch (...) {
            std::lock_guard<std::mutex> lock(error_mutex);
            if (!worker_error) worker_error = std::current_exception();
            next_column.store(grid_.coarse_size());
        }
    };
    if (thread_count_ == 1) {
        worker();
    } else {
        std::vector<std::thread> workers;
        workers.reserve(static_cast<std::size_t>(thread_count_));
        for (int index = 0; index < thread_count_; ++index) {
            workers.emplace_back(worker);
        }
        for (auto& thread : workers) thread.join();
    }
    if (worker_error) std::rethrow_exception(worker_error);
    steps_ = target_steps;
    solve_ms_ += std::chrono::duration<double, std::milli>(
        Clock::now() - solve_begin).count();
}

inline SparseMatrix GlobalEnergyPcgPath::prolongation(
    double drop_tolerance) {
    std::vector<Triplet> entries;
    std::size_t entry_count = static_cast<std::size_t>(grid_.coarse_size());
    for (const ColumnState& state : columns_) {
        entry_count += static_cast<std::size_t>(std::count_if(
            state.solution.begin(), state.solution.end(),
            [drop_tolerance](double value) {
                return std::abs(value) > drop_tolerance;
            }));
    }
    entries.reserve(entry_count);
    for (int coarse = 0; coarse < grid_.coarse_size(); ++coarse) {
        entries.push_back({grid_.coarse_fine_id(coarse), coarse, 1.0});
        const Vector& weights =
            columns_[static_cast<std::size_t>(coarse)].solution;
        for (std::size_t local = 0; local < weights.size(); ++local) {
            if (std::abs(weights[local]) > drop_tolerance) {
                entries.push_back(
                    {system_.f_nodes[local], coarse, weights[local]});
            }
        }
    }
    SparseMatrix result(
        grid_.fine_size(), grid_.coarse_size(), entries, 0.0);
    return result;
}

inline GlobalPcgPathReport GlobalEnergyPcgPath::report() const {
    GlobalPcgPathReport result;
    result.steps = steps_;
    result.total_ms = system_.assembly_ms + solve_ms_;
    double rz_sum = 0.0;
    double initial_rz_sum = 0.0;
    for (const ColumnState& state : columns_) {
        rz_sum += std::max(state.rz, 0.0);
        initial_rz_sum += state.initial_rz;
    }
    result.relative_preconditioned_residual = std::sqrt(
        rz_sum / initial_rz_sum);
    return result;
}

}

namespace tgi {

struct AdaptiveGlobalPcgOptions {
    int minimum_steps = 12;
    int maximum_steps = 60;
    int step_quantum = 2;
    int screening_stride = 8;
    int pilot_iterations = 64;
    int confirmation_iterations = 384;
    int confirmation_candidates = 3;
    int tail_window = 24;
    double early_accept_forecast_multiple = 2.0;
    double acceptable_cycle_slack = 0.03;
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
    bool pilot_converged = false;
    bool confirmed = false;
    bool selected = false;
};

struct AdaptiveGlobalPcgReport {
    int selected_steps = 0;
    int estimated_selected_cycles = 0;
    double selection_wall_ms = 0.0;
    std::vector<AdaptiveGlobalPcgCheckpoint> history;
};

struct AdaptiveGlobalPcgResult {
    std::shared_ptr<SparseMatrix> prolongation;
    std::unique_ptr<TwoGridCycle> cycle;
    AdaptiveGlobalPcgReport report;
};

namespace adaptive_global_pcg_detail {

using Clock = std::chrono::steady_clock;

inline double milliseconds(Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

inline double median(std::vector<double> values) {
    if (values.empty()) {
        throw std::invalid_argument("median requires at least one value");
    }
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
    if (completed <= 0) return {};
    if (completed == 1) {
        const double current = std::max(norms.back(), 1.0e-300);
        const double previous = std::max(norms.front(), 1.0e-300);
        TailModel model;
        model.rho = std::exp(std::min(
            -1.0e-10, std::log(current / previous)));
        model.relative_uncertainty = 1.0;
        return model;
    }
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
    const int first = std::max(0, completed - tail_window);
    const double long_rate = (
        std::log(std::max(norms.back(), 1.0e-300)) -
        std::log(std::max(
            norms[static_cast<std::size_t>(first)], 1.0e-300))) /
        static_cast<double>(completed - first);
    const double blended = 0.65 * center + 0.35 * long_rate;
    const double conservative = blended + 0.25 * scale;
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
    Vector solution;
    Vector residual;
    TwoGridCycle::Workspace workspace;
    std::vector<double> norms;
    double initial_norm = 1.0;
};

inline std::unique_ptr<Candidate> make_candidate(
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

    candidate->cycle = std::make_unique<TwoGridCycle>(
        a, *candidate->prolongation,
        options.smoothing_steps, options.thread_count);
    candidate->solution.assign(rhs.size(), 0.0);
    candidate->residual = rhs;
    candidate->initial_norm = std::max(norm2(rhs), 1.0e-300);
    candidate->norms.reserve(
        static_cast<std::size_t>(options.confirmation_iterations) + 1U);
    candidate->norms.push_back(candidate->initial_norm);
    return candidate;
}

inline void advance_candidate(
    Candidate& candidate, const Vector& rhs, int target_iterations,
    const AdaptiveGlobalPcgOptions& options, bool confirmation) {
    const auto begin = Clock::now();
    while (candidate.checkpoint.pilot_iterations < target_iterations &&
           !candidate.checkpoint.pilot_converged) {
        const double squared = candidate.cycle->iterate(
            rhs, candidate.solution, candidate.residual,
            candidate.workspace);
        const double current = std::max(std::sqrt(squared), 1.0e-300);
        candidate.norms.push_back(current);
        ++candidate.checkpoint.pilot_iterations;
        candidate.checkpoint.pilot_converged =
            current / candidate.initial_norm <= options.solve_tolerance;
    }
    candidate.checkpoint.pilot_ms += milliseconds(begin, Clock::now());
    candidate.checkpoint.confirmed =
        candidate.checkpoint.confirmed || confirmation;
    candidate.checkpoint.pilot_relative_residual =
        candidate.norms.back() / candidate.initial_norm;
    const TailModel model = fit_tail(
        candidate.norms,
        std::min(options.tail_window,
                 candidate.checkpoint.pilot_iterations));
    candidate.checkpoint.rho_rhs_pilot = model.rho;
    candidate.checkpoint.forecast_relative_uncertainty =
        model.relative_uncertainty;
    candidate.checkpoint.predicted_cycles = forecast_cycles(
        candidate.checkpoint.pilot_iterations,
        candidate.checkpoint.pilot_relative_residual,
        candidate.checkpoint.rho_rhs_pilot,
        options.solve_tolerance, options.maximum_cycles);
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

inline std::vector<int> screening_steps(
    const AdaptiveGlobalPcgOptions& options) {
    std::vector<int> steps;
    for (int value = options.minimum_steps;
         value < options.maximum_steps;
         value += options.screening_stride) {
        steps.push_back(quantize(
            value, options.step_quantum,
            options.minimum_steps, options.maximum_steps));
    }
    if (steps.empty() || steps.back() != options.maximum_steps) {
        steps.push_back(options.maximum_steps);
    }
    steps.erase(std::unique(steps.begin(), steps.end()), steps.end());
    return steps;
}

}

inline AdaptiveGlobalPcgResult build_adaptive_global_pcg_interpolation(
    const StructuredGrid& grid, const SparseMatrix& a,
    const SparseMatrix& initial_prolongation,
    const AdaptiveGlobalPcgOptions& options,
    const Vector* representative_rhs) {
    using namespace adaptive_global_pcg_detail;
    if (representative_rhs == nullptr) {
        throw std::invalid_argument(
            "adaptive global PCG requires a representative right-hand side");
    }
    if (representative_rhs->size() !=
        static_cast<std::size_t>(a.rows())) {
        throw std::invalid_argument(
            "adaptive global PCG right-hand side has the wrong size");
    }
    if (options.minimum_steps < 0 ||
        options.maximum_steps <= options.minimum_steps ||
        options.step_quantum <= 0 || options.screening_stride <= 0 ||
        options.pilot_iterations < 2 ||
        options.confirmation_iterations < options.pilot_iterations ||
        options.confirmation_candidates <= 0 ||
        options.tail_window < 2 || options.maximum_cycles <= 0 ||
        options.smoothing_steps <= 0 ||
        !(options.solve_tolerance > 0.0) ||
        !(options.solve_tolerance < 1.0) ||
        options.drop_tolerance < 0.0 ||
        options.acceptable_cycle_slack < 0.0) {
        throw std::invalid_argument(
            "adaptive global PCG received invalid options");
    }
    const auto begin = Clock::now();
    const Vector& rhs = *representative_rhs;
    std::vector<std::unique_ptr<Candidate>> candidates;
    std::vector<AdaptiveGlobalPcgCheckpoint> archived_checkpoints;
    GlobalPcgPathReport initial_report;
    initial_report.relative_preconditioned_residual = 1.0;
    candidates.push_back(make_candidate(
        a, rhs, initial_prolongation, 0, "initial",
        options, initial_report));
    advance_candidate(
        *candidates.back(), rhs, options.pilot_iterations,
        options, false);

    AdaptiveGlobalPcgResult result;
    const double easy_limit = options.early_accept_forecast_multiple *
        static_cast<double>(options.pilot_iterations);
    if (candidates.front()->checkpoint.predicted_cycles > easy_limit) {
        GlobalEnergyPcgPath path(
            grid, a, initial_prolongation, options.thread_count);
        for (int steps : screening_steps(options)) {
            path.advance_to(steps);
            candidates.push_back(make_candidate(
                a, rhs, path.prolongation(options.drop_tolerance),
                steps, "screen", options, path.report()));
            advance_candidate(
                *candidates.back(), rhs, options.pilot_iterations,
                options, false);
        }

        std::vector<std::size_t> ranking(candidates.size());
        for (std::size_t index = 0; index < ranking.size(); ++index) {
            ranking[index] = index;
        }
        std::stable_sort(
            ranking.begin(), ranking.end(),
            [&](std::size_t lhs, std::size_t rhs_index) {
                const auto& left = candidates[lhs]->checkpoint;
                const auto& right = candidates[rhs_index]->checkpoint;
                if (left.predicted_cycles != right.predicted_cycles) {
                    return left.predicted_cycles < right.predicted_cycles;
                }
                return left.preconditioned_pcg_residual <
                    right.preconditioned_pcg_residual;
            });
        const int ranked_confirmations = std::min(
            options.confirmation_candidates,
            static_cast<int>(ranking.size()));
        std::vector<std::size_t> confirmation_indices;
        confirmation_indices.reserve(
            static_cast<std::size_t>(ranked_confirmations) + 1U);
        for (int index = 0; index < ranked_confirmations; ++index) {
            confirmation_indices.push_back(
                ranking[static_cast<std::size_t>(index)]);
        }
        const std::size_t parsimonious_screen_winner =
            select_candidate(candidates, options.acceptable_cycle_slack);
        if (std::find(
                confirmation_indices.begin(), confirmation_indices.end(),
                parsimonious_screen_winner) ==
            confirmation_indices.end()) {
            confirmation_indices.push_back(parsimonious_screen_winner);
        }
        const auto anchor = std::find_if(
            candidates.begin(), candidates.end(),
            [&](const auto& candidate) {
                return candidate->checkpoint.steps == options.minimum_steps;
            });
        if (anchor != candidates.end()) {
            const std::size_t anchor_index = static_cast<std::size_t>(
                anchor - candidates.begin());
            if (std::find(
                    confirmation_indices.begin(), confirmation_indices.end(),
                    anchor_index) == confirmation_indices.end()) {
                confirmation_indices.push_back(anchor_index);
            }
        }
        std::vector<bool> confirmed(candidates.size(), false);
        for (std::size_t candidate_index : confirmation_indices) {
            confirmed[candidate_index] = true;
            advance_candidate(
                *candidates[candidate_index],
                rhs, options.confirmation_iterations, options, true);
        }

        std::vector<std::unique_ptr<Candidate>> finalists;
        finalists.reserve(
            confirmation_indices.size() + 3U);
        for (std::size_t index = 0; index < candidates.size(); ++index) {
            if (confirmed[index]) {
                finalists.push_back(std::move(candidates[index]));
            } else {
                archived_checkpoints.push_back(candidates[index]->checkpoint);
            }
        }
        candidates = std::move(finalists);

        const std::size_t provisional = select_candidate(candidates, 0.0);
        std::size_t runner_up = provisional;
        for (std::size_t index = 0; index < candidates.size(); ++index) {
            if (index == provisional) continue;
            if (runner_up == provisional ||
                candidates[index]->checkpoint.predicted_cycles <
                    candidates[runner_up]->checkpoint.predicted_cycles) {
                runner_up = index;
            }
        }
        const int best_steps = candidates[provisional]->checkpoint.steps;
        std::vector<int> refinement_steps;
        for (int value : {
                 best_steps - options.step_quantum,
                 best_steps + options.step_quantum}) {
            if (value > options.minimum_steps &&
                value < options.maximum_steps) {
                refinement_steps.push_back(value);
            }
        }
        if (runner_up != provisional) {
            const int other_steps =
                candidates[runner_up]->checkpoint.steps;
            if (std::abs(best_steps - other_steps) >
                2 * options.step_quantum) {
                refinement_steps.push_back(quantize(
                    (best_steps + other_steps) / 2,
                    options.step_quantum,
                    options.minimum_steps, options.maximum_steps));
            }
        }
        std::sort(refinement_steps.begin(), refinement_steps.end());
        refinement_steps.erase(
            std::unique(refinement_steps.begin(), refinement_steps.end()),
            refinement_steps.end());
        refinement_steps.erase(
            std::remove_if(
                refinement_steps.begin(), refinement_steps.end(),
                [&](int steps) {
                    return std::any_of(
                        candidates.begin(), candidates.end(),
                        [&](const auto& candidate) {
                            return candidate->checkpoint.steps == steps;
                        });
                }),
            refinement_steps.end());
        if (!refinement_steps.empty()) {
            GlobalEnergyPcgPath refinement_path(
                grid, a, initial_prolongation, options.thread_count);
            for (int steps : refinement_steps) {
                refinement_path.advance_to(steps);
                candidates.push_back(make_candidate(
                    a, rhs,
                    refinement_path.prolongation(options.drop_tolerance),
                    steps, "refine", options,
                    refinement_path.report()));
                advance_candidate(
                    *candidates.back(), rhs,
                    options.confirmation_iterations, options, true);
            }
        }
    }

    const std::size_t selected = select_candidate(
        candidates, options.acceptable_cycle_slack);
    candidates[selected]->checkpoint.selected = true;
    result.prolongation = std::move(candidates[selected]->prolongation);
    result.cycle = std::move(candidates[selected]->cycle);
    result.report.selected_steps = candidates[selected]->checkpoint.steps;
    result.report.estimated_selected_cycles =
        candidates[selected]->checkpoint.predicted_cycles;
    for (const auto& checkpoint : archived_checkpoints) {
        result.report.history.push_back(checkpoint);
    }
    for (const auto& candidate : candidates) {
        result.report.history.push_back(candidate->checkpoint);
    }
    result.report.selection_wall_ms = milliseconds(begin, Clock::now());
    return result;
}

}
