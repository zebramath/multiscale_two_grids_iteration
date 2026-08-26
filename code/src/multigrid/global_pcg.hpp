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
    double expected_rhs_count = 1.0;
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
    double cycle_setup_ms = 0.0;
    double pilot_ms = 0.0;
    double predicted_solve_ms = 0.0;
    double amortized_cost_ms = 0.0;
    bool pilot_converged = false;
    bool selected = false;
};

struct AdaptiveGlobalPcgReport {
    int selected_steps = 0;
    int estimated_selected_cycles = 0;
    int checkpoint_stride = 0;
    int maximum_sampled_steps = 0;
    int pilot_iterations = 0;
    double expected_rhs_count = 1.0;
    bool used_local_refinement = false;
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

inline int quantize(int value, int quantum) {
    return std::max(quantum, ((value + quantum / 2) / quantum) * quantum);
}

struct SelectionProfile {
    int stride = 8;
    int maximum_steps = 60;
    int pilot_iterations = 160;
    double near_optimal_slack = 0.02;
    bool refine = false;
};

inline SelectionProfile selection_profile(
    const AdaptiveGlobalPcgOptions& options) {
    const double detail = std::min(
        1.0, std::max(0.0, std::log2(options.expected_rhs_count) / 8.0));
    SelectionProfile profile;
    profile.stride = quantize(static_cast<int>(
        std::lround(20.0 - 12.0 * detail)), options.step_quantum);
    const int removable_tail = std::min(
        8, options.maximum_steps - options.minimum_steps);
    profile.maximum_steps = std::max(
        options.minimum_steps,
        quantize(static_cast<int>(std::lround(
            static_cast<double>(options.maximum_steps - removable_tail) +
            static_cast<double>(removable_tail) * detail)),
            options.step_quantum));
    profile.maximum_steps = std::min(
        profile.maximum_steps, options.maximum_steps);
    profile.pilot_iterations = static_cast<int>(
        std::lround(16.0 + 144.0 * detail));
    profile.near_optimal_slack = 0.10 - 0.08 * detail;
    profile.refine = detail >= 0.5;
    return profile;
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

inline TailModel fit_tail(const std::vector<double>& norms) {
    const int completed = static_cast<int>(norms.size()) - 1;
    if (completed <= 0) return {};
    const int begin = std::max(1, completed / 2);
    std::vector<double> rates;
    rates.reserve(static_cast<std::size_t>(completed - begin + 1));
    for (int iteration = begin; iteration <= completed; ++iteration) {
        const double current = std::max(
            norms[static_cast<std::size_t>(iteration)], 1.0e-300);
        const double previous = std::max(norms[
            static_cast<std::size_t>(iteration - 1)], 1.0e-300);
        rates.push_back(std::log(current / previous));
    }
    const double center = median(rates);
    std::vector<double> deviations;
    deviations.reserve(rates.size());
    for (double rate : rates) deviations.push_back(std::abs(rate - center));
    const double scale = 1.4826 * median(std::move(deviations));
    TailModel model;
    model.rho = std::exp(std::min(-1.0e-10, center + 0.25 * scale));
    model.relative_uncertainty = std::min(
        2.0, scale / std::max(std::abs(center), 1.0e-10));
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
    checkpoint.cycle_setup_ms =
        candidate->cycle->setup_report().total_ms;
    candidate->solution.assign(rhs.size(), 0.0);
    candidate->residual = rhs;
    candidate->initial_norm = std::max(norm2(rhs), 1.0e-300);
    candidate->norms.push_back(candidate->initial_norm);
    return candidate;
}

inline void probe_candidate(
    Candidate& candidate, const Vector& rhs, int pilot_iterations,
    const AdaptiveGlobalPcgOptions& options) {
    const auto begin = Clock::now();
    while (candidate.checkpoint.pilot_iterations < pilot_iterations &&
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
    candidate.checkpoint.pilot_relative_residual =
        candidate.norms.back() / candidate.initial_norm;
    const TailModel model = fit_tail(candidate.norms);
    candidate.checkpoint.rho_rhs_pilot = model.rho;
    candidate.checkpoint.forecast_relative_uncertainty =
        model.relative_uncertainty;
    candidate.checkpoint.predicted_cycles = forecast_cycles(
        candidate.checkpoint.pilot_iterations,
        candidate.checkpoint.pilot_relative_residual,
        candidate.checkpoint.rho_rhs_pilot,
        options.solve_tolerance, options.maximum_cycles);
    const double cycle_ms = candidate.checkpoint.pilot_ms /
        static_cast<double>(
            std::max(1, candidate.checkpoint.pilot_iterations));
    candidate.checkpoint.predicted_solve_ms = cycle_ms *
        static_cast<double>(candidate.checkpoint.predicted_cycles);
    candidate.checkpoint.amortized_cost_ms =
        candidate.checkpoint.cycle_setup_ms /
            options.expected_rhs_count +
        candidate.checkpoint.predicted_solve_ms;
}

inline std::size_t select_candidate(
    const std::vector<std::unique_ptr<Candidate>>& candidates,
    double slack) {
    int best = std::numeric_limits<int>::max();
    for (const auto& candidate : candidates) {
        best = std::min(best, candidate->checkpoint.predicted_cycles);
    }
    const double limit = (1.0 + slack) * static_cast<double>(best);
    std::size_t selected = 0;
    bool found = false;
    for (std::size_t index = 0; index < candidates.size(); ++index) {
        const auto& checkpoint = candidates[index]->checkpoint;
        if (static_cast<double>(checkpoint.predicted_cycles) <= limit &&
            (!found || candidates[index]->checkpoint.steps <
                candidates[selected]->checkpoint.steps)) {
            selected = index;
            found = true;
        }
    }
    return selected;
}

inline std::vector<int> sampled_steps(
    const AdaptiveGlobalPcgOptions& options,
    const SelectionProfile& profile) {
    std::vector<int> steps;
    for (int value = options.minimum_steps;
         value < profile.maximum_steps;
         value += profile.stride) {
        steps.push_back(value);
    }
    if (steps.empty() || steps.back() != profile.maximum_steps) {
        steps.push_back(profile.maximum_steps);
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
        options.step_quantum <= 0 || options.maximum_cycles <= 0 ||
        options.smoothing_steps <= 0 ||
        !(options.expected_rhs_count >= 1.0) ||
        !std::isfinite(options.expected_rhs_count) ||
        !(options.solve_tolerance > 0.0) ||
        !(options.solve_tolerance < 1.0) ||
        options.drop_tolerance < 0.0) {
        throw std::invalid_argument(
            "adaptive global PCG received invalid options");
    }
    const auto begin = Clock::now();
    const Vector& rhs = *representative_rhs;
    const SelectionProfile profile = selection_profile(options);
    std::vector<std::unique_ptr<Candidate>> candidates;
    GlobalPcgPathReport initial_report;
    initial_report.relative_preconditioned_residual = 1.0;
    candidates.push_back(make_candidate(
        a, rhs, initial_prolongation, 0, "sample",
        options, initial_report));
    probe_candidate(
        *candidates.back(), rhs, profile.pilot_iterations, options);

    AdaptiveGlobalPcgResult result;
    GlobalEnergyPcgPath path(
        grid, a, initial_prolongation, options.thread_count);
    for (int steps : sampled_steps(options, profile)) {
        path.advance_to(steps);
        candidates.push_back(make_candidate(
            a, rhs, path.prolongation(options.drop_tolerance),
            steps, "sample", options, path.report()));
        probe_candidate(
            *candidates.back(), rhs, profile.pilot_iterations, options);
    }

    if (profile.refine) {
        const std::size_t provisional = select_candidate(candidates, 0.0);
        const int center = candidates[provisional]->checkpoint.steps;
        std::vector<int> refinement_steps;
        for (int steps : {center - options.step_quantum,
                          center + options.step_quantum}) {
            if (steps >= options.minimum_steps &&
                steps <= profile.maximum_steps &&
                std::none_of(
                    candidates.begin(), candidates.end(),
                    [&](const auto& candidate) {
                        return candidate->checkpoint.steps == steps;
                    })) {
                refinement_steps.push_back(steps);
            }
        }
        std::sort(refinement_steps.begin(), refinement_steps.end());
        if (!refinement_steps.empty()) {
            GlobalEnergyPcgPath refinement_path(
                grid, a, initial_prolongation, options.thread_count);
            for (int steps : refinement_steps) {
                refinement_path.advance_to(steps);
                candidates.push_back(make_candidate(
                    a, rhs,
                    refinement_path.prolongation(options.drop_tolerance),
                    steps, "refine", options, refinement_path.report()));
                probe_candidate(
                    *candidates.back(), rhs,
                    profile.pilot_iterations, options);
            }
            result.report.used_local_refinement = true;
        }
    }

    const std::size_t selected = select_candidate(
        candidates, profile.near_optimal_slack);
    candidates[selected]->checkpoint.selected = true;
    result.prolongation = std::move(candidates[selected]->prolongation);
    result.cycle = std::move(candidates[selected]->cycle);
    result.report.selected_steps = candidates[selected]->checkpoint.steps;
    result.report.estimated_selected_cycles =
        candidates[selected]->checkpoint.predicted_cycles;
    result.report.checkpoint_stride = profile.stride;
    result.report.maximum_sampled_steps = profile.maximum_steps;
    result.report.pilot_iterations = profile.pilot_iterations;
    result.report.expected_rhs_count = options.expected_rhs_count;
    for (const auto& candidate : candidates) {
        result.report.history.push_back(candidate->checkpoint);
    }
    result.report.selection_wall_ms = milliseconds(begin, Clock::now());
    return result;
}

}
