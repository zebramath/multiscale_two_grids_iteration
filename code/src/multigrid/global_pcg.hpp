#pragma once

#include "multigrid/energy_interpolation.hpp"
#include "multigrid/two_grid_solver.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

namespace tgi {

class GlobalEnergyPcgPath {
public:
    GlobalEnergyPcgPath(
        const StructuredGrid& grid, const SparseMatrix& a,
        const SparseMatrix& initial_prolongation, int thread_count = 1);

    void advance_to(int target_steps);
    SparseMatrix prolongation(double drop_tolerance = 0.0);

private:
    struct ColumnState {
        Vector solution;
        Vector residual;
        Vector direction;
        double rz = 0.0;
        int iterations = 0;
        bool active = true;
    };

    const StructuredGrid& grid_;
    energy_interpolation_detail::GlobalFSystem system_;
    std::vector<ColumnState> columns_;
    int thread_count_ = 1;
    int steps_ = 0;
};

inline GlobalEnergyPcgPath::GlobalEnergyPcgPath(
    const StructuredGrid& grid, const SparseMatrix& a,
    const SparseMatrix& initial_prolongation, int thread_count)
    : grid_(grid) {
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
        const double threshold = std::numeric_limits<double>::epsilon() *
            residual_scale;
        state.active = std::isfinite(state.rz) && state.rz > 0.0 &&
            norm2(state.residual) > threshold;
    }
}

inline void GlobalEnergyPcgPath::advance_to(int target_steps) {
    if (target_steps < 0) {
        throw std::invalid_argument(
            "GlobalEnergyPcgPath target steps must be nonnegative");
    }
    if (target_steps < steps_) {
        throw std::invalid_argument(
            "GlobalEnergyPcgPath checkpoints must be nondecreasing");
    }
    if (target_steps == steps_) return;
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
    return SparseMatrix(
        grid_.fine_size(), grid_.coarse_size(), entries, 0.0);
}

enum class AdaptiveGlobalPcgPolicy {
    Fast,
    Reuse
};

struct AdaptiveGlobalPcgOptions {
    AdaptiveGlobalPcgPolicy policy = AdaptiveGlobalPcgPolicy::Fast;
    int maximum_cycles = 40000;
    int smoothing_steps = 1;
    int thread_count = 1;
    double solve_tolerance = 1.0e-6;
    double drop_tolerance = 0.0;
};

struct AdaptiveGlobalPcgReport {
    int selected_steps = 0;
    int candidate_count = 0;
    int maximum_sampled_steps = 0;
    int pilot_cycles = 0;
};

struct AdaptiveGlobalPcgResult {
    std::shared_ptr<SparseMatrix> prolongation;
    std::unique_ptr<TwoGridCycle> cycle;
    AdaptiveGlobalPcgReport report;
};

namespace adaptive_global_pcg_detail {

struct SelectionProfile {
    std::vector<int> checkpoints;
    int pilot_cycles = 0;
};

inline int scaled_checkpoint(
    int resolution, int numerator, int denominator) {
    return (numerator * resolution + denominator / 2) / denominator;
}

inline SelectionProfile selection_profile(
    const StructuredGrid& grid, const SparseMatrix& a,
    AdaptiveGlobalPcgPolicy policy) {
    const int resolution = grid.intervals();
    SelectionProfile profile;
    if (policy == AdaptiveGlobalPcgPolicy::Fast) {
        const int coarse_resolution = resolution / grid.ratio();
        if (coarse_resolution <= 8) {
            profile.checkpoints.push_back(
                scaled_checkpoint(resolution, 1, 8));
            return profile;
        }
        const Vector diagonal = a.diagonal();
        const auto extrema = std::minmax_element(
            diagonal.begin(), diagonal.end());
        const double stiffness_ratio = *extrema.second / *extrema.first;
        if (stiffness_ratio < 1.0e3) {
            profile.checkpoints.push_back(
                scaled_checkpoint(resolution, 1, 4));
        } else if (stiffness_ratio >= 1.0e5) {
            profile.checkpoints.push_back(
                scaled_checkpoint(resolution, 1, 2));
        } else {
            profile.checkpoints.push_back(
                scaled_checkpoint(resolution, 1, 3));
        }
        return profile;
    }
    profile.checkpoints = {
        scaled_checkpoint(resolution, 1, 8),
        scaled_checkpoint(resolution, 3, 16),
        scaled_checkpoint(resolution, 1, 4),
        scaled_checkpoint(resolution, 1, 3),
        scaled_checkpoint(resolution, 1, 2)};
    profile.checkpoints.erase(
        std::unique(
            profile.checkpoints.begin(), profile.checkpoints.end()),
        profile.checkpoints.end());
    profile.pilot_cycles = std::max(
        1, scaled_checkpoint(resolution, 1, 2));
    return profile;
}

inline double fit_tail(const std::vector<double>& norms) {
    const int completed = static_cast<int>(norms.size()) - 1;
    if (completed <= 0) return 1.0;
    const int begin = completed / 2;
    const int span = completed - begin;
    const double first = std::max(
        norms[static_cast<std::size_t>(begin)], 1.0e-300);
    const double last = std::max(norms.back(), 1.0e-300);
    return std::exp(std::log(last / first) /
                    static_cast<double>(span));
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
    int steps = 0;
    int pilot_iterations = 0;
    int predicted_cycles = 0;
    bool pilot_converged = false;
    Vector solution;
    Vector residual;
    TwoGridCycle::Workspace workspace;
    std::vector<double> norms;
    double initial_norm = 1.0;
};

inline std::unique_ptr<Candidate> make_candidate(
    const SparseMatrix& a, const Vector& rhs, SparseMatrix prolongation,
    int steps, const AdaptiveGlobalPcgOptions& options) {
    auto candidate = std::make_unique<Candidate>();
    candidate->prolongation = std::make_shared<SparseMatrix>(
        std::move(prolongation));
    candidate->steps = steps;
    candidate->cycle = std::make_unique<TwoGridCycle>(
        a, *candidate->prolongation,
        options.smoothing_steps, options.thread_count);
    candidate->solution.assign(rhs.size(), 0.0);
    candidate->residual = rhs;
    candidate->initial_norm = std::max(norm2(rhs), 1.0e-300);
    candidate->norms.push_back(candidate->initial_norm);
    return candidate;
}

inline void probe_candidate(
    Candidate& candidate, const Vector& rhs, int pilot_limit,
    const AdaptiveGlobalPcgOptions& options) {
    while (candidate.pilot_iterations < pilot_limit &&
           !candidate.pilot_converged) {
        const double squared = candidate.cycle->iterate(
            rhs, candidate.solution, candidate.residual,
            candidate.workspace);
        const double current = std::max(std::sqrt(squared), 1.0e-300);
        candidate.norms.push_back(current);
        ++candidate.pilot_iterations;
        candidate.pilot_converged =
            current / candidate.initial_norm <= options.solve_tolerance;
    }
    const double relative_residual =
        candidate.norms.back() / candidate.initial_norm;
    const double pilot_rho = fit_tail(candidate.norms);
    candidate.predicted_cycles = forecast_cycles(
        candidate.pilot_iterations, relative_residual,
        pilot_rho,
        options.solve_tolerance, options.maximum_cycles);
}

inline bool is_better_candidate(
    const Candidate& candidate, const Candidate& incumbent) {
    return candidate.predicted_cycles < incumbent.predicted_cycles ||
        (candidate.predicted_cycles == incumbent.predicted_cycles &&
         candidate.steps < incumbent.steps);
}

}

inline AdaptiveGlobalPcgResult build_adaptive_global_pcg_interpolation(
    const StructuredGrid& grid, const SparseMatrix& a,
    const SparseMatrix& initial_prolongation,
    const AdaptiveGlobalPcgOptions& options,
    const Vector& representative_rhs) {
    using namespace adaptive_global_pcg_detail;
    if (options.maximum_cycles <= 0 ||
        options.smoothing_steps <= 0 ||
        !(options.solve_tolerance > 0.0) ||
        !(options.solve_tolerance < 1.0) ||
        options.drop_tolerance < 0.0) {
        throw std::invalid_argument(
            "adaptive global PCG received invalid options");
    }
    const SelectionProfile profile = selection_profile(
        grid, a, options.policy);
    GlobalEnergyPcgPath path(
        grid, a, initial_prolongation, options.thread_count);
    std::unique_ptr<Candidate> best;
    for (int steps : profile.checkpoints) {
        std::unique_ptr<Candidate> candidate;
        if (steps > 0) {
            path.advance_to(steps);
            candidate = make_candidate(
                a, representative_rhs,
                path.prolongation(options.drop_tolerance),
                steps, options);
        } else {
            candidate = make_candidate(
                a, representative_rhs, initial_prolongation,
                steps, options);
        }
        if (profile.pilot_cycles > 0) {
            probe_candidate(
                *candidate, representative_rhs,
                profile.pilot_cycles, options);
        }
        if (!best || is_better_candidate(*candidate, *best)) {
            best = std::move(candidate);
        }
    }

    AdaptiveGlobalPcgResult result;
    result.prolongation = std::move(best->prolongation);
    result.cycle = std::move(best->cycle);
    result.report.selected_steps = best->steps;
    result.report.candidate_count =
        static_cast<int>(profile.checkpoints.size());
    result.report.maximum_sampled_steps = profile.checkpoints.back();
    result.report.pilot_cycles = profile.pilot_cycles;
    return result;
}

}
