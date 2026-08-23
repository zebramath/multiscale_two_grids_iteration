#pragma once

#include "multigrid/energy_interpolation.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

namespace tgi {

struct GlobalPcgPathReport {
    int steps = 0;
    int active_columns = 0;
    long long total_column_iterations = 0;
    double maximum_relative_residual = 0.0;
    double rms_relative_residual = 0.0;
    double relative_preconditioned_residual = 0.0;
    double assembly_ms = 0.0;
    double solve_ms = 0.0;
    double finalize_ms = 0.0;
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
    int steps() const { return steps_; }

private:
    using Clock = std::chrono::steady_clock;

    struct ColumnState {
        Vector solution;
        Vector residual;
        Vector direction;
        double residual_scale = 1.0;
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
    long long total_column_iterations_ = 0;
    double solve_ms_ = 0.0;
    double finalize_ms_ = 0.0;
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
        state.residual_scale = std::max(norm2(rhs), 1.0e-30);
        state.direction.resize(n);
        for (std::size_t index = 0; index < n; ++index) {
            state.direction[index] = system_.inverse_diagonal[index] *
                state.residual[index];
        }
        state.rz = dot(state.residual, state.direction);
        state.initial_rz = std::max(state.rz, 1.0e-300);
        const double threshold = std::numeric_limits<double>::epsilon() *
            state.residual_scale;
        state.active = std::isfinite(state.rz) && state.rz > 0.0 &&
            norm2(state.residual) > threshold;
    }
    system_.assembly_ms = std::chrono::duration<double, std::milli>(
        Clock::now() - begin_).count();
}

inline void GlobalEnergyPcgPath::advance_to(int target_steps) {
    if (target_steps == steps_) return;
    const auto solve_begin = Clock::now();
    std::atomic<int> next_column{0};
    std::atomic<long long> performed{0};
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
                    performed.fetch_add(1, std::memory_order_relaxed);
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
    total_column_iterations_ += performed.load(std::memory_order_relaxed);
    solve_ms_ += std::chrono::duration<double, std::milli>(
        Clock::now() - solve_begin).count();
}

inline SparseMatrix GlobalEnergyPcgPath::prolongation(
    double drop_tolerance) {
    const auto finalize_begin = Clock::now();
    std::vector<Triplet> entries;
    entries.reserve(
        static_cast<std::size_t>(grid_.coarse_size()) *
        (system_.f_nodes.size() + 1U));
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
    finalize_ms_ += std::chrono::duration<double, std::milli>(
        Clock::now() - finalize_begin).count();
    return result;
}

inline GlobalPcgPathReport GlobalEnergyPcgPath::report() const {
    GlobalPcgPathReport result;
    result.steps = steps_;
    result.total_column_iterations = total_column_iterations_;
    result.assembly_ms = system_.assembly_ms;
    result.solve_ms = solve_ms_;
    result.finalize_ms = finalize_ms_;
    result.total_ms = result.assembly_ms + result.solve_ms;
    double squared_sum = 0.0;
    double rz_sum = 0.0;
    double initial_rz_sum = 0.0;
    for (const ColumnState& state : columns_) {
        if (state.active) ++result.active_columns;
        const double relative = norm2(state.residual) / state.residual_scale;
        result.maximum_relative_residual = std::max(
            result.maximum_relative_residual, relative);
        squared_sum += relative * relative;
        rz_sum += std::max(state.rz, 0.0);
        initial_rz_sum += state.initial_rz;
    }
    result.rms_relative_residual = std::sqrt(
        squared_sum / static_cast<double>(columns_.size()));
    result.relative_preconditioned_residual = std::sqrt(
        rz_sum / initial_rz_sum);
    return result;
}

}
