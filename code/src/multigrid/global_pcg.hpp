#pragma once

#include "multigrid/energy_interpolation.hpp"
#include "multigrid/two_grid_solver.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

namespace tgi {

struct GlobalPcgPathReport {
    int systems = 0;
    int total_iterations = 0;
    int minimum_iterations = 0;
    int maximum_iterations = 0;
    int failed_systems = 0;
    double maximum_relative_residual = 0.0;
};

struct ColumnPropagationReport {
    int coarse_column = 0;
    int iterations = 0;
    int initial_residual_support = 0;
    int correction_support = 0;
    int reachable_support = 0;
    int f_unknowns = 0;
    int maximum_correction_distance = -1;
    int theoretical_distance_bound = -1;
    int bound_violations = 0;
};

class GlobalEnergyPcgPath {
public:
    GlobalEnergyPcgPath(
        const StructuredGrid& grid, const SparseMatrix& a,
        const SparseMatrix& initial_prolongation, int thread_count = 1);

    void advance_to(int target_steps);
    GlobalPcgPathReport advance_until_relative_residual(
        double tolerance, int maximum_steps = 40000);
    GlobalPcgPathReport report(double tolerance = 0.0) const;
    ColumnPropagationReport column_propagation_report(
        int coarse_column, double zero_tolerance = 0.0) const;
    SparseMatrix prolongation();

private:
    struct ColumnState {
        Vector solution;
        Vector residual;
        Vector direction;
        double rz = 0.0;
        double initial_residual_norm = 1.0;
        std::vector<std::pair<int, double>> initial_nonzeros;
        std::vector<int> initial_residual_support;
        int iterations = 0;
        bool active = true;
    };

    bool advance_one_iteration(
        ColumnState& state, Vector& product, Vector& z) const;

    const StructuredGrid& grid_;
    energy_interpolation_detail::GlobalFSystem system_;
    std::vector<ColumnState> columns_;
    int thread_count_ = 1;
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
            state.solution[static_cast<std::size_t>(local)] =
                initial_transpose.values()[static_cast<std::size_t>(position)];
            state.initial_nonzeros.push_back({
                local,
                initial_transpose.values()[static_cast<std::size_t>(position)]});
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
        for (std::size_t local = 0; local < state.residual.size(); ++local) {
            if (state.residual[local] != 0.0) {
                state.initial_residual_support.push_back(
                    static_cast<int>(local));
            }
        }
        const double residual_scale = std::max(norm2(rhs), 1.0e-30);
        state.initial_residual_norm = residual_scale;
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

inline bool GlobalEnergyPcgPath::advance_one_iteration(
    ColumnState& state, Vector& product, Vector& z) const {
    system_.matrix.multiply(state.direction, product);
    const double denominator = dot(state.direction, product);
    if (!(denominator > 0.0) || !std::isfinite(denominator) ||
        !(state.rz > 0.0) || !std::isfinite(state.rz)) {
        throw std::runtime_error("global PCG path lost positive curvature");
    }
    const double alpha = state.rz / denominator;
    axpy(alpha, state.direction, state.solution);
    axpy(-alpha, product, state.residual);
    ++state.iterations;
    for (std::size_t index = 0; index < z.size(); ++index) {
        z[index] = system_.inverse_diagonal[index] * state.residual[index];
    }
    const double rz_new = dot(state.residual, z);
    if (!(rz_new > 0.0) || !std::isfinite(rz_new)) {
        state.rz = rz_new;
        state.active = false;
        const double threshold = std::numeric_limits<double>::epsilon() *
            state.initial_residual_norm;
        if (!std::isfinite(rz_new) || norm2(state.residual) > threshold) {
            throw std::runtime_error(
                "global PCG path broke down before convergence");
        }
        return false;
    }
    const double beta = rz_new / state.rz;
    for (std::size_t index = 0; index < state.direction.size(); ++index) {
        state.direction[index] = z[index] + beta * state.direction[index];
    }
    state.rz = rz_new;
    return true;
}

inline void GlobalEnergyPcgPath::advance_to(int target_steps) {
    if (target_steps < 0) {
        throw std::invalid_argument("global PCG step count must be nonnegative");
    }
    for (const ColumnState& state : columns_) {
        if (state.iterations > target_steps) {
            throw std::invalid_argument("global PCG path cannot move backward");
        }
    }
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
                    if (!advance_one_iteration(state, product, z)) break;
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
}

inline GlobalPcgPathReport GlobalEnergyPcgPath::report(
    double tolerance) const {
    GlobalPcgPathReport value;
    value.systems = static_cast<int>(columns_.size());
    value.minimum_iterations = std::numeric_limits<int>::max();
    for (const ColumnState& state : columns_) {
        value.total_iterations += state.iterations;
        value.minimum_iterations = std::min(
            value.minimum_iterations, state.iterations);
        value.maximum_iterations = std::max(
            value.maximum_iterations, state.iterations);
        const double relative =
            norm2(state.residual) / state.initial_residual_norm;
        value.maximum_relative_residual = std::max(
            value.maximum_relative_residual, relative);
        if (tolerance > 0.0 && relative > tolerance) {
            ++value.failed_systems;
        }
    }
    return value;
}

inline GlobalPcgPathReport
GlobalEnergyPcgPath::advance_until_relative_residual(
    double tolerance, int maximum_steps) {
    if (!(tolerance > 0.0) || !std::isfinite(tolerance)) {
        throw std::invalid_argument(
            "global PCG tolerance must be finite and positive");
    }
    if (maximum_steps < 0) {
        throw std::invalid_argument(
            "global PCG step limit must be nonnegative");
    }
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
                const double target = tolerance * state.initial_residual_norm;
                while (state.active && state.iterations < maximum_steps &&
                       norm2(state.residual) > target) {
                    if (!advance_one_iteration(state, product, z)) break;
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
    const GlobalPcgPathReport value = report(tolerance);
    if (value.failed_systems != 0) {
        throw std::runtime_error(
            "one or more PCG columns missed the residual tolerance");
    }
    return value;
}

inline ColumnPropagationReport
GlobalEnergyPcgPath::column_propagation_report(
    int coarse_column, double zero_tolerance) const {
    if (coarse_column < 0 || coarse_column >= grid_.coarse_size()) {
        throw std::out_of_range("coarse column is outside the interpolation");
    }
    if (!(zero_tolerance >= 0.0) || !std::isfinite(zero_tolerance)) {
        throw std::invalid_argument(
            "support zero tolerance must be finite and nonnegative");
    }
    const ColumnState& state =
        columns_[static_cast<std::size_t>(coarse_column)];
    const int f_size = static_cast<int>(state.solution.size());
    std::vector<double> correction = state.solution;
    for (const auto& [local, value] : state.initial_nonzeros) {
        correction[static_cast<std::size_t>(local)] -= value;
    }

    std::vector<int> distance(static_cast<std::size_t>(f_size), -1);
    std::vector<int> queue;
    queue.reserve(static_cast<std::size_t>(f_size));
    for (int local : state.initial_residual_support) {
        if (distance[static_cast<std::size_t>(local)] == -1) {
            distance[static_cast<std::size_t>(local)] = 0;
            queue.push_back(local);
        }
    }
    for (std::size_t head = 0; head < queue.size(); ++head) {
        const int row = queue[head];
        const int next_distance =
            distance[static_cast<std::size_t>(row)] + 1;
        for (int position =
                 system_.matrix.row_ptr()[static_cast<std::size_t>(row)];
             position < system_.matrix.row_ptr()[
                 static_cast<std::size_t>(row) + 1U]; ++position) {
            const int column = system_.matrix.col_idx()[
                static_cast<std::size_t>(position)];
            if (column == row ||
                distance[static_cast<std::size_t>(column)] != -1) {
                continue;
            }
            distance[static_cast<std::size_t>(column)] = next_distance;
            queue.push_back(column);
        }
    }

    ColumnPropagationReport report;
    report.coarse_column = coarse_column;
    report.iterations = state.iterations;
    report.initial_residual_support =
        static_cast<int>(state.initial_residual_support.size());
    report.f_unknowns = f_size;
    report.theoretical_distance_bound = state.iterations - 1;
    for (int local = 0; local < f_size; ++local) {
        const int graph_distance = distance[static_cast<std::size_t>(local)];
        if (graph_distance >= 0 &&
            graph_distance <= report.theoretical_distance_bound) {
            ++report.reachable_support;
        }
        if (std::abs(correction[static_cast<std::size_t>(local)]) <=
            zero_tolerance) {
            continue;
        }
        ++report.correction_support;
        report.maximum_correction_distance = std::max(
            report.maximum_correction_distance, graph_distance);
        if (graph_distance < 0 ||
            graph_distance > report.theoretical_distance_bound) {
            ++report.bound_violations;
        }
    }
    return report;
}

inline SparseMatrix GlobalEnergyPcgPath::prolongation() {
    std::vector<Triplet> entries;
    std::size_t entry_count = static_cast<std::size_t>(grid_.coarse_size());
    for (const ColumnState& state : columns_) {
        entry_count += static_cast<std::size_t>(std::count_if(
            state.solution.begin(), state.solution.end(),
            [](double value) { return value != 0.0; }));
    }
    entries.reserve(entry_count);
    for (int coarse = 0; coarse < grid_.coarse_size(); ++coarse) {
        entries.push_back({grid_.coarse_fine_id(coarse), coarse, 1.0});
        const Vector& weights =
            columns_[static_cast<std::size_t>(coarse)].solution;
        for (std::size_t local = 0; local < weights.size(); ++local) {
            if (weights[local] != 0.0) {
                entries.push_back(
                    {system_.f_nodes[local], coarse, weights[local]});
            }
        }
    }
    return SparseMatrix(
        grid_.fine_size(), grid_.coarse_size(), entries, 0.0);
}

}
