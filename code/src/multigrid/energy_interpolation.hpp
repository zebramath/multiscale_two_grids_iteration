#pragma once

#include "pde/diffusion_problem.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <exception>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

namespace tgi {

struct GlobalEnergyOptions {
    double tolerance = 1.0e-10;
    int maximum_iterations = 40000;
    int thread_count = 1;
    double drop_tolerance = 0.0;
    bool require_convergence = true;
};

struct GlobalSolveStats {
    int systems = 0;
    int total_iterations = 0;
    int maximum_iterations = 0;
    int failed_systems = 0;
    double maximum_relative_residual = 0.0;
};

struct InterpolationReport {
    GlobalSolveStats column_solves;
    int threads_used = 1;
};

struct InterpolationResult {
    SparseMatrix prolongation;
    InterpolationReport report;
};

namespace energy_interpolation_detail {

struct ConjugateGradientStats {
    int iterations = 0;
    bool converged = false;
    double relative_residual = 0.0;
};

inline Vector solve_pcg(
    const SparseMatrix& matrix, const Vector& rhs, double tolerance,
    int maximum_iterations, ConjugateGradientStats& stats,
    const Vector* initial = nullptr,
    const Vector* inverse_diagonal_override = nullptr) {
    const int n = matrix.rows();
    Vector x = initial != nullptr
        ? *initial
        : Vector(static_cast<std::size_t>(n), 0.0);
    if (n == 0) {
        stats.converged = true;
        return x;
    }

    Vector inverse_diagonal_storage;
    const Vector* inverse_diagonal = inverse_diagonal_override;
    if (inverse_diagonal == nullptr) {
        const Vector diagonal = matrix.diagonal();
        inverse_diagonal_storage.resize(diagonal.size());
        for (std::size_t index = 0; index < diagonal.size(); ++index) {
            inverse_diagonal_storage[index] = 1.0 / diagonal[index];
        }
        inverse_diagonal = &inverse_diagonal_storage;
    }

    Vector residual = rhs;
    Vector product;
    if (initial != nullptr) {
        matrix.multiply(x, product);
        axpy(-1.0, product, residual);
    }
    const double residual_scale = std::max(norm2(rhs), 1.0e-30);
    const bool fixed_budget = tolerance == 0.0;
    const double target_squared = fixed_budget
        ? -1.0
        : tolerance * tolerance * residual_scale * residual_scale;
    double residual_squared = dot(residual, residual);
    if (residual_squared <= target_squared) {
        stats.converged = true;
        stats.relative_residual =
            std::sqrt(residual_squared) / residual_scale;
        return x;
    }

    Vector z(residual.size());
    for (std::size_t index = 0; index < z.size(); ++index) {
        z[index] = (*inverse_diagonal)[index] * residual[index];
    }
    Vector direction = z;
    Vector action;
    double rz = dot(residual, z);
    if (!(rz > 0.0) || !std::isfinite(rz)) {
        stats.relative_residual =
            std::sqrt(residual_squared) / residual_scale;
        return x;
    }
    for (int iteration = 0; iteration < maximum_iterations; ++iteration) {
        matrix.multiply(direction, action);
        const double denominator = dot(direction, action);
        if (!(denominator > 0.0) || !std::isfinite(denominator)) break;
        const double alpha = rz / denominator;
        axpy(alpha, direction, x);
        axpy(-alpha, action, residual);
        stats.iterations = iteration + 1;
        residual_squared = dot(residual, residual);
        if (residual_squared <= target_squared) {
            matrix.multiply(x, product);
            residual = rhs;
            axpy(-1.0, product, residual);
            residual_squared = dot(residual, residual);
            if (residual_squared <= target_squared) {
                stats.converged = true;
                break;
            }
            for (std::size_t index = 0; index < z.size(); ++index) {
                z[index] = (*inverse_diagonal)[index] * residual[index];
            }
            rz = dot(residual, z);
            if (!(rz > 0.0) || !std::isfinite(rz)) break;
            direction = z;
            continue;
        }
        for (std::size_t index = 0; index < z.size(); ++index) {
            z[index] = (*inverse_diagonal)[index] * residual[index];
        }
        const double rz_new = dot(residual, z);
        if (!(rz_new > 0.0) || !std::isfinite(rz_new)) break;
        const double beta = rz_new / rz;
        for (std::size_t index = 0; index < direction.size(); ++index) {
            direction[index] = z[index] + beta * direction[index];
        }
        rz = rz_new;
    }
    matrix.multiply(x, product);
    residual = rhs;
    axpy(-1.0, product, residual);
    residual_squared = dot(residual, residual);
    if (!fixed_budget && residual_squared <= target_squared) {
        stats.converged = true;
    }
    stats.relative_residual =
        std::sqrt(residual_squared) / residual_scale;
    return x;
}

struct ColumnResult {
    std::vector<Triplet> triplets;
    ConjugateGradientStats solve;
};

struct GlobalFSystem {
    SparseMatrix matrix;
    std::vector<int> f_nodes;
    std::vector<int> local_index;
    Vector inverse_diagonal;
    std::vector<std::vector<std::pair<int, double>>> rhs_entries;
};

inline int coarse_id_from_fine_node(
    const StructuredGrid& grid, int fine) {
    const auto [ix, iy] = grid.fine_coords(fine);
    return grid.coarse_id(
        (ix + 1) / grid.ratio() - 1,
        (iy + 1) / grid.ratio() - 1);
}

inline GlobalFSystem assemble_global_f_system(
    const StructuredGrid& grid, const SparseMatrix& matrix) {
    GlobalFSystem system;
    system.f_nodes = grid.all_f_nodes();
    system.rhs_entries.resize(
        static_cast<std::size_t>(grid.coarse_size()));
    system.local_index.assign(
        static_cast<std::size_t>(grid.fine_size()), -1);
    for (std::size_t local = 0; local < system.f_nodes.size(); ++local) {
        system.local_index[static_cast<std::size_t>(system.f_nodes[local])] =
            static_cast<int>(local);
    }

    std::vector<int> row_ptr(system.f_nodes.size() + 1U, 0);
    std::vector<int> col_idx;
    Vector values;
    col_idx.reserve(system.f_nodes.size() * 5U);
    values.reserve(system.f_nodes.size() * 5U);
    for (std::size_t local_row = 0;
         local_row < system.f_nodes.size(); ++local_row) {
        const int global_row = system.f_nodes[local_row];
        for (int position =
                 matrix.row_ptr()[static_cast<std::size_t>(global_row)];
             position < matrix.row_ptr()[
                 static_cast<std::size_t>(global_row) + 1U];
             ++position) {
            const int global_col =
                matrix.col_idx()[static_cast<std::size_t>(position)];
            const double value =
                matrix.values()[static_cast<std::size_t>(position)];
            const int local_col =
                system.local_index[static_cast<std::size_t>(global_col)];
            if (local_col >= 0) {
                col_idx.push_back(local_col);
                values.push_back(value);
            } else if (grid.is_coarse_node(global_col)) {
                const int coarse =
                    coarse_id_from_fine_node(grid, global_col);
                system.rhs_entries[static_cast<std::size_t>(coarse)]
                    .push_back({static_cast<int>(local_row), -value});
            }
        }
        row_ptr[local_row + 1U] = static_cast<int>(values.size());
    }

    system.matrix = SparseMatrix(
        static_cast<int>(system.f_nodes.size()),
        static_cast<int>(system.f_nodes.size()), std::move(row_ptr),
        std::move(col_idx), std::move(values));
    const Vector diagonal = system.matrix.diagonal();
    system.inverse_diagonal.resize(diagonal.size());
    for (std::size_t index = 0; index < diagonal.size(); ++index) {
        system.inverse_diagonal[index] = 1.0 / diagonal[index];
    }
    return system;
}

inline SparseMatrix assemble_prolongation(
    const StructuredGrid& grid, const std::vector<ColumnResult>& columns) {
    std::size_t entry_count = 0;
    for (const ColumnResult& column : columns) {
        entry_count += column.triplets.size();
    }
    std::vector<int> row_ptr(
        static_cast<std::size_t>(grid.fine_size()) + 1U, 0);
    for (const ColumnResult& column : columns) {
        for (const Triplet& entry : column.triplets) {
            ++row_ptr[static_cast<std::size_t>(entry.row) + 1U];
        }
    }
    for (int row = 0; row < grid.fine_size(); ++row) {
        row_ptr[static_cast<std::size_t>(row) + 1U] +=
            row_ptr[static_cast<std::size_t>(row)];
    }
    std::vector<int> next = row_ptr;
    std::vector<int> col_idx(entry_count);
    Vector values(entry_count);
    for (const ColumnResult& column : columns) {
        for (const Triplet& entry : column.triplets) {
            const int target = next[static_cast<std::size_t>(entry.row)]++;
            col_idx[static_cast<std::size_t>(target)] = entry.column;
            values[static_cast<std::size_t>(target)] = entry.value;
        }
    }
    return SparseMatrix(
        grid.fine_size(), grid.coarse_size(), std::move(row_ptr),
        std::move(col_idx), std::move(values));
}

inline InterpolationResult solve_global_energy_columns(
    const StructuredGrid& grid, const SparseMatrix& matrix,
    const GlobalEnergyOptions& options,
    const SparseMatrix* initial_transpose) {
    const GlobalFSystem system = assemble_global_f_system(grid, matrix);
    unsigned int requested = options.thread_count > 0
        ? static_cast<unsigned int>(options.thread_count)
        : std::thread::hardware_concurrency();
    if (requested == 0U) requested = 1U;
    const int thread_count = std::max(
        1, std::min(grid.coarse_size(), static_cast<int>(requested)));

    std::vector<ColumnResult> columns(
        static_cast<std::size_t>(grid.coarse_size()));
    std::atomic<int> next_coarse{0};
    std::exception_ptr worker_error;
    std::mutex error_mutex;
    auto worker = [&]() {
        Vector rhs(system.f_nodes.size(), 0.0);
        Vector initial(system.f_nodes.size(), 0.0);
        try {
            while (true) {
                const int coarse = next_coarse.fetch_add(
                    1, std::memory_order_relaxed);
                if (coarse >= grid.coarse_size()) break;
                std::fill(rhs.begin(), rhs.end(), 0.0);
                for (const auto& [row, value] :
                     system.rhs_entries[static_cast<std::size_t>(coarse)]) {
                    rhs[static_cast<std::size_t>(row)] += value;
                }

                const Vector* initial_pointer = nullptr;
                if (initial_transpose != nullptr) {
                    std::fill(initial.begin(), initial.end(), 0.0);
                    for (int position = initial_transpose->row_ptr()[
                             static_cast<std::size_t>(coarse)];
                         position < initial_transpose->row_ptr()[
                             static_cast<std::size_t>(coarse) + 1U];
                         ++position) {
                        const int fine = initial_transpose->col_idx()[
                            static_cast<std::size_t>(position)];
                        if (grid.is_coarse_node(fine)) continue;
                        const int local = system.local_index[
                            static_cast<std::size_t>(fine)];
                        if (local >= 0) {
                            initial[static_cast<std::size_t>(local)] =
                                initial_transpose->values()[
                                    static_cast<std::size_t>(position)];
                        }
                    }
                    initial_pointer = &initial;
                }

                ColumnResult result;
                const Vector weights = solve_pcg(
                    system.matrix, rhs, options.tolerance,
                    options.maximum_iterations, result.solve,
                    initial_pointer, &system.inverse_diagonal);
                result.triplets.reserve(weights.size() + 1U);
                result.triplets.push_back(
                    {grid.coarse_fine_id(coarse), coarse, 1.0});
                for (std::size_t local = 0; local < weights.size(); ++local) {
                    if (std::abs(weights[local]) > options.drop_tolerance) {
                        result.triplets.push_back(
                            {system.f_nodes[local], coarse, weights[local]});
                    }
                }
                columns[static_cast<std::size_t>(coarse)] =
                    std::move(result);
            }
        } catch (...) {
            std::lock_guard<std::mutex> lock(error_mutex);
            if (!worker_error) worker_error = std::current_exception();
            next_coarse.store(grid.coarse_size());
        }
    };

    if (thread_count == 1) {
        worker();
    } else {
        std::vector<std::thread> workers;
        workers.reserve(static_cast<std::size_t>(thread_count));
        for (int index = 0; index < thread_count; ++index) {
            workers.emplace_back(worker);
        }
        for (auto& thread : workers) thread.join();
    }
    if (worker_error) std::rethrow_exception(worker_error);
    InterpolationReport report;
    report.threads_used = thread_count;
    for (const ColumnResult& column : columns) {
        ++report.column_solves.systems;
        report.column_solves.total_iterations += column.solve.iterations;
        report.column_solves.maximum_iterations = std::max(
            report.column_solves.maximum_iterations,
            column.solve.iterations);
        report.column_solves.maximum_relative_residual = std::max(
            report.column_solves.maximum_relative_residual,
            column.solve.relative_residual);
        if (!column.solve.converged) {
            ++report.column_solves.failed_systems;
        }
    }
    if (options.require_convergence &&
        report.column_solves.failed_systems != 0) {
        throw std::runtime_error(
            "one or more global energy solves did not converge");
    }

    SparseMatrix prolongation = assemble_prolongation(grid, columns);
    return {std::move(prolongation), report};
}

}

inline InterpolationResult build_geometric_interpolation(
    const StructuredGrid& grid) {
    std::vector<int> row_ptr(
        static_cast<std::size_t>(grid.fine_size()) + 1U, 0);
    std::vector<int> col_idx;
    Vector values;
    col_idx.reserve(static_cast<std::size_t>(4 * grid.fine_size()));
    values.reserve(static_cast<std::size_t>(4 * grid.fine_size()));
    const int coarse_intervals = grid.intervals() / grid.ratio();

    for (int fine = 0; fine < grid.fine_size(); ++fine) {
        const auto [ix, iy] = grid.fine_coords(fine);
        const int lattice_x = ix + 1;
        const int lattice_y = iy + 1;
        const int left_x = lattice_x / grid.ratio();
        const int left_y = lattice_y / grid.ratio();
        const double tx = static_cast<double>(lattice_x % grid.ratio()) /
            static_cast<double>(grid.ratio());
        const double ty = static_cast<double>(lattice_y % grid.ratio()) /
            static_cast<double>(grid.ratio());
        const int qx[2] = {left_x, left_x + 1};
        const int qy[2] = {left_y, left_y + 1};
        const double wx[2] = {1.0 - tx, tx};
        const double wy[2] = {1.0 - ty, ty};

        for (int ay = 0; ay < 2; ++ay) {
            for (int ax = 0; ax < 2; ++ax) {
                if (wx[ax] == 0.0 || wy[ay] == 0.0 ||
                    qx[ax] <= 0 || qx[ax] >= coarse_intervals ||
                    qy[ay] <= 0 || qy[ay] >= coarse_intervals) {
                    continue;
                }
                col_idx.push_back(
                    grid.coarse_id(qx[ax] - 1, qy[ay] - 1));
                values.push_back(wx[ax] * wy[ay]);
            }
        }
        row_ptr[static_cast<std::size_t>(fine) + 1U] =
            static_cast<int>(values.size());
    }

    SparseMatrix prolongation(
        grid.fine_size(), grid.coarse_size(), std::move(row_ptr),
        std::move(col_idx), std::move(values));
    return {std::move(prolongation), {}};
}

inline InterpolationResult build_global_energy_interpolation(
    const StructuredGrid& grid, const SparseMatrix& matrix,
    const GlobalEnergyOptions& options = {}) {
    return energy_interpolation_detail::solve_global_energy_columns(
        grid, matrix, options, nullptr);
}

inline InterpolationResult refine_global_energy_interpolation(
    const StructuredGrid& grid, const SparseMatrix& matrix,
    const SparseMatrix& initial, const GlobalEnergyOptions& options) {
    const SparseMatrix initial_transpose =
        initial.transpose(options.thread_count);
    return energy_interpolation_detail::solve_global_energy_columns(
        grid, matrix, options, &initial_transpose);
}

}
