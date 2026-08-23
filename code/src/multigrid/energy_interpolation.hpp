#pragma once

#include "pde/diffusion_problem.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <exception>
#include <iterator>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

namespace tgi {

enum class InterpolationStrategy {
    GeometricBilinear,
    GlobalEnergyMinimum,
    LocalEnergyMinimum
};

struct InterpolationOptions {
    InterpolationStrategy strategy = InterpolationStrategy::LocalEnergyMinimum;
    int patch_layers = 3;
    double local_tolerance = 1e-3;
    int local_max_iterations = 4000;
    int thread_count = 1;
    double drop_tolerance = 1e-14;
    bool require_convergence = true;
};

struct LocalSolveStats {
    int systems = 0;
    int total_iterations = 0;
    int max_iterations = 0;
    int failed_systems = 0;
    double sum_relative_residual = 0.0;
    double max_relative_residual = 0.0;
};

struct InterpolationTiming {

    double patch_assembly_work_ms = 0.0;
    double local_solve_work_ms = 0.0;
    double basis_scatter_work_ms = 0.0;
    double parallel_wall_ms = 0.0;
    double finalize_ms = 0.0;
    double total_ms = 0.0;
};

struct InterpolationReport {
    LocalSolveStats local_solves;
    InterpolationTiming timing;
    int threads_used = 1;
};

struct InterpolationResult {
    SparseMatrix prolongation;
    InterpolationReport report;
};

inline InterpolationResult build_interpolation(const StructuredGrid& grid,
                                        const SparseMatrix& a,
                                        const InterpolationOptions& options);
inline InterpolationResult build_energy_interpolation_on_supports(
    const StructuredGrid& grid, const SparseMatrix& a,
    const std::vector<std::vector<int>>& supports,
    const InterpolationOptions& options);
inline InterpolationResult refine_energy_interpolation_on_supports(
    const StructuredGrid& grid, const SparseMatrix& a,
    const std::vector<std::vector<int>>& supports,
    const SparseMatrix& initial,
    const InterpolationOptions& options);
inline InterpolationResult refine_selected_energy_interpolation_on_supports(
    const StructuredGrid& grid, const SparseMatrix& a,
    const std::vector<std::vector<int>>& supports,
    const std::vector<unsigned char>& refine_column,
    const SparseMatrix& initial,
    const InterpolationOptions& options);
inline InterpolationResult refine_global_energy_interpolation(
    const StructuredGrid& grid, const SparseMatrix& a,
    const SparseMatrix& initial, const InterpolationOptions& options);


namespace energy_interpolation_detail {

using Clock = std::chrono::steady_clock;

inline double milliseconds(Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

struct ConjugateGradientStats {
    int iterations = 0;
    bool converged = false;
    double relative_residual = 0.0;
};

inline Vector solve_local_cg(const SparseMatrix& matrix, const Vector& rhs,
                      double tolerance, int max_iterations, ConjugateGradientStats& stats,
                      const Vector* initial = nullptr,
                      const Vector* inverse_diagonal_override = nullptr) {
    const int n = matrix.rows();
    Vector x = initial != nullptr
        ? *initial
        : Vector(static_cast<std::size_t>(n), 0.0);
    if (x.size() != static_cast<std::size_t>(n)) {
        throw std::invalid_argument(
            "solve_local_cg: initial vector has wrong size");
    }
    const double rhs_norm = norm2(rhs);
    if (n == 0) {
        stats.converged = true;
        return x;
    }

    Vector inverse_diagonal_storage;
    const Vector* inverse_diagonal = inverse_diagonal_override;
    if (inverse_diagonal == nullptr) {
        const Vector diagonal = matrix.diagonal();
        inverse_diagonal_storage.resize(diagonal.size());
        for (std::size_t i = 0; i < diagonal.size(); ++i) {
            inverse_diagonal_storage[i] = 1.0 / diagonal[i];
        }
        inverse_diagonal = &inverse_diagonal_storage;
    } else if (inverse_diagonal->size() != static_cast<std::size_t>(n)) {
        throw std::invalid_argument(
            "solve_local_cg: inverse diagonal has wrong size");
    }

    Vector residual = rhs;
    Vector product;
    if (initial != nullptr) {
        matrix.multiply(x, product);
        axpy(-1.0, product, residual);
    }
    const double residual_scale = std::max(rhs_norm, 1.0e-30);
    const double target = tolerance * residual_scale;
    const double target_squared = target * target;
    double residual_squared = dot(residual, residual);
    if (residual_squared <= target_squared) {
        stats.converged = true;
        stats.relative_residual =
            std::sqrt(residual_squared) / residual_scale;
        return x;
    }
    Vector z(residual.size());
    for (std::size_t i = 0; i < z.size(); ++i) {
        z[i] = (*inverse_diagonal)[i] * residual[i];
    }
    Vector direction = z;
    Vector ad;
    double rz = dot(residual, z);
    for (int iteration = 1; iteration <= max_iterations; ++iteration) {
        matrix.multiply(direction, ad);
        const double denominator = dot(direction, ad);
        if (!(denominator > 0.0) || !std::isfinite(denominator)) break;
        const double alpha = rz / denominator;
        axpy(alpha, direction, x);
        axpy(-alpha, ad, residual);
        stats.iterations = iteration;
        residual_squared = dot(residual, residual);
        if (residual_squared <= target_squared) {
            stats.converged = true;
            break;
        }
        for (std::size_t i = 0; i < z.size(); ++i) {
            z[i] = (*inverse_diagonal)[i] * residual[i];
        }
        const double rz_new = dot(residual, z);
        const double beta = rz_new / rz;
        for (std::size_t i = 0; i < direction.size(); ++i) {
            direction[i] = z[i] + beta * direction[i];
        }
        rz = rz_new;
    }
    stats.relative_residual =
        std::sqrt(residual_squared) / residual_scale;
    return x;
}

struct LocalBasisResult {
    std::vector<Triplet> triplets;
    ConjugateGradientStats solve;
    double patch_assembly_ms = 0.0;
    double local_solve_ms = 0.0;
    double scatter_ms = 0.0;
};

struct GlobalFSystem {
    SparseMatrix matrix;
    std::vector<int> f_nodes;
    std::vector<int> local_index;
    Vector inverse_diagonal;
    std::vector<std::vector<std::pair<int, double>>> rhs_entries;
    double assembly_ms = 0.0;
};

inline int coarse_id_from_fine_node(
    const StructuredGrid& grid, int fine) {
    const auto [ix, iy] = grid.fine_coords(fine);
    return grid.coarse_id(
        (ix + 1) / grid.ratio() - 1,
        (iy + 1) / grid.ratio() - 1);
}

inline GlobalFSystem assemble_global_f_system(
    const StructuredGrid& grid, const SparseMatrix& a) {
    const auto begin = Clock::now();
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
                 a.row_ptr()[static_cast<std::size_t>(global_row)];
             position <
                 a.row_ptr()[static_cast<std::size_t>(global_row) + 1U];
             ++position) {
            const int global_col =
                a.col_idx()[static_cast<std::size_t>(position)];
            const double value =
                a.values()[static_cast<std::size_t>(position)];
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
    system.assembly_ms = milliseconds(begin, Clock::now());
    return system;
}

inline InterpolationResult build_global_energy_interpolation(
    const StructuredGrid& grid, const SparseMatrix& a,
    const InterpolationOptions& options,
    const SparseMatrix* initial_transpose = nullptr) {
    const auto total_begin = Clock::now();
    const GlobalFSystem system = assemble_global_f_system(grid, a);
    unsigned int requested = options.thread_count > 0
        ? static_cast<unsigned int>(options.thread_count)
        : std::thread::hardware_concurrency();
    if (requested == 0U) requested = 1U;
    const int thread_count = std::max(
        1, std::min(grid.coarse_size(), static_cast<int>(requested)));
    std::vector<LocalBasisResult> bases(
        static_cast<std::size_t>(grid.coarse_size()));
    std::atomic<int> next_coarse{0};
    std::exception_ptr worker_error;
    std::mutex error_mutex;

    const auto solve_begin = Clock::now();
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
                    for (int position =
                             initial_transpose->row_ptr()[
                                 static_cast<std::size_t>(coarse)];
                         position <
                             initial_transpose->row_ptr()[
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
                LocalBasisResult basis;
                const Vector weights = solve_local_cg(
                    system.matrix, rhs, options.local_tolerance,
                    options.local_max_iterations, basis.solve,
                    initial_pointer, &system.inverse_diagonal);
                basis.triplets.reserve(weights.size() + 1U);
                basis.triplets.push_back(
                    {grid.coarse_fine_id(coarse), coarse, 1.0});
                for (std::size_t local = 0; local < weights.size(); ++local) {
                    if (std::abs(weights[local]) > options.drop_tolerance) {
                        basis.triplets.push_back(
                            {system.f_nodes[local], coarse, weights[local]});
                    }
                }
                bases[static_cast<std::size_t>(coarse)] = std::move(basis);
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
    const auto solve_end = Clock::now();

    InterpolationReport report;
    report.threads_used = thread_count;
    report.timing.patch_assembly_work_ms = system.assembly_ms;
    report.timing.parallel_wall_ms = milliseconds(solve_begin, solve_end);
    std::size_t total_entries = 0;
    for (const LocalBasisResult& basis : bases) {
        ++report.local_solves.systems;
        report.local_solves.total_iterations += basis.solve.iterations;
        report.local_solves.max_iterations = std::max(
            report.local_solves.max_iterations, basis.solve.iterations);
        report.local_solves.sum_relative_residual +=
            basis.solve.relative_residual;
        report.local_solves.max_relative_residual = std::max(
            report.local_solves.max_relative_residual,
            basis.solve.relative_residual);
        if (!basis.solve.converged) ++report.local_solves.failed_systems;
        total_entries += basis.triplets.size();
    }
    if (options.require_convergence &&
        report.local_solves.failed_systems != 0) {
        throw std::runtime_error(
            "one or more global energy solves did not converge");
    }

    const auto finalize_begin = Clock::now();
    std::vector<int> row_ptr(
        static_cast<std::size_t>(grid.fine_size()) + 1U, 0);
    for (const LocalBasisResult& basis : bases) {
        for (const Triplet& entry : basis.triplets) {
            ++row_ptr[static_cast<std::size_t>(entry.row) + 1U];
        }
    }
    for (int row = 0; row < grid.fine_size(); ++row) {
        row_ptr[static_cast<std::size_t>(row) + 1U] +=
            row_ptr[static_cast<std::size_t>(row)];
    }
    std::vector<int> next = row_ptr;
    std::vector<int> col_idx(total_entries);
    Vector values(total_entries);
    for (const LocalBasisResult& basis : bases) {
        for (const Triplet& entry : basis.triplets) {
            const int target = next[static_cast<std::size_t>(entry.row)]++;
            col_idx[static_cast<std::size_t>(target)] = entry.column;
            values[static_cast<std::size_t>(target)] = entry.value;
        }
    }
    SparseMatrix prolongation(
        grid.fine_size(), grid.coarse_size(), std::move(row_ptr),
        std::move(col_idx), std::move(values));
    report.timing.finalize_ms = milliseconds(finalize_begin, Clock::now());
    report.timing.total_ms = milliseconds(total_begin, Clock::now());
    return {std::move(prolongation), report};
}

inline LocalBasisResult build_basis_on_nodes(
    const StructuredGrid& grid, const SparseMatrix& a, int coarse,
    const InterpolationOptions& options, std::vector<int>& local_index,
    const std::vector<int>& local_nodes,
    const SparseMatrix* initial_transpose = nullptr) {
    LocalBasisResult result;
    const auto patch_begin = Clock::now();
    const int coarse_fine = grid.coarse_fine_id(coarse);
    for (std::size_t local = 0; local < local_nodes.size(); ++local) {
        local_index[static_cast<std::size_t>(local_nodes[local])] =
            static_cast<int>(local);
    }

    std::vector<int> local_row_ptr(local_nodes.size() + 1U, 0);
    std::vector<int> local_col_idx;
    Vector local_values;
    local_col_idx.reserve(local_nodes.size() * 5U);
    local_values.reserve(local_nodes.size() * 5U);
    Vector rhs(local_nodes.size(), 0.0);
    for (std::size_t local_row = 0; local_row < local_nodes.size(); ++local_row) {
        const int global_row = local_nodes[local_row];
        for (int pos = a.row_ptr()[static_cast<std::size_t>(global_row)];
             pos < a.row_ptr()[static_cast<std::size_t>(global_row) + 1U]; ++pos) {
            const int global_col = a.col_idx()[static_cast<std::size_t>(pos)];
            const double value = a.values()[static_cast<std::size_t>(pos)];
            const int local_col = local_index[static_cast<std::size_t>(global_col)];
            if (local_col >= 0) {
                local_col_idx.push_back(local_col);
                local_values.push_back(value);
            } else if (global_col == coarse_fine) {
                rhs[local_row] -= value;
            }
        }
        local_row_ptr[local_row + 1U] =
            static_cast<int>(local_values.size());
    }
    const SparseMatrix local_matrix(static_cast<int>(local_nodes.size()),
                                    static_cast<int>(local_nodes.size()),
                                    std::move(local_row_ptr),
                                    std::move(local_col_idx),
                                    std::move(local_values));
    result.patch_assembly_ms = milliseconds(patch_begin, Clock::now());

    const auto solve_begin = Clock::now();
    Vector initial_weights;
    const Vector* initial_pointer = nullptr;
    if (initial_transpose != nullptr) {
        initial_weights.assign(local_nodes.size(), 0.0);
        for (int position =
                 initial_transpose->row_ptr()[
                     static_cast<std::size_t>(coarse)];
             position <
                 initial_transpose->row_ptr()[
                     static_cast<std::size_t>(coarse) + 1U];
             ++position) {
            const int fine = initial_transpose->col_idx()[
                static_cast<std::size_t>(position)];
            const int local =
                local_index[static_cast<std::size_t>(fine)];
            if (local >= 0) {
                initial_weights[static_cast<std::size_t>(local)] =
                    initial_transpose->values()[
                        static_cast<std::size_t>(position)];
            }
        }
        initial_pointer = &initial_weights;
    }
    const Vector weights = solve_local_cg(
        local_matrix, rhs, options.local_tolerance,
        options.local_max_iterations, result.solve,
        initial_pointer);
    result.local_solve_ms = milliseconds(solve_begin, Clock::now());

    const auto scatter_begin = Clock::now();
    result.triplets.reserve(local_nodes.size() + 1U);
    result.triplets.push_back({coarse_fine, coarse, 1.0});
    for (std::size_t local = 0; local < local_nodes.size(); ++local) {
        if (std::abs(weights[local]) > options.drop_tolerance) {
            result.triplets.push_back(
                {local_nodes[local], coarse, weights[local]});
        }
        local_index[static_cast<std::size_t>(local_nodes[local])] = -1;
    }
    result.scatter_ms = milliseconds(scatter_begin, Clock::now());
    return result;
}

inline LocalBasisResult build_local_basis(const StructuredGrid& grid, const SparseMatrix& a,
                              int coarse, const InterpolationOptions& options,
                              std::vector<int>& local_index,
                              std::vector<int>& patch_nodes,
                              const std::vector<int>& global_f_nodes) {
    const std::vector<int>* local_nodes = &global_f_nodes;
    if (options.patch_layers != 0) {
        grid.patch_f_nodes(coarse, options.patch_layers, patch_nodes);
        local_nodes = &patch_nodes;
    }
    return build_basis_on_nodes(
        grid, a, coarse, options, local_index, *local_nodes);
}

inline InterpolationResult build_local_energy_interpolation(
    const StructuredGrid& grid, const SparseMatrix& a,
    const InterpolationOptions& options) {
    const auto total_begin = Clock::now();
    unsigned int requested = options.thread_count > 0
        ? static_cast<unsigned int>(options.thread_count)
        : std::thread::hardware_concurrency();
    if (requested == 0U) requested = 1U;
    const int thread_count = std::min(
        grid.coarse_size(), static_cast<int>(requested));

    std::vector<LocalBasisResult> bases(
        static_cast<std::size_t>(grid.coarse_size()));
    const std::vector<int> global_f_nodes = options.patch_layers == 0
        ? grid.all_f_nodes()
        : std::vector<int>{};
    std::atomic<int> next_coarse{0};
    std::exception_ptr worker_error;
    std::mutex error_mutex;

    const auto parallel_begin = Clock::now();
    auto worker = [&]() {
        std::vector<int> local_index(
            static_cast<std::size_t>(grid.fine_size()), -1);
        std::vector<int> patch_nodes;
        try {
            while (true) {
                const int coarse = next_coarse.fetch_add(
                    1, std::memory_order_relaxed);
                if (coarse >= grid.coarse_size()) break;
                bases[static_cast<std::size_t>(coarse)] =
                    build_local_basis(
                        grid, a, coarse, options, local_index,
                        patch_nodes, global_f_nodes);
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
    const auto parallel_end = Clock::now();
    if (worker_error) std::rethrow_exception(worker_error);

    InterpolationReport report;
    report.threads_used = thread_count;
    report.timing.parallel_wall_ms =
        milliseconds(parallel_begin, parallel_end);
    std::size_t total_entries = 0;
    for (const LocalBasisResult& basis : bases) {
        ++report.local_solves.systems;
        report.local_solves.total_iterations += basis.solve.iterations;
        report.local_solves.max_iterations = std::max(
            report.local_solves.max_iterations, basis.solve.iterations);
        report.local_solves.sum_relative_residual +=
            basis.solve.relative_residual;
        report.local_solves.max_relative_residual = std::max(
            report.local_solves.max_relative_residual,
            basis.solve.relative_residual);
        if (!basis.solve.converged) ++report.local_solves.failed_systems;
        report.timing.patch_assembly_work_ms += basis.patch_assembly_ms;
        report.timing.local_solve_work_ms += basis.local_solve_ms;
        report.timing.basis_scatter_work_ms += basis.scatter_ms;
        total_entries += basis.triplets.size();
    }
    if (report.local_solves.failed_systems != 0) {
        throw std::runtime_error("one or more local energy solves did not converge");
    }

    const auto finalize_begin = Clock::now();
    std::vector<int> row_ptr(
        static_cast<std::size_t>(grid.fine_size()) + 1U, 0);
    for (const LocalBasisResult& basis : bases) {
        for (const Triplet& entry : basis.triplets) {
            ++row_ptr[static_cast<std::size_t>(entry.row) + 1U];
        }
    }
    for (int row = 0; row < grid.fine_size(); ++row) {
        row_ptr[static_cast<std::size_t>(row) + 1U] +=
            row_ptr[static_cast<std::size_t>(row)];
    }
    std::vector<int> next = row_ptr;
    std::vector<int> col_idx(total_entries);
    std::vector<double> values(total_entries);

    for (const LocalBasisResult& basis : bases) {
        for (const Triplet& entry : basis.triplets) {
            const int target = next[static_cast<std::size_t>(entry.row)]++;
            col_idx[static_cast<std::size_t>(target)] = entry.column;
            values[static_cast<std::size_t>(target)] = entry.value;
        }
    }
    SparseMatrix prolongation(
        grid.fine_size(), grid.coarse_size(), std::move(row_ptr),
        std::move(col_idx), std::move(values));
    report.timing.finalize_ms = milliseconds(finalize_begin, Clock::now());
    report.timing.total_ms = milliseconds(total_begin, Clock::now());
    return {std::move(prolongation), report};
}

inline InterpolationResult build_supported_energy_interpolation(
    const StructuredGrid& grid, const SparseMatrix& a,
    const std::vector<std::vector<int>>& supports,
    const InterpolationOptions& options,
    const SparseMatrix* initial_transpose = nullptr,
    const std::vector<unsigned char>* refine_column = nullptr) {
    if (supports.size() !=
        static_cast<std::size_t>(grid.coarse_size())) {
        throw std::invalid_argument(
            "one support is required for every coarse basis");
    }
    if (refine_column != nullptr &&
        (initial_transpose == nullptr ||
         refine_column->size() !=
             static_cast<std::size_t>(grid.coarse_size()))) {
        throw std::invalid_argument(
            "selective refinement requires one mask value per coarse basis");
    }
    const auto total_begin = Clock::now();
    unsigned int requested = options.thread_count > 0
        ? static_cast<unsigned int>(options.thread_count)
        : std::thread::hardware_concurrency();
    if (requested == 0U) requested = 1U;
    const int thread_count = std::min(
        grid.coarse_size(), static_cast<int>(requested));

    std::vector<LocalBasisResult> bases(
        static_cast<std::size_t>(grid.coarse_size()));
    std::atomic<int> next_coarse{0};
    std::exception_ptr worker_error;
    std::mutex error_mutex;

    const auto parallel_begin = Clock::now();
    auto worker = [&]() {
        std::vector<int> local_index(
            static_cast<std::size_t>(grid.fine_size()), -1);
        try {
            while (true) {
                const int coarse = next_coarse.fetch_add(
                    1, std::memory_order_relaxed);
                if (coarse >= grid.coarse_size()) break;
                const auto& support =
                    supports[static_cast<std::size_t>(coarse)];
                for (int node : support) {
                    if (node < 0 || node >= grid.fine_size() ||
                        grid.is_coarse_node(node)) {
                        throw std::invalid_argument(
                            "energy support contains an invalid F node");
                    }
                }
                if (refine_column != nullptr &&
                    (*refine_column)[static_cast<std::size_t>(coarse)] == 0U) {
                    LocalBasisResult copied;
                    copied.solve.converged = true;
                    for (int position =
                             initial_transpose->row_ptr()[
                                 static_cast<std::size_t>(coarse)];
                         position <
                             initial_transpose->row_ptr()[
                                 static_cast<std::size_t>(coarse) + 1U];
                         ++position) {
                        copied.triplets.push_back({
                            initial_transpose->col_idx()[
                                static_cast<std::size_t>(position)],
                            coarse,
                            initial_transpose->values()[
                                static_cast<std::size_t>(position)]
                        });
                    }
                    bases[static_cast<std::size_t>(coarse)] =
                        std::move(copied);
                } else {
                    bases[static_cast<std::size_t>(coarse)] =
                        build_basis_on_nodes(
                            grid, a, coarse, options, local_index, support,
                            initial_transpose);
                }
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
    const auto parallel_end = Clock::now();
    if (worker_error) std::rethrow_exception(worker_error);

    InterpolationReport report;
    report.threads_used = thread_count;
    report.timing.parallel_wall_ms =
        milliseconds(parallel_begin, parallel_end);
    std::size_t total_entries = 0;
    for (const LocalBasisResult& basis : bases) {
        ++report.local_solves.systems;
        report.local_solves.total_iterations += basis.solve.iterations;
        report.local_solves.max_iterations = std::max(
            report.local_solves.max_iterations, basis.solve.iterations);
        report.local_solves.sum_relative_residual +=
            basis.solve.relative_residual;
        report.local_solves.max_relative_residual = std::max(
            report.local_solves.max_relative_residual,
            basis.solve.relative_residual);
        if (!basis.solve.converged) ++report.local_solves.failed_systems;
        report.timing.patch_assembly_work_ms += basis.patch_assembly_ms;
        report.timing.local_solve_work_ms += basis.local_solve_ms;
        report.timing.basis_scatter_work_ms += basis.scatter_ms;
        total_entries += basis.triplets.size();
    }
    if (options.require_convergence &&
        report.local_solves.failed_systems != 0) {
        throw std::runtime_error(
            "one or more supported energy solves did not converge");
    }

    const auto finalize_begin = Clock::now();
    std::vector<int> row_ptr(
        static_cast<std::size_t>(grid.fine_size()) + 1U, 0);
    for (const LocalBasisResult& basis : bases) {
        for (const Triplet& entry : basis.triplets) {
            ++row_ptr[static_cast<std::size_t>(entry.row) + 1U];
        }
    }
    for (int row = 0; row < grid.fine_size(); ++row) {
        row_ptr[static_cast<std::size_t>(row) + 1U] +=
            row_ptr[static_cast<std::size_t>(row)];
    }
    std::vector<int> next = row_ptr;
    std::vector<int> col_idx(total_entries);
    std::vector<double> values(total_entries);
    for (const LocalBasisResult& basis : bases) {
        for (const Triplet& entry : basis.triplets) {
            const int target = next[static_cast<std::size_t>(entry.row)]++;
            col_idx[static_cast<std::size_t>(target)] = entry.column;
            values[static_cast<std::size_t>(target)] = entry.value;
        }
    }
    SparseMatrix prolongation(
        grid.fine_size(), grid.coarse_size(), std::move(row_ptr),
        std::move(col_idx), std::move(values));
    report.timing.finalize_ms = milliseconds(finalize_begin, Clock::now());
    report.timing.total_ms = milliseconds(total_begin, Clock::now());
    return {std::move(prolongation), report};
}

inline InterpolationResult build_geometric_interpolation(
    const StructuredGrid& grid, const InterpolationOptions& options) {
    const auto begin = Clock::now();
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
                const int coarse =
                    grid.coarse_id(qx[ax] - 1, qy[ay] - 1);
                const double value = wx[ax] * wy[ay];
                if (std::abs(value) > options.drop_tolerance) {
                    col_idx.push_back(coarse);
                    values.push_back(value);
                }
            }
        }
        row_ptr[static_cast<std::size_t>(fine) + 1U] =
            static_cast<int>(values.size());
    }

    InterpolationReport report;
    const auto finalize_begin = Clock::now();
    SparseMatrix prolongation(
        grid.fine_size(), grid.coarse_size(), std::move(row_ptr),
        std::move(col_idx), std::move(values));
    report.timing.finalize_ms = milliseconds(finalize_begin, Clock::now());
    report.timing.total_ms = milliseconds(begin, Clock::now());
    report.threads_used = 1;
    return {std::move(prolongation), report};
}

}

inline InterpolationResult build_interpolation(const StructuredGrid& grid,
                                        const SparseMatrix& a,
                                        const InterpolationOptions& options) {
    if (options.strategy == InterpolationStrategy::GeometricBilinear) {
        return energy_interpolation_detail::build_geometric_interpolation(
            grid, options);
    }
    if (options.strategy == InterpolationStrategy::GlobalEnergyMinimum) {
        InterpolationOptions global_options = options;
        global_options.patch_layers = 0;
        return energy_interpolation_detail::build_global_energy_interpolation(
            grid, a, global_options);
    }
    return energy_interpolation_detail::build_local_energy_interpolation(
        grid, a, options);
}

inline InterpolationResult build_energy_interpolation_on_supports(
    const StructuredGrid& grid, const SparseMatrix& a,
    const std::vector<std::vector<int>>& supports,
    const InterpolationOptions& options) {
    return energy_interpolation_detail::build_supported_energy_interpolation(
        grid, a, supports, options);
}

inline InterpolationResult refine_energy_interpolation_on_supports(
    const StructuredGrid& grid, const SparseMatrix& a,
    const std::vector<std::vector<int>>& supports,
    const SparseMatrix& initial,
    const InterpolationOptions& options) {
    if (initial.rows() != grid.fine_size() ||
        initial.cols() != grid.coarse_size()) {
        throw std::invalid_argument(
            "refine_energy_interpolation_on_supports: "
            "incompatible initial interpolation");
    }
    const SparseMatrix initial_transpose =
        initial.transpose(options.thread_count);
    return energy_interpolation_detail::build_supported_energy_interpolation(
        grid, a, supports, options, &initial_transpose);
}

inline InterpolationResult refine_selected_energy_interpolation_on_supports(
    const StructuredGrid& grid, const SparseMatrix& a,
    const std::vector<std::vector<int>>& supports,
    const std::vector<unsigned char>& refine_column,
    const SparseMatrix& initial,
    const InterpolationOptions& options) {
    if (initial.rows() != grid.fine_size() ||
        initial.cols() != grid.coarse_size()) {
        throw std::invalid_argument(
            "refine_selected_energy_interpolation_on_supports: "
            "incompatible initial interpolation");
    }
    const SparseMatrix initial_transpose =
        initial.transpose(options.thread_count);
    return energy_interpolation_detail::build_supported_energy_interpolation(
        grid, a, supports, options, &initial_transpose, &refine_column);
}

inline InterpolationResult refine_global_energy_interpolation(
    const StructuredGrid& grid, const SparseMatrix& a,
    const SparseMatrix& initial, const InterpolationOptions& options) {
    if (initial.rows() != grid.fine_size() ||
        initial.cols() != grid.coarse_size()) {
        throw std::invalid_argument(
            "refine_global_energy_interpolation: incompatible initial interpolation");
    }
    const SparseMatrix initial_transpose =
        initial.transpose(options.thread_count);
    InterpolationOptions global_options = options;
    global_options.strategy = InterpolationStrategy::GlobalEnergyMinimum;
    global_options.patch_layers = 0;
    return energy_interpolation_detail::build_global_energy_interpolation(
        grid, a, global_options, &initial_transpose);
}

}
