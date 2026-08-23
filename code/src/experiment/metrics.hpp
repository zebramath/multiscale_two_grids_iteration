#pragma once

#include "multigrid/two_grid_solver.hpp"
#include "pde/diffusion_problem.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace experiment_support {

struct GlobalProlongationMetrics {
    double aggregate_relative_l2_error = 0.0;
    double aggregate_relative_energy_error = 0.0;
};

namespace metrics_detail {

inline void load_column(
    const tgi::SparseMatrix& transpose, int column,
    tgi::Vector& dense, std::vector<int>& active) {
    active.clear();
    for (int position =
             transpose.row_ptr()[static_cast<std::size_t>(column)];
         position <
             transpose.row_ptr()[static_cast<std::size_t>(column) + 1U];
         ++position) {
        const int row =
            transpose.col_idx()[static_cast<std::size_t>(position)];
        dense[static_cast<std::size_t>(row)] =
            transpose.values()[static_cast<std::size_t>(position)];
        active.push_back(row);
    }
}

inline void clear_column(
    tgi::Vector& dense, const std::vector<int>& active) {
    for (int row : active) dense[static_cast<std::size_t>(row)] = 0.0;
}

} // namespace metrics_detail

inline GlobalProlongationMetrics compare_prolongations_global(
    const tgi::StructuredGrid& grid, const tgi::SparseMatrix& a,
    const tgi::SparseMatrix& reference,
    const tgi::SparseMatrix& candidate) {
    if (a.rows() != grid.fine_size() || a.cols() != grid.fine_size() ||
        reference.rows() != grid.fine_size() ||
        candidate.rows() != grid.fine_size() ||
        reference.cols() != grid.coarse_size() ||
        candidate.cols() != grid.coarse_size()) {
        throw std::invalid_argument(
            "compare_prolongations_global: incompatible dimensions");
    }
    const tgi::SparseMatrix reference_transpose = reference.transpose();
    const tgi::SparseMatrix candidate_transpose = candidate.transpose();
    tgi::Vector reference_dense(
        static_cast<std::size_t>(grid.fine_size()), 0.0);
    tgi::Vector candidate_dense(
        static_cast<std::size_t>(grid.fine_size()), 0.0);
    tgi::Vector difference(
        static_cast<std::size_t>(grid.fine_size()), 0.0);
    std::vector<int> reference_active;
    std::vector<int> candidate_active;
    double error_l2 = 0.0;
    double reference_l2 = 0.0;
    double error_energy = 0.0;
    double reference_energy = 0.0;

    for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
        metrics_detail::load_column(
            reference_transpose, coarse, reference_dense,
            reference_active);
        metrics_detail::load_column(
            candidate_transpose, coarse, candidate_dense,
            candidate_active);
        for (int row = 0; row < grid.fine_size(); ++row) {
            const double exact =
                reference_dense[static_cast<std::size_t>(row)];
            const double delta =
                candidate_dense[static_cast<std::size_t>(row)] - exact;
            difference[static_cast<std::size_t>(row)] = delta;
            error_l2 += delta * delta;
            reference_l2 += exact * exact;
        }
        const tgi::Vector a_difference = a.multiply(difference);
        const tgi::Vector a_reference = a.multiply(reference_dense);
        error_energy += std::max(
            0.0, tgi::dot(difference, a_difference));
        reference_energy += std::max(
            0.0, tgi::dot(reference_dense, a_reference));
        std::fill(difference.begin(), difference.end(), 0.0);
        metrics_detail::clear_column(reference_dense, reference_active);
        metrics_detail::clear_column(candidate_dense, candidate_active);
    }

    GlobalProlongationMetrics metrics;
    metrics.aggregate_relative_l2_error = reference_l2 > 0.0
        ? std::sqrt(error_l2 / reference_l2)
        : 0.0;
    metrics.aggregate_relative_energy_error = reference_energy > 0.0
        ? std::sqrt(error_energy / reference_energy)
        : 0.0;
    return metrics;
}

struct CycleMetrics {
    int cycles = 0;
    double relative_residual = 0.0;
    bool converged = false;
    double convergence_factor = 0.0;
    double coarse_setup_ms = 0.0;
    double solve_ms = 0.0;
    double total_ms = 0.0;
    std::size_t coarse_nnz = 0;
    std::size_t factor_nnz = 0;
};

inline CycleMetrics evaluate_two_grid(
    const tgi::SparseMatrix& a, const tgi::Vector& rhs,
    const tgi::SparseMatrix& prolongation, int setup_threads,
    double outer_tolerance, int max_cycles,
    int spectral_iterations = 0) {
    const auto begin = std::chrono::steady_clock::now();
    const tgi::TwoGridCycle cycle(a, prolongation, 1, setup_threads);
    const auto solve_begin = std::chrono::steady_clock::now();
    const tgi::TwoGridIterationResult result = tgi::solve_two_grid(
        a, rhs, cycle, outer_tolerance, max_cycles);
    const auto end = std::chrono::steady_clock::now();
    CycleMetrics metrics;
    metrics.cycles = result.cycles;
    metrics.relative_residual = result.relative_residual;
    metrics.converged = result.converged;
    if (spectral_iterations > 0) {
        metrics.convergence_factor =
            cycle.estimate_convergence_factor(spectral_iterations, 271828U);
    }
    metrics.coarse_setup_ms = cycle.setup_report().total_ms;
    metrics.solve_ms =
        std::chrono::duration<double, std::milli>(end - solve_begin).count();
    metrics.total_ms =
        std::chrono::duration<double, std::milli>(end - begin).count();
    metrics.coarse_nnz = cycle.setup_report().coarse_nnz;
    metrics.factor_nnz = cycle.setup_report().factor_nnz;
    return metrics;
}

} // namespace experiment_support
