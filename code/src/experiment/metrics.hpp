#pragma once

#include "experiment/reporting.hpp"
#include "multigrid/energy_interpolation.hpp"
#include "multigrid/two_grid_solver.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace experiment_support {

struct Statistics {
    double mean = 0.0;
    double p50 = 0.0;
    double p95 = 0.0;
    double p99 = 0.0;
};

inline Statistics summarize(std::vector<double> values) {
    if (values.empty()) return {};
    Statistics result;
    result.mean = std::accumulate(values.begin(), values.end(), 0.0) /
                  static_cast<double>(values.size());
    std::sort(values.begin(), values.end());
    const auto quantile = [&](double probability) {
        const double index =
            probability * static_cast<double>(values.size() - 1U);
        const std::size_t lower = static_cast<std::size_t>(index);
        const std::size_t upper =
            std::min(lower + 1U, values.size() - 1U);
        const double fraction = index - static_cast<double>(lower);
        return values[lower] * (1.0 - fraction) +
               values[upper] * fraction;
    };
    result.p50 = quantile(0.50);
    result.p95 = quantile(0.95);
    result.p99 = quantile(0.99);
    return result;
}

struct BasisMetrics {
    std::vector<double> relative_residual;
    std::vector<double> relative_l2_error;
    std::vector<double> relative_energy_error;
};

struct GlobalProlongationMetrics {
    BasisMetrics basis;
    double aggregate_relative_l2_error = 0.0;
    double aggregate_relative_energy_error = 0.0;
};

inline void append(BasisMetrics& destination, const BasisMetrics& source) {
    destination.relative_residual.insert(
        destination.relative_residual.end(),
        source.relative_residual.begin(), source.relative_residual.end());
    destination.relative_l2_error.insert(
        destination.relative_l2_error.end(),
        source.relative_l2_error.begin(), source.relative_l2_error.end());
    destination.relative_energy_error.insert(
        destination.relative_energy_error.end(),
        source.relative_energy_error.begin(),
        source.relative_energy_error.end());
}

inline void load_column(const tgi::SparseMatrix& transpose, int column,
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

inline void clear_column(tgi::Vector& dense,
                         const std::vector<int>& active) {
    for (int row : active) dense[static_cast<std::size_t>(row)] = 0.0;
}

inline double quadratic_energy(const tgi::SparseMatrix& a,
                               const tgi::Vector& vector,
                               const std::vector<int>& active) {
    double energy = 0.0;
    for (int row : active) {
        const double row_value = vector[static_cast<std::size_t>(row)];
        if (row_value == 0.0) continue;
        double product = 0.0;
        for (int position =
                 a.row_ptr()[static_cast<std::size_t>(row)];
             position <
                 a.row_ptr()[static_cast<std::size_t>(row) + 1U];
             ++position) {
            product += a.values()[static_cast<std::size_t>(position)] *
                vector[static_cast<std::size_t>(
                    a.col_idx()[static_cast<std::size_t>(position)])];
        }
        energy += row_value * product;
    }
    return std::max(0.0, energy);
}

inline BasisMetrics compare_prolongations(
    const tgi::StructuredGrid& grid, const tgi::SparseMatrix& a,
    const tgi::SparseMatrix& reference,
    const tgi::SparseMatrix& candidate, int patch_layers) {
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
    std::vector<int> energy_active;

    BasisMetrics metrics;
    metrics.relative_residual.reserve(
        static_cast<std::size_t>(grid.coarse_size()));
    metrics.relative_l2_error.reserve(
        static_cast<std::size_t>(grid.coarse_size()));
    metrics.relative_energy_error.reserve(
        static_cast<std::size_t>(grid.coarse_size()));

    for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
        load_column(reference_transpose, coarse, reference_dense,
                    reference_active);
        load_column(candidate_transpose, coarse, candidate_dense,
                    candidate_active);
        const std::vector<int> local_nodes =
            grid.patch_f_nodes(coarse, patch_layers);
        const int coarse_fine = grid.coarse_fine_id(coarse);

        double error_l2_squared = 0.0;
        double reference_l2_squared = 0.0;
        energy_active.clear();
        energy_active.reserve(local_nodes.size() + 1U);
        for (int node : local_nodes) {
            const double exact =
                reference_dense[static_cast<std::size_t>(node)];
            const double error =
                candidate_dense[static_cast<std::size_t>(node)] - exact;
            difference[static_cast<std::size_t>(node)] = error;
            error_l2_squared += error * error;
            reference_l2_squared += exact * exact;
            energy_active.push_back(node);
        }
        energy_active.push_back(coarse_fine);
        const double error_energy =
            quadratic_energy(a, difference, local_nodes);
        const double reference_energy =
            quadratic_energy(a, reference_dense, energy_active);

        double residual_squared = 0.0;
        double rhs_squared = 0.0;
        for (int node : local_nodes) {
            double ap = 0.0;
            double rhs_entry = 0.0;
            for (int position =
                     a.row_ptr()[static_cast<std::size_t>(node)];
                 position <
                     a.row_ptr()[static_cast<std::size_t>(node) + 1U];
                 ++position) {
                const int col =
                    a.col_idx()[static_cast<std::size_t>(position)];
                const double value =
                    a.values()[static_cast<std::size_t>(position)];
                ap += value *
                    candidate_dense[static_cast<std::size_t>(col)];
                if (col == coarse_fine) rhs_entry = -value;
            }
            residual_squared += ap * ap;
            rhs_squared += rhs_entry * rhs_entry;
        }

        metrics.relative_l2_error.push_back(
            reference_l2_squared > 0.0
                ? std::sqrt(error_l2_squared / reference_l2_squared)
                : 0.0);
        metrics.relative_energy_error.push_back(
            reference_energy > 0.0
                ? std::sqrt(error_energy / reference_energy)
                : 0.0);
        metrics.relative_residual.push_back(
            rhs_squared > 0.0
                ? std::sqrt(residual_squared / rhs_squared)
                : 0.0);

        for (int node : local_nodes) {
            difference[static_cast<std::size_t>(node)] = 0.0;
        }
        clear_column(reference_dense, reference_active);
        clear_column(candidate_dense, candidate_active);
    }
    return metrics;
}

inline GlobalProlongationMetrics compare_prolongations_global(
    const tgi::StructuredGrid& grid, const tgi::SparseMatrix& a,
    const tgi::SparseMatrix& reference,
    const tgi::SparseMatrix& candidate) {
    if (reference.rows() != grid.fine_size() ||
        candidate.rows() != grid.fine_size() ||
        reference.cols() != grid.coarse_size() ||
        candidate.cols() != grid.coarse_size()) {
        throw std::invalid_argument(
            "compare_prolongations_global: incompatible dimensions");
    }
    const tgi::SparseMatrix reference_transpose = reference.transpose();
    const tgi::SparseMatrix candidate_transpose = candidate.transpose();
    const tgi::Vector diagonal = a.diagonal();
    tgi::Vector reference_dense(
        static_cast<std::size_t>(grid.fine_size()), 0.0);
    tgi::Vector candidate_dense(
        static_cast<std::size_t>(grid.fine_size()), 0.0);
    tgi::Vector difference(
        static_cast<std::size_t>(grid.fine_size()), 0.0);
    std::vector<int> reference_active;
    std::vector<int> candidate_active;

    GlobalProlongationMetrics metrics;
    metrics.basis.relative_residual.reserve(
        static_cast<std::size_t>(grid.coarse_size()));
    metrics.basis.relative_l2_error.reserve(
        static_cast<std::size_t>(grid.coarse_size()));
    metrics.basis.relative_energy_error.reserve(
        static_cast<std::size_t>(grid.coarse_size()));
    double aggregate_error_l2 = 0.0;
    double aggregate_reference_l2 = 0.0;
    double aggregate_error_energy = 0.0;
    double aggregate_reference_energy = 0.0;

    for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
        load_column(
            reference_transpose, coarse, reference_dense,
            reference_active);
        load_column(
            candidate_transpose, coarse, candidate_dense,
            candidate_active);
        double error_l2 = 0.0;
        double reference_l2 = 0.0;
        double candidate_l2 = 0.0;
        for (int row = 0; row < grid.fine_size(); ++row) {
            const double exact =
                reference_dense[static_cast<std::size_t>(row)];
            const double approximate =
                candidate_dense[static_cast<std::size_t>(row)];
            const double error = approximate - exact;
            difference[static_cast<std::size_t>(row)] = error;
            error_l2 += error * error;
            reference_l2 += exact * exact;
            candidate_l2 += approximate * approximate;
        }
        const tgi::Vector a_difference = a.multiply(difference);
        const tgi::Vector a_reference = a.multiply(reference_dense);
        const tgi::Vector a_candidate = a.multiply(candidate_dense);
        const double error_energy = std::max(
            0.0, tgi::dot(difference, a_difference));
        const double reference_energy = std::max(
            0.0, tgi::dot(reference_dense, a_reference));
        double scaled_f_residual = 0.0;
        for (int row = 0; row < grid.fine_size(); ++row) {
            if (grid.is_coarse_node(row)) continue;
            const double value =
                a_candidate[static_cast<std::size_t>(row)] /
                diagonal[static_cast<std::size_t>(row)];
            scaled_f_residual += value * value;
        }

        metrics.basis.relative_l2_error.push_back(
            reference_l2 > 0.0
                ? std::sqrt(error_l2 / reference_l2)
                : 0.0);
        metrics.basis.relative_energy_error.push_back(
            reference_energy > 0.0
                ? std::sqrt(error_energy / reference_energy)
                : 0.0);
        metrics.basis.relative_residual.push_back(
            candidate_l2 > 0.0
                ? std::sqrt(scaled_f_residual / candidate_l2)
                : 0.0);
        aggregate_error_l2 += error_l2;
        aggregate_reference_l2 += reference_l2;
        aggregate_error_energy += error_energy;
        aggregate_reference_energy += reference_energy;

        std::fill(difference.begin(), difference.end(), 0.0);
        clear_column(reference_dense, reference_active);
        clear_column(candidate_dense, candidate_active);
    }
    metrics.aggregate_relative_l2_error =
        aggregate_reference_l2 > 0.0
            ? std::sqrt(aggregate_error_l2 / aggregate_reference_l2)
            : 0.0;
    metrics.aggregate_relative_energy_error =
        aggregate_reference_energy > 0.0
            ? std::sqrt(
                aggregate_error_energy / aggregate_reference_energy)
            : 0.0;
    return metrics;
}

struct CycleMetrics {
    int cycles = 0;
    double relative_residual = 0.0;
    bool converged = false;
    double coarse_setup_ms = 0.0;
    double solve_ms = 0.0;
    double total_ms = 0.0;
    std::size_t coarse_nnz = 0;
    std::size_t factor_nnz = 0;
};

inline CycleMetrics evaluate_two_grid(
    const tgi::SparseMatrix& a, const tgi::Vector& rhs,
    const tgi::SparseMatrix& prolongation, int setup_threads,
    double outer_tolerance, int max_cycles) {
    const auto begin = Clock::now();
    const tgi::TwoGridCycle cycle(a, prolongation, 1, setup_threads);
    const auto solve_begin = Clock::now();
    const tgi::TwoGridIterationResult result =
        tgi::solve_two_grid(
            a, rhs, cycle, outer_tolerance, max_cycles);
    const auto end = Clock::now();
    CycleMetrics metrics;
    metrics.cycles = result.cycles;
    metrics.relative_residual = result.relative_residual;
    metrics.converged = result.converged;
    metrics.coarse_setup_ms = cycle.setup_report().total_ms;
    metrics.solve_ms = milliseconds(solve_begin, end);
    metrics.total_ms = milliseconds(begin, end);
    metrics.coarse_nnz = cycle.setup_report().coarse_nnz;
    metrics.factor_nnz = cycle.setup_report().factor_nnz;
    return metrics;
}

inline tgi::SparseMatrix drop_prolongation(
    const tgi::SparseMatrix& matrix, double threshold) {
    std::vector<int> row_ptr(
        static_cast<std::size_t>(matrix.rows()) + 1U, 0);
    std::vector<int> col_idx;
    std::vector<double> values;
    col_idx.reserve(matrix.nnz());
    values.reserve(matrix.nnz());
    for (int row = 0; row < matrix.rows(); ++row) {
        for (int position =
                 matrix.row_ptr()[static_cast<std::size_t>(row)];
             position <
                 matrix.row_ptr()[static_cast<std::size_t>(row) + 1U];
             ++position) {
            const double value =
                matrix.values()[static_cast<std::size_t>(position)];
            if (std::abs(value) > threshold) {
                col_idx.push_back(
                    matrix.col_idx()[static_cast<std::size_t>(position)]);
                values.push_back(value);
            }
        }
        row_ptr[static_cast<std::size_t>(row) + 1U] =
            static_cast<int>(values.size());
    }
    return tgi::SparseMatrix(
        matrix.rows(), matrix.cols(), std::move(row_ptr),
        std::move(col_idx), std::move(values));
}

}
