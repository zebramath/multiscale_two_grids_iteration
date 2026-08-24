#pragma once

#include "multigrid/two_grid_solver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace tgi {

struct MultilevelSetupReport {
    int levels = 0;
    std::size_t operator_nnz = 0;
    std::size_t interpolation_nnz = 0;
    double operator_complexity = 0.0;
    double interpolation_complexity = 0.0;
    double setup_ms = 0.0;
};

struct MultilevelSolveResult {
    Vector solution;
    int iterations = 0;
    double relative_residual = 0.0;
    bool converged = false;
};

class MultilevelCycle {
public:
    MultilevelCycle(
        std::vector<SparseMatrix> matrices,
        std::vector<SparseMatrix> prolongations,
        int smoothing_steps = 1, int application_threads = 1);

    Vector apply(const Vector& residual) const;
    double iterate(
        const Vector& rhs, Vector& solution, Vector& residual) const;

    const SparseMatrix& fine_matrix() const { return levels_.front().matrix; }
    const MultilevelSetupReport& setup_report() const { return report_; }

private:
    struct Level {
        SparseMatrix matrix;
        SparseMatrix prolongation;
        SparseMatrix restriction;
        Vector inverse_diagonal;
        std::vector<int> diagonal_position;
    };

    void correction_cycle(
        std::size_t level, Vector residual, Vector& correction) const;
    void triangular_correction(
        const Level& level, const Vector& residual,
        bool forward, Vector& correction) const;

    std::vector<Level> levels_;
    SparseCholesky coarsest_solver_;
    int smoothing_steps_ = 1;
    int application_threads_ = 1;
    MultilevelSetupReport report_;
};

inline MultilevelCycle::MultilevelCycle(
    std::vector<SparseMatrix> matrices,
    std::vector<SparseMatrix> prolongations,
    int smoothing_steps, int application_threads)
    : smoothing_steps_(smoothing_steps),
      application_threads_(std::max(1, application_threads)) {
    const auto begin = std::chrono::steady_clock::now();
    if (matrices.empty() || prolongations.size() + 1U != matrices.size()) {
        throw std::invalid_argument(
            "multilevel hierarchy requires one prolongation per transition");
    }
    if (smoothing_steps_ <= 0) {
        throw std::invalid_argument(
            "multilevel hierarchy requires positive smoothing steps");
    }

    levels_.reserve(matrices.size());
    for (std::size_t index = 0; index < matrices.size(); ++index) {
        Level level;
        level.matrix = std::move(matrices[index]);
        if (level.matrix.rows() != level.matrix.cols() ||
            level.matrix.rows() <= 0) {
            throw std::invalid_argument(
                "multilevel hierarchy matrices must be nonempty and square");
        }
        report_.operator_nnz += level.matrix.nnz();
        if (index < prolongations.size()) {
            level.prolongation = std::move(prolongations[index]);
            if (level.prolongation.rows() != level.matrix.rows() ||
                level.prolongation.cols() != matrices[index + 1U].rows()) {
                throw std::invalid_argument(
                    "multilevel prolongation has incompatible dimensions");
            }
            report_.interpolation_nnz += level.prolongation.nnz();
            level.restriction = level.prolongation.transpose(
                application_threads_);
        }

        level.inverse_diagonal.resize(
            static_cast<std::size_t>(level.matrix.rows()));
        level.diagonal_position.resize(
            static_cast<std::size_t>(level.matrix.rows()));
        for (int row = 0; row < level.matrix.rows(); ++row) {
            const std::size_t offset = static_cast<std::size_t>(row);
            int position = level.matrix.row_ptr()[offset];
            const int end = level.matrix.row_ptr()[offset + 1U];
            while (position < end && level.matrix.col_idx()[
                       static_cast<std::size_t>(position)] < row) {
                ++position;
            }
            if (position == end || level.matrix.col_idx()[
                    static_cast<std::size_t>(position)] != row) {
                throw std::runtime_error(
                    "multilevel matrix is missing a diagonal entry");
            }
            const double diagonal = level.matrix.values()[
                static_cast<std::size_t>(position)];
            if (!(diagonal > 0.0)) {
                throw std::runtime_error(
                    "multilevel matrix has nonpositive diagonal");
            }
            level.inverse_diagonal[offset] = 1.0 / diagonal;
            level.diagonal_position[offset] = position;
        }
        levels_.push_back(std::move(level));
    }

    coarsest_solver_.factorize(levels_.back().matrix);
    const double fine_nnz = static_cast<double>(levels_.front().matrix.nnz());
    report_.levels = static_cast<int>(levels_.size());
    report_.operator_complexity =
        static_cast<double>(report_.operator_nnz) / fine_nnz;
    report_.interpolation_complexity =
        static_cast<double>(report_.interpolation_nnz) / fine_nnz;
    report_.setup_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - begin).count();
}

inline void MultilevelCycle::triangular_correction(
    const Level& level, const Vector& residual,
    bool forward, Vector& correction) const {
    correction.assign(residual.size(), 0.0);
    const int begin = forward ? 0 : level.matrix.rows() - 1;
    const int end = forward ? level.matrix.rows() : -1;
    const int increment = forward ? 1 : -1;
    for (int row = begin; row != end; row += increment) {
        double value = residual[static_cast<std::size_t>(row)];
        for (int position = level.matrix.row_ptr()[
                 static_cast<std::size_t>(row)];
             position < level.matrix.row_ptr()[
                 static_cast<std::size_t>(row) + 1U]; ++position) {
            const int column = level.matrix.col_idx()[
                static_cast<std::size_t>(position)];
            const bool triangular = forward ? column < row : column > row;
            if (triangular) {
                value -= level.matrix.values()[
                             static_cast<std::size_t>(position)] *
                    correction[static_cast<std::size_t>(column)];
            }
        }
        correction[static_cast<std::size_t>(row)] =
            value * level.inverse_diagonal[static_cast<std::size_t>(row)];
    }
}

inline void MultilevelCycle::correction_cycle(
    std::size_t level_index, Vector residual, Vector& correction) const {
    if (level_index + 1U == levels_.size()) {
        Vector work;
        coarsest_solver_.solve(residual, correction, work);
        return;
    }

    const Level& level = levels_[level_index];
    correction.assign(residual.size(), 0.0);
    Vector sweep;
    Vector product;
    for (int step = 0; step < smoothing_steps_; ++step) {
        triangular_correction(level, residual, true, sweep);
        axpy(1.0, sweep, correction);
        level.matrix.multiply(sweep, product, application_threads_);
        axpy(-1.0, product, residual);
    }

    Vector coarse_residual;
    level.restriction.multiply(
        residual, coarse_residual, application_threads_);
    Vector coarse_correction;
    correction_cycle(level_index + 1U, coarse_residual, coarse_correction);
    Vector fine_correction;
    level.prolongation.multiply(
        coarse_correction, fine_correction, application_threads_);
    axpy(1.0, fine_correction, correction);
    level.matrix.multiply(
        fine_correction, product, application_threads_);
    axpy(-1.0, product, residual);

    for (int step = 0; step < smoothing_steps_; ++step) {
        triangular_correction(level, residual, false, sweep);
        axpy(1.0, sweep, correction);
        level.matrix.multiply(sweep, product, application_threads_);
        axpy(-1.0, product, residual);
    }
}

inline Vector MultilevelCycle::apply(const Vector& residual) const {
    if (residual.size() !=
        static_cast<std::size_t>(levels_.front().matrix.rows())) {
        throw std::invalid_argument(
            "multilevel residual has the wrong size");
    }
    Vector correction;
    correction_cycle(0U, residual, correction);
    return correction;
}

inline double MultilevelCycle::iterate(
    const Vector& rhs, Vector& solution, Vector& residual) const {
    levels_.front().matrix.residual_squared(
        solution, rhs, residual, application_threads_);
    const Vector correction = apply(residual);
    axpy(1.0, correction, solution);
    return levels_.front().matrix.residual_squared(
        solution, rhs, residual, application_threads_);
}

inline MultilevelSolveResult solve_multilevel(
    const Vector& rhs, const MultilevelCycle& cycle,
    double relative_tolerance = 1.0e-8, int maximum_iterations = 40000) {
    if (!(relative_tolerance > 0.0) || maximum_iterations <= 0) {
        throw std::invalid_argument("invalid multilevel solve options");
    }
    MultilevelSolveResult result;
    result.solution.assign(rhs.size(), 0.0);
    Vector residual = rhs;
    const double initial = norm2(rhs);
    if (initial == 0.0) {
        result.converged = true;
        return result;
    }
    for (int iteration = 0; iteration < maximum_iterations; ++iteration) {
        const double squared = cycle.iterate(
            rhs, result.solution, residual);
        result.iterations = iteration + 1;
        result.relative_residual = std::sqrt(squared) / initial;
        if (result.relative_residual <= relative_tolerance) {
            result.converged = true;
            break;
        }
    }
    return result;
}

inline MultilevelSolveResult solve_preconditioned_cg(
    const SparseMatrix& matrix, const Vector& rhs,
    const MultilevelCycle& preconditioner,
    double relative_tolerance = 1.0e-8, int maximum_iterations = 40000) {
    if (matrix.rows() != matrix.cols() ||
        rhs.size() != static_cast<std::size_t>(matrix.rows()) ||
        !(relative_tolerance > 0.0) || maximum_iterations <= 0) {
        throw std::invalid_argument("invalid preconditioned CG inputs");
    }
    MultilevelSolveResult result;
    result.solution.assign(rhs.size(), 0.0);
    Vector residual = rhs;
    const double initial = norm2(rhs);
    if (initial == 0.0) {
        result.converged = true;
        return result;
    }
    Vector preconditioned = preconditioner.apply(residual);
    Vector direction = preconditioned;
    double rz = dot(residual, preconditioned);
    if (!(rz > 0.0) || !std::isfinite(rz)) {
        throw std::runtime_error(
            "multilevel preconditioner is not positive definite");
    }
    for (int iteration = 0; iteration < maximum_iterations; ++iteration) {
        const Vector product = matrix.multiply(direction);
        const double denominator = dot(direction, product);
        if (!(denominator > 0.0) || !std::isfinite(denominator)) {
            throw std::runtime_error(
                "preconditioned CG encountered nonpositive curvature");
        }
        const double alpha = rz / denominator;
        axpy(alpha, direction, result.solution);
        axpy(-alpha, product, residual);
        result.iterations = iteration + 1;
        result.relative_residual = norm2(residual) / initial;
        if (result.relative_residual <= relative_tolerance) {
            result.converged = true;
            break;
        }
        preconditioned = preconditioner.apply(residual);
        const double next_rz = dot(residual, preconditioned);
        if (!(next_rz > 0.0) || !std::isfinite(next_rz)) {
            throw std::runtime_error(
                "multilevel preconditioner lost positive definiteness");
        }
        const double beta = next_rz / rz;
        for (std::size_t index = 0; index < direction.size(); ++index) {
            direction[index] = preconditioned[index] +
                beta * direction[index];
        }
        rz = next_rz;
    }
    return result;
}

}
