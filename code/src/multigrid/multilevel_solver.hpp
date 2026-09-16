#pragma once
#include "multigrid/two_grid_solver.hpp"
#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>
namespace tgi {
class MultilevelVCycle {
public:
    struct Workspace {
        std::vector<Vector> residual;
        std::vector<Vector> coarse_rhs;
        std::vector<Vector> correction;
        Vector direct_work;
    };
    MultilevelVCycle(std::vector<SparseMatrix> level_matrices,
                     std::vector<SparseMatrix> prolongations,
                     int smoothing_steps = 1, int setup_threads = 1)
        : matrices_(std::move(level_matrices)),
          prolongations_(std::move(prolongations)),
          smoothing_steps_(smoothing_steps) {
        if (matrices_.size() != prolongations_.size() + 1U ||
            matrices_.size() < 2U) {
            throw std::invalid_argument("invalid multilevel hierarchy");
        }
        restrictions_.reserve(prolongations_.size());
        inverse_diagonal_.resize(prolongations_.size());
        diagonal_position_.resize(prolongations_.size());
        for (std::size_t level = 0; level < prolongations_.size(); ++level) {
            const SparseMatrix& matrix = matrices_[level];
            restrictions_.push_back(
                prolongations_[level].transpose(setup_threads));
            Vector& inverse = inverse_diagonal_[level];
            std::vector<int>& positions = diagonal_position_[level];
            inverse.resize(static_cast<std::size_t>(matrix.rows()));
            positions.resize(static_cast<std::size_t>(matrix.rows()));
            for (int row = 0; row < matrix.rows(); ++row) {
                const std::size_t index = static_cast<std::size_t>(row);
                int position = matrix.row_ptr()[index];
                const int end = matrix.row_ptr()[index + 1U];
                while (position < end &&
                       matrix.col_idx()[static_cast<std::size_t>(position)] <
                           row) {
                    ++position;
                }
                if (position == end ||
                    matrix.col_idx()[static_cast<std::size_t>(position)] !=
                        row) {
                    throw std::runtime_error("missing diagonal entry");
                }
                inverse[index] =
                    1.0 / matrix.values()[static_cast<std::size_t>(position)];
                positions[index] = position;
            }
        }
        const SparseMatrix& coarsest = matrices_.back();
        coarse_solver_.factorize(
            coarsest, two_grid_solver_detail::coarse_ordering(coarsest));
    }
    double iterate(const Vector& rhs, Vector& solution, Vector& residual,
                   Workspace& workspace) const {
        prepare_workspace(workspace);
        apply_level(0, rhs, solution, workspace);
        return matrices_.front().residual_squared(solution, rhs, residual);
    }
    const SparseMatrix& prolongation(int level) const {
        return prolongations_.at(static_cast<std::size_t>(level));
    }
    double operator_complexity() const {
        double entries = 0.0;
        for (const SparseMatrix& matrix : matrices_) {
            entries += static_cast<double>(matrix.nnz());
        }
        return entries / static_cast<double>(matrices_.front().nnz());
    }
    double interpolation_complexity() const {
        double entries = 0.0;
        for (const SparseMatrix& interpolation : prolongations_) {
            entries += static_cast<double>(interpolation.nnz());
        }
        return entries / static_cast<double>(matrices_.front().nnz());
    }
private:
    void prepare_workspace(Workspace& workspace) const {
        const std::size_t transitions = prolongations_.size();
        workspace.residual.resize(transitions);
        workspace.coarse_rhs.resize(transitions);
        workspace.correction.resize(transitions);
        for (std::size_t level = 0; level < transitions; ++level) {
            workspace.residual[level].resize(
                static_cast<std::size_t>(matrices_[level].rows()));
            workspace.coarse_rhs[level].resize(
                static_cast<std::size_t>(matrices_[level + 1U].rows()));
            workspace.correction[level].resize(
                static_cast<std::size_t>(matrices_[level + 1U].rows()));
        }
    }
    void smooth(int level, const Vector& rhs, bool forward,
                Vector& solution) const {
        const std::size_t index = static_cast<std::size_t>(level);
        two_grid_solver_detail::gauss_seidel_sweep(
            matrices_[index], inverse_diagonal_[index],
            diagonal_position_[index], rhs, forward, solution);
    }
    void apply_level(int level, const Vector& rhs, Vector& solution,
                     Workspace& workspace) const {
        const std::size_t index = static_cast<std::size_t>(level);
        for (int step = 0; step < smoothing_steps_; ++step) {
            smooth(level, rhs, true, solution);
        }
        matrices_[index].residual_squared(
            solution, rhs, workspace.residual[index]);
        restrictions_[index].multiply(
            workspace.residual[index], workspace.coarse_rhs[index]);
        Vector& correction = workspace.correction[index];
        std::fill(correction.begin(), correction.end(), 0.0);
        if (index + 1U == matrices_.size() - 1U) {
            coarse_solver_.solve(
                workspace.coarse_rhs[index], correction,
                workspace.direct_work);
        } else {
            apply_level(
                level + 1, workspace.coarse_rhs[index], correction,
                workspace);
        }
        prolongations_[index].multiply_add(1.0, correction, solution);
        for (int step = 0; step < smoothing_steps_; ++step) {
            smooth(level, rhs, false, solution);
        }
    }
    std::vector<SparseMatrix> matrices_;
    std::vector<SparseMatrix> prolongations_;
    std::vector<SparseMatrix> restrictions_;
    std::vector<Vector> inverse_diagonal_;
    std::vector<std::vector<int>> diagonal_position_;
    int smoothing_steps_ = 1;
    SparseCholesky coarse_solver_;
};
inline StationaryIterationResult solve_multilevel(
    const Vector& rhs, const MultilevelVCycle& cycle,
    double relative_tolerance = 1e-8, int max_cycles = 40000) {
    return solve_stationary_cycles(
        rhs, cycle, relative_tolerance, max_cycles);
}
}
