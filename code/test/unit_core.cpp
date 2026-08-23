#include "experiment/test_problem.hpp"
#include "multigrid/energy_interpolation.hpp"
#include "multigrid/residual_budget_support.hpp"
#include "multigrid/two_grid_solver.hpp"
#include "pde/diffusion_problem.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

double half_trace(const tgi::SparseMatrix& matrix) {
    double trace = 0.0;
    for (int row = 0; row < matrix.rows(); ++row) {
        for (int position = matrix.row_ptr()[static_cast<std::size_t>(row)];
             position < matrix.row_ptr()[static_cast<std::size_t>(row) + 1U];
             ++position) {
            if (matrix.col_idx()[static_cast<std::size_t>(position)] == row) {
                trace += matrix.values()[static_cast<std::size_t>(position)];
                break;
            }
        }
    }
    return 0.5 * trace;
}

tgi::InterpolationOptions local_options(int layers) {
    tgi::InterpolationOptions options;
    options.strategy = tgi::InterpolationStrategy::LocalEnergyMinimum;
    options.patch_layers = layers;
    options.local_tolerance = 1.0e-8;
    options.local_max_iterations = 4000;
    options.thread_count = 2;
    return options;
}

} // namespace

int main() {
    const tgi::StructuredGrid grid(15, 4);
    tgi::CoefficientOptions coefficient_options;
    coefficient_options.distribution =
        tgi::CoefficientDistribution::ChannelizedBinary;
    coefficient_options.contrast = 1.0e3;
    coefficient_options.seed = 3;
    coefficient_options.channel_background_block_size = 4;
    coefficient_options.channel_width_fine_cells = 1;
    const auto coefficient = tgi::make_coefficient(
        grid, coefficient_options);
    const auto a = tgi::assemble_diffusion(grid, coefficient.values);

    const auto local1 = tgi::build_interpolation(
        grid, a, local_options(1));
    const auto local2 = tgi::build_interpolation(
        grid, a, local_options(2));
    require(local1.prolongation.rows() == grid.fine_size(),
            "local interpolation has wrong row count");
    require(local1.prolongation.cols() == grid.coarse_size(),
            "local interpolation has wrong column count");

    const tgi::TwoGridCycle cycle1(a, local1.prolongation, 1, 2);
    const tgi::TwoGridCycle cycle2(a, local2.prolongation, 1, 2);
    require(half_trace(cycle2.coarse_matrix()) <=
                half_trace(cycle1.coarse_matrix()) * (1.0 + 1.0e-9),
            "energy did not decrease when support grew");

    tgi::ResidualBudgetSupportOptions support;
    support.base_patch_layers = 1;
    support.maximum_rounds = 3;
    support.maximum_extra_nodes_per_column = 12;
    support.maximum_nodes_per_round = 4;
    support.target_residual_ratio = 0.5;
    support.strength_scaling = tgi::StrengthScaling::RowMaximum;
    support.strong_edge_fraction = 0.25;
    support.thread_count = 2;
    const auto adaptive = tgi::build_residual_budget_interpolation(
        grid, a, local1.prolongation, local_options(1), support);
    require(adaptive.prolongation.rows() == grid.fine_size(),
            "adaptive interpolation has wrong row count");
    require(adaptive.report.total_extra_nodes <=
                grid.coarse_size() * support.maximum_extra_nodes_per_column,
            "adaptive support exceeded its budget");
    require(adaptive.report.final_mean_scaled_residual <=
                adaptive.report.initial_mean_scaled_residual *
                    (1.0 + 1.0e-8),
            "adaptive refinement increased mean scaled residual");

    const auto rhs = a.multiply(
        experiment_support::manufactured_solution(grid));
    const tgi::TwoGridCycle adaptive_cycle(
        a, adaptive.prolongation, 1, 2);
    const auto solved = tgi::solve_two_grid(
        a, rhs, adaptive_cycle, 1.0e-6, 4000);
    require(solved.converged, "small adaptive two-grid solve did not converge");

    const double rho = adaptive_cycle.estimate_convergence_factor(20, 9);
    require(std::isfinite(rho) && rho >= 0.0 && rho <= 1.0,
            "spectral proxy is outside [0,1]");
    return 0;
}
