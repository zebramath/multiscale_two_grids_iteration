#include "experiment/test_problem.hpp"
#include "multigrid/adaptive_support.hpp"
#include "multigrid/algebraic_interpolation.hpp"
#include "multigrid/energy_interpolation.hpp"
#include "multigrid/reference_pruning.hpp"
#include "multigrid/residual_budget_support.hpp"
#include "multigrid/two_grid_solver.hpp"
#include "pde/diffusion_problem.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

double row_sum(const tgi::SparseMatrix& matrix, int row) {
    double sum = 0.0;
    for (int position = matrix.row_ptr()[static_cast<std::size_t>(row)];
         position < matrix.row_ptr()[static_cast<std::size_t>(row) + 1U];
         ++position) {
        sum += matrix.values()[static_cast<std::size_t>(position)];
    }
    return sum;
}

} // namespace

int main() {
    const tgi::StructuredGrid grid(15, 4);
    tgi::CoefficientOptions coefficient_options;
    coefficient_options.distribution =
        tgi::CoefficientDistribution::ChannelizedBinary;
    coefficient_options.contrast = 1.0e4;
    coefficient_options.channel_background_block_size = 4;
    const auto coefficient = tgi::make_coefficient(grid, coefficient_options);
    const tgi::SparseMatrix a = tgi::assemble_diffusion(
        grid, coefficient.values);
    require(a.rows() == grid.fine_size(), "wrong matrix dimension");
    for (double diagonal : a.diagonal()) {
        require(diagonal > 0.0, "diffusion diagonal is not positive");
    }

    tgi::InterpolationOptions geometric_options;
    geometric_options.strategy =
        tgi::InterpolationStrategy::GeometricBilinear;
    const auto geometric = tgi::build_interpolation(
        grid, a, geometric_options);
    tgi::InterpolationOptions local_options;
    local_options.strategy =
        tgi::InterpolationStrategy::LocalEnergyMinimum;
    local_options.patch_layers = 1;
    local_options.local_tolerance = 1.0e-8;
    local_options.local_max_iterations = 5000;
    local_options.thread_count = 2;
    const auto local = tgi::build_interpolation(grid, a, local_options);
    require(local.prolongation.cols() == grid.coarse_size(),
            "wrong local interpolation dimension");

    tgi::JacobiInterpolationOptions jacobi_options;
    jacobi_options.steps = 1;
    jacobi_options.maximum_entries_per_row = 0;
    jacobi_options.relative_drop_tolerance = 0.0;
    jacobi_options.thread_count = 2;
    const auto jacobi1 = tgi::build_jacobi_interpolation(
        grid, a, geometric.prolongation, jacobi_options);
    jacobi_options.steps = 4;
    const auto jacobi4 = tgi::build_jacobi_interpolation(
        grid, a, geometric.prolongation, jacobi_options);
    require(jacobi4.report.final_f_residual <=
                jacobi1.report.final_f_residual * (1.0 + 1.0e-10),
            "Jacobi F residual did not decrease");
    for (int fine = 0; fine < grid.fine_size(); ++fine) {
        if (!grid.is_coarse_node(fine)) continue;
        require(std::abs(row_sum(jacobi4.prolongation, fine) - 1.0) <
                    1.0e-14,
                "Jacobi interpolation lost C-point injection");
    }

    tgi::StrengthDistanceOptions distance_options;
    distance_options.coarse_candidates_per_row = 3;
    distance_options.thread_count = 2;
    const auto distance = tgi::build_strength_distance_interpolation(
        grid, a, distance_options);
    for (int fine = 0; fine < grid.fine_size(); ++fine) {
        if (!grid.is_coarse_node(fine)) continue;
        require(std::abs(row_sum(distance.prolongation, fine) - 1.0) <
                    1.0e-12,
                "strength-distance interpolation lost C-point injection");
    }

    tgi::InterpolationOptions global_options = local_options;
    global_options.strategy =
        tgi::InterpolationStrategy::GlobalEnergyMinimum;
    global_options.patch_layers = 0;
    global_options.local_tolerance = 1.0e-10;
    const auto global = tgi::build_interpolation(grid, a, global_options);
    const auto matched = tgi::build_budget_matched_reference(
        grid, a, global.prolongation, local.prolongation,
        local_options);
    require(
        tgi::interpolation_f_entries_per_column(
            grid, matched.prolongation, 2) ==
        tgi::interpolation_f_entries_per_column(
            grid, local.prolongation, 2),
        "budget-matched reference changed per-column support counts");

    tgi::ResidualStrongSupportOptions strong_options;
    strong_options.base_patch_layers = 1;
    strong_options.maximum_extra_nodes_per_column = 8;
    strong_options.thread_count = 2;
    const auto strong = tgi::build_residual_strong_supports(
        grid, a, local.prolongation, strong_options);
    require(strong.supports.size() ==
                static_cast<std::size_t>(grid.coarse_size()),
            "strong support count is wrong");

    tgi::ResidualBudgetSupportOptions budget_options;
    budget_options.base_patch_layers = 1;
    budget_options.maximum_rounds = 2;
    budget_options.maximum_extra_nodes_per_column = 8;
    budget_options.maximum_nodes_per_round = 4;
    budget_options.strength_scaling =
        tgi::StrengthScaling::SymmetricDiagonal;
    budget_options.thread_count = 2;
    const auto budget = tgi::build_residual_budget_interpolation(
        grid, a, local.prolongation, local_options, budget_options);
    require(budget.prolongation.rows() == grid.fine_size(),
            "residual-budget interpolation dimension is wrong");

    const tgi::Vector exact =
        experiment_support::manufactured_solution(grid);
    const tgi::Vector rhs = a.multiply(exact);
    const tgi::TwoGridCycle cycle(a, global.prolongation, 1, 2);
    const auto solved = tgi::solve_two_grid(a, rhs, cycle, 1.0e-6, 1000);
    require(solved.converged, "two-grid solve did not converge");
    return 0;
}
