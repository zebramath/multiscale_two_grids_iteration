#include "multigrid/energy_interpolation.hpp"
#include "multigrid/multilevel_solver.hpp"
#include "multigrid/two_grid_solver.hpp"
#include "pde/diffusion_problem.hpp"
#include "version.hpp"

#include <algorithm>
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

}

int main() {
    require(tgi::version == "6.4.0", "wrong package version");
    require(tgi::version_major == 6 && tgi::version_minor == 4 &&
                tgi::version_patch == 0,
            "inconsistent numeric package version");

    const tgi::StructuredGrid grid(15, 4);
    tgi::CoefficientOptions coefficient_options;
    coefficient_options.contrast = 1.0e4;
    coefficient_options.channel_background_block_size = 4;
    const auto coefficient = tgi::make_coefficient(grid, coefficient_options);
    const tgi::SparseMatrix matrix =
        tgi::assemble_diffusion(grid, coefficient.values);
    require(matrix.rows() == grid.fine_size(), "wrong matrix dimension");
    for (double diagonal : matrix.diagonal()) {
        require(diagonal > 0.0, "diffusion diagonal is not positive");
    }

    std::vector<tgi::Vector> complex_topologies;
    for (const auto distribution : {
             tgi::CoefficientDistribution::BranchingChannelsBinary,
             tgi::CoefficientDistribution::WindingRingBinary}) {
        coefficient_options.distribution = distribution;
        const auto field = tgi::make_coefficient(grid, coefficient_options);
        require(field.actual_contrast == coefficient_options.contrast,
                "channel topology lost the requested contrast");
        complex_topologies.push_back(field.values);
    }
    require(tgi::norm2(tgi::subtract(
                complex_topologies[0], complex_topologies[1])) > 0.0,
            "channel topologies are indistinguishable");

    const auto geometric = tgi::build_geometric_interpolation(grid);
    for (int fine = 0; fine < grid.fine_size(); ++fine) {
        if (!grid.is_coarse_node(fine)) continue;
        require(std::abs(row_sum(geometric.prolongation, fine) - 1.0) <
                    1.0e-12,
                "geometric interpolation lost C-point injection");
    }

    tgi::GlobalEnergyOptions global_options;
    global_options.thread_count = 2;
    const auto global = tgi::build_global_energy_interpolation(
        grid, matrix, global_options);
    require(global.report.column_solves.failed_systems == 0,
            "global energy interpolation did not converge");
    require(global.report.column_solves.maximum_relative_residual <=
                1.01 * global_options.tolerance,
            "global energy interpolation missed its residual tolerance");
    require(global.report.threads_used == 2,
            "global energy interpolation ignored its thread budget");

    const tgi::TwoGridCycle cycle(matrix, global.prolongation, 1, 2);
    const auto& setup = cycle.setup_report();
    double coarse_trace = 0.0;
    for (double value : cycle.coarse_matrix().diagonal()) {
        coarse_trace += value;
    }
    require(std::abs(coarse_trace - setup.interpolation_energy) <=
                1.0e-12 * coarse_trace,
            "interpolation energy does not match trace(P^T A P)");

    tgi::Vector coarse_rhs(
        static_cast<std::size_t>(cycle.coarse_matrix().rows()), 1.0);
    tgi::Vector coarse_solution;
    tgi::Vector coarse_work;
    cycle.solve_coarse_system(coarse_rhs, coarse_solution, coarse_work);
    const tgi::Vector coarse_residual = tgi::subtract(
        coarse_rhs, cycle.coarse_matrix().multiply(coarse_solution));
    require(tgi::norm2(coarse_residual) / tgi::norm2(coarse_rhs) < 1.0e-10,
            "coarse direct solve is not numerically accurate");

    const tgi::Vector rhs(
        static_cast<std::size_t>(grid.fine_size()), 1.0);
    const auto solved = tgi::solve_two_grid(rhs, cycle, 1.0e-6, 1000);
    require(solved.converged, "two-grid solve did not converge");
    require(solved.status == tgi::TwoGridIterationStatus::Converged,
            "converged two-grid solve reported the wrong status");
    const auto limited = tgi::solve_two_grid(rhs, cycle, 1.0e-30, 1);
    require(!limited.converged &&
                limited.status == tgi::TwoGridIterationStatus::SlowAtLimit,
            "cycle-limited contraction was not classified as slow");
    require(std::isfinite(limited.tail_factor) &&
                limited.tail_factor > 0.0,
            "cycle-limited solve reported an invalid tail factor");
    require(std::string(tgi::two_grid_status_name(
                static_cast<tgi::TwoGridIterationStatus>(99))) == "unknown",
            "invalid two-grid status was misclassified");

    const tgi::StructuredGrid middle_grid(3, 2);
    const auto middle_geometric =
        tgi::build_geometric_interpolation(middle_grid);
    const tgi::SparseMatrix geometric_coarse =
        tgi::galerkin_coarse_operator(
            matrix, geometric.prolongation, 2);
    const tgi::TwoGridCycle geometric_two_grid(
        matrix, geometric.prolongation, 1, 2);
    tgi::MultilevelVCycle two_level(
        {matrix, geometric_coarse}, {geometric.prolongation}, 1, 2);
    tgi::Vector two_grid_solution(rhs.size(), 0.0);
    tgi::Vector two_grid_residual;
    tgi::TwoGridCycle::Workspace two_grid_workspace;
    (void)geometric_two_grid.iterate(
        rhs, two_grid_solution, two_grid_residual, two_grid_workspace);
    tgi::Vector two_level_solution(rhs.size(), 0.0);
    tgi::Vector two_level_residual;
    tgi::MultilevelVCycle::Workspace two_level_workspace;
    (void)two_level.iterate(
        rhs, two_level_solution, two_level_residual, two_level_workspace);
    require(tgi::norm2(tgi::subtract(
                two_grid_solution, two_level_solution)) <=
                1.0e-11 * std::max(1.0, tgi::norm2(two_grid_solution)),
            "two-level V-cycle did not reduce to exact two-grid");

    const tgi::SparseMatrix coarsest = tgi::galerkin_coarse_operator(
        geometric_coarse, middle_geometric.prolongation, 2);
    tgi::MultilevelVCycle multilevel(
        {matrix, geometric_coarse, coarsest},
        {geometric.prolongation, middle_geometric.prolongation}, 1, 2);
    require(multilevel.levels() == 3 &&
                multilevel.operator_complexity() > 1.0,
            "multilevel hierarchy metadata is invalid");
    const auto multilevel_solved = tgi::solve_multilevel(
        rhs, multilevel, 1.0e-6, 2000);
    require(multilevel_solved.converged,
            "three-level geometric V-cycle did not converge");
    return 0;
}
