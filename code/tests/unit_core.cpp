#include "multigrid/energy_interpolation.hpp"
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
    require(tgi::version == "5.5.0", "wrong package version");
    require(tgi::version_major == 5 && tgi::version_minor == 5 &&
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

    tgi::GlobalEnergyOptions fixed_step_options = global_options;
    fixed_step_options.tolerance = 0.0;
    fixed_step_options.maximum_iterations = 2;
    fixed_step_options.require_convergence = false;
    const auto pcg2 = tgi::refine_global_energy_interpolation(
        grid, matrix, geometric.prolongation, fixed_step_options);
    require(pcg2.report.column_solves.maximum_iterations == 2,
            "fixed-step PCG did not honor the iteration budget");
    require(pcg2.report.column_solves.total_iterations ==
                2 * pcg2.report.column_solves.systems,
            "fixed-step PCG stopped before its explicit budget");

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
    return 0;
}
