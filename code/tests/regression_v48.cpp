#include "experiment/problem.hpp"
#include "multigrid/global_pcg.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

template <class Action>
void require_throws(Action&& action, const std::string& message) {
    try {
        action();
    } catch (const std::exception&) {
        return;
    }
    throw std::runtime_error(message);
}

double relative_action_difference(
    const tgi::SparseMatrix& lhs, const tgi::SparseMatrix& rhs) {
    tgi::Vector input(static_cast<std::size_t>(lhs.cols()));
    for (int index = 0; index < lhs.cols(); ++index) {
        input[static_cast<std::size_t>(index)] =
            std::sin(0.31 * static_cast<double>(index + 1));
    }
    const tgi::Vector lhs_value = lhs.multiply(input);
    const tgi::Vector rhs_value = rhs.multiply(input);
    return tgi::norm2(tgi::subtract(lhs_value, rhs_value)) /
        std::max(tgi::norm2(rhs_value), 1.0e-30);
}

}

int main() {
    const tgi::StructuredGrid grid(15, 4);
    tgi::CoefficientOptions coefficient_options;
    coefficient_options.distribution =
        tgi::CoefficientDistribution::MeanderingChannelBinary;
    coefficient_options.contrast = 1.0e4;
    coefficient_options.seed = 7;
    coefficient_options.channel_background_block_size = 4;
    const auto coefficient = tgi::make_coefficient(grid, coefficient_options);
    require(coefficient.actual_contrast >= 0.99e4,
            "new channel topology lost the requested contrast");
    const tgi::SparseMatrix a =
        tgi::assemble_diffusion(grid, coefficient.values);

    const auto geometric = tgi::build_geometric_interpolation(grid);

    tgi::GlobalEnergyPcgPath path(
        grid, a, geometric.prolongation, 2);
    path.advance_to(2);
    const tgi::SparseMatrix checkpoint_two = path.prolongation(0.0);
    path.advance_to(5);
    const tgi::SparseMatrix continued = path.prolongation(0.0);
    tgi::GlobalEnergyOptions fixed_options;
    fixed_options.tolerance = 0.0;
    fixed_options.maximum_iterations = 5;
    fixed_options.thread_count = 2;
    fixed_options.drop_tolerance = 0.0;
    fixed_options.require_convergence = false;
    const auto restarted = tgi::refine_global_energy_interpolation(
        grid, a, geometric.prolongation, fixed_options);
    require(relative_action_difference(
                continued, restarted.prolongation) < 1.0e-12,
            "continued PCG path differs from the fixed-budget reference");
    require(path.steps() == 5,
            "continued PCG path reported the wrong checkpoint");
    const tgi::TwoGridCycle cycle_two(a, checkpoint_two, 1, 2);
    const tgi::TwoGridCycle cycle_five(a, continued, 1, 2);
    require(cycle_five.setup_report().interpolation_energy <=
                cycle_two.setup_report().interpolation_energy,
            "global PCG checkpoint increased interpolation energy");
    require_throws(
        [&]() { path.advance_to(4); },
        "continued PCG path accepted a decreasing checkpoint");

    tgi::AdaptiveGlobalPcgOptions adaptive_options;
    adaptive_options.policy = tgi::AdaptiveGlobalPcgPolicy::Fast;
    adaptive_options.maximum_cycles = 1000;
    adaptive_options.thread_count = 2;
    const tgi::Vector adaptive_rhs(
        static_cast<std::size_t>(grid.fine_size()), 1.0);
    const auto adaptive = tgi::build_adaptive_global_pcg_interpolation(
        grid, a, geometric.prolongation, adaptive_options, adaptive_rhs);
    require(adaptive.report.selected_steps == 5,
            "fast PCG selected the wrong scaled checkpoint");
    require(adaptive.report.candidate_count == 1 &&
                adaptive.report.pilot_cycles == 0,
            "fast PCG performed avoidable candidate probes");
    require(adaptive.report.maximum_sampled_steps == 5,
            "fast PCG reported the wrong scale-aware path budget");
    require(adaptive.prolongation->rows() == grid.fine_size(),
            "adaptive PCG returned an invalid prolongation");
    require(adaptive.cycle->coarse_size() == grid.coarse_size(),
            "adaptive PCG returned an invalid reusable hierarchy");
    const auto adaptive_solve = tgi::solve_two_grid(
        adaptive_rhs, *adaptive.cycle, 1.0e-6, 1000);
    require(adaptive_solve.converged,
            "single-RHS adaptive hierarchy did not converge");
    tgi::AdaptiveGlobalPcgOptions detailed_options = adaptive_options;
    detailed_options.policy = tgi::AdaptiveGlobalPcgPolicy::Reuse;
    const auto detailed = tgi::build_adaptive_global_pcg_interpolation(
        grid, a, geometric.prolongation, detailed_options, adaptive_rhs);
    require(detailed.report.pilot_cycles == 8,
            "reuse-aware selection reported the wrong scaled pilot budget");
    require(detailed.report.maximum_sampled_steps == 8,
            "reuse-aware selection reported the wrong path budget");
    require(detailed.report.candidate_count == 5,
            "reuse-aware PCG did not use its fixed candidate budget");
    const auto detailed_solve = tgi::solve_two_grid(
        adaptive_rhs, *detailed.cycle, 1.0e-6, 1000);
    require(detailed_solve.converged,
            "reuse-aware adaptive hierarchy did not converge");

    tgi::AdaptiveGlobalPcgOptions invalid_options = adaptive_options;
    invalid_options.maximum_cycles = 0;
    require_throws(
        [&]() {
            (void)tgi::build_adaptive_global_pcg_interpolation(
                grid, a, geometric.prolongation, invalid_options,
                adaptive_rhs);
        },
        "adaptive PCG accepted a nonpositive cycle limit");

    return 0;
}
