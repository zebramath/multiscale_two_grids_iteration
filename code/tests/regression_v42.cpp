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

    tgi::InterpolationOptions geometric_options;
    geometric_options.strategy =
        tgi::InterpolationStrategy::GeometricBilinear;
    const auto geometric =
        tgi::build_interpolation(grid, a, geometric_options);

    tgi::GlobalEnergyPcgPath path(
        grid, a, geometric.prolongation, 2);
    path.advance_to(2);
    const tgi::SparseMatrix checkpoint_two = path.prolongation(0.0);
    path.advance_to(5);
    const tgi::SparseMatrix continued = path.prolongation(0.0);
    tgi::InterpolationOptions fixed_options;
    fixed_options.strategy =
        tgi::InterpolationStrategy::GlobalEnergyMinimum;
    fixed_options.local_tolerance = 0.0;
    fixed_options.local_max_iterations = 5;
    fixed_options.thread_count = 2;
    fixed_options.drop_tolerance = 0.0;
    fixed_options.require_convergence = false;
    const auto restarted = tgi::refine_global_energy_interpolation(
        grid, a, geometric.prolongation, fixed_options);
    require(relative_action_difference(
                continued, restarted.prolongation) < 1.0e-12,
            "continued PCG path differs from the fixed-budget reference");
    require(path.report().steps == 5,
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
    adaptive_options.minimum_steps = 2;
    adaptive_options.maximum_steps = 8;
    adaptive_options.step_quantum = 2;
    adaptive_options.expected_rhs_count = 1.0;
    adaptive_options.maximum_cycles = 1000;
    adaptive_options.thread_count = 2;
    const tgi::Vector adaptive_rhs(
        static_cast<std::size_t>(grid.fine_size()), 1.0);
    const auto adaptive = tgi::build_adaptive_global_pcg_interpolation(
        grid, a, geometric.prolongation, adaptive_options, &adaptive_rhs);
    require(adaptive.report.selected_steps >= 0 &&
                adaptive.report.selected_steps <= 8,
            "adaptive PCG selected an invalid checkpoint");
    require(!adaptive.report.history.empty(),
            "adaptive PCG did not record its decisions");
    require(adaptive.report.history.size() <= 3U,
            "economy PCG exceeded its uniform candidate budget");
    require(!adaptive.report.used_local_refinement,
            "single-right-hand-side selection unexpectedly refined");
    require(adaptive.prolongation->rows() == grid.fine_size(),
            "adaptive PCG returned an invalid prolongation");
    require(adaptive.cycle->coarse_size() == grid.coarse_size(),
            "adaptive PCG returned an invalid reusable hierarchy");
    require_throws(
        [&]() {
            (void)tgi::build_adaptive_global_pcg_interpolation(
                grid, a, geometric.prolongation, adaptive_options, nullptr);
        },
        "adaptive PCG accepted a missing representative right-hand side");
    tgi::AdaptiveGlobalPcgOptions detailed_options = adaptive_options;
    detailed_options.expected_rhs_count = 128.0;
    const auto detailed = tgi::build_adaptive_global_pcg_interpolation(
        grid, a, geometric.prolongation, detailed_options, &adaptive_rhs);
    require(detailed.report.pilot_iterations >
                adaptive.report.pilot_iterations,
            "reuse-aware selection did not increase its evidence budget");
    require(detailed.report.history.size() <= 5U,
            "reuse-aware PCG exceeded its candidate budget");

    tgi::AdaptiveGlobalPcgOptions invalid_options = adaptive_options;
    invalid_options.expected_rhs_count = 0.5;
    require_throws(
        [&]() {
            (void)tgi::build_adaptive_global_pcg_interpolation(
                grid, a, geometric.prolongation, invalid_options,
                &adaptive_rhs);
        },
        "adaptive PCG accepted fewer than one expected right-hand side");

    return 0;
}
