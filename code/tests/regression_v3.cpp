#include "experiment/problem.hpp"
#include "multigrid/adaptive_global_pcg.hpp"
#include "multigrid/global_pcg_path.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
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

} // namespace

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

    tgi::AdaptiveGlobalPcgOptions adaptive_options;
    adaptive_options.minimum_steps = 2;
    adaptive_options.maximum_steps = 8;
    adaptive_options.maximum_screening_steps = 8;
    adaptive_options.screening_increment = 2;
    adaptive_options.screening_pilot_iterations = 4;
    adaptive_options.screening_tail_window = 2;
    adaptive_options.minimum_screened_positive_candidates = 2;
    adaptive_options.refinement_backtrack_steps = 2;
    adaptive_options.refinement_stop_before_anchor_steps = 0;
    adaptive_options.refinement_increment = 2;
    adaptive_options.refinement_pilot_iterations = 4;
    adaptive_options.refinement_tail_window = 2;
    adaptive_options.confirmation_candidates = 2;
    adaptive_options.easy_accept_cycles = 1000;
    adaptive_options.medium_accept_cycles = 1000;
    adaptive_options.initial_safety_forecast_cycles = 1000;
    adaptive_options.maximum_confirmation_cycles = 1000;
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
    require(adaptive.report.history.size() <= 4U,
            "budget adaptive PCG exceeded its finite candidate budget");
    require(!adaptive.report.selected_cycles_confirmed,
            "budget adaptive PCG unexpectedly ran a full confirmation solve");
    require(adaptive.prolongation.rows() == grid.fine_size(),
            "adaptive PCG returned an invalid prolongation");

    tgi::AdaptiveGlobalPcgOptions invalid_options = adaptive_options;
    invalid_options.maximum_confirmation_cycles = 0;
    bool rejected_invalid_budget = false;
    try {
        (void)tgi::build_adaptive_global_pcg_interpolation(
            grid, a, geometric.prolongation,
            invalid_options, &adaptive_rhs);
    } catch (const std::invalid_argument&) {
        rejected_invalid_budget = true;
    }
    require(rejected_invalid_budget,
            "adaptive PCG accepted a nonpositive confirmation cap");

    adaptive_options.cost_aware_mode = false;
    const auto staged = tgi::build_adaptive_global_pcg_interpolation(
        grid, a, geometric.prolongation, adaptive_options, &adaptive_rhs);
    require(staged.report.selected_cycles > 0,
            "staged adaptive PCG did not confirm a candidate");
    require(staged.report.selected_cycles_confirmed,
            "staged adaptive PCG lost its confirmation status");
    require(staged.prolongation.rows() == grid.fine_size(),
            "staged adaptive PCG returned an invalid prolongation");

    return 0;
}
