#include "experiment/problem.hpp"
#include "multigrid/adaptive_global_pcg.hpp"
#include "multigrid/frontier_gain_support.hpp"
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
    adaptive_options.step_increment = 2;
    adaptive_options.patience = 2;
    adaptive_options.probe_count = 2;
    adaptive_options.power_iterations = 4;
    adaptive_options.thread_count = 2;
    adaptive_options.expected_rhs = 2;
    const auto adaptive = tgi::build_adaptive_global_pcg_interpolation(
        grid, a, geometric.prolongation, adaptive_options);
    require(adaptive.report.selected_steps >= 0 &&
                adaptive.report.selected_steps <= 8,
            "adaptive PCG selected an invalid checkpoint");
    require(!adaptive.report.history.empty(),
            "adaptive PCG did not record its decisions");
    require(adaptive.prolongation.rows() == grid.fine_size(),
            "adaptive PCG returned an invalid prolongation");

    tgi::InterpolationOptions local_options;
    local_options.strategy =
        tgi::InterpolationStrategy::LocalEnergyMinimum;
    local_options.patch_layers = 1;
    local_options.local_tolerance = 1.0e-7;
    local_options.local_max_iterations = 1000;
    local_options.thread_count = 2;
    const auto local = tgi::build_interpolation(grid, a, local_options);
    tgi::FrontierGainSupportOptions support_options;
    support_options.base_patch_layers = 1;
    support_options.maximum_rounds = 2;
    support_options.maximum_extra_nodes_per_column = 8;
    support_options.maximum_total_nodes_per_round = 12;
    support_options.maximum_nodes_per_column_per_round = 3;
    support_options.thread_count = 2;
    const auto frontier = tgi::build_frontier_gain_interpolation(
        grid, a, local.prolongation, local_options, support_options);
    require(frontier.supports.size() ==
                static_cast<std::size_t>(grid.coarse_size()),
            "frontier-gain returned the wrong number of supports");
    for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
        const auto& support =
            frontier.supports[static_cast<std::size_t>(coarse)];
        require(std::all_of(support.begin(), support.end(),
                    [&](int node) { return !grid.is_coarse_node(node); }),
                "frontier-gain inserted a coarse node into a local support");
        const int extra = static_cast<int>(support.size()) -
            static_cast<int>(grid.patch_f_nodes(coarse, 1).size());
        require(extra <= support_options.maximum_extra_nodes_per_column,
                "frontier-gain exceeded the per-column budget");
    }
    return 0;
}
