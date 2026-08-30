#include "multigrid/global_pcg.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

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
    using tgi::adaptive_global_pcg_detail::forecast_cycles;
    require(forecast_cycles(10, 1.0e-2, 0.5, 1.0e-6, 1000) == 24,
            "contractive pilot produced the wrong cycle forecast");
    require(forecast_cycles(10, 1.0e-7, 0.5, 1.0e-6, 1000) == 10,
            "converged pilot did not preserve its cycle count");
    require(forecast_cycles(
                10, std::numeric_limits<double>::infinity(),
                0.9, 1.0e-6, 1000) == 1000 &&
                forecast_cycles(
                    10, 1.0e-2,
                    std::numeric_limits<double>::quiet_NaN(),
                    1.0e-6, 1000) == 1000 &&
                forecast_cycles(
                    10, -1.0, 0.5, 1.0e-6, 1000) == 1000 &&
                forecast_cycles(
                    10, 1.0e-2, 0.0, 1.0e-6, 1000) == 1000 &&
                forecast_cycles(
                    10, 1.0e-2, std::nextafter(1.0, 0.0),
                    1.0e-300, 1000) == 1000 &&
                forecast_cycles(1000, 1.0e-2, 0.5,
                                1.0e-6, 1000) == 1000,
            "invalid or unbounded pilot forecast escaped the cycle cap");

    const tgi::StructuredGrid grid(15, 4);
    tgi::CoefficientOptions coefficient_options;
    coefficient_options.distribution =
        tgi::CoefficientDistribution::MeanderingChannelBinary;
    coefficient_options.contrast = 1.0e4;
    coefficient_options.seed = 7;
    coefficient_options.channel_background_block_size = 4;
    const auto coefficient = tgi::make_coefficient(grid, coefficient_options);
    require(coefficient.actual_contrast >= 0.99e4,
            "channel topology lost the requested contrast");
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
    const tgi::TwoGridCycle cycle_two(a, checkpoint_two, 1, 2);
    const tgi::TwoGridCycle cycle_five(a, continued, 1, 2);
    require(cycle_five.setup_report().interpolation_energy <=
                cycle_two.setup_report().interpolation_energy,
            "global PCG checkpoint increased interpolation energy");
    require_throws(
        [&]() { path.advance_to(4); },
        "continued PCG path accepted a decreasing checkpoint");
    require_throws(
        [&]() { path.advance_to(-1); },
        "continued PCG path accepted a negative checkpoint");

    tgi::AdaptiveGlobalPcgOptions adaptive_options;
    adaptive_options.policy = tgi::AdaptiveGlobalPcgPolicy::Fast;
    adaptive_options.maximum_cycles = 1000;
    adaptive_options.thread_count = 2;
    const tgi::Vector adaptive_rhs(
        static_cast<std::size_t>(grid.fine_size()), 1.0);
    const auto adaptive = tgi::build_adaptive_global_pcg_interpolation(
        grid, a, geometric.prolongation, adaptive_options, adaptive_rhs);
    require(adaptive.report.selected_steps == 2,
            "fast PCG selected the wrong scaled checkpoint");
    require(adaptive.report.candidate_count == 1 &&
                adaptive.report.pilot_cycles == 0,
            "fast PCG performed avoidable candidate probes");
    require(adaptive.report.maximum_sampled_steps == 2,
            "fast PCG reported the wrong scale-aware path budget");
    require(adaptive.prolongation->rows() == grid.fine_size(),
            "adaptive PCG returned an invalid prolongation");
    require(adaptive.cycle->coarse_matrix().rows() == grid.coarse_size(),
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

    const tgi::StructuredGrid doubled_grid(31, 4);
    coefficient_options.contrast = 1.0e4;
    const auto doubled_field = tgi::make_coefficient(
        doubled_grid, coefficient_options);
    const auto doubled_matrix = tgi::assemble_diffusion(
        doubled_grid, doubled_field.values);
    const auto doubled_fast =
        tgi::adaptive_global_pcg_detail::selection_profile(
            doubled_grid, doubled_matrix,
            tgi::AdaptiveGlobalPcgPolicy::Fast);
    const auto doubled_reuse =
        tgi::adaptive_global_pcg_detail::selection_profile(
            doubled_grid, doubled_matrix,
            tgi::AdaptiveGlobalPcgPolicy::Reuse);
    require(doubled_fast.checkpoints == std::vector<int>{4},
            "fast checkpoint ignored the coarse-resolution regime");
    require(doubled_reuse.checkpoints ==
                std::vector<int>({4, 6, 8, 11, 16}) &&
                doubled_reuse.pilot_cycles == 16,
            "reuse checkpoint set did not preserve normalized positions");

    for (int resolution = 8; resolution <= 257; ++resolution) {
        auto within_rounding_bound = [resolution](
                                         int numerator, int denominator) {
            const int checkpoint =
                tgi::adaptive_global_pcg_detail::scaled_checkpoint(
                    resolution, numerator, denominator);
            const double realized = static_cast<double>(checkpoint) /
                static_cast<double>(resolution);
            const double requested = static_cast<double>(numerator) /
                static_cast<double>(denominator);
            return std::abs(realized - requested) <=
                0.5 / static_cast<double>(resolution) + 1.0e-15;
        };
        require(within_rounding_bound(1, 8) &&
                    within_rounding_bound(3, 16) &&
                    within_rounding_bound(1, 4) &&
                    within_rounding_bound(1, 3) &&
                    within_rounding_bound(1, 2),
                "normalized checkpoint exceeded the h/2 rounding bound");
    }

    coefficient_options.contrast = 1.0e6;
    const auto high_contrast_field = tgi::make_coefficient(
        doubled_grid, coefficient_options);
    const auto high_contrast_matrix = tgi::assemble_diffusion(
        doubled_grid, high_contrast_field.values);
    const auto high_contrast_fast =
        tgi::adaptive_global_pcg_detail::selection_profile(
            doubled_grid, high_contrast_matrix,
            tgi::AdaptiveGlobalPcgPolicy::Fast);
    require(high_contrast_fast.checkpoints == std::vector<int>{4},
            "fast coarse-resolution regime depended on contrast");

    const tgi::StructuredGrid banded_grid(35, 3);
    auto diagonal_matrix = [&](double ratio) {
        std::vector<tgi::Triplet> entries;
        entries.reserve(static_cast<std::size_t>(banded_grid.fine_size()));
        for (int row = 0; row < banded_grid.fine_size(); ++row) {
            entries.push_back({row, row, row == 0 ? ratio : 1.0});
        }
        return tgi::SparseMatrix(
            banded_grid.fine_size(), banded_grid.fine_size(), entries);
    };
    const auto low_profile =
        tgi::adaptive_global_pcg_detail::selection_profile(
            banded_grid, diagonal_matrix(1.0e2),
            tgi::AdaptiveGlobalPcgPolicy::Fast);
    const auto medium_profile =
        tgi::adaptive_global_pcg_detail::selection_profile(
            banded_grid, diagonal_matrix(1.0e4),
            tgi::AdaptiveGlobalPcgPolicy::Fast);
    const auto high_profile =
        tgi::adaptive_global_pcg_detail::selection_profile(
            banded_grid, diagonal_matrix(1.0e6),
            tgi::AdaptiveGlobalPcgPolicy::Fast);
    require(low_profile.checkpoints == std::vector<int>{9} &&
                medium_profile.checkpoints == std::vector<int>{12} &&
                high_profile.checkpoints == std::vector<int>{18},
            "fast stiffness bands selected the wrong normalized checkpoint");
    const auto detailed_solve = tgi::solve_two_grid(
        adaptive_rhs, *detailed.cycle, 1.0e-6, 1000);
    require(detailed_solve.converged,
            "reuse-aware adaptive hierarchy did not converge");
    require_throws(
        [&]() {
            (void)tgi::solve_two_grid(
                adaptive_rhs, *detailed.cycle, 0.0, 1000);
        },
        "two-grid solve accepted a nonpositive tolerance");
    require_throws(
        [&]() {
            (void)tgi::solve_two_grid(
                adaptive_rhs, *detailed.cycle, 1.0, 1000);
        },
        "two-grid solve accepted a noncontractive tolerance");
    require_throws(
        [&]() {
            (void)tgi::solve_two_grid(
                adaptive_rhs, *detailed.cycle, 1.0e-6, 0);
        },
        "two-grid solve accepted a nonpositive cycle limit");

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
