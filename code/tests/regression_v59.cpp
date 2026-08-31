#include "multigrid/global_pcg.hpp"

#include <cmath>
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
            "channel topology lost the requested contrast");
    const tgi::SparseMatrix a =
        tgi::assemble_diffusion(grid, coefficient.values);

    const auto geometric = tgi::build_geometric_interpolation(grid);

    tgi::GlobalEnergyPcgPath path(
        grid, a, geometric.prolongation, 2);
    path.advance_to(2);
    const tgi::SparseMatrix checkpoint_two = path.prolongation();
    path.advance_to(5);
    const tgi::SparseMatrix continued = path.prolongation();
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

    const tgi::Vector adaptive_rhs(
        static_cast<std::size_t>(grid.fine_size()), 1.0);
    const auto adaptive = tgi::build_adaptive_global_pcg_interpolation(
        grid, a, geometric.prolongation, 2);
    require(adaptive.report.selected_steps == 2,
            "adaptive PCG selected the wrong scaled checkpoint");
    require(adaptive.prolongation->rows() == grid.fine_size(),
            "adaptive PCG returned an invalid prolongation");
    require(adaptive.cycle->coarse_matrix().rows() == grid.coarse_size(),
            "adaptive PCG returned an invalid reusable hierarchy");
    const auto adaptive_solve = tgi::solve_two_grid(
        adaptive_rhs, *adaptive.cycle, 1.0e-6, 1000);
    require(adaptive_solve.converged,
            "adaptive hierarchy did not converge");

    const tgi::StructuredGrid doubled_grid(31, 4);
    coefficient_options.contrast = 1.0e4;
    const auto doubled_field = tgi::make_coefficient(
        doubled_grid, coefficient_options);
    const auto doubled_matrix = tgi::assemble_diffusion(
        doubled_grid, doubled_field.values);
    const int doubled_steps =
        tgi::adaptive_global_pcg_detail::select_steps(
            doubled_grid, doubled_matrix);
    require(doubled_steps == 4,
            "adaptive checkpoint ignored the coarse-resolution regime");

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
    const int high_contrast_steps =
        tgi::adaptive_global_pcg_detail::select_steps(
            doubled_grid, high_contrast_matrix);
    require(high_contrast_steps == 4,
            "adaptive coarse-resolution regime depended on contrast");

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
    const int low_steps = tgi::adaptive_global_pcg_detail::select_steps(
        banded_grid, diagonal_matrix(1.0e2));
    const int medium_steps = tgi::adaptive_global_pcg_detail::select_steps(
        banded_grid, diagonal_matrix(1.0e4));
    const int high_steps = tgi::adaptive_global_pcg_detail::select_steps(
        banded_grid, diagonal_matrix(1.0e6));
    require(low_steps == 9 && medium_steps == 12 && high_steps == 18,
            "adaptive stiffness bands selected the wrong checkpoint");
    require_throws(
        [&]() {
            (void)tgi::solve_two_grid(
                adaptive_rhs, *adaptive.cycle, 0.0, 1000);
        },
        "two-grid solve accepted a nonpositive tolerance");
    require_throws(
        [&]() {
            (void)tgi::solve_two_grid(
                adaptive_rhs, *adaptive.cycle, 1.0, 1000);
        },
        "two-grid solve accepted a noncontractive tolerance");
    require_throws(
        [&]() {
            (void)tgi::solve_two_grid(
                adaptive_rhs, *adaptive.cycle, 1.0e-6, 0);
        },
        "two-grid solve accepted a nonpositive cycle limit");

    return 0;
}
