#include "experiment/reporting.hpp"
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
    require_throws(
        []() {
            (void)tgi::dot(tgi::Vector{1.0}, tgi::Vector{1.0, 2.0});
        },
        "dot product accepted mismatched dimensions");
    require_throws(
        []() {
            (void)tgi::SparseMatrix(
                2, 2, std::vector<int>{0, 2, 2},
                std::vector<int>{1, 0}, tgi::Vector{1.0, 1.0});
        },
        "CSR constructor accepted unsorted row structure");
    require_throws(
        []() {
            const tgi::SparseMatrix identity(
                2, 2, std::vector<tgi::Triplet>{{0, 0, 1.0}, {1, 1, 1.0}});
            tgi::SparseCholesky factor;
            factor.factorize(identity, {0, 0});
        },
        "SparseCholesky accepted a non-bijective permutation");
    require_throws(
        []() {
            tgi::SparseCholesky factor;
            tgi::Vector result;
            tgi::Vector work;
            factor.solve(tgi::Vector{1.0}, result, work);
        },
        "unfactorized SparseCholesky accepted a solve");
    require_throws(
        []() {
            experiment_support::Report report("schema test");
            report.add_table(
                "bad table", {"A", "B"}, {4, 4}, {{"only-one"}});
        },
        "report accepted a row with the wrong column count");

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
            "adaptive PCG returned an invalid two-grid cycle");
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

    tgi::FixedPhysicalCoefficientOptions physical_options;
    physical_options.distribution =
        tgi::CoefficientDistribution::ChannelizedBinary;
    physical_options.contrast = 1.0e4;
    physical_options.seed = 7;
    physical_options.background_blocks_per_direction = 8;
    physical_options.channel_width = 1.0 / 16.0;
    const auto physical_coarse = tgi::make_fixed_physical_coefficient(
        grid, physical_options);
    const auto physical_fine = tgi::make_fixed_physical_coefficient(
        doubled_grid, physical_options);
    for (int iy = 0; iy < grid.fine_n(); ++iy) {
        for (int ix = 0; ix < grid.fine_n(); ++ix) {
            const int coarse_id = grid.fine_id(ix, iy);
            const int fine_id = doubled_grid.fine_id(
                2 * (ix + 1) - 1, 2 * (iy + 1) - 1);
            require(
                physical_coarse.values[
                    static_cast<std::size_t>(coarse_id)] ==
                    physical_fine.values[
                        static_cast<std::size_t>(fine_id)],
                "fixed-physical coefficient changed at a shared node");
        }
    }

    tgi::GlobalEnergyPcgPath residual_path(
        grid, a, geometric.prolongation, 2);
    require_throws(
        [&]() { residual_path.advance_until_relative_residual(1.0); },
        "residual stopping accepted a noncontractive tolerance");
    require_throws(
        [&]() {
            residual_path.advance_until_relative_residual(1.0e-2, 0);
        },
        "residual stopping accepted a nonpositive iteration limit");
    const auto residual_report =
        residual_path.advance_until_relative_residual(1.0e-2, 1000);
    require(residual_report.failed_systems == 0 &&
                residual_report.maximum_relative_residual <= 1.0e-2,
            "columnwise residual stopping missed its tolerance");
    require_throws(
        [&]() { residual_path.advance_to(10); },
        "residual-stopped PCG path was incorrectly continued");

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
