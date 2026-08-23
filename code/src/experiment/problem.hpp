#pragma once

#include "experiment/config.hpp"

#include <cmath>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace experiment_support {

struct FieldCase {
    std::string name;
    tgi::CoefficientDistribution distribution;
};

inline const std::vector<FieldCase>& standard_fields() {
    static const std::vector<FieldCase> fields{
        {"continuous", tgi::CoefficientDistribution::RandomContinuous},
        {"channel", tgi::CoefficientDistribution::ChannelizedBinary},
        {"checker", tgi::CoefficientDistribution::RandomBinaryCheckerboard}
    };
    return fields;
}

inline tgi::Vector manufactured_solution(
    const tgi::StructuredGrid& grid) {
    constexpr double pi = 3.141592653589793238462643383279502884;
    tgi::Vector exact(static_cast<std::size_t>(grid.fine_size()));
    for (int id = 0; id < grid.fine_size(); ++id) {
        const auto [ix, iy] = grid.fine_coords(id);
        const double x = static_cast<double>(ix + 1) * grid.h();
        const double y = static_cast<double>(iy + 1) * grid.h();
        exact[static_cast<std::size_t>(id)] =
            std::sin(pi * x) * std::sin(pi * y) +
            0.23 * std::sin(3.0 * pi * x) *
                std::sin(2.0 * pi * y) +
            0.11 * std::sin(7.0 * pi * x) *
                std::sin(5.0 * pi * y);
    }
    return exact;
}

struct ExperimentProblem {
    std::string field_name;
    tgi::SparseMatrix matrix;
    tgi::Vector rhs;
};

inline ExperimentProblem make_problem(
    const tgi::StructuredGrid& grid, const FieldCase& field,
    const BasicConfig& config, std::uint64_t seed = 1) {
    tgi::CoefficientOptions options;
    options.distribution = field.distribution;
    options.contrast = config.contrast;
    options.seed = seed;
    options.checkerboard_block_size = grid.ratio();
    options.channel_background_block_size = grid.ratio();
    options.channel_width_fine_cells = 2;
    auto coefficient = tgi::make_coefficient(grid, options);
    auto matrix = tgi::assemble_diffusion(grid, coefficient.values);
    auto rhs = matrix.multiply(manufactured_solution(grid));
    return {field.name, std::move(matrix), std::move(rhs)};
}

inline tgi::InterpolationResult geometric_interpolation(
    const tgi::StructuredGrid& grid, const tgi::SparseMatrix& matrix) {
    tgi::InterpolationOptions options;
    options.strategy = tgi::InterpolationStrategy::GeometricBilinear;
    return tgi::build_interpolation(grid, matrix, options);
}

inline tgi::InterpolationResult build_global_reference(
    const tgi::StructuredGrid& grid, const tgi::SparseMatrix& matrix,
    int threads) {
    auto options = energy_options(0, threads, 1.0e-10);
    options.drop_tolerance = 0.0;
    return tgi::build_interpolation(grid, matrix, options);
}

} // namespace experiment_support
