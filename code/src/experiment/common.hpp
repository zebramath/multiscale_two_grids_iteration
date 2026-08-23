#pragma once

#include "experiment/metrics.hpp"
#include "experiment/reporting.hpp"
#include "experiment/test_problem.hpp"
#include "multigrid/energy_interpolation.hpp"
#include "pde/diffusion_problem.hpp"

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace experiment_support {

struct FieldCase {
    std::string name;
    tgi::CoefficientDistribution distribution;
};

inline std::vector<FieldCase> standard_fields() {
    return {
        {"continuous", tgi::CoefficientDistribution::RandomContinuous},
        {"channel", tgi::CoefficientDistribution::ChannelizedBinary},
        {"checker", tgi::CoefficientDistribution::RandomBinaryCheckerboard}
    };
}

inline tgi::CoefficientField make_field(
    const tgi::StructuredGrid& grid, const FieldCase& field,
    double contrast, std::uint64_t seed = 1) {
    tgi::CoefficientOptions options;
    options.distribution = field.distribution;
    options.contrast = contrast;
    options.seed = seed;
    options.checkerboard_block_size = grid.ratio();
    options.channel_background_block_size = grid.ratio();
    options.channel_width_fine_cells = 2;
    return tgi::make_coefficient(grid, options);
}

inline tgi::InterpolationOptions energy_options(
    int layers, int threads, double tolerance = 1.0e-6) {
    tgi::InterpolationOptions options;
    options.strategy = layers == 0
        ? tgi::InterpolationStrategy::GlobalEnergyMinimum
        : tgi::InterpolationStrategy::LocalEnergyMinimum;
    options.patch_layers = layers;
    options.local_tolerance = tolerance;
    options.local_max_iterations = 20000;
    options.thread_count = threads;
    return options;
}

inline tgi::InterpolationResult geometric_interpolation(
    const tgi::StructuredGrid& grid, const tgi::SparseMatrix& a) {
    tgi::InterpolationOptions options;
    options.strategy = tgi::InterpolationStrategy::GeometricBilinear;
    return tgi::build_interpolation(grid, a, options);
}

struct BasicConfig {
    // The v2.5 baseline is deliberately fixed at H/h = 8 while doubling
    // both grids relative to v2.4.  This makes the localization trade-offs
    // visible without changing the coefficient-field or solver settings.
    int fine_intervals = 128;
    int coarse_intervals = 16;
    int threads = 4;
    int max_cycles = 4000;
    double contrast = 1.0e4;
};

inline void parse_basic_argument(
    BasicConfig& config, const std::string& argument) {
    if (argument.rfind("--fine=", 0) == 0) {
        config.fine_intervals = std::stoi(argument.substr(7));
    } else if (argument.rfind("--coarse=", 0) == 0) {
        config.coarse_intervals = std::stoi(argument.substr(9));
    } else if (argument.rfind("--threads=", 0) == 0) {
        config.threads = std::stoi(argument.substr(10));
    } else if (argument.rfind("--contrast=", 0) == 0) {
        config.contrast = std::stod(argument.substr(11));
    } else if (argument.rfind("--max-cycles=", 0) == 0) {
        config.max_cycles = std::stoi(argument.substr(13));
    } else {
        throw std::invalid_argument("unknown argument: " + argument);
    }
}

inline tgi::StructuredGrid make_grid(const BasicConfig& config) {
    if (config.fine_intervals <= 0 || config.coarse_intervals <= 0 ||
        config.fine_intervals % config.coarse_intervals != 0 ||
        config.threads <= 0 || config.max_cycles <= 0 ||
        !(config.contrast >= 1.0)) {
        throw std::invalid_argument("invalid experiment configuration");
    }
    return tgi::StructuredGrid(
        config.fine_intervals - 1,
        config.fine_intervals / config.coarse_intervals);
}

} // namespace experiment_support
