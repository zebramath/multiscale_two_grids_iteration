#pragma once

#include "multigrid/energy_interpolation.hpp"
#include "pde/diffusion_problem.hpp"

#include <stdexcept>
#include <string>

namespace experiment_support {

struct BasicConfig {
    int fine_intervals = 128;
    int coarse_intervals = 16;
    int threads = 4;
    int max_cycles = 40000;
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

inline BasicConfig parse_config(int argc, char** argv) {
    BasicConfig config;
    for (int index = 1; index < argc; ++index) {
        parse_basic_argument(config, argv[index]);
    }
    if (config.fine_intervals <= 0 || config.coarse_intervals <= 0 ||
        config.fine_intervals % config.coarse_intervals != 0 ||
        config.threads <= 0 || config.max_cycles <= 0 ||
        !(config.contrast >= 1.0)) {
        throw std::invalid_argument("invalid experiment configuration");
    }
    return config;
}

inline tgi::StructuredGrid make_grid(const BasicConfig& config) {
    return tgi::StructuredGrid(
        config.fine_intervals - 1,
        config.fine_intervals / config.coarse_intervals);
}

inline tgi::InterpolationOptions energy_options(
    int layers, int threads, double tolerance = 1.0e-6) {
    tgi::InterpolationOptions options;
    options.strategy = layers == 0
        ? tgi::InterpolationStrategy::GlobalEnergyMinimum
        : tgi::InterpolationStrategy::LocalEnergyMinimum;
    options.patch_layers = layers;
    options.local_tolerance = tolerance;
    options.local_max_iterations = 40000;
    options.thread_count = threads;
    return options;
}

} // namespace experiment_support
