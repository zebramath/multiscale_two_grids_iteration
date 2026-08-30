#pragma once

#include "multigrid/energy_interpolation.hpp"
#include "pde/diffusion_problem.hpp"

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace experiment_support {

inline constexpr int maximum_two_grid_cycles = 20000;
inline constexpr int maximum_geometric_cycles = 30000;

struct BasicConfig {
    int fine_intervals = 128;
    int coarse_intervals = 16;
    int threads = 4;
    double contrast = 1.0e4;
};

inline tgi::StructuredGrid make_grid(const BasicConfig& config) {
    return tgi::StructuredGrid(
        config.fine_intervals - 1,
        config.fine_intervals / config.coarse_intervals);
}

struct FieldCase {
    std::string name;
    tgi::CoefficientDistribution distribution;
};

inline const std::vector<FieldCase>& channel_topologies() {
    static const std::vector<FieldCase> fields{
        {"cross-channel", tgi::CoefficientDistribution::ChannelizedBinary},
        {"meandering-channel",
         tgi::CoefficientDistribution::MeanderingChannelBinary},
        {"diagonal-channels",
         tgi::CoefficientDistribution::DiagonalChannelsBinary},
        {"parallel-channels",
         tgi::CoefficientDistribution::ParallelChannelsBinary},
        {"branching-channels",
         tgi::CoefficientDistribution::BranchingChannelsBinary},
        {"winding-ring",
         tgi::CoefficientDistribution::WindingRingBinary}
    };
    return fields;
}

struct ExperimentProblem {
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
    options.channel_background_block_size = grid.ratio();
    options.channel_width_fine_cells = 2;
    auto coefficient = tgi::make_coefficient(grid, options);
    auto matrix = tgi::assemble_diffusion(grid, coefficient.values);
    tgi::Vector rhs(static_cast<std::size_t>(grid.fine_size()), 1.0);
    return {std::move(matrix), std::move(rhs)};
}

inline tgi::InterpolationResult build_global_reference(
    const tgi::StructuredGrid& grid, const tgi::SparseMatrix& matrix,
    int threads) {
    tgi::GlobalEnergyOptions options;
    options.tolerance = 1.0e-10;
    options.thread_count = threads;
    options.drop_tolerance = 0.0;
    return tgi::build_global_energy_interpolation(grid, matrix, options);
}

}
