#pragma once

#include "experiment/problem.hpp"

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace experiment_support {

struct ComparisonCase {
    int fine = 32;
    int coarse = 8;
    double contrast = 1.0e4;
    std::uint64_t seed = 1;
    FieldCase field;
};

inline std::vector<ComparisonCase> comparison_cases(bool quick) {
    const std::array<std::pair<int, int>, 4> grids{
        std::pair<int, int>{32, 8}, {64, 8}, {64, 16}, {128, 16}};
    const std::array<double, 3> contrasts{1.0e2, 1.0e4, 1.0e6};
    const auto& topologies = channel_topologies();
    std::vector<ComparisonCase> cases;
    int index = 0;
    for (const auto& grid : grids) {
        for (double contrast : contrasts) {
            cases.push_back({
                grid.first, grid.second, contrast,
                index % 2 == 0 ? 1U : 17U,
                topologies[static_cast<std::size_t>(index) %
                           topologies.size()]});
            ++index;
        }
    }
    if (!quick) return cases;
    return {cases[0], cases[4], cases[8], cases[10]};
}

inline BasicConfig comparison_config(
    const ComparisonCase& item, int threads, int maximum_cycles) {
    BasicConfig config;
    config.fine_intervals = item.fine;
    config.coarse_intervals = item.coarse;
    config.contrast = item.contrast;
    config.threads = threads;
    config.max_cycles = maximum_cycles;
    return config;
}

}
