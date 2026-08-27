#pragma once

#include "experiment/problem.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace experiment_support {

struct ComparisonCase {
    std::string axis;
    int fine = 32;
    int coarse = 8;
    double contrast = 1.0e4;
    std::uint64_t seed = 1;
    FieldCase field;
};

inline std::vector<ComparisonCase> comparison_cases(bool quick) {
    const auto& topologies = channel_topologies();
    const std::vector<ComparisonCase> cases{
        {"size", 32, 8, 1.0e4, 1, topologies[0]},
        {"size", 64, 8, 1.0e4, 1, topologies[0]},
        {"size", 64, 16, 1.0e4, 1, topologies[0]},
        {"center", 128, 16, 1.0e4, 1, topologies[0]},
        {"large-scale", 256, 16, 1.0e4, 1, topologies[0]},
        {"large-scale", 256, 16, 1.0e4, 1, topologies[5]},
        {"contrast", 128, 16, 1.0e2, 1, topologies[0]},
        {"contrast", 128, 16, 1.0e6, 1, topologies[0]},
        {"topology", 128, 16, 1.0e4, 1, topologies[1]},
        {"topology", 128, 16, 1.0e4, 1, topologies[2]},
        {"topology", 128, 16, 1.0e4, 1, topologies[3]},
        {"topology", 128, 16, 1.0e4, 1, topologies[4]},
        {"topology", 128, 16, 1.0e4, 1, topologies[5]}};
    if (!quick) return cases;
    return {cases[0], cases[3], cases[12]};
}

inline BasicConfig comparison_config(
    const ComparisonCase& item, int threads) {
    BasicConfig config;
    config.fine_intervals = item.fine;
    config.coarse_intervals = item.coarse;
    config.contrast = item.contrast;
    config.threads = threads;
    return config;
}

}
