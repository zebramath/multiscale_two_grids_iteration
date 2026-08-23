#pragma once

#include "multigrid/energy_interpolation.hpp"

#include <cstddef>
#include <string>
#include <utility>

namespace experiment_support {

struct StudyCandidate {
    std::string method;
    std::string parameter;
    tgi::SparseMatrix prolongation;
    double build_ms = 0.0;
    double mean_construction_iterations = 0.0;
};

inline StudyCandidate make_candidate(
    std::string method, std::string parameter,
    tgi::InterpolationResult interpolation) {
    const int systems = interpolation.report.local_solves.systems;
    const double mean_iterations = systems > 0
        ? static_cast<double>(
              interpolation.report.local_solves.total_iterations) /
              static_cast<double>(systems)
        : 0.0;
    return {
        std::move(method), std::move(parameter),
        std::move(interpolation.prolongation),
        interpolation.report.timing.total_ms, mean_iterations};
}

inline double interpolation_density_percent(
    const tgi::SparseMatrix& prolongation) {
    const double entries = static_cast<double>(prolongation.rows()) *
        static_cast<double>(prolongation.cols());
    return entries > 0.0
        ? 100.0 * static_cast<double>(prolongation.nnz()) / entries
        : 0.0;
}

} // namespace experiment_support
