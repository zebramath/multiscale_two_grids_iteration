#pragma once

#include "experiment/problem.hpp"
#include "experiment/reporting.hpp"
#include "multigrid/two_grid_solver.hpp"
#include "version.hpp"

#include <chrono>
#include <cstddef>
#include <string>
#include <utility>

namespace experiment_support {

struct StudyCandidate {
    std::string method;
    std::string parameter;
    tgi::SparseMatrix prolongation;
    double build_ms = 0.0;
};

inline StudyCandidate make_candidate(
    std::string method, std::string parameter,
    tgi::InterpolationResult interpolation) {
    return {
        std::move(method), std::move(parameter),
        std::move(interpolation.prolongation),
        interpolation.report.timing.total_ms};
}

inline double interpolation_density_percent(
    const tgi::SparseMatrix& prolongation) {
    const double entries = static_cast<double>(prolongation.rows()) *
        static_cast<double>(prolongation.cols());
    return entries > 0.0
        ? 100.0 * static_cast<double>(prolongation.nnz()) / entries
        : 0.0;
}

struct CycleMetrics {
    int cycles = 0;
    bool converged = false;
    double coarse_setup_ms = 0.0;
    double total_ms = 0.0;
};

inline CycleMetrics evaluate_two_grid(
    const tgi::SparseMatrix& matrix, const tgi::Vector& rhs,
    const tgi::SparseMatrix& prolongation, int setup_threads,
    double solve_tolerance, int max_cycles) {
    const auto begin = std::chrono::steady_clock::now();
    const tgi::TwoGridCycle cycle(
        matrix, prolongation, 1, setup_threads);
    const auto result = tgi::solve_two_grid(
        matrix, rhs, cycle, solve_tolerance, max_cycles);
    const auto end = std::chrono::steady_clock::now();
    return {
        result.cycles,
        result.converged,
        cycle.setup_report().total_ms,
        std::chrono::duration<double, std::milli>(end - begin).count()};
}

inline const Row& study_headers() {
    static const Row headers{
        "Field", "Method", "Parameter", "P density %",
        "Setup ms", "Total ms", "Cycles"};
    return headers;
}

inline std::vector<int> study_widths() {
    return {11, 18, 20, 11, 10, 10, 12};
}

inline Row evaluate_candidate(
    const std::string& field, const tgi::SparseMatrix& matrix,
    const tgi::Vector& rhs, const StudyCandidate& candidate,
    const BasicConfig& config) {
    const CycleMetrics cycle = evaluate_two_grid(
        matrix, rhs, candidate.prolongation, config.threads,
        1.0e-6, config.max_cycles);
    const double setup_ms = candidate.build_ms + cycle.coarse_setup_ms;
    const double total_ms = candidate.build_ms + cycle.total_ms;
    return {
        field,
        candidate.method,
        candidate.parameter,
        fixed(interpolation_density_percent(candidate.prolongation), 4),
        fixed(setup_ms),
        fixed(total_ms),
        cycle.converged
            ? std::to_string(cycle.cycles)
            : "failed@" + std::to_string(cycle.cycles)};
}

inline Summary fixed_study_summary(
    const BasicConfig& config, const std::string& extra_label = {},
    const std::string& extra_value = {}) {
    Summary summary{
        {"Version", std::string(tgi::version)},
        {"Fine grid", "h=1/" + std::to_string(config.fine_intervals)},
        {"Coarse grid", "H=1/" + std::to_string(config.coarse_intervals)},
        {"Contrast", scientific(config.contrast, 0)},
        {"Smoother", "one forward + one backward Gauss-Seidel"},
        {"Solve tolerance", "1e-6"},
        {"Timing", "one end-to-end run per candidate"}};
    if (!extra_label.empty()) summary.push_back({extra_label, extra_value});
    return summary;
}

} // namespace experiment_support
