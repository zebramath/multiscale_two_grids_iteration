#pragma once

#include "experiment/candidate.hpp"
#include "experiment/config.hpp"
#include "experiment/reporting.hpp"
#include "multigrid/two_grid_solver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>

namespace experiment_support {

struct CycleMetrics {
    int cycles = 0;
    bool converged = false;
    double average_convergence_factor = 0.0;
    double coarse_setup_ms = 0.0;
    double solve_ms = 0.0;
    double total_ms = 0.0;
    std::size_t coarse_nnz = 0;
    std::size_t factor_nnz = 0;
};

inline CycleMetrics evaluate_two_grid(
    const tgi::SparseMatrix& matrix, const tgi::Vector& rhs,
    const tgi::SparseMatrix& prolongation, int setup_threads,
    double solve_tolerance, int max_cycles) {
    const auto begin = std::chrono::steady_clock::now();
    const tgi::TwoGridCycle cycle(
        matrix, prolongation, 1, setup_threads);
    const auto solve_begin = std::chrono::steady_clock::now();
    const auto result = tgi::solve_two_grid(
        matrix, rhs, cycle, solve_tolerance, max_cycles);
    const auto end = std::chrono::steady_clock::now();

    CycleMetrics metrics;
    metrics.cycles = result.cycles;
    metrics.converged = result.converged;
    if (result.cycles > 0) {
        const double ratio = std::max(
            result.relative_residual,
            std::numeric_limits<double>::min());
        metrics.average_convergence_factor = std::pow(
            ratio, 1.0 / static_cast<double>(result.cycles));
    }
    metrics.coarse_setup_ms = cycle.setup_report().total_ms;
    metrics.solve_ms =
        std::chrono::duration<double, std::milli>(end - solve_begin).count();
    metrics.total_ms =
        std::chrono::duration<double, std::milli>(end - begin).count();
    metrics.coarse_nnz = cycle.setup_report().coarse_nnz;
    metrics.factor_nnz = cycle.setup_report().factor_nnz;
    return metrics;
}

inline const Row& study_headers() {
    static const Row headers{
        "Field", "Method", "Parameter", "P nnz", "P density %",
        "Ac nnz", "L nnz", "Build ms", "TG setup ms", "Setup ms",
        "Solve ms", "Total ms", "Cycles", "Rho avg", "Build iters",
        "Converged"};
    return headers;
}

inline std::vector<int> study_widths() {
    return {11, 18, 20, 9, 11, 9, 9, 10,
            11, 10, 10, 10, 8, 8, 11, 9};
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
        std::to_string(candidate.prolongation.nnz()),
        fixed(interpolation_density_percent(candidate.prolongation), 4),
        std::to_string(cycle.coarse_nnz),
        std::to_string(cycle.factor_nnz),
        fixed(candidate.build_ms),
        fixed(cycle.coarse_setup_ms),
        fixed(setup_ms),
        fixed(cycle.solve_ms),
        fixed(total_ms),
        std::to_string(cycle.cycles),
        fixed(cycle.average_convergence_factor, 4),
        fixed(candidate.mean_construction_iterations, 1),
        cycle.converged ? "yes" : "no"};
}

inline Summary fixed_study_summary(
    const BasicConfig& config, const std::string& extra_label = {},
    const std::string& extra_value = {}) {
    Summary summary{
        {"Version", "2.7.0"},
        {"Fine grid", "h=1/" + std::to_string(config.fine_intervals)},
        {"Coarse grid", "H=1/" + std::to_string(config.coarse_intervals)},
        {"Contrast", scientific(config.contrast, 0)},
        {"Smoother", "one forward + one backward Gauss-Seidel"},
        {"Solve tolerance", "1e-6"},
        {"Rho avg", "geometric mean over the actual solve"},
        {"Timing", "one end-to-end run per candidate"}};
    if (!extra_label.empty()) {
        summary.push_back({extra_label, extra_value});
    }
    return summary;
}

} // namespace experiment_support
