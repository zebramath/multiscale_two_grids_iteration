#pragma once

#include "experiment/common.hpp"
#include "multigrid/algebraic_interpolation.hpp"

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

inline const Row& study_headers() {
    static const Row headers{
        "Field", "Method", "Parameter", "P nnz", "Density %", "Ac nnz",
        "F residual", "Energy error", "Build ms", "TG setup ms",
        "Setup ms", "Solve ms", "Total ms", "Cycles", "Rho",
        "Construct iters"};
    return headers;
}

inline std::vector<int> study_widths() {
    return {11, 18, 13, 9, 10, 9, 12, 12,
            10, 11, 10, 10, 10, 9, 8, 15};
}

inline Row evaluate_candidate(
    const std::string& field, const tgi::StructuredGrid& grid,
    const tgi::SparseMatrix& a, const tgi::Vector& rhs,
    const tgi::SparseMatrix& global_reference,
    const StudyCandidate& candidate, const BasicConfig& config,
    int spectral_iterations = 80) {
    const auto error = compare_prolongations_global(
        grid, a, global_reference, candidate.prolongation);
    const double f_residual =
        tgi::algebraic_interpolation_detail::scaled_f_residual(
            grid, a, candidate.prolongation, config.threads);
    // One timed run is intentional in v2.5: the table reports a reproducible
    // workload estimate without hiding setup/solve costs behind repeated
    // measurements.  The first run still performs the optional spectral
    // estimate, which is explicitly excluded by evaluate_two_grid.
    const CycleMetrics cycle = evaluate_two_grid(
        a, rhs, candidate.prolongation, config.threads,
        1.0e-6, config.max_cycles, spectral_iterations);
    const double coarse_setup_ms = cycle.coarse_setup_ms;
    const double solve_ms = cycle.solve_ms;
    const double cycle_total_ms = cycle.total_ms;
    const double setup_ms = candidate.build_ms + coarse_setup_ms;
    const double total_ms = candidate.build_ms + cycle_total_ms;
    return {
        field,
        candidate.method,
        candidate.parameter,
        std::to_string(candidate.prolongation.nnz()),
        fixed(interpolation_density_percent(candidate.prolongation), 4),
        std::to_string(cycle.coarse_nnz),
        scientific(f_residual),
        scientific(error.aggregate_relative_energy_error),
        fixed(candidate.build_ms),
        fixed(coarse_setup_ms),
        fixed(setup_ms),
        fixed(solve_ms),
        fixed(total_ms),
        cycle.converged
            ? std::to_string(cycle.cycles)
            : "failed@" + std::to_string(cycle.cycles),
        fixed(cycle.convergence_factor, 4),
        fixed(candidate.mean_construction_iterations, 1)};
}

inline Summary fixed_study_summary(
    const BasicConfig& config, const std::string& extra_label = {},
    const std::string& extra_value = {}) {
    Summary summary{
        {"Fine grid", "h=1/" + std::to_string(config.fine_intervals)},
        {"Coarse grid", "H=1/" + std::to_string(config.coarse_intervals)},
        {"Contrast", scientific(config.contrast, 0)},
        {"Smoother", "one forward + one backward Gauss-Seidel"},
        {"Outer tolerance", "1e-6"},
        {"Spectral iterations", "80 (excluded from timing)"}};
    summary.push_back({"Timing", "single two-grid run (no median/repetition)"});
    if (!extra_label.empty()) {
        summary.push_back({extra_label, extra_value});
    }
    return summary;
}

inline tgi::InterpolationResult build_global_reference(
    const tgi::StructuredGrid& grid, const tgi::SparseMatrix& a,
    int threads) {
    auto options = energy_options(0, threads, 1.0e-10);
    options.drop_tolerance = 0.0;
    return tgi::build_interpolation(grid, a, options);
}

} // namespace experiment_support
