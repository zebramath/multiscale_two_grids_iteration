#include "experiment/study.hpp"
#include "multigrid/global_pcg.hpp"

#include <chrono>
#include <string>
#include <utility>
#include <vector>

namespace {

struct DiagnosticMeasurement {
    std::string method;
    std::string parameter;
    double density_percent = 0.0;
    double build_ms = 0.0;
    double hierarchy_and_solve_ms = 0.0;
    double rho_estimate = 0.0;
    double pcg_energy_residual = 0.0;
    bool has_pcg_energy_residual = false;
    int cycles = 0;
    bool converged = false;
    tgi::CoarseSetupReport setup;
};

DiagnosticMeasurement measure_candidate(
    const tgi::SparseMatrix& a, const tgi::Vector& rhs,
    const tgi::SparseMatrix& p, std::string method,
    std::string parameter, double build_ms, int threads,
    int maximum_cycles, double pcg_energy_residual = 0.0,
    bool has_pcg_energy_residual = false) {
    DiagnosticMeasurement measurement;
    measurement.method = std::move(method);
    measurement.parameter = std::move(parameter);
    measurement.density_percent =
        experiment_support::interpolation_density_percent(p);
    measurement.build_ms = build_ms;
    measurement.pcg_energy_residual = pcg_energy_residual;
    measurement.has_pcg_energy_residual = has_pcg_energy_residual;

    const auto begin = std::chrono::steady_clock::now();
    const tgi::TwoGridCycle cycle(a, p, 1, threads);
    const auto solved = tgi::solve_two_grid(
        rhs, cycle, 1.0e-6, maximum_cycles);
    measurement.hierarchy_and_solve_ms =
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - begin).count();
    measurement.setup = cycle.setup_report();
    measurement.cycles = solved.cycles;
    measurement.converged = solved.converged;
    measurement.rho_estimate = cycle.estimate_convergence_factor(60);
    return measurement;
}

experiment_support::Row standard_row(
    const std::string& field, const DiagnosticMeasurement& item) {
    return {
        field,
        item.method,
        item.parameter,
        experiment_support::fixed(item.density_percent, 4),
        experiment_support::fixed(
            item.build_ms + item.setup.total_ms),
        experiment_support::fixed(
            item.build_ms + item.hierarchy_and_solve_ms),
        item.converged
            ? std::to_string(item.cycles)
            : "failed@" + std::to_string(item.cycles)};
}

experiment_support::Row diagnostic_row(
    const DiagnosticMeasurement& item, double exact_energy,
    double initial_energy_excess) {
    const double normalized_excess = initial_energy_excess > 0.0
        ? (item.setup.interpolation_energy - exact_energy) /
            initial_energy_excess
        : 0.0;
    return {
        item.method,
        item.parameter,
        item.has_pcg_energy_residual
            ? experiment_support::scientific(
                  item.pcg_energy_residual, 3)
            : "-",
        experiment_support::scientific(normalized_excess, 3),
        experiment_support::fixed(item.density_percent, 4),
        experiment_support::fixed(
            item.setup.interpolation_complexity, 4),
        std::to_string(item.setup.coarse_nnz),
        std::to_string(item.setup.factor_nnz),
        experiment_support::fixed(
            item.setup.two_grid_operator_complexity, 4),
        experiment_support::fixed(item.setup.factor_fill_ratio, 4),
        experiment_support::fixed(item.rho_estimate, 6),
        item.converged
            ? std::to_string(item.cycles)
            : "failed@" + std::to_string(item.cycles)};
}

}

int run_pcg_dense_scan(int argc, char** argv) {
    const auto config = experiment_support::parse_config(argc, argv);
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);
    const auto& channel = experiment_support::standard_fields().at(1);
    const auto problem = experiment_support::make_problem(
        grid, channel, config);
    const auto& a = problem.matrix;
    const auto& rhs = problem.rhs;

    auto geometric = experiment_support::geometric_interpolation(grid, a);
    const double geometric_ms = geometric.report.timing.total_ms;
    auto global = experiment_support::build_global_reference(
        grid, a, config.threads);
    const double exact_build_ms = global.report.timing.total_ms;

    const DiagnosticMeasurement geometric_measurement = measure_candidate(
        a, rhs, geometric.prolongation, "geometric", "P_G",
        geometric_ms, config.threads, config.max_cycles, 1.0, true);
    const DiagnosticMeasurement exact_measurement = measure_candidate(
        a, rhs, global.prolongation, "global-exact", "tol=1e-10",
        exact_build_ms, config.threads, config.max_cycles);
    const double exact_energy =
        exact_measurement.setup.interpolation_energy;
    const double initial_energy_excess =
        geometric_measurement.setup.interpolation_energy - exact_energy;

    experiment_support::Rows rows;
    std::vector<DiagnosticMeasurement> measurements;
    measurements.push_back(geometric_measurement);
    tgi::GlobalEnergyPcgPath path(
        grid, a, geometric.prolongation, config.threads);
    for (int steps = 16; steps <= 64; steps += 2) {
        path.advance_to(steps);
        const auto path_report = path.report();
        DiagnosticMeasurement measurement = measure_candidate(
            a, rhs, path.prolongation(0.0), "PCG-global",
            "m=" + std::to_string(steps),
            geometric_ms + path_report.total_ms,
            config.threads, config.max_cycles,
            path_report.relative_preconditioned_residual, true);
        rows.push_back(standard_row(channel.name, measurement));
        measurements.push_back(std::move(measurement));
    }
    rows.push_back(standard_row(channel.name, exact_measurement));
    measurements.push_back(exact_measurement);

    experiment_support::Rows diagnostic_rows;
    for (const auto& measurement : measurements) {
        diagnostic_rows.push_back(diagnostic_row(
            measurement, exact_energy, initial_energy_excess));
    }

    experiment_support::Report report(
        "Channel finite PCG-step scan");
    report.add_summary(experiment_support::fixed_study_summary(
        config, "PCG scan", "m=16,18,...,64"));
    report.add_note(
        "Every finite checkpoint reuses one incremental global PCG path. "
        "Each interpolation is still evaluated by an independent Galerkin "
        "hierarchy and an independent solve. Setup ms is the cost of reaching "
        "that checkpoint plus its production hierarchy.");
    report.add_table(
        "Channel finite-step PCG scan", experiment_support::study_headers(),
        experiment_support::study_widths(), rows);
    report.add_note(
        "The normalized energy excess is "
        "[trace(P^T A P)-trace(P_*^T A P_*)] divided by the same difference "
        "for P_G. C_P=nnz(P)/nnz(A), C_op=(nnz(A)+nnz(A_c))/nnz(A), and "
        "factor fill=nnz(L+L^T-diag(L))/nnz(A_c). rho_est uses 60 power "
        "iterations of "
        "the implemented symmetric two-grid error propagation. Diagnostic "
        "power iterations are excluded from Setup and Total timing.");
    report.add_table(
        "Energy, hierarchy and spectral diagnostics",
        {"Method", "Parameter", "PCG eta", "Energy excess", "P density %",
         "C_P", "Ac nnz", "L nnz", "C_op", "L/Ac", "rho_est", "Cycles"},
        {12, 10, 11, 13, 11, 8, 9, 9, 8, 8, 10, 12},
        diagnostic_rows, true);
    report.save("pcg_dense_scan");
    return 0;
}
