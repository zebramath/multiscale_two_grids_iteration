#include "experiment/problem.hpp"
#include "experiment/reporting.hpp"
#include "multigrid/two_grid_solver.hpp"

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace {

struct Measurement {
    std::string method;
    std::string construction;
    std::uint64_t cycles = 0;
    double relative_residual = 1.0;
    double effective_factor = 1.0;
    double density = 0.0;
    double energy = 0.0;
    std::size_t coarse_nnz = 0;
};

Measurement measure_until_converged(
    const std::string& method, const std::string& construction,
    const tgi::SparseMatrix& prolongation, const tgi::TwoGridCycle& cycle,
    const tgi::Vector& rhs, double tolerance) {
    tgi::Vector solution(rhs.size(), 0.0);
    tgi::Vector residual = rhs;
    tgi::TwoGridCycle::Workspace workspace;
    const double initial_norm = tgi::norm2(residual);

    Measurement result;
    result.method = method;
    result.construction = construction;
    result.density =
        experiment_support::interpolation_density_percent(prolongation);
    result.energy = cycle.setup_report().interpolation_energy;
    result.coarse_nnz = cycle.coarse_matrix().nnz();

    while (result.relative_residual > tolerance) {
        const double residual_squared =
            cycle.iterate(rhs, solution, residual, workspace);
        ++result.cycles;
        if (!(residual_squared >= 0.0) ||
            !std::isfinite(residual_squared)) {
            throw std::runtime_error(
                method + " produced a nonfinite residual");
        }
        result.relative_residual =
            std::sqrt(residual_squared) / initial_norm;
    }
    if (result.relative_residual > 0.0) {
        result.effective_factor = std::exp(
            std::log(result.relative_residual) /
            static_cast<double>(result.cycles));
    } else {
        result.effective_factor = 0.0;
    }
    return result;
}

experiment_support::Row measurement_row(const Measurement& value) {
    return {
        value.method,
        value.construction,
        std::to_string(value.cycles),
        experiment_support::scientific(value.relative_residual, 8),
        experiment_support::fixed(value.effective_factor, 9),
        experiment_support::fixed(value.density, 6),
        experiment_support::scientific(value.energy, 8),
        std::to_string(value.coarse_nnz)};
}

}

int main(int argc, char** argv) {
    int threads = 4;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument.rfind("--threads=", 0) == 0) {
            threads = std::stoi(argument.substr(10));
        }
    }

    experiment_support::BasicConfig config;
    config.fine_intervals = 128;
    config.coarse_intervals = 16;
    config.contrast = 1.0e4;
    config.threads = threads;
    const experiment_support::FieldCase field{
        "cross-channel", tgi::CoefficientDistribution::ChannelizedBinary};
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);
    const auto problem = experiment_support::make_problem(
        grid, field, config, 1);

    experiment_support::progress("geometric interpolation comparison");
    const auto geometric = tgi::build_geometric_interpolation(grid);
    const tgi::TwoGridCycle geometric_cycle(
        problem.matrix, geometric.prolongation, 1, threads);
    const Measurement geometric_measurement = measure_until_converged(
        "geometric", "bilinear", geometric.prolongation,
        geometric_cycle, problem.rhs, 1.0e-6);

    experiment_support::progress("energy-minimizing interpolation comparison");
    const auto energy = experiment_support::build_global_reference(
        grid, problem.matrix, threads);
    const tgi::TwoGridCycle energy_cycle(
        problem.matrix, energy.prolongation, 1, threads);
    const Measurement energy_measurement = measure_until_converged(
        "energy-minimizing", "column relres=1e-10",
        energy.prolongation, energy_cycle, problem.rhs, 1.0e-6);

    const double cycle_ratio =
        static_cast<double>(geometric_measurement.cycles) /
        static_cast<double>(energy_measurement.cycles);
    experiment_support::Report report(
        "Geometric and energy-minimizing interpolation on the center problem");
    report.add_summary({
        {"Grid 1/h, 1/H", "128, 16"},
        {"Contrast", "1e4"},
        {"Topology", "cross-channel"},
        {"Coefficient seed", "1"},
        {"RHS", "constant"},
        {"Threads", std::to_string(threads)},
        {"Solve tolerance", "1e-6"},
        {"Measurements", "one deterministic solve per method"},
        {"Geometric/energy cycle ratio",
         experiment_support::fixed(cycle_ratio, 3)}});
    report.add_note(
        "The matrix, coarse nodes, Galerkin construction, one forward and "
        "one backward Gauss--Seidel sweep, zero initial guess, RHS and "
        "stopping tolerance are identical. Only the interpolation matrix "
        "changes. Each solve continues until the tolerance is attained.");
    report.add_table(
        "Single-run comparison",
        {"Interpolation", "Construction", "Cycles", "Final relres",
         "Effective factor", "P density %", "J(P)", "Ac nnz"},
        {19, 22, 10, 16, 16, 12, 16, 10},
        {measurement_row(geometric_measurement),
         measurement_row(energy_measurement)});
    report.save("experiment6_endpoint_comparison");
    return 0;
}
