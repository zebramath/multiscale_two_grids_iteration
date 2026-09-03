#include "experiment/problem.hpp"
#include "experiment/reporting.hpp"
#include "multigrid/global_pcg.hpp"
#include "version.hpp"

#include <limits>
#include <string>

namespace {

struct Measurement {
    int steps;
    double density;
    tgi::StationaryIterationResult solved;
};

Measurement measure(
    int steps, const tgi::SparseMatrix& matrix, const tgi::Vector& rhs,
    const tgi::SparseMatrix& prolongation, int threads, int maximum_cycles) {
    const tgi::TwoGridCycle cycle(matrix, prolongation, 1, threads);
    return {
        steps,
        experiment_support::interpolation_density_percent(prolongation),
        tgi::solve_two_grid(rhs, cycle, 1.0e-6, maximum_cycles)};
}

experiment_support::Row result_row(const Measurement& value) {
    return {
        std::to_string(value.steps),
        experiment_support::fixed(value.density, 4),
        std::to_string(value.solved.cycles),
        tgi::stationary_status_name(value.solved.status),
        experiment_support::scientific(value.solved.relative_residual, 3),
        experiment_support::fixed(value.solved.effective_factor, 7),
        experiment_support::fixed(value.solved.tail_factor, 7)};
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
    constexpr int maximum_cycles = 12000;
    const auto& field = experiment_support::channel_topologies().front();
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);
    const auto problem = experiment_support::make_problem(grid, field, config);
    const auto initial = tgi::build_geometric_interpolation(grid);
    tgi::GlobalEnergyPcgPath path(
        grid, problem.matrix, initial.prolongation, threads);

    experiment_support::Rows rows;
    int best_steps = 0;
    int best_cycles = std::numeric_limits<int>::max();
    double best_factor = 1.0;
    for (int steps = 1; steps <= config.fine_intervals / 2; ++steps) {
        experiment_support::progress(
            "central step scan " + std::to_string(steps) + "/" +
            std::to_string(config.fine_intervals / 2));
        path.advance_to(steps);
        const Measurement measured = measure(
            steps, problem.matrix, problem.rhs, path.prolongation(),
            threads, maximum_cycles);
        rows.push_back(result_row(measured));
        if (measured.solved.converged &&
            measured.solved.cycles < best_cycles) {
            best_steps = steps;
            best_cycles = measured.solved.cycles;
            best_factor = measured.solved.effective_factor;
        }
    }

    const auto reference = experiment_support::build_global_reference(
        grid, problem.matrix, threads);
    const Measurement reference_measurement = measure(
        0, problem.matrix, problem.rhs, reference.prolongation,
        threads, experiment_support::maximum_two_grid_cycles);

    const experiment_support::Row headers{
        "m", "P density %", "Cycles", "Status", "Final relres",
        "Effective factor", "Tail factor"};
    experiment_support::save_csv(
        "experiment2_central_step_scan", headers, rows);

    experiment_support::Report report("Central finite-PCG step scan");
    report.add_summary({
        {"Version", std::string(tgi::version)},
        {"Problem", "128/16 cross-channel, contrast 1e4"},
        {"Scanned steps", "m=1,...,64"},
        {"Threads", std::to_string(threads)},
        {"Solve tolerance", "1e-6"},
        {"Cycle cap", std::to_string(maximum_cycles)},
        {"Best converged m", std::to_string(best_steps)},
        {"Best cycles", std::to_string(best_cycles)},
        {"Best effective factor", experiment_support::fixed(best_factor, 7)},
        {"Global-reference cycles",
         std::to_string(reference_measurement.solved.cycles)},
        {"Global-reference effective factor",
         experiment_support::fixed(
             reference_measurement.solved.effective_factor, 7)}});
    report.add_note(
        "Every integer finite-PCG step count in the stated interval is "
        "evaluated on one fixed central problem. Cycle counts at the cap are "
        "reported as slow-limit. The complete scan, rather than sparse path "
        "checkpoints, resolves the nonmonotone dependence of two-grid "
        "convergence on interpolation work.");
    report.add_table(
        "Complete step scan", headers,
        {5, 11, 8, 10, 13, 16, 12}, rows);
    report.save("experiment2_step_scan");
    return 0;
}
