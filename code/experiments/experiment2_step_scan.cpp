#include "experiment/problem.hpp"
#include "experiment/reporting.hpp"
#include "multigrid/global_pcg.hpp"
#include "version.hpp"

#include <limits>
#include <string>

namespace {

double convergence_factor(
    const tgi::SparseMatrix& matrix, const tgi::Vector& rhs,
    const tgi::SparseMatrix& prolongation, int threads,
    int observation_limit) {
    const tgi::TwoGridCycle cycle(matrix, prolongation, 1, threads);
    return tgi::solve_two_grid(
        rhs, cycle, 1.0e-6, observation_limit).effective_factor;
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
    constexpr int observation_limit = 12000;
    const auto& field = experiment_support::channel_topologies().front();
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);
    const auto problem = experiment_support::make_problem(grid, field, config);
    const auto initial = tgi::build_geometric_interpolation(grid);
    tgi::GlobalEnergyPcgPath path(
        grid, problem.matrix, initial.prolongation, threads);

    const int adaptive_steps =
        tgi::adaptive_global_pcg_detail::select_steps(grid, problem.matrix);
    int best_steps = 0;
    double best_factor = std::numeric_limits<double>::infinity();
    double adaptive_factor = 1.0;
    experiment_support::Rows rows;
    for (int steps = 1; steps <= config.fine_intervals; ++steps) {
        experiment_support::progress(
            "convergence-factor scan " + std::to_string(steps) + "/" +
            std::to_string(config.fine_intervals));
        path.advance_to(steps);
        const double factor = convergence_factor(
            problem.matrix, problem.rhs, path.prolongation(),
            threads, observation_limit);
        rows.push_back({
            std::to_string(steps), experiment_support::fixed(factor, 7)});
        if (factor < best_factor) {
            best_steps = steps;
            best_factor = factor;
        }
        if (steps == adaptive_steps) adaptive_factor = factor;
    }

    const auto reference = experiment_support::build_global_reference(
        grid, problem.matrix, threads);
    const double reference_factor = convergence_factor(
        problem.matrix, problem.rhs, reference.prolongation,
        threads, experiment_support::maximum_two_grid_cycles);

    experiment_support::save_csv(
        "experiment2_central_step_scan",
        {"m", "Effective factor"}, rows);

    experiment_support::Report report(
        "Central finite-PCG convergence-factor scan");
    report.add_summary({
        {"Version", std::string(tgi::version)},
        {"Problem", "128/16 cross-channel, contrast 1e4"},
        {"Scanned steps", "m=1,...,128"},
        {"Threads", std::to_string(threads)},
        {"Solve tolerance", "1e-6"},
        {"Minimum-factor m", std::to_string(best_steps)},
        {"Minimum effective factor",
         experiment_support::fixed(best_factor, 7)},
        {"Adaptive m", std::to_string(adaptive_steps)},
        {"Adaptive effective factor",
         experiment_support::fixed(adaptive_factor, 7)},
        {"Global-reference effective factor",
         experiment_support::fixed(reference_factor, 7)}});
    report.add_note(
        "Every integer finite-PCG step count is evaluated on the same "
        "central problem. The CSV and figure contain only the effective "
        "convergence factor, exposing its nonmonotone dependence on m.");
    report.save("experiment2_step_scan");
    return 0;
}
