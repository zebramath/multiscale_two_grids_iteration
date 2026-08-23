#include "experiment/study.hpp"
#include "multigrid/global_pcg_path.hpp"

#include <cmath>
#include <string>

int main(int argc, char** argv) {
    const auto config = experiment_support::parse_config(argc, argv);
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);
    const auto problem = experiment_support::make_problem(
        grid, experiment_support::channel_topologies().front(), config, 1);
    const auto geometric = experiment_support::geometric_interpolation(
        grid, problem.matrix);
    tgi::GlobalEnergyPcgPath path(
        grid, problem.matrix, geometric.prolongation, config.threads);
    const experiment_support::Row headers{
        "m", "r32", "r64", "r96", "r128", "r160", "cycles"};
    experiment_support::Rows rows;
    for (int steps = 28; steps <= 48; steps += 2) {
        path.advance_to(steps);
        const tgi::SparseMatrix p = path.prolongation(0.0);
        const tgi::TwoGridCycle cycle(
            problem.matrix, p, 1, config.threads);
        tgi::Vector solution(problem.rhs.size(), 0.0);
        tgi::Vector residual = problem.rhs;
        tgi::TwoGridCycle::Workspace workspace;
        const double initial = tgi::norm2(residual);
        double values[5]{};
        int next = 0;
        const int targets[5]{32, 64, 96, 128, 160};
        for (int iteration = 1; iteration <= 160; ++iteration) {
            const double squared = cycle.iterate(
                problem.rhs, solution, residual, workspace);
            if (iteration == targets[next]) {
                values[next] = std::sqrt(squared) / initial;
                ++next;
                if (next == 5) break;
            }
        }
        const auto solved = tgi::solve_two_grid(
            problem.matrix, problem.rhs, cycle, 1.0e-6,
            config.max_cycles);
        rows.push_back({
            std::to_string(steps),
            experiment_support::scientific(values[0], 4),
            experiment_support::scientific(values[1], 4),
            experiment_support::scientific(values[2], 4),
            experiment_support::scientific(values[3], 4),
            experiment_support::scientific(values[4], 4),
            solved.converged ? std::to_string(solved.cycles)
                : "failed@" + std::to_string(solved.cycles)});
    }
    experiment_support::Report report(
        "Pilot residual ranking diagnostic");
    report.add_summary(experiment_support::fixed_study_summary(config));
    report.add_table(
        "Finite PCG path", headers,
        {5, 13, 13, 13, 13, 13, 10}, rows);
    report.save("experiment10");
    experiment_support::write_csv("experiment10", headers, rows);
    return 0;
}
