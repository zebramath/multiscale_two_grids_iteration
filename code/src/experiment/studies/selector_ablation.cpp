#include "experiment/study.hpp"
#include "multigrid/global_pcg.hpp"

#include <array>
#include <cmath>
#include <string>

namespace {

int measured_cycles(
    const tgi::Vector& rhs, const tgi::TwoGridCycle& cycle,
    int maximum_cycles) {
    const auto solved = tgi::solve_two_grid(
        rhs, cycle, 1.0e-6, maximum_cycles);
    return solved.converged ? solved.cycles : maximum_cycles + 1;
}

int measured_cycles(
    const tgi::SparseMatrix& a, const tgi::Vector& rhs,
    const tgi::SparseMatrix& p, int threads, int maximum_cycles) {
    const tgi::TwoGridCycle cycle(a, p, 1, threads);
    return measured_cycles(rhs, cycle, maximum_cycles);
}

}

int run_selector_ablation(int argc, char** argv) {
    const auto config = experiment_support::parse_config(argc, argv);
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);
    const auto problem = experiment_support::make_problem(
        grid, experiment_support::channel_topologies().front(), config, 1);
    const auto geometric = experiment_support::geometric_interpolation(
        grid, problem.matrix);

    int oracle_steps = 0;
    int oracle_cycles = measured_cycles(
        problem.matrix, problem.rhs, geometric.prolongation,
        config.threads, config.max_cycles);
    tgi::GlobalEnergyPcgPath oracle_path(
        grid, problem.matrix, geometric.prolongation, config.threads);
    for (int steps = 16; steps <= 64; steps += 2) {
        oracle_path.advance_to(steps);
        const int cycles = measured_cycles(
            problem.matrix, problem.rhs, oracle_path.prolongation(0.0),
            config.threads, config.max_cycles);
        if (cycles < oracle_cycles) {
            oracle_steps = steps;
            oracle_cycles = cycles;
        }
    }

    constexpr std::array<int, 3> pilot_lengths{12, 24, 32};
    constexpr std::array<double, 3> slacks{0.0, 0.13, 0.20};
    const experiment_support::Row headers{
        "Pilot", "Slack %", "Selected m", "Candidates", "Forecast",
        "Cycles", "Oracle m", "Oracle cycles", "Gap %",
        "Selection ms", "P density %"};
    experiment_support::Rows rows;
    for (int pilot : pilot_lengths) {
        for (double slack : slacks) {
            tgi::AdaptiveGlobalPcgOptions options;
            options.minimum_steps = 16;
            options.maximum_steps = 64;
            options.pilot_iterations = pilot;
            options.acceptable_cycle_slack = slack;
            options.maximum_cycles = config.max_cycles;
            options.thread_count = config.threads;
            const auto adaptive = tgi::build_adaptive_global_pcg_interpolation(
                grid, problem.matrix, geometric.prolongation,
                options, &problem.rhs);
            const int cycles = measured_cycles(
                problem.rhs, *adaptive.cycle, config.max_cycles);
            const double gap = oracle_cycles > 0
                ? 100.0 * static_cast<double>(cycles - oracle_cycles) /
                    static_cast<double>(oracle_cycles)
                : 0.0;
            rows.push_back({
                std::to_string(pilot),
                experiment_support::fixed(100.0 * slack, 0),
                std::to_string(adaptive.report.selected_steps),
                std::to_string(adaptive.report.history.size()),
                std::to_string(adaptive.report.estimated_selected_cycles),
                cycles <= config.max_cycles
                    ? std::to_string(cycles)
                    : "failed@" + std::to_string(config.max_cycles),
                std::to_string(oracle_steps),
                std::to_string(oracle_cycles),
                experiment_support::fixed(gap, 2),
                experiment_support::fixed(
                    adaptive.report.selection_wall_ms),
                experiment_support::fixed(
                    experiment_support::interpolation_density_percent(
                        *adaptive.prolongation), 4)});
        }
    }

    experiment_support::Report report(
        "Initial pilot-length and selection-slack ablation");
    report.add_summary(experiment_support::fixed_study_summary(
        config, "Oracle candidates", "m=0 and m=16,18,...,64"));
    report.add_note(
        "This v3.7 table is a deliberately small main-problem ablation. "
        "It varies only pilot length and final selection slack while keeping "
        "the candidate interval, tail window and energy-residual direction "
        "rule fixed. It is an initial sensitivity check, not a frozen "
        "training/validation study.");
    report.add_table(
        "Pilot and slack sensitivity", headers,
        {7, 8, 10, 10, 10, 10, 9, 13, 8, 12, 11}, rows);
    report.save("selector_ablation");
    return 0;
}
