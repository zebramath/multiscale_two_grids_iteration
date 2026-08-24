#include "experiment/study.hpp"
#include "multigrid/global_pcg.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace {

struct OracleCase {
    int fine;
    int coarse;
    double contrast;
    std::uint64_t seed;
    experiment_support::FieldCase field;
};

int cycles_for(
    const tgi::SparseMatrix& a, const tgi::Vector& rhs,
    const tgi::SparseMatrix& p, int threads, int maximum) {
    const tgi::TwoGridCycle cycle(a, p, 1, threads);
    const auto solved = tgi::solve_two_grid(
        rhs, cycle, 1.0e-6, maximum);
    return solved.converged ? solved.cycles : maximum + 1;
}

int cycles_for(
    const tgi::Vector& rhs, const tgi::TwoGridCycle& cycle, int maximum) {
    const auto solved = tgi::solve_two_grid(
        rhs, cycle, 1.0e-6, maximum);
    return solved.converged ? solved.cycles : maximum + 1;
}

}

int run_oracle_quality(int argc, char** argv) {
    int threads = 4;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument.rfind("--threads=", 0) == 0)
            threads = std::stoi(argument.substr(10));
    }
    const auto& topology = experiment_support::channel_topologies();
    const std::array<OracleCase, 7> cases{{
        {32, 8, 1.0e4, 1, topology[0]},
        {64, 8, 1.0e4, 1, topology[0]},
        {64, 16, 1.0e4, 1, topology[0]},
        {128, 16, 1.0e4, 1, topology[0]},
        {128, 16, 1.0e2, 1, topology[0]},
        {128, 16, 1.0e6, 1, topology[0]},
        {128, 16, 1.0e4, 1, topology[5]}
    }};
    constexpr int maximum_cycles = 6000;
    const experiment_support::Row headers{
        "1/h", "1/H", "Contrast", "Topology", "Selected m",
        "Adaptive cycles", "Oracle m", "Oracle cycles", "Gap %",
        "Selection ms", "P density %"};
    experiment_support::Rows rows;

    int case_index = 0;
    for (const OracleCase& item : cases) {
        ++case_index;
        experiment_support::progress(
            "oracle case " + std::to_string(case_index) + "/" +
            std::to_string(cases.size()) + ": " + item.field.name);
        experiment_support::BasicConfig config;
        config.fine_intervals = item.fine;
        config.coarse_intervals = item.coarse;
        config.contrast = item.contrast;
        config.threads = threads;
        config.max_cycles = maximum_cycles;
        const tgi::StructuredGrid grid = experiment_support::make_grid(config);
        const auto problem = experiment_support::make_problem(
            grid, item.field, config, item.seed);
        const auto geometric = experiment_support::geometric_interpolation(
            grid, problem.matrix);

        int oracle_steps = 0;
        int oracle_cycles = cycles_for(
            problem.matrix, problem.rhs, geometric.prolongation,
            threads, maximum_cycles);
        tgi::GlobalEnergyPcgPath path(
            grid, problem.matrix, geometric.prolongation, threads);
        for (int steps = 12; steps <= 60; steps += 2) {
            path.advance_to(steps);
            const tgi::SparseMatrix candidate = path.prolongation(0.0);
            const int cycles = cycles_for(
                problem.matrix, problem.rhs, candidate,
                threads, maximum_cycles);
            if (cycles < oracle_cycles) {
                oracle_cycles = cycles;
                oracle_steps = steps;
            }
        }

        tgi::AdaptiveGlobalPcgOptions options;
        options.minimum_steps = 12;
        options.maximum_steps = 60;
        options.maximum_cycles = maximum_cycles;
        options.thread_count = threads;
        const auto adaptive = tgi::build_adaptive_global_pcg_interpolation(
            grid, problem.matrix, geometric.prolongation,
            options, &problem.rhs);
        const int adaptive_cycles = cycles_for(
            problem.rhs, *adaptive.cycle, maximum_cycles);
        const double gap = oracle_cycles > 0
            ? 100.0 * static_cast<double>(
                  adaptive_cycles - oracle_cycles) /
                  static_cast<double>(oracle_cycles)
            : 0.0;
        rows.push_back({
            std::to_string(item.fine), std::to_string(item.coarse),
            experiment_support::scientific(item.contrast, 0),
            item.field.name,
            std::to_string(adaptive.report.selected_steps),
            std::to_string(adaptive_cycles),
            std::to_string(oracle_steps), std::to_string(oracle_cycles),
            experiment_support::fixed(gap, 2),
            experiment_support::fixed(adaptive.report.selection_wall_ms),
            experiment_support::fixed(
                experiment_support::interpolation_density_percent(
                    *adaptive.prolongation), 4)});
    }

    experiment_support::Report report(
        "Adaptive PCG versus an offline step-two oracle");
    report.add_summary({
        {"Version", std::string(tgi::version)},
        {"Oracle candidates", "m=0 and m=12,14,...,60"},
        {"Solve tolerance", "1e-6"}});
    report.add_note(
        "The step-two oracle is evaluation-only. The practical selector "
        "screens m=0,12,20,...,60 for 64 cycles, confirms the three best "
        "forecasts plus the minimum-step anchor and parsimonious screen "
        "winner, and locally refines the best interval. It never reads "
        "coefficient labels.");
    report.add_table(
        "Representative oracle gaps", headers,
        {5, 5, 10, 20, 10, 16, 9, 14, 8, 13, 11}, rows);
    report.save("oracle_quality");
    return 0;
}
