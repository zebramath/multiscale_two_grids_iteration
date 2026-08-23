#include "experiment/study.hpp"
#include "multigrid/global_pcg.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

struct ValidationCase {
    int fine;
    int coarse;
    double contrast;
    std::uint64_t seed;
    experiment_support::FieldCase field;
};

int reported_cycles(
    const experiment_support::CycleMetrics& metric, int maximum) {
    return metric.converged ? metric.cycles : maximum + 1;
}

}

int run_diagnostic_cases(int argc, char** argv) {
    int threads = 4;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument.rfind("--threads=", 0) == 0)
            threads = std::stoi(argument.substr(10));
    }

    const auto& channels = experiment_support::channel_topologies();
    const auto& standard = experiment_support::standard_fields();
    const std::array<ValidationCase, 8> cases{{
        {32, 8, 1.0e4, 101, channels[3]},
        {32, 8, 1.0e6, 313, channels[1]},
        {64, 8, 1.0e2, 101, channels[0]},
        {64, 8, 1.0e4, 313, channels[1]},
        {64, 16, 1.0e4, 101, channels[2]},
        {64, 16, 1.0e6, 313, channels[3]},
        {64, 8, 1.0e4, 701, standard[0]},
        {64, 16, 1.0e6, 907, standard[2]}
    }};
    constexpr int maximum_cycles = 4000;
    const experiment_support::Row headers{
        "1/h", "1/H", "Contrast", "Seed", "Topology",
        "Selected m", "Adaptive cycles", "Adaptive setup ms",
        "Oracle m", "Oracle cycles", "Gap %", "Geometric cycles",
        "Fixed-40 cycles", "Forecast"};
    experiment_support::Rows rows;

    int case_index = 0;
    for (const ValidationCase& item : cases) {
        ++case_index;
        experiment_support::progress(
            "diagnostic case " + std::to_string(case_index) + "/" +
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
        const auto geometric_metric = experiment_support::evaluate_two_grid(
            problem.matrix, problem.rhs, geometric.prolongation,
            threads, 1.0e-6, maximum_cycles);
        const int geometric_cycles = reported_cycles(
            geometric_metric, maximum_cycles);

        int oracle_steps = 0;
        int oracle_cycles = geometric_cycles;
        int fixed_cycles = maximum_cycles + 1;
        tgi::GlobalEnergyPcgPath oracle_path(
            grid, problem.matrix, geometric.prolongation, threads);
        for (int steps = 12; steps <= 48; steps += 4) {
            oracle_path.advance_to(steps);
            const tgi::SparseMatrix candidate = oracle_path.prolongation(0.0);
            const auto metric = experiment_support::evaluate_two_grid(
                problem.matrix, problem.rhs, candidate,
                threads, 1.0e-6, maximum_cycles);
            const int cycles = reported_cycles(metric, maximum_cycles);
            if (steps == 40) fixed_cycles = cycles;
            if (cycles < oracle_cycles) {
                oracle_cycles = cycles;
                oracle_steps = steps;
            }
        }

        tgi::AdaptiveGlobalPcgOptions options;
        options.minimum_steps = 12;
        options.maximum_steps = 56;
        options.maximum_cycles = maximum_cycles;
        options.thread_count = threads;
        const auto adaptive = tgi::build_adaptive_global_pcg_interpolation(
            grid, problem.matrix, geometric.prolongation,
            options, &problem.rhs);
        const auto adaptive_metric = experiment_support::evaluate_two_grid(
            problem.matrix, problem.rhs, adaptive.prolongation,
            threads, 1.0e-6, maximum_cycles);
        const int adaptive_cycles = reported_cycles(
            adaptive_metric, maximum_cycles);
        const double adaptive_setup = geometric.report.timing.total_ms +
            adaptive.report.selection_wall_ms +
            adaptive_metric.coarse_setup_ms;
        const double gap = oracle_cycles > 0
            ? 100.0 * static_cast<double>(adaptive_cycles - oracle_cycles) /
                static_cast<double>(oracle_cycles)
            : 0.0;
        rows.push_back({
            std::to_string(item.fine), std::to_string(item.coarse),
            experiment_support::scientific(item.contrast, 0),
            std::to_string(item.seed), item.field.name,
            std::to_string(adaptive.report.selected_steps),
            adaptive_cycles <= maximum_cycles
                ? std::to_string(adaptive_cycles) : "failed",
            experiment_support::fixed(adaptive_setup),
            std::to_string(oracle_steps), std::to_string(oracle_cycles),
            experiment_support::fixed(gap, 2),
            geometric_cycles <= maximum_cycles
                ? std::to_string(geometric_cycles) : "failed",
            fixed_cycles <= maximum_cycles
                ? std::to_string(fixed_cycles) : "failed",
            std::to_string(adaptive.report.estimated_selected_cycles)});
    }

    experiment_support::Report report(
        "Diagnostic extension on additional seeds and coefficient fields");
    report.add_summary({
        {"Version", std::string(tgi::version)},
        {"Cases", std::to_string(cases.size())},
        {"Seeds", "101, 313, 701, 907"},
        {"Oracle", "m=0 and m=12,16,...,48"},
        {"Solve tolerance", "1e-6"}});
    report.add_note(
        "This diagnostic set is rerun with the midpoint adaptive selector. "
        "Candidate locations use the common search interval and pilot cycle "
        "forecasts rather than case-specific step values. The oracle is "
        "evaluation-only.");
    report.add_table(
        "Post-diagnosis validation", headers,
        {5, 5, 10, 6, 20, 10, 16, 18, 9, 14, 8, 16, 15, 9}, rows);
    report.save("diagnostic_cases");
    return 0;
}
