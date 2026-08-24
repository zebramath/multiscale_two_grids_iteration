#include "experiment/comparison_cases.hpp"
#include "experiment/study.hpp"
#include "multigrid/global_pcg.hpp"

#include <chrono>
#include <map>
#include <string>

namespace {

struct Aggregate {
    int cases = 0;
    int converged = 0;
    long long cycles = 0;
    double density = 0.0;
    double setup_ms = 0.0;
    double solve_ms = 0.0;
};

experiment_support::Row measurement_row(
    const experiment_support::ComparisonCase& item,
    const std::string& method, const std::string& parameter,
    const tgi::SparseMatrix& prolongation, const tgi::TwoGridCycle& cycle,
    const tgi::Vector& rhs, double setup_ms, int maximum_cycles,
    Aggregate& aggregate) {
    const auto begin = std::chrono::steady_clock::now();
    const auto solved = tgi::solve_two_grid(
        rhs, cycle, 1.0e-6, maximum_cycles);
    const double solve_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - begin).count();
    const double density =
        experiment_support::interpolation_density_percent(prolongation);
    ++aggregate.cases;
    aggregate.converged += solved.converged ? 1 : 0;
    aggregate.cycles += solved.cycles;
    aggregate.density += density;
    aggregate.setup_ms += setup_ms;
    aggregate.solve_ms += solve_ms;
    return {
        std::to_string(item.fine), std::to_string(item.coarse),
        experiment_support::scientific(item.contrast, 0),
        std::to_string(item.seed), item.field.name, method, parameter,
        experiment_support::fixed(density, 4),
        std::to_string(cycle.setup_report().coarse_nnz),
        experiment_support::fixed(setup_ms),
        experiment_support::fixed(solve_ms),
        experiment_support::fixed(setup_ms + solve_ms),
        solved.converged ? std::to_string(solved.cycles)
            : "failed@" + std::to_string(solved.cycles)};
}

}

int run_two_grid_comparison(int argc, char** argv) {
    bool quick = false;
    int threads = 4;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--quick") quick = true;
        if (argument.rfind("--threads=", 0) == 0) {
            threads = std::stoi(argument.substr(10));
        }
    }
    constexpr int maximum_cycles = 6000;
    const auto cases = experiment_support::comparison_cases(quick);
    const experiment_support::Row headers{
        "1/h", "1/H", "Contrast", "Seed", "Topology", "Method",
        "Parameter", "P density %", "Ac nnz", "Setup ms", "Solve ms",
        "Total ms", "Cycles"};
    experiment_support::Rows rows;
    std::map<std::string, Aggregate> aggregates;

    int case_number = 0;
    for (const auto& item : cases) {
        ++case_number;
        experiment_support::progress(
            "two-grid comparison " + std::to_string(case_number) + "/" +
            std::to_string(cases.size()) + ": " + item.field.name);
        const auto config = experiment_support::comparison_config(
            item, threads, maximum_cycles);
        const tgi::StructuredGrid grid = experiment_support::make_grid(config);
        const auto problem = experiment_support::make_problem(
            grid, item.field, config, item.seed);

        const auto geometric = experiment_support::geometric_interpolation(
            grid, problem.matrix);
        const double geometric_ms = geometric.report.timing.total_ms;

        tgi::AdaptiveGlobalPcgOptions adaptive_options;
        adaptive_options.minimum_steps = 12;
        adaptive_options.maximum_steps = 56;
        adaptive_options.maximum_cycles = maximum_cycles;
        adaptive_options.thread_count = threads;
        const auto adaptive = tgi::build_adaptive_global_pcg_interpolation(
            grid, problem.matrix, geometric.prolongation, adaptive_options,
            &problem.rhs);
        rows.push_back(measurement_row(
            item, "adaptive-PCG",
            "m=" + std::to_string(adaptive.report.selected_steps),
            *adaptive.prolongation, *adaptive.cycle, problem.rhs,
            geometric_ms + adaptive.report.selection_wall_ms,
            maximum_cycles, aggregates["adaptive-PCG"]));

        const auto exact = experiment_support::build_global_reference(
            grid, problem.matrix, threads);
        const tgi::TwoGridCycle exact_cycle(
            problem.matrix, exact.prolongation, 1, threads);
        rows.push_back(measurement_row(
            item, "global-exact", "tol=1e-10", exact.prolongation,
            exact_cycle, problem.rhs,
            exact.report.timing.total_ms +
                exact_cycle.setup_report().total_ms,
            maximum_cycles, aggregates["global-exact"]));

        const tgi::TwoGridCycle geometric_cycle(
            problem.matrix, geometric.prolongation, 1, threads);
        rows.push_back(measurement_row(
            item, "geometric", "P_G", geometric.prolongation,
            geometric_cycle, problem.rhs,
            geometric_ms + geometric_cycle.setup_report().total_ms,
            maximum_cycles, aggregates["geometric"]));
    }

    experiment_support::Rows summary_rows;
    for (const std::string method :
         {"adaptive-PCG", "global-exact", "geometric"}) {
        const Aggregate& value = aggregates[method];
        summary_rows.push_back({
            method,
            std::to_string(value.converged) + "/" +
                std::to_string(value.cases),
            std::to_string(value.cycles),
            experiment_support::fixed(
                value.density / static_cast<double>(value.cases), 4),
            experiment_support::fixed(value.setup_ms),
            experiment_support::fixed(value.solve_ms),
            experiment_support::fixed(value.setup_ms + value.solve_ms)});
    }

    experiment_support::Report report(
        "Cross-problem two-grid interpolation comparison");
    report.add_summary({
        {"Version", std::string(tgi::version)},
        {"Cases", std::to_string(cases.size())},
        {"Mode", quick ? "quick" : "full"},
        {"Solve tolerance", "1e-6"},
        {"Maximum cycles", std::to_string(maximum_cycles)}});
    report.add_note(
        "The primary v4.0 experiment varies fine/coarse size, contrast, seed "
        "and all six channel topologies. Every case uses the same coarse "
        "nodes and smoother. The comparison therefore isolates interpolation "
        "choice rather than comparing unrelated AMG implementations.");
    report.add_table(
        "All two-grid cases", headers,
        {5, 5, 10, 6, 20, 15, 12, 11, 9, 10, 10, 10, 12}, rows, true);
    report.add_table(
        "Aggregate comparison",
        {"Method", "Converged", "Cycle sum", "Mean density %",
         "Setup sum ms", "Solve sum ms", "Total sum ms"},
        {15, 10, 11, 14, 13, 12, 12}, summary_rows);
    report.save("two_grid_comparison");
    return 0;
}
