#include "experiment/comparison_cases.hpp"
#include "experiment/reporting.hpp"
#include "multigrid/global_pcg.hpp"
#include "version.hpp"

#include <algorithm>
#include <chrono>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace {

struct Measurement {
    std::string method;
    std::string parameter;
    double density = 0.0;
    std::size_t coarse_nnz = 0;
    double setup_ms = 0.0;
    double solve_ms = 0.0;
    tgi::TwoGridIterationResult solved;
};

struct Aggregate {
    int cases = 0;
    int converged = 0;
    int slow = 0;
    int diverged = 0;
    long long converged_cycles = 0;
    double density = 0.0;
    double setup_ms = 0.0;
    int comparable_cases = 0;
    double comparable_setup_ms = 0.0;
    double comparable_solve_ms = 0.0;
};

Measurement measure(
    const std::string& method, const std::string& parameter,
    const tgi::SparseMatrix& prolongation, const tgi::TwoGridCycle& cycle,
    const tgi::Vector& rhs, double setup_ms) {
    const auto begin = std::chrono::steady_clock::now();
    auto solved = tgi::solve_two_grid(
        rhs, cycle, 1.0e-6,
        experiment_support::maximum_two_grid_cycles);
    const double solve_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - begin).count();
    return {
        method, parameter,
        experiment_support::interpolation_density_percent(prolongation),
        cycle.setup_report().coarse_nnz, setup_ms, solve_ms,
        std::move(solved)};
}

experiment_support::Row measurement_row(
    const experiment_support::ComparisonCase& item,
    const Measurement& measurement) {
    const bool timing_valid = measurement.solved.converged;
    return {
        item.axis,
        std::to_string(item.fine), std::to_string(item.coarse),
        experiment_support::scientific(item.contrast, 0),
        std::to_string(item.seed), item.field.name,
        measurement.method, measurement.parameter,
        experiment_support::fixed(measurement.density, 4),
        std::to_string(measurement.coarse_nnz),
        experiment_support::fixed(measurement.setup_ms),
        timing_valid ? experiment_support::fixed(measurement.solve_ms) : "--",
        timing_valid
            ? experiment_support::fixed(
                measurement.setup_ms + measurement.solve_ms)
            : "--",
        std::to_string(measurement.solved.cycles),
        tgi::two_grid_status_name(measurement.solved.status),
        experiment_support::scientific(
            measurement.solved.relative_residual, 2),
        experiment_support::fixed(measurement.solved.tail_factor, 6)};
}

void accumulate(
    Aggregate& aggregate, const Measurement& measurement,
    bool jointly_converged) {
    ++aggregate.cases;
    aggregate.density += measurement.density;
    aggregate.setup_ms += measurement.setup_ms;
    if (measurement.solved.converged) {
        ++aggregate.converged;
        aggregate.converged_cycles += measurement.solved.cycles;
    } else if (measurement.solved.status ==
               tgi::TwoGridIterationStatus::Diverged) {
        ++aggregate.diverged;
    } else {
        ++aggregate.slow;
    }
    if (jointly_converged) {
        ++aggregate.comparable_cases;
        aggregate.comparable_setup_ms += measurement.setup_ms;
        aggregate.comparable_solve_ms += measurement.solve_ms;
    }
}

}

int main(int argc, char** argv) {
    bool quick = false;
    int threads = 4;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--quick") quick = true;
        if (argument.rfind("--threads=", 0) == 0) {
            threads = std::stoi(argument.substr(10));
        }
    }
    const auto cases = experiment_support::comparison_cases(quick);
    experiment_support::Rows rows;
    std::map<std::string, Aggregate> aggregates;

    int case_number = 0;
    for (const auto& item : cases) {
        ++case_number;
        experiment_support::progress(
            "two-grid comparison " + std::to_string(case_number) + "/" +
            std::to_string(cases.size()) + ": " + item.field.name);
        const auto config = experiment_support::comparison_config(
            item, threads);
        const tgi::StructuredGrid grid = experiment_support::make_grid(config);
        const auto problem = experiment_support::make_problem(
            grid, item.field, config, item.seed);

        const auto geometric =
            experiment_support::geometric_interpolation(grid);
        const double geometric_ms = geometric.report.timing.total_ms;
        std::vector<Measurement> case_measurements;

        for (const auto& policy : {
                 std::pair<const char*, double>{"adaptive-fast", 1.0},
                 std::pair<const char*, double>{"adaptive-reuse", 256.0}}) {
            tgi::AdaptiveGlobalPcgOptions adaptive_options;
            adaptive_options.minimum_steps = 12;
            adaptive_options.maximum_steps = 60;
            adaptive_options.expected_rhs_count = policy.second;
            adaptive_options.maximum_cycles =
                experiment_support::maximum_two_grid_cycles;
            adaptive_options.thread_count = threads;
            const auto adaptive =
                tgi::build_adaptive_global_pcg_interpolation(
                    grid, problem.matrix, geometric.prolongation,
                    adaptive_options, problem.rhs);
            case_measurements.push_back(measure(
                policy.first,
                "R=" + experiment_support::fixed(policy.second, 0) +
                    ",m=" +
                    std::to_string(adaptive.report.selected_steps),
                *adaptive.prolongation, *adaptive.cycle, problem.rhs,
                geometric_ms + adaptive.report.selection_wall_ms));
        }

        const auto reference = experiment_support::build_global_reference(
            grid, problem.matrix, threads);
        const tgi::TwoGridCycle reference_cycle(
            problem.matrix, reference.prolongation, 1, threads);
        case_measurements.push_back(measure(
            "global-reference", "tol=1e-10", reference.prolongation,
            reference_cycle, problem.rhs,
            reference.report.timing.total_ms +
                reference_cycle.setup_report().total_ms));

        const tgi::TwoGridCycle geometric_cycle(
            problem.matrix, geometric.prolongation, 1, threads);
        case_measurements.push_back(measure(
            "geometric", "P_G", geometric.prolongation,
            geometric_cycle, problem.rhs,
            geometric_ms + geometric_cycle.setup_report().total_ms));

        const bool jointly_converged = std::all_of(
            case_measurements.begin(), case_measurements.end(),
            [](const Measurement& value) { return value.solved.converged; });
        for (const auto& measurement : case_measurements) {
            rows.push_back(measurement_row(item, measurement));
            accumulate(
                aggregates[measurement.method], measurement,
                jointly_converged);
        }
    }

    experiment_support::Rows convergence_rows;
    experiment_support::Rows timing_rows;
    for (const std::string method :
         {"adaptive-fast", "adaptive-reuse",
          "global-reference", "geometric"}) {
        const Aggregate& value = aggregates[method];
        convergence_rows.push_back({
            method,
            std::to_string(value.converged) + "/" +
                std::to_string(value.slow) + "/" +
                std::to_string(value.diverged),
            std::to_string(value.converged_cycles),
            experiment_support::fixed(
                value.density / static_cast<double>(value.cases), 4),
            experiment_support::fixed(value.setup_ms)});
        timing_rows.push_back({
            method,
            std::to_string(value.comparable_cases),
            experiment_support::fixed(value.comparable_setup_ms),
            experiment_support::fixed(value.comparable_solve_ms),
            experiment_support::fixed(
                value.comparable_setup_ms + value.comparable_solve_ms),
            experiment_support::fixed(
                value.comparable_setup_ms +
                256.0 * value.comparable_solve_ms)});
    }

    experiment_support::Report report(
        "Cross-problem two-grid interpolation comparison");
    report.add_summary({
        {"Version", std::string(tgi::version)},
        {"Cases", std::to_string(cases.size())},
        {"Mode", quick ? "quick" : "full"},
        {"Threads", std::to_string(threads)},
        {"Solve tolerance", "1e-6"},
        {"Maximum cycles",
         std::to_string(experiment_support::maximum_two_grid_cycles)}});
    report.add_note(
        "The v4.7 matrix varies fine/coarse scale, contrast and six channel "
        "topologies. Two fixed-H scale extensions, 160/16 and 192/16, test "
        "larger fine-grid separation. Fast evaluates m=0,12,32,52 with a "
        "16-cycle pilot cap; reuse evaluates m=0,12,20,...,60 with a "
        "160-cycle cap and at most two step-two refinements. A slow-limit "
        "row still contracts at the 12000-cycle cap; a diverged row has "
        "nonfinite residual or sustained growing tail. Solve and total times "
        "are shown only for converged rows.");
    report.add_table(
        "All two-grid cases",
        {"Axis", "1/h", "1/H", "Contrast", "Seed", "Topology", "Method",
         "Parameter", "P density %", "Ac nnz", "Setup ms", "Solve ms",
         "Total ms", "Cycles", "Status", "Final relres", "Tail factor"},
        {12, 5, 5, 10, 6, 20, 16, 12, 11, 9, 10, 10, 10, 8, 10, 12, 11},
        rows, true);
    report.add_table(
        "Convergence and setup summary",
        {"Method", "Conv/slow/div", "Converged cycle sum",
         "Mean density %", "Setup sum ms"},
        {16, 13, 19, 14, 13}, convergence_rows);
    report.add_table(
        "Timing on the jointly converged subset",
        {"Method", "Cases", "Setup sum ms", "Solve sum ms",
         "Total R=1 ms", "Total R=256 ms"},
        {16, 7, 13, 12, 13, 14}, timing_rows);
    report.save("experiment1_two_grid_comparison");
    return 0;
}
