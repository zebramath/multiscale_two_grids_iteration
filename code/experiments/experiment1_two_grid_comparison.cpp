#include "experiment/comparison_cases.hpp"
#include "experiment/reporting.hpp"
#include "multigrid/global_pcg.hpp"
#include "version.hpp"

#include <algorithm>
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
    tgi::StationaryIterationResult solved;
};

struct Aggregate {
    int cases = 0;
    int converged = 0;
    int slow = 0;
    int diverged = 0;
    long long converged_cycles = 0;
    double density = 0.0;
};

Measurement measure(
    const std::string& method, const std::string& parameter,
    const tgi::SparseMatrix& prolongation, const tgi::TwoGridCycle& cycle,
    const tgi::Vector& rhs, int maximum_cycles) {
    auto solved = tgi::solve_two_grid(
        rhs, cycle, 1.0e-6, maximum_cycles);
    return {
        method, parameter,
        experiment_support::interpolation_density_percent(prolongation),
        cycle.coarse_matrix().nnz(),
        std::move(solved)};
}

experiment_support::Row measurement_row(
    const experiment_support::ComparisonCase& item,
    const Measurement& measurement) {
    return {
        item.axis,
        std::to_string(item.fine), std::to_string(item.coarse),
        experiment_support::scientific(item.contrast, 0),
        item.field.name,
        measurement.method, measurement.parameter,
        experiment_support::fixed(measurement.density, 4),
        std::to_string(measurement.coarse_nnz),
        std::to_string(measurement.solved.cycles),
        tgi::stationary_status_name(measurement.solved.status),
        experiment_support::scientific(
            measurement.solved.relative_residual, 2),
        experiment_support::fixed(
            measurement.solved.effective_factor, 6),
        experiment_support::fixed(measurement.solved.tail_factor, 6)};
}

void accumulate(
    Aggregate& aggregate, const Measurement& measurement) {
    ++aggregate.cases;
    aggregate.density += measurement.density;
    if (measurement.solved.converged) {
        ++aggregate.converged;
        aggregate.converged_cycles += measurement.solved.cycles;
    } else if (measurement.solved.status ==
               tgi::StationaryIterationStatus::Diverged) {
        ++aggregate.diverged;
    } else {
        ++aggregate.slow;
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
            grid, item.field, config);

        const auto geometric = tgi::build_geometric_interpolation(grid);
        std::vector<Measurement> case_measurements;

        const auto adaptive = tgi::build_adaptive_global_pcg_interpolation(
            grid, problem.matrix, geometric.prolongation, threads);
        case_measurements.push_back(measure(
            "adaptive", "m=" + std::to_string(adaptive.report.selected_steps),
            *adaptive.prolongation, *adaptive.cycle, problem.rhs,
            experiment_support::maximum_two_grid_cycles));

        const auto reference = experiment_support::build_global_reference(
            grid, problem.matrix, threads);
        const tgi::TwoGridCycle reference_cycle(
            problem.matrix, reference.prolongation, 1, threads);
        case_measurements.push_back(measure(
            "global-reference", "tol=1e-10", reference.prolongation,
            reference_cycle, problem.rhs,
            experiment_support::maximum_two_grid_cycles));

        const tgi::TwoGridCycle geometric_cycle(
            problem.matrix, geometric.prolongation, 1, threads);
        case_measurements.push_back(measure(
            "geometric", "P_G", geometric.prolongation,
            geometric_cycle, problem.rhs,
            experiment_support::maximum_geometric_cycles));

        for (const auto& measurement : case_measurements) {
            rows.push_back(measurement_row(item, measurement));
            accumulate(aggregates[measurement.method], measurement);
        }
    }

    experiment_support::Rows convergence_rows;
    for (const std::string method :
         {"adaptive", "global-reference", "geometric"}) {
        const Aggregate& value = aggregates[method];
        convergence_rows.push_back({
            method,
            std::to_string(value.converged) + "/" +
                std::to_string(value.slow) + "/" +
                std::to_string(value.diverged),
            std::to_string(value.converged_cycles),
            experiment_support::fixed(
                value.density / static_cast<double>(value.cases), 4)});
    }

    experiment_support::Report report(
        "Cross-problem two-grid interpolation comparison");
    report.add_summary({
        {"Version", std::string(tgi::version)},
        {"Cases", std::to_string(cases.size())},
        {"Mode", quick ? "quick" : "full"},
        {"Threads", std::to_string(threads)},
        {"Solve tolerance", "1e-6"},
        {"Maximum cycles", "adaptive/reference 20000; geometric 30000"}});
    report.add_note(
        "The matrix varies fine/coarse scale, contrast and six channel "
        "topologies. The two 256/16 extensions test cross-channel and "
        "winding-ring coefficients at the largest scale. Adaptive uses (1/h)/8 "
        "when 1/H<=8; otherwise it uses (1/h)/4, (1/h)/3 or (1/h)/2 in "
        "the low, intermediate or high diagonal-ratio band. Effective factor "
        "is the whole-solve geometric mean contraction; tail factor uses the "
        "last 32 cycles.");
    report.add_table(
        "All two-grid cases",
        {"Axis", "1/h", "1/H", "Contrast", "Topology", "Method",
         "Parameter", "P density %", "Ac nnz", "Cycles", "Status",
         "Final relres", "Effective factor", "Tail factor"},
        {12, 5, 5, 10, 20, 16, 12, 11, 9, 8, 10, 12, 16, 11},
        rows, true);
    report.add_table(
        "Convergence and setup summary",
        {"Method", "Conv/slow/div", "Converged cycle sum",
         "Mean density %"},
        {16, 13, 19, 14}, convergence_rows);
    report.save("experiment1_two_grid_comparison");
    return 0;
}
