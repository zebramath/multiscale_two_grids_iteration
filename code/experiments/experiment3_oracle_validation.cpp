#include "experiment/problem.hpp"
#include "experiment/reporting.hpp"
#include "multigrid/global_pcg.hpp"
#include "version.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace {

struct OracleCase {
    int fine;
    int coarse;
    double contrast;
    std::uint64_t seed;
    experiment_support::FieldCase field;
};

struct CycleMeasurement {
    int cycles = 0;
    tgi::TwoGridIterationStatus status =
        tgi::TwoGridIterationStatus::SlowAtLimit;
    bool converged = false;
};

CycleMeasurement cycles_for(
    const tgi::SparseMatrix& a, const tgi::Vector& rhs,
    const tgi::SparseMatrix& p, int threads, int maximum) {
    const tgi::TwoGridCycle cycle(a, p, 1, threads);
    const auto solved = tgi::solve_two_grid(
        rhs, cycle, 1.0e-6, maximum);
    return {solved.cycles, solved.status, solved.converged};
}

CycleMeasurement cycles_for(
    const tgi::Vector& rhs, const tgi::TwoGridCycle& cycle, int maximum) {
    const auto solved = tgi::solve_two_grid(
        rhs, cycle, 1.0e-6, maximum);
    return {solved.cycles, solved.status, solved.converged};
}

}

int main(int argc, char** argv) {
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
    constexpr int maximum_cycles =
        experiment_support::maximum_two_grid_cycles;
    const experiment_support::Row headers{
        "1/h", "1/H", "Contrast", "Topology", "Policy", "R",
        "Selected m", "Cycles", "Status", "Oracle m", "Oracle cycles",
        "Oracle status", "Gap %",
        "Selection ms", "Candidates", "Pilot cap", "Stride", "Max m",
        "Refine"};
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
        const tgi::StructuredGrid grid = experiment_support::make_grid(config);
        const auto problem = experiment_support::make_problem(
            grid, item.field, config, item.seed);
        const auto geometric =
            experiment_support::geometric_interpolation(grid);

        int oracle_steps = 0;
        CycleMeasurement oracle = cycles_for(
            problem.matrix, problem.rhs, geometric.prolongation,
            threads, maximum_cycles);
        tgi::GlobalEnergyPcgPath path(
            grid, problem.matrix, geometric.prolongation, threads);
        for (int steps = 12; steps <= 60; steps += 2) {
            path.advance_to(steps);
            const tgi::SparseMatrix candidate = path.prolongation(0.0);
            const CycleMeasurement measured = cycles_for(
                problem.matrix, problem.rhs, candidate,
                threads, maximum_cycles);
            if (measured.converged &&
                (!oracle.converged || measured.cycles < oracle.cycles)) {
                oracle = measured;
                oracle_steps = steps;
            }
        }

        for (const auto& policy : {
                 std::pair<const char*, double>{"fast", 1.0},
                 std::pair<const char*, double>{"reuse", 256.0}}) {
            tgi::AdaptiveGlobalPcgOptions options;
            options.minimum_steps = 12;
            options.maximum_steps = 60;
            options.expected_rhs_count = policy.second;
            options.maximum_cycles = maximum_cycles;
            options.thread_count = threads;
            const auto adaptive =
                tgi::build_adaptive_global_pcg_interpolation(
                    grid, problem.matrix, geometric.prolongation,
                    options, problem.rhs);
            const CycleMeasurement adaptive_cycles = cycles_for(
                problem.rhs, *adaptive.cycle, maximum_cycles);
            const double gap = adaptive_cycles.converged && oracle.converged
                ? 100.0 * static_cast<double>(
                      adaptive_cycles.cycles - oracle.cycles) /
                      static_cast<double>(oracle.cycles)
                : 0.0;
            rows.push_back({
                std::to_string(item.fine), std::to_string(item.coarse),
                experiment_support::scientific(item.contrast, 0),
                item.field.name, policy.first,
                experiment_support::fixed(policy.second, 0),
                std::to_string(adaptive.report.selected_steps),
                std::to_string(adaptive_cycles.cycles),
                tgi::two_grid_status_name(adaptive_cycles.status),
                std::to_string(oracle_steps),
                std::to_string(oracle.cycles),
                tgi::two_grid_status_name(oracle.status),
                adaptive_cycles.converged && oracle.converged
                    ? experiment_support::fixed(gap, 2)
                    : "--",
                experiment_support::fixed(
                    adaptive.report.selection_wall_ms),
                std::to_string(adaptive.report.candidate_count),
                std::to_string(adaptive.report.pilot_limit),
                std::to_string(adaptive.report.checkpoint_stride),
                std::to_string(adaptive.report.maximum_sampled_steps),
                adaptive.report.used_local_refinement ? "yes" : "no"});
        }
    }

    experiment_support::Report report(
        "Adaptive PCG versus an offline step-two oracle");
    report.add_summary({
        {"Version", std::string(tgi::version)},
        {"Threads", std::to_string(threads)},
        {"Oracle candidates", "m=0 and m=12,14,...,60"},
        {"Solve tolerance", "1e-6"},
        {"Maximum cycles", std::to_string(maximum_cycles)}});
    report.add_note(
        "The step-two oracle is evaluation-only. Fast samples "
        "m=0,12,32,52 with a 16-cycle pilot cap, a 10% near-optimality slack "
        "and no refinement. Reuse samples m=0,12,20,...,60 with a "
        "160-cycle pilot cap, a 2% slack and at most two neighboring step-two "
        "refinements. "
        "Neither policy reads scale, contrast or topology labels.");
    report.add_table(
        "Representative oracle gaps", headers,
        {5, 5, 10, 20, 7, 4, 10, 8, 10, 9, 14, 13, 8, 13, 10, 9, 7, 7, 7}, rows,
        true);
    report.save("experiment3_oracle_validation");
    return 0;
}
