#include "experiment/problem.hpp"
#include "experiment/reporting.hpp"
#include "multigrid/global_pcg.hpp"
#include "version.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
constexpr double pi = 3.141592653589793238462643383279502884;

struct SolveResult {
    int cycles = 0;
    double milliseconds = 0.0;
    tgi::TwoGridIterationStatus status =
        tgi::TwoGridIterationStatus::SlowAtLimit;
    bool converged = false;
};

struct CycleAggregate {
    int cases = 0;
    int converged = 0;
    int slow = 0;
    int diverged = 0;
    long long cycles = 0;
    int maximum_cycles = 0;
};

struct TimingSample {
    double setup_ms = 0.0;
    double solve_ms = 0.0;
    int cycles = 0;
};

struct RhsCase {
    std::string name;
    tgi::Vector values;
};

double elapsed_ms(Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

SolveResult measure_solve(
    const tgi::Vector& rhs, const tgi::TwoGridCycle& cycle,
    int maximum_cycles) {
    const auto begin = Clock::now();
    const auto solved = tgi::solve_two_grid(
        rhs, cycle, 1.0e-6, maximum_cycles);
    return {
        solved.cycles, elapsed_ms(begin, Clock::now()),
        solved.status, solved.converged};
}

tgi::AdaptiveGlobalPcgResult build_adaptive(
    const tgi::StructuredGrid& grid, const tgi::SparseMatrix& matrix,
    const tgi::SparseMatrix& geometric, const tgi::Vector& training_rhs,
    double expected_rhs_count, int threads, int maximum_cycles) {
    tgi::AdaptiveGlobalPcgOptions options;
    options.minimum_steps = 12;
    options.maximum_steps = 60;
    options.expected_rhs_count = expected_rhs_count;
    options.maximum_cycles = maximum_cycles;
    options.thread_count = threads;
    return tgi::build_adaptive_global_pcg_interpolation(
        grid, matrix, geometric, options, training_rhs);
}

void add_aggregate(
    CycleAggregate& aggregate, const SolveResult& measurement) {
    ++aggregate.cases;
    if (measurement.converged) {
        ++aggregate.converged;
        aggregate.cycles += measurement.cycles;
        aggregate.maximum_cycles = std::max(
            aggregate.maximum_cycles, measurement.cycles);
    } else if (measurement.status ==
               tgi::TwoGridIterationStatus::Diverged) {
        ++aggregate.diverged;
    } else {
        ++aggregate.slow;
    }
}

tgi::Vector normalized(tgi::Vector values) {
    const double length = tgi::norm2(values);
    for (double& value : values) value /= length;
    return values;
}

std::vector<RhsCase> make_rhs_cases(const tgi::StructuredGrid& grid) {
    const std::size_t size = static_cast<std::size_t>(grid.fine_size());
    std::vector<RhsCase> cases;
    cases.push_back({"constant", normalized(tgi::Vector(size, 1.0))});

    for (const auto& specification : {
             std::pair<std::string, int>{"smooth-low", 0},
             {"smooth-mixed", 1},
             {"localized-left", 2},
             {"localized-right", 3},
             {"oscillatory", 4}}) {
        tgi::Vector rhs(size, 0.0);
        for (int node = 0; node < grid.fine_size(); ++node) {
            const auto [ix, iy] = grid.fine_coords(node);
            const double x = static_cast<double>(ix + 1) * grid.h();
            const double y = static_cast<double>(iy + 1) * grid.h();
            double value = 0.0;
            if (specification.second == 0) {
                value = std::sin(pi * x) * std::sin(pi * y);
            } else if (specification.second == 1) {
                value = std::sin(2.0 * pi * x) *
                        std::sin(3.0 * pi * y) +
                    0.35 *
                        std::sin(5.0 * pi * x) *
                        std::sin(2.0 * pi * y);
            } else if (specification.second == 2) {
                const double dx = x - 0.28;
                const double dy = y - 0.72;
                value = std::exp(-90.0 * (dx * dx + dy * dy));
            } else if (specification.second == 3) {
                const double dx = x - 0.74;
                const double dy = y - 0.31;
                value = std::exp(-75.0 * (dx * dx + dy * dy));
            } else {
                value = std::sin(7.0 * pi * x) *
                    std::sin(6.0 * pi * y);
            }
            rhs[static_cast<std::size_t>(node)] = value;
        }
        cases.push_back({specification.first, normalized(std::move(rhs))});
    }
    return cases;
}

double quantile(std::vector<double> values, double probability) {
    std::sort(values.begin(), values.end());
    const double position = probability *
        static_cast<double>(values.size() - 1U);
    const std::size_t lower = static_cast<std::size_t>(std::floor(position));
    const std::size_t upper = static_cast<std::size_t>(std::ceil(position));
    const double weight = position - static_cast<double>(lower);
    return (1.0 - weight) * values[lower] + weight * values[upper];
}

TimingSample timing_sample(
    const tgi::StructuredGrid& grid, const tgi::SparseMatrix& matrix,
    const tgi::Vector& rhs, double expected_rhs_count, int threads,
    int maximum_cycles) {
    const auto setup_begin = Clock::now();
    const auto geometric =
        experiment_support::geometric_interpolation(grid);
    const auto adaptive = build_adaptive(
        grid, matrix, geometric.prolongation, rhs, expected_rhs_count,
        threads, maximum_cycles);
    const double setup_ms = elapsed_ms(setup_begin, Clock::now());
    const SolveResult solved = measure_solve(
        rhs, *adaptive.cycle, maximum_cycles);
    return {setup_ms, solved.milliseconds, solved.cycles};
}

experiment_support::Row timing_row(
    const std::string& policy, const std::vector<TimingSample>& samples) {
    std::vector<double> setup;
    std::vector<double> solve;
    std::vector<double> total;
    std::vector<double> cycles;
    for (const auto& sample : samples) {
        setup.push_back(sample.setup_ms);
        solve.push_back(sample.solve_ms);
        total.push_back(sample.setup_ms + sample.solve_ms);
        cycles.push_back(static_cast<double>(sample.cycles));
    }
    return {
        policy, std::to_string(samples.size()),
        experiment_support::fixed(quantile(setup, 0.50)),
        experiment_support::fixed(quantile(setup, 0.25)),
        experiment_support::fixed(quantile(setup, 0.75)),
        experiment_support::fixed(quantile(solve, 0.50)),
        experiment_support::fixed(quantile(solve, 0.25)),
        experiment_support::fixed(quantile(solve, 0.75)),
        experiment_support::fixed(quantile(total, 0.50)),
        experiment_support::fixed(quantile(total, 0.25)),
        experiment_support::fixed(quantile(total, 0.75)),
        experiment_support::fixed(quantile(cycles, 0.50), 0)};
}

}  // namespace

int main(int argc, char** argv) {
    int threads = 4;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument.rfind("--threads=", 0) == 0) {
            threads = std::stoi(argument.substr(10));
        }
    }
    constexpr int maximum_cycles =
        experiment_support::maximum_two_grid_cycles;
    constexpr std::array<std::uint64_t, 5> seeds{1, 3, 7, 11, 19};
    const auto& cross = experiment_support::channel_topologies().front();
    experiment_support::BasicConfig config;
    config.fine_intervals = 64;
    config.coarse_intervals = 16;
    config.contrast = 1.0e4;
    config.threads = threads;
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);

    experiment_support::Rows seed_rows;
    std::map<std::string, CycleAggregate> seed_aggregates;
    for (std::uint64_t seed : seeds) {
        experiment_support::progress(
            "submission seed " + std::to_string(seed));
        const auto problem = experiment_support::make_problem(
            grid, cross, config, seed);
        const auto geometric =
            experiment_support::geometric_interpolation(grid);
        for (const auto& policy : {
                 std::pair<const char*, double>{"fast", 1.0},
                 std::pair<const char*, double>{"reuse", 256.0}}) {
            const auto adaptive = build_adaptive(
                grid, problem.matrix, geometric.prolongation, problem.rhs,
                policy.second, threads, maximum_cycles);
            const SolveResult solved = measure_solve(
                problem.rhs, *adaptive.cycle, maximum_cycles);
            add_aggregate(seed_aggregates[policy.first], solved);
            seed_rows.push_back({
                std::to_string(seed), policy.first,
                std::to_string(adaptive.report.selected_steps),
                std::to_string(solved.cycles),
                experiment_support::fixed(
                    experiment_support::interpolation_density_percent(
                        *adaptive.prolongation), 4),
                tgi::two_grid_status_name(solved.status)});
        }
        const auto exact = experiment_support::build_global_reference(
            grid, problem.matrix, threads);
        const tgi::TwoGridCycle exact_cycle(
            problem.matrix, exact.prolongation, 1, threads);
        const SolveResult solved = measure_solve(
            problem.rhs, exact_cycle, maximum_cycles);
        add_aggregate(seed_aggregates["global-reference"], solved);
        seed_rows.push_back({
            std::to_string(seed), "global-reference", "exact",
            std::to_string(solved.cycles),
            experiment_support::fixed(
                experiment_support::interpolation_density_percent(
                    exact.prolongation), 4),
            tgi::two_grid_status_name(solved.status)});
    }
    experiment_support::Rows seed_summary;
    for (const std::string method : {
             "fast", "reuse", "global-reference"}) {
        const auto& aggregate = seed_aggregates[method];
        seed_summary.push_back({
            method,
            std::to_string(aggregate.converged) + "/" +
                std::to_string(aggregate.slow) + "/" +
                std::to_string(aggregate.diverged),
            std::to_string(aggregate.cycles),
            aggregate.converged > 0
                ? experiment_support::fixed(
                    static_cast<double>(aggregate.cycles) /
                    static_cast<double>(aggregate.converged), 2)
                : "--",
            aggregate.converged > 0
                ? std::to_string(aggregate.maximum_cycles)
                : "--"});
    }

    experiment_support::progress("submission RHS transfer");
    const auto transfer_problem = experiment_support::make_problem(
        grid, cross, config, 1);
    const auto transfer_geometric =
        experiment_support::geometric_interpolation(grid);
    const auto fast = build_adaptive(
        grid, transfer_problem.matrix, transfer_geometric.prolongation,
        transfer_problem.rhs, 1.0, threads, maximum_cycles);
    const auto reuse = build_adaptive(
        grid, transfer_problem.matrix, transfer_geometric.prolongation,
        transfer_problem.rhs, 256.0, threads, maximum_cycles);
    const auto exact = experiment_support::build_global_reference(
        grid, transfer_problem.matrix, threads);
    const tgi::TwoGridCycle geometric_cycle(
        transfer_problem.matrix, transfer_geometric.prolongation, 1, threads);
    const tgi::TwoGridCycle exact_cycle(
        transfer_problem.matrix, exact.prolongation, 1, threads);
    const std::vector<RhsCase> rhs_cases = make_rhs_cases(grid);
    experiment_support::Rows rhs_rows;
    std::map<std::string, CycleAggregate> rhs_aggregates;
    for (const auto& rhs_case : rhs_cases) {
        for (const auto& method : {
                 std::pair<const char*, const tgi::TwoGridCycle*>{
                     "adaptive-fast", fast.cycle.get()},
                 {"adaptive-reuse", reuse.cycle.get()},
                 {"global-reference", &exact_cycle},
                 {"geometric", &geometric_cycle}}) {
            const SolveResult solved = measure_solve(
                rhs_case.values, *method.second, maximum_cycles);
            add_aggregate(rhs_aggregates[method.first], solved);
            rhs_rows.push_back({
                rhs_case.name, method.first,
                method.first == std::string("adaptive-fast")
                    ? std::to_string(fast.report.selected_steps)
                    : method.first == std::string("adaptive-reuse")
                        ? std::to_string(reuse.report.selected_steps)
                        : "-",
                std::to_string(solved.cycles),
                tgi::two_grid_status_name(solved.status)});
        }
    }
    experiment_support::Rows rhs_summary;
    for (const std::string method : {
             "adaptive-fast", "adaptive-reuse",
             "global-reference", "geometric"}) {
        const auto& aggregate = rhs_aggregates[method];
        rhs_summary.push_back({
            method,
            std::to_string(aggregate.converged) + "/" +
                std::to_string(aggregate.slow) + "/" +
                std::to_string(aggregate.diverged),
            std::to_string(aggregate.cycles),
            aggregate.converged > 0
                ? experiment_support::fixed(
                    static_cast<double>(aggregate.cycles) /
                    static_cast<double>(aggregate.converged), 2)
                : "--",
            aggregate.converged > 0
                ? std::to_string(aggregate.maximum_cycles)
                : "--"});
    }

    experiment_support::progress("submission repeated timing");
    (void)timing_sample(
        grid, transfer_problem.matrix, transfer_problem.rhs,
        1.0, threads, maximum_cycles);
    (void)timing_sample(
        grid, transfer_problem.matrix, transfer_problem.rhs,
        256.0, threads, maximum_cycles);
    std::vector<TimingSample> fast_timings;
    std::vector<TimingSample> reuse_timings;
    constexpr int repetitions = 5;
    for (int repetition = 0; repetition < repetitions; ++repetition) {
        experiment_support::progress(
            "timing repetition " + std::to_string(repetition + 1) + "/" +
            std::to_string(repetitions));
        if (repetition % 2 == 0) {
            fast_timings.push_back(timing_sample(
                grid, transfer_problem.matrix, transfer_problem.rhs,
                1.0, threads, maximum_cycles));
            reuse_timings.push_back(timing_sample(
                grid, transfer_problem.matrix, transfer_problem.rhs,
                256.0, threads, maximum_cycles));
        } else {
            reuse_timings.push_back(timing_sample(
                grid, transfer_problem.matrix, transfer_problem.rhs,
                256.0, threads, maximum_cycles));
            fast_timings.push_back(timing_sample(
                grid, transfer_problem.matrix, transfer_problem.rhs,
                1.0, threads, maximum_cycles));
        }
    }
    experiment_support::Rows timing_rows{
        timing_row("adaptive-fast", fast_timings),
        timing_row("adaptive-reuse", reuse_timings)};
    std::vector<double> fast_setup;
    std::vector<double> fast_solve;
    std::vector<double> reuse_setup;
    std::vector<double> reuse_solve;
    for (const auto& sample : fast_timings) {
        fast_setup.push_back(sample.setup_ms);
        fast_solve.push_back(sample.solve_ms);
    }
    for (const auto& sample : reuse_timings) {
        reuse_setup.push_back(sample.setup_ms);
        reuse_solve.push_back(sample.solve_ms);
    }
    const double setup_difference =
        quantile(reuse_setup, 0.50) - quantile(fast_setup, 0.50);
    const double solve_difference =
        quantile(fast_solve, 0.50) - quantile(reuse_solve, 0.50);
    const double break_even = solve_difference > 0.0
        ? setup_difference / solve_difference
        : std::numeric_limits<double>::infinity();

    experiment_support::Report report(
        "Submission-focused robustness checks");
    report.add_summary({
        {"Version", std::string(tgi::version)},
        {"Grid", "64/16"},
        {"Contrast", "1e4"},
        {"Topology", "cross-channel"},
        {"Threads", std::to_string(threads)},
        {"Solve tolerance", "1e-6"},
        {"Maximum cycles", std::to_string(maximum_cycles)}});
    report.add_note(
        "This single compact experiment adds only the three checks needed "
        "before submission: five coefficient seeds, six right-hand sides "
        "evaluated with interpolation trained only on the constant RHS, and "
        "five post-warmup timing repetitions. Timing order alternates by "
        "repetition; Q1, median and Q3 are reported instead of a single run.");
    report.add_table(
        "Coefficient-seed stability",
        {"Seed", "Method", "m", "Cycles", "P density %", "Status"},
        {7, 18, 7, 10, 12, 10}, seed_rows, true);
    report.add_table(
        "Seed aggregate",
        {"Method", "Conv/slow/div", "Converged cycle sum",
         "Mean converged", "Worst converged"},
        {18, 13, 19, 14, 15}, seed_summary);
    report.add_table(
        "Right-hand-side transfer (constant-RHS training)",
        {"Evaluation RHS", "Method", "m", "Cycles", "Status"},
        {18, 18, 7, 8, 10}, rhs_rows, true);
    report.add_table(
        "RHS aggregate",
        {"Method", "Conv/slow/div", "Converged cycle sum",
         "Mean converged", "Worst converged"},
        {18, 13, 19, 14, 15}, rhs_summary);
    report.add_table(
        "Repeated wall-clock summary",
        {"Policy", "Runs", "Setup med ms", "Setup Q1", "Setup Q3",
         "Solve med ms", "Solve Q1", "Solve Q3", "Total med ms", "Total Q1",
         "Total Q3", "Cycles med"},
        {16, 6, 11, 10, 10, 10, 9, 9, 10, 9, 9, 11}, timing_rows);
    report.add_note(
        "Using policy-wise medians in S+R*T, the measured fast/reuse "
        "break-even is R=" + experiment_support::fixed(break_even, 2) +
        " right-hand sides on this case; it is a workload measurement, not "
        "a universal threshold.");
    report.save("experiment4_submission_robustness");
    return 0;
}
