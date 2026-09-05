#include "experiment/problem.hpp"
#include "experiment/reporting.hpp"
#include "multigrid/global_pcg.hpp"
#include "version.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
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
    double effective_factor = 1.0;
    tgi::StationaryIterationStatus status =
        tgi::StationaryIterationStatus::SlowAtLimit;
    bool converged = false;
};

struct CycleAggregate {
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

struct TimingStatistics {
    double mean = 0.0;
    double standard_deviation = 0.0;
};

enum class TimingMethod {
    Adaptive,
    GlobalReference
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
        solved.effective_factor,
        solved.status, solved.converged};
}

void add_aggregate(
    CycleAggregate& aggregate, const SolveResult& measurement) {
    if (measurement.converged) {
        ++aggregate.converged;
        aggregate.cycles += measurement.cycles;
        aggregate.maximum_cycles = std::max(
            aggregate.maximum_cycles, measurement.cycles);
    } else if (measurement.status ==
               tgi::StationaryIterationStatus::Diverged) {
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

TimingSample timing_sample(
    const tgi::StructuredGrid& grid, const tgi::SparseMatrix& matrix,
    const tgi::Vector& rhs, TimingMethod method, int threads,
    int maximum_cycles) {
    const auto setup_begin = Clock::now();
    const auto initial = tgi::build_geometric_interpolation(grid);
    if (method == TimingMethod::GlobalReference) {
        tgi::GlobalEnergyPcgPath path(
            grid, matrix, initial.prolongation, threads);
        path.advance_until_relative_residual(1.0e-10);
        const tgi::SparseMatrix reference = path.prolongation();
        const tgi::TwoGridCycle cycle(matrix, reference, 1, threads);
        const double setup_ms = elapsed_ms(setup_begin, Clock::now());
        const SolveResult solved = measure_solve(
            rhs, cycle, maximum_cycles);
        return {setup_ms, solved.milliseconds, solved.cycles};
    }
    const auto adaptive = tgi::build_adaptive_global_pcg_interpolation(
        grid, matrix, initial.prolongation, threads);
    const double setup_ms = elapsed_ms(setup_begin, Clock::now());
    const SolveResult solved = measure_solve(
        rhs, *adaptive.cycle, maximum_cycles);
    return {setup_ms, solved.milliseconds, solved.cycles};
}

TimingStatistics statistics(const std::vector<double>& values) {
    double sum = 0.0;
    for (double value : values) sum += value;
    const double count = static_cast<double>(values.size());
    const double mean = sum / count;
    double squared_deviations = 0.0;
    for (double value : values) {
        const double difference = value - mean;
        squared_deviations += difference * difference;
    }
    return {
        mean,
        std::sqrt(squared_deviations /
                  static_cast<double>(values.size() - 1U))};
}

experiment_support::Row timing_row(
    const std::string& policy, const std::vector<TimingSample>& samples) {
    std::vector<double> setup;
    std::vector<double> solve;
    std::vector<double> total;
    double cycles = 0.0;
    setup.reserve(samples.size());
    solve.reserve(samples.size());
    total.reserve(samples.size());
    for (const auto& sample : samples) {
        setup.push_back(sample.setup_ms);
        solve.push_back(sample.solve_ms);
        total.push_back(sample.setup_ms + sample.solve_ms);
        cycles += static_cast<double>(sample.cycles);
    }
    const double count = static_cast<double>(samples.size());
    const TimingStatistics setup_statistics = statistics(setup);
    const TimingStatistics solve_statistics = statistics(solve);
    const TimingStatistics total_statistics = statistics(total);
    return {
        policy, std::to_string(samples.size()),
        experiment_support::fixed(setup_statistics.mean),
        experiment_support::fixed(setup_statistics.standard_deviation),
        experiment_support::fixed(solve_statistics.mean),
        experiment_support::fixed(solve_statistics.standard_deviation),
        experiment_support::fixed(total_statistics.mean),
        experiment_support::fixed(total_statistics.standard_deviation),
        experiment_support::fixed(cycles / count, 0)};
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
        const auto initial = tgi::build_geometric_interpolation(grid);
        const auto adaptive = tgi::build_adaptive_global_pcg_interpolation(
            grid, problem.matrix, initial.prolongation, threads);
        const SolveResult adaptive_solved = measure_solve(
            problem.rhs, *adaptive.cycle, maximum_cycles);
        add_aggregate(seed_aggregates["adaptive"], adaptive_solved);
        seed_rows.push_back({
            std::to_string(seed), "adaptive",
            std::to_string(adaptive.report.selected_steps),
            std::to_string(adaptive_solved.cycles),
            experiment_support::fixed(
                adaptive_solved.effective_factor, 6),
            experiment_support::fixed(
                experiment_support::interpolation_density_percent(
                    *adaptive.prolongation), 4),
            tgi::stationary_status_name(adaptive_solved.status)});
        const auto reference = experiment_support::build_global_reference(
            grid, problem.matrix, threads);
        const tgi::TwoGridCycle reference_cycle(
            problem.matrix, reference.prolongation, 1, threads);
        const SolveResult solved = measure_solve(
            problem.rhs, reference_cycle, maximum_cycles);
        add_aggregate(seed_aggregates["global-reference"], solved);
        seed_rows.push_back({
            std::to_string(seed), "global-reference", "tol=1e-10",
            std::to_string(solved.cycles),
            experiment_support::fixed(solved.effective_factor, 6),
            experiment_support::fixed(
                experiment_support::interpolation_density_percent(
                    reference.prolongation), 4),
            tgi::stationary_status_name(solved.status)});
    }
    experiment_support::Rows seed_summary;
    for (const std::string method : {"adaptive", "global-reference"}) {
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

    experiment_support::progress("RHS robustness");
    const auto transfer_problem = experiment_support::make_problem(
        grid, cross, config, 1);
    const auto transfer_initial = tgi::build_geometric_interpolation(grid);
    const auto adaptive = tgi::build_adaptive_global_pcg_interpolation(
        grid, transfer_problem.matrix, transfer_initial.prolongation,
        threads);
    const auto reference = experiment_support::build_global_reference(
        grid, transfer_problem.matrix, threads);
    const tgi::TwoGridCycle reference_cycle(
        transfer_problem.matrix, reference.prolongation, 1, threads);
    const std::vector<RhsCase> rhs_cases = make_rhs_cases(grid);
    experiment_support::Rows rhs_rows;
    std::map<std::string, CycleAggregate> rhs_aggregates;
    for (const auto& rhs_case : rhs_cases) {
        for (const auto& method : {
                 std::pair<const char*, const tgi::TwoGridCycle*>{
                     "adaptive", adaptive.cycle.get()},
                 {"global-reference", &reference_cycle}}) {
            const SolveResult solved = measure_solve(
                rhs_case.values, *method.second, maximum_cycles);
            add_aggregate(rhs_aggregates[method.first], solved);
            rhs_rows.push_back({
                rhs_case.name, method.first,
                method.first == std::string("adaptive")
                    ? std::to_string(adaptive.report.selected_steps)
                    : "-",
                std::to_string(solved.cycles),
                experiment_support::fixed(solved.effective_factor, 6),
                tgi::stationary_status_name(solved.status)});
        }
    }
    experiment_support::Rows rhs_summary;
    for (const std::string method : {"adaptive", "global-reference"}) {
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
    experiment_support::BasicConfig timing_config = config;
    timing_config.fine_intervals = 128;
    const tgi::StructuredGrid timing_grid =
        experiment_support::make_grid(timing_config);
    const auto timing_problem = experiment_support::make_problem(
        timing_grid, cross, timing_config, 1);
    (void)timing_sample(
        timing_grid, timing_problem.matrix, timing_problem.rhs,
        TimingMethod::Adaptive, threads, maximum_cycles);
    (void)timing_sample(
        timing_grid, timing_problem.matrix, timing_problem.rhs,
        TimingMethod::GlobalReference, threads, maximum_cycles);
    std::vector<TimingSample> adaptive_timings;
    std::vector<TimingSample> reference_timings;
    constexpr int repetitions = 5;
    for (int repetition = 0; repetition < repetitions; ++repetition) {
        experiment_support::progress(
            "timing repetition " + std::to_string(repetition + 1) + "/" +
            std::to_string(repetitions));
        const std::array<TimingMethod, 2> order = repetition % 2 == 0
            ? std::array<TimingMethod, 2>{
                  TimingMethod::Adaptive, TimingMethod::GlobalReference}
            : std::array<TimingMethod, 2>{
                  TimingMethod::GlobalReference, TimingMethod::Adaptive};
        for (TimingMethod method : order) {
            TimingSample sample = timing_sample(
                timing_grid, timing_problem.matrix, timing_problem.rhs,
                method, threads, maximum_cycles);
            if (method == TimingMethod::Adaptive) {
                adaptive_timings.push_back(sample);
            } else {
                reference_timings.push_back(sample);
            }
        }
    }
    experiment_support::Rows timing_rows{
        timing_row("adaptive", adaptive_timings),
        timing_row("global-reference", reference_timings)};
    experiment_support::Report report(
        "Robustness and repeated-timing checks");
    report.add_summary({
        {"Version", std::string(tgi::version)},
        {"Grid", "64/16"},
        {"Timing grid", "128/16"},
        {"Contrast", "1e4"},
        {"Topology", "cross-channel"},
        {"Threads", std::to_string(threads)},
        {"Solve tolerance", "1e-6"},
        {"Maximum cycles", "20000"}});
    report.add_note(
        "Five coefficient seeds and six right-hand sides test robustness; "
        "the latter reuse each method's matrix-dependent interpolation. "
        "Timing is restricted to the central 128/16 cross-channel case. "
        "Adaptive and global-reference are measured in five post-warmup "
        "repetitions with alternating order. Both start from the geometric "
        "interpolation and use the same Jacobi-PCG path; global-reference "
        "continues until every column reaches relative residual 1e-10. "
        "Setup includes interpolation, Galerkin assembly and cycle "
        "construction; solve starts from zero. Means and sample standard "
        "deviations summarize the post-warmup measurements.");
    report.add_table(
        "Coefficient-seed stability",
        {"Seed", "Method", "Parameter", "Cycles", "Effective factor",
         "P density %", "Status"},
        {7, 18, 10, 10, 16, 12, 10}, seed_rows, true);
    report.add_table(
        "Seed aggregate",
        {"Method", "Conv/slow/div", "Converged cycle sum",
         "Mean converged", "Worst converged"},
        {18, 13, 19, 14, 15}, seed_summary);
    report.add_table(
        "Right-hand-side robustness (fixed interpolation)",
        {"Evaluation RHS", "Method", "m", "Cycles", "Effective factor",
         "Status"},
        {18, 18, 7, 8, 16, 10}, rhs_rows, true);
    report.add_table(
        "RHS aggregate",
        {"Method", "Conv/slow/div", "Converged cycle sum",
         "Mean converged", "Worst converged"},
        {18, 13, 19, 14, 15}, rhs_summary);
    report.add_table(
        "Central 128/16 repeated wall-clock comparison",
        {"Policy", "Runs", "Setup mean", "Setup SD", "Solve mean",
         "Solve SD", "Total mean", "Total SD", "Mean cycles"},
        {16, 6, 11, 9, 10, 8, 10, 8, 11}, timing_rows);
    report.save("experiment4_submission_robustness");
    return 0;
}
