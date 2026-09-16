#include "experiment/problem.hpp"
#include "experiment/reporting.hpp"
#include "multigrid/global_pcg.hpp"
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
struct RhsCase {
    std::string name;
    tgi::Vector values;
};
enum class Method {
    FiniteOneThird,
    EnergyEndpoint
};
double elapsed_ms(Clock::time_point begin) {
    return std::chrono::duration<double, std::milli>(
               Clock::now() - begin)
        .count();
}
SolveResult measure_solve(
    const tgi::Vector& rhs, const tgi::TwoGridCycle& cycle,
    int maximum_cycles) {
    const auto begin = Clock::now();
    const auto solved = tgi::solve_two_grid(
        rhs, cycle, 1.0e-6, maximum_cycles);
    return {
        solved.cycles, elapsed_ms(begin), solved.effective_factor,
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
    constexpr double pi = 3.141592653589793238462643383279502884;
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
                    0.35 * std::sin(5.0 * pi * x) *
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
tgi::SparseMatrix make_transfer(
    const tgi::StructuredGrid& grid, const tgi::SparseMatrix& matrix,
    const tgi::SparseMatrix& initial, Method method, int threads) {
    tgi::GlobalEnergyPcgPath path(grid, matrix, initial, threads);
    if (method == Method::FiniteOneThird) {
        path.advance_to(experiment_support::fixed_path_steps(
            grid.intervals(), 1, 3));
    } else {
        path.advance_until_relative_residual(1.0e-10);
    }
    return path.prolongation();
}
TimingSample timing_sample(
    const tgi::StructuredGrid& grid, const tgi::SparseMatrix& matrix,
    const tgi::Vector& rhs, Method method, int threads,
    int maximum_cycles) {
    const auto begin = Clock::now();
    const auto initial = tgi::build_geometric_interpolation(grid);
    const auto transfer = make_transfer(
        grid, matrix, initial, method, threads);
    const tgi::TwoGridCycle cycle(matrix, transfer, 1, threads);
    const double setup_ms = elapsed_ms(begin);
    const SolveResult solved = measure_solve(rhs, cycle, maximum_cycles);
    return {setup_ms, solved.milliseconds, solved.cycles};
}
experiment_support::Row timing_row(
    const std::string& method, const std::vector<TimingSample>& samples) {
    double setup = 0.0;
    double solve = 0.0;
    double total = 0.0;
    double cycles = 0.0;
    for (const auto& sample : samples) {
        setup += sample.setup_ms;
        solve += sample.solve_ms;
        total += sample.setup_ms + sample.solve_ms;
        cycles += static_cast<double>(sample.cycles);
    }
    const double count = static_cast<double>(samples.size());
    return {
        method, std::to_string(samples.size()),
        experiment_support::fixed(setup / count),
        experiment_support::fixed(solve / count),
        experiment_support::fixed(total / count),
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
    const auto& field = experiment_support::channel_topologies().front();
    experiment_support::BasicConfig config;
    config.fine_intervals = 64;
    config.coarse_intervals = 16;
    config.contrast = 1.0e4;
    config.threads = threads;
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);
    experiment_support::Rows seed_rows;
    std::map<std::string, CycleAggregate> aggregates;
    for (std::uint64_t seed : seeds) {
        experiment_support::progress(
            "coefficient seed " + std::to_string(seed));
        const auto problem = experiment_support::make_problem(
            grid, field, config, seed);
        const auto initial = tgi::build_geometric_interpolation(grid);
        for (const auto& item : {
                 std::pair<std::string, Method>{
                     "finite-1/3", Method::FiniteOneThird},
                 {"energy-endpoint", Method::EnergyEndpoint}}) {
            const auto transfer = make_transfer(
                grid, problem.matrix, initial, item.second, threads);
            const tgi::TwoGridCycle cycle(
                problem.matrix, transfer, 1, threads);
            const SolveResult solved = measure_solve(
                problem.rhs, cycle, maximum_cycles);
            add_aggregate(aggregates[item.first], solved);
            seed_rows.push_back({
                std::to_string(seed), item.first,
                item.second == Method::FiniteOneThird
                    ? "m=" + std::to_string(
                          experiment_support::fixed_path_steps(
                              config.fine_intervals, 1, 3))
                    : "relres<=1e-10",
                std::to_string(solved.cycles),
                experiment_support::fixed(solved.effective_factor, 6),
                experiment_support::fixed(
                    experiment_support::interpolation_density_percent(
                        transfer), 4),
                tgi::stationary_status_name(solved.status)});
        }
    }
    experiment_support::Rows aggregate_rows;
    for (const std::string method : {"finite-1/3", "energy-endpoint"}) {
        const auto& value = aggregates[method];
        aggregate_rows.push_back({
            method,
            std::to_string(value.converged) + "/" +
                std::to_string(value.slow) + "/" +
                std::to_string(value.diverged),
            std::to_string(value.cycles),
            experiment_support::fixed(
                static_cast<double>(value.cycles) /
                    static_cast<double>(value.converged),
                2),
            std::to_string(value.maximum_cycles)});
    }
    experiment_support::progress("right-hand-side robustness");
    const auto rhs_problem = experiment_support::make_problem(
        grid, field, config, 1);
    const auto rhs_initial = tgi::build_geometric_interpolation(grid);
    const auto finite_transfer = make_transfer(
        grid, rhs_problem.matrix, rhs_initial,
        Method::FiniteOneThird, threads);
    const auto endpoint_transfer = make_transfer(
        grid, rhs_problem.matrix, rhs_initial,
        Method::EnergyEndpoint, threads);
    const tgi::TwoGridCycle finite_cycle(
        rhs_problem.matrix, finite_transfer, 1, threads);
    const tgi::TwoGridCycle endpoint_cycle(
        rhs_problem.matrix, endpoint_transfer, 1, threads);
    experiment_support::Rows rhs_rows;
    std::map<std::string, CycleAggregate> rhs_aggregates;
    for (const auto& rhs_case : make_rhs_cases(grid)) {
        for (const auto& item : {
                 std::pair<std::string, const tgi::TwoGridCycle*>{
                     "finite-1/3", &finite_cycle},
                 {"energy-endpoint", &endpoint_cycle}}) {
            const SolveResult solved = measure_solve(
                rhs_case.values, *item.second, maximum_cycles);
            add_aggregate(rhs_aggregates[item.first], solved);
            rhs_rows.push_back({
                rhs_case.name, item.first,
                item.first == "finite-1/3"
                    ? "m=" + std::to_string(
                          experiment_support::fixed_path_steps(
                              config.fine_intervals, 1, 3))
                    : "relres<=1e-10",
                std::to_string(solved.cycles),
                experiment_support::fixed(solved.effective_factor, 6),
                tgi::stationary_status_name(solved.status)});
        }
    }
    experiment_support::Rows rhs_aggregate_rows;
    for (const std::string method : {"finite-1/3", "energy-endpoint"}) {
        const auto& value = rhs_aggregates[method];
        rhs_aggregate_rows.push_back({
            method,
            std::to_string(value.converged) + "/" +
                std::to_string(value.slow) + "/" +
                std::to_string(value.diverged),
            std::to_string(value.cycles),
            experiment_support::fixed(
                static_cast<double>(value.cycles) /
                    static_cast<double>(value.converged),
                2),
            std::to_string(value.maximum_cycles)});
    }
    experiment_support::BasicConfig timing_config = config;
    timing_config.fine_intervals = 128;
    const tgi::StructuredGrid timing_grid =
        experiment_support::make_grid(timing_config);
    const auto timing_problem = experiment_support::make_problem(
        timing_grid, field, timing_config, 1);
    for (Method method : {Method::FiniteOneThird, Method::EnergyEndpoint}) {
        (void)timing_sample(
            timing_grid, timing_problem.matrix, timing_problem.rhs,
            method, threads, maximum_cycles);
    }
    std::vector<TimingSample> finite_timings;
    std::vector<TimingSample> endpoint_timings;
    constexpr int repetitions = 5;
    for (int repetition = 0; repetition < repetitions; ++repetition) {
        experiment_support::progress(
            "timing repetition " + std::to_string(repetition + 1) + "/" +
            std::to_string(repetitions));
        const std::array<Method, 2> order = repetition % 2 == 0
            ? std::array<Method, 2>{
                  Method::FiniteOneThird, Method::EnergyEndpoint}
            : std::array<Method, 2>{
                  Method::EnergyEndpoint, Method::FiniteOneThird};
        for (Method method : order) {
            auto sample = timing_sample(
                timing_grid, timing_problem.matrix, timing_problem.rhs,
                method, threads, maximum_cycles);
            (method == Method::FiniteOneThird
                 ? finite_timings : endpoint_timings)
                .push_back(sample);
        }
    }
    experiment_support::Report report(
        "Coefficient and right-hand-side robustness with repeated timing");
    report.add_summary({
        {"Seed grid", "64/16"},
        {"Timing grid", "128/16"},
        {"Contrast", "1e4"},
        {"Topology", "cross-channel"},
        {"Threads", std::to_string(threads)},
        {"Solve tolerance", "1e-6"}});
    report.add_table(
        "Coefficient-seed stability",
        {"Seed", "Method", "Parameter", "Cycles", "Effective factor",
         "P density %", "Status"},
        {7, 18, 16, 10, 16, 12, 10}, seed_rows, true);
    report.add_table(
        "Seed aggregate",
        {"Method", "Conv/slow/div", "Converged cycle sum",
         "Mean converged", "Worst converged"},
        {18, 13, 19, 14, 15}, aggregate_rows);
    report.add_table(
        "Right-hand-side robustness with fixed interpolation",
        {"Evaluation RHS", "Method", "Parameter", "Cycles",
         "Effective factor", "Status"},
        {18, 18, 16, 10, 16, 10}, rhs_rows, true);
    report.add_table(
        "RHS aggregate",
        {"Method", "Conv/slow/div", "Converged cycle sum",
         "Mean converged", "Worst converged"},
        {18, 13, 19, 14, 15}, rhs_aggregate_rows);
    report.add_table(
        "Central 128/16 repeated wall-clock comparison",
        {"Method", "Runs", "Setup mean", "Solve mean", "Total mean",
         "Mean cycles"},
        {18, 6, 11, 10, 10, 11},
        {timing_row("finite-1/3", finite_timings),
         timing_row("energy-endpoint", endpoint_timings)});
    report.save("experiment5_robustness");
    return 0;
}
