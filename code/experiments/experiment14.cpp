#include "experiment/study.hpp"
#include "multigrid/adaptive_global_pcg.hpp"
#include "multigrid/global_pcg_path.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct WorkloadCase {
    int fine;
    int coarse;
    double contrast;
    std::uint64_t seed;
    experiment_support::FieldCase field;
};

std::vector<tgi::Vector> make_rhs_set(const tgi::StructuredGrid& grid) {
    constexpr double pi = 3.14159265358979323846;
    std::vector<tgi::Vector> result(
        5, tgi::Vector(static_cast<std::size_t>(grid.fine_size()), 0.0));
    for (int id = 0; id < grid.fine_size(); ++id) {
        const auto [ix, iy] = grid.fine_coords(id);
        const double x = static_cast<double>(ix + 1) * grid.h();
        const double y = static_cast<double>(iy + 1) * grid.h();
        result[0][static_cast<std::size_t>(id)] = 1.0;
        result[1][static_cast<std::size_t>(id)] =
            std::sin(pi * x) * std::sin(2.0 * pi * y) + 0.15;
        const double dx = x - 0.31;
        const double dy = y - 0.69;
        result[2][static_cast<std::size_t>(id)] =
            std::exp(-80.0 * (dx * dx + dy * dy));
        result[3][static_cast<std::size_t>(id)] =
            std::cos(2.0 * pi * x) * std::sin(pi * y) +
            0.25 * std::sin(3.0 * pi * x) * std::sin(2.0 * pi * y);
        result[4][static_cast<std::size_t>(id)] =
            std::sin(3.0 * pi * x + 0.2) *
                std::sin(4.0 * pi * y + 0.7) +
            0.35 * std::sin(7.0 * pi * x) * std::sin(2.0 * pi * y);
    }
    return result;
}

struct WorkloadMetric {
    double hierarchy_setup_ms = 0.0;
    double average_solve_ms = 0.0;
    double average_cycles = 0.0;
    int maximum_cycles = 0;
    bool all_converged = true;
};

WorkloadMetric evaluate_workload(
    const tgi::SparseMatrix& a, const tgi::SparseMatrix& p,
    const std::vector<tgi::Vector>& rhs_set, int threads,
    int maximum_cycles) {
    WorkloadMetric result;
    const tgi::TwoGridCycle cycle(a, p, 1, threads);
    result.hierarchy_setup_ms = cycle.setup_report().total_ms;
    // Warm allocations with two cycles.  A previous full-solve warmup could
    // spend another maximum_cycles iterations on a failing hierarchy and
    // looked like a hang even though it was bounded.
    tgi::Vector warm_solution(rhs_set.front().size(), 0.0);
    tgi::Vector warm_residual = rhs_set.front();
    tgi::TwoGridCycle::Workspace warm_workspace;
    for (int step = 0; step < 2; ++step) {
        (void)cycle.iterate(
            rhs_set.front(), warm_solution, warm_residual, warm_workspace);
    }
    double solve_ms = 0.0;
    long long cycles = 0;
    int repeats_completed = 0;
    constexpr int requested_repeats = 3;
    for (int repeat = 0; repeat < requested_repeats; ++repeat) {
        bool repeat_converged = true;
        for (const tgi::Vector& rhs : rhs_set) {
            const auto begin = std::chrono::steady_clock::now();
            const auto solved = tgi::solve_two_grid(
                a, rhs, cycle, 1.0e-6, maximum_cycles);
            solve_ms += std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - begin).count();
            cycles += solved.cycles;
            result.maximum_cycles = std::max(
                result.maximum_cycles, solved.cycles);
            repeat_converged = repeat_converged && solved.converged;
        }
        ++repeats_completed;
        result.all_converged = result.all_converged && repeat_converged;
        if (!repeat_converged) break;
    }
    const double solves = static_cast<double>(
        rhs_set.size() * static_cast<std::size_t>(repeats_completed));
    result.average_solve_ms = solve_ms / solves;
    result.average_cycles = static_cast<double>(cycles) /
        solves;
    return result;
}

std::string break_even_rhs(
    double adaptive_setup, double adaptive_solve,
    double baseline_setup, double baseline_solve) {
    if (adaptive_setup <= baseline_setup) return "0";
    const double saving = baseline_solve - adaptive_solve;
    if (!(saving > 0.0)) return "never";
    return std::to_string(static_cast<int>(std::ceil(
        (adaptive_setup - baseline_setup) / saving)));
}

tgi::AdaptiveGlobalPcgOptions quality_options(
    const tgi::StructuredGrid& grid, int threads, int maximum_cycles) {
    tgi::AdaptiveGlobalPcgOptions options;
    options.cost_aware_mode = false;
    options.minimum_steps = grid.intervals() >= 128 ? 16 : 12;
    options.maximum_steps = grid.intervals() >= 128 ? 64 : 56;
    options.maximum_screening_steps = grid.intervals() >= 128 ? 48 : 44;
    options.screening_increment = 16;
    options.screening_pilot_iterations = grid.ratio() >= 8 ? 32 : 24;
    options.screening_tail_window = grid.ratio() >= 8 ? 8 : 6;
    options.minimum_screened_positive_candidates =
        grid.intervals() >= 128 ? 3 : 4;
    options.refinement_backtrack_steps = 10;
    options.refinement_stop_before_anchor_steps = 4;
    options.refinement_increment = 2;
    options.refinement_pilot_iterations = grid.ratio() >= 8 ? 64 : 48;
    options.refinement_tail_window = grid.ratio() >= 8 ? 16 : 12;
    options.confirmation_candidates = 2;
    options.maximum_confirmation_cycles = maximum_cycles;
    options.thread_count = threads;
    return options;
}

} // namespace

int main(int argc, char** argv) {
    int threads = 4;
    if (argc == 2) {
        const std::string argument = argv[1];
        if (argument.rfind("--threads=", 0) != 0) {
            throw std::invalid_argument("unknown argument: " + argument);
        }
        threads = std::stoi(argument.substr(10));
    } else if (argc > 2) {
        throw std::invalid_argument("experiment14 accepts only --threads");
    }
    const auto& topology = experiment_support::channel_topologies();
    const std::array<WorkloadCase, 4> cases{{
        {64, 8, 1.0e4, 17, topology[1]},
        {64, 16, 1.0e4, 1, topology[2]},
        {64, 16, 1.0e6, 1, topology[0]},
        {128, 16, 1.0e4, 1, topology[0]}
    }};
    constexpr int maximum_cycles = 4000;
    const experiment_support::Row headers{
        "1/h", "1/H", "Contrast", "Topology", "Method", "m",
        "One-time setup ms", "Avg cycles", "Max cycles", "Avg solve ms",
        "Break-even vs G", "Break-even vs m40"};
    experiment_support::Rows rows;

    int case_index = 0;
    for (const WorkloadCase& item : cases) {
        ++case_index;
        experiment_support::progress(
            "experiment14 case " + std::to_string(case_index) + "/" +
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
        const auto rhs_set = make_rhs_set(grid);
        const auto geometric = experiment_support::geometric_interpolation(
            grid, problem.matrix);
        const WorkloadMetric geometric_metric = evaluate_workload(
            problem.matrix, geometric.prolongation, rhs_set,
            threads, maximum_cycles);
        const double geometric_setup = geometric.report.timing.total_ms +
            geometric_metric.hierarchy_setup_ms;

        tgi::GlobalEnergyPcgPath fixed_path(
            grid, problem.matrix, geometric.prolongation, threads);
        fixed_path.advance_to(40);
        const auto fixed_report = fixed_path.report();
        const tgi::SparseMatrix fixed_p = fixed_path.prolongation(0.0);
        const WorkloadMetric fixed_metric = evaluate_workload(
            problem.matrix, fixed_p, rhs_set, threads, maximum_cycles);
        const double fixed_setup = geometric.report.timing.total_ms +
            fixed_report.total_ms + fixed_metric.hierarchy_setup_ms;

        tgi::AdaptiveGlobalPcgOptions budget_options;
        budget_options.minimum_steps = grid.intervals() >= 128 ? 16 : 12;
        budget_options.maximum_steps = 64;
        budget_options.maximum_confirmation_cycles = maximum_cycles;
        budget_options.thread_count = threads;
        const auto budget = tgi::build_adaptive_global_pcg_interpolation(
            grid, problem.matrix, geometric.prolongation,
            budget_options, &problem.rhs);
        const WorkloadMetric budget_metric = evaluate_workload(
            problem.matrix, budget.prolongation, rhs_set,
            threads, maximum_cycles);
        const double budget_setup = geometric.report.timing.total_ms +
            budget.report.selection_wall_ms +
            budget_metric.hierarchy_setup_ms;

        auto staged_options = quality_options(
            grid, threads, maximum_cycles);
        const auto staged = tgi::build_adaptive_global_pcg_interpolation(
            grid, problem.matrix, geometric.prolongation,
            staged_options, &problem.rhs);
        const WorkloadMetric staged_metric = evaluate_workload(
            problem.matrix, staged.prolongation, rhs_set,
            threads, maximum_cycles);
        const double staged_setup = geometric.report.timing.total_ms +
            staged.report.selection_wall_ms +
            staged_metric.hierarchy_setup_ms;

        const auto append = [&](const std::string& method, int steps,
                                double setup, const WorkloadMetric& metric) {
            const std::string versus_fixed =
                method == "budget-v3.2" && steps == 40
                    ? "never"
                    : break_even_rhs(
                          setup, metric.average_solve_ms,
                          fixed_setup, fixed_metric.average_solve_ms);
            rows.push_back({
                std::to_string(item.fine), std::to_string(item.coarse),
                experiment_support::scientific(item.contrast, 0),
                item.field.name, method, std::to_string(steps),
                experiment_support::fixed(setup),
                experiment_support::fixed(metric.average_cycles, 1),
                metric.all_converged
                    ? std::to_string(metric.maximum_cycles) : "failed",
                experiment_support::fixed(metric.average_solve_ms),
                break_even_rhs(
                    setup, metric.average_solve_ms,
                    geometric_setup, geometric_metric.average_solve_ms),
                versus_fixed});
        };
        append("geometric", 0, geometric_setup, geometric_metric);
        append("fixed", 40, fixed_setup, fixed_metric);
        append("budget-v3.2", budget.report.selected_steps,
               budget_setup, budget_metric);
        append("staged-v3.1", staged.report.selected_steps,
               staged_setup, staged_metric);
    }

    experiment_support::Report report(
        "Five-right-hand-side workload and setup break-even study");
    report.add_summary({
        {"Version", std::string(tgi::version)},
        {"Problems", std::to_string(cases.size())},
        {"Right-hand sides per problem", "5"},
        {"Solve tolerance", "1e-6"}});
    report.add_note(
        "Each hierarchy is constructed once and reused for five distinct "
        "smooth/localized right-hand sides. Break-even is the smallest RHS "
        "count at which setup plus mean solve time favors the method over the "
        "named baseline; 'never' means the measured per-RHS time is not lower. "
        "Timings are single-run wall-clock measurements and should be read as "
        "order-of-magnitude evidence, not microbenchmarks.");
    report.add_table(
        "Amortized workload", headers,
        {5, 5, 10, 20, 13, 5, 18, 11, 10, 13, 15, 17}, rows, true);
    report.save("experiment14");
    experiment_support::write_csv("experiment14", headers, rows);
    return 0;
}
