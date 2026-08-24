#include "experiment/comparison_cases.hpp"
#include "experiment/reporting.hpp"
#include "multigrid/multilevel_hierarchy.hpp"
#include "version.hpp"

#include <chrono>
#include <map>
#include <string>
#include <vector>

namespace {

std::string level_steps(
    const std::vector<int>& steps,
    tgi::MultilevelInterpolationPolicy policy) {
    if (policy == tgi::MultilevelInterpolationPolicy::Geometric) {
        return "geometric";
    }
    if (policy == tgi::MultilevelInterpolationPolicy::ExactGlobalEnergy) {
        return "exact/all";
    }
    std::string text;
    for (int value : steps) {
        if (!text.empty()) text += "/";
        text += std::to_string(value);
    }
    return text;
}

struct PolicyCase {
    std::string name;
    tgi::MultilevelInterpolationPolicy policy;
};

struct Aggregate {
    int cases = 0;
    int stationary_converged = 0;
    int pcg_converged = 0;
    long long stationary_iterations = 0;
    long long pcg_iterations = 0;
    double operator_complexity = 0.0;
    double interpolation_complexity = 0.0;
    double setup_ms = 0.0;
    double stationary_ms = 0.0;
    double pcg_ms = 0.0;
};

}

int run_multilevel_comparison(int argc, char** argv) {
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
    const std::vector<PolicyCase> policies{
        {"adaptive-PCG", tgi::MultilevelInterpolationPolicy::AdaptiveGlobalPcg},
        {"global-exact", tgi::MultilevelInterpolationPolicy::ExactGlobalEnergy},
        {"geometric", tgi::MultilevelInterpolationPolicy::Geometric}};
    constexpr int stationary_limit = 3000;
    constexpr int pcg_limit = 1000;
    experiment_support::Rows rows;
    std::map<std::string, Aggregate> aggregates;
    int case_number = 0;

    for (const auto& item : cases) {
        ++case_number;
        experiment_support::progress(
            "multilevel comparison " + std::to_string(case_number) + "/" +
            std::to_string(cases.size()) + ": " + item.field.name);
        const auto config = experiment_support::comparison_config(
            item, threads, stationary_limit);
        const tgi::StructuredGrid grid = experiment_support::make_grid(config);
        const auto problem = experiment_support::make_problem(
            grid, item.field, config, item.seed);

        for (const auto& policy : policies) {
            tgi::MultilevelHierarchyOptions options;
            options.first_coarse_intervals = item.coarse;
            options.thread_count = threads;
            options.maximum_pilot_cycles = stationary_limit;
            options.first_level_adaptive.minimum_steps = 12;
            options.first_level_adaptive.maximum_steps = 56;
            options.first_level_adaptive.maximum_cycles = stationary_limit;
            options.first_level_adaptive.thread_count = threads;
            const auto hierarchy = tgi::build_multilevel_hierarchy(
                item.fine, problem.matrix, problem.rhs,
                policy.policy, options);

            const auto stationary_begin =
                std::chrono::steady_clock::now();
            const auto stationary = tgi::solve_multilevel(
                problem.rhs, *hierarchy.cycle, 1.0e-6,
                stationary_limit);
            const double stationary_ms =
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() -
                    stationary_begin).count();

            const auto pcg_begin = std::chrono::steady_clock::now();
            const auto pcg = tgi::solve_preconditioned_cg(
                problem.matrix, problem.rhs, *hierarchy.cycle,
                1.0e-6, pcg_limit);
            const double pcg_ms =
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - pcg_begin).count();
            const auto& setup = hierarchy.cycle->setup_report();
            Aggregate& aggregate = aggregates[policy.name];
            ++aggregate.cases;
            aggregate.stationary_converged += stationary.converged ? 1 : 0;
            aggregate.pcg_converged += pcg.converged ? 1 : 0;
            aggregate.stationary_iterations += stationary.iterations;
            aggregate.pcg_iterations += pcg.iterations;
            aggregate.operator_complexity += setup.operator_complexity;
            aggregate.interpolation_complexity +=
                setup.interpolation_complexity;
            aggregate.setup_ms += hierarchy.hierarchy.total_setup_ms;
            aggregate.stationary_ms += stationary_ms;
            aggregate.pcg_ms += pcg_ms;
            rows.push_back({
                std::to_string(item.fine), std::to_string(item.coarse),
                experiment_support::scientific(item.contrast, 0),
                item.field.name, policy.name,
                std::to_string(setup.levels),
                level_steps(
                    hierarchy.hierarchy.selected_steps, policy.policy),
                experiment_support::fixed(setup.operator_complexity, 4),
                experiment_support::fixed(
                    setup.interpolation_complexity, 4),
                experiment_support::fixed(
                    hierarchy.hierarchy.total_setup_ms),
                stationary.converged
                    ? std::to_string(stationary.iterations)
                    : "failed@" + std::to_string(stationary.iterations),
                experiment_support::fixed(stationary_ms),
                pcg.converged ? std::to_string(pcg.iterations)
                    : "failed@" + std::to_string(pcg.iterations),
                experiment_support::fixed(pcg_ms)});
        }
    }

    experiment_support::Rows aggregate_rows;
    for (const auto& policy : policies) {
        const Aggregate& value = aggregates[policy.name];
        const double count = static_cast<double>(value.cases);
        aggregate_rows.push_back({
            policy.name,
            std::to_string(value.stationary_converged) + "/" +
                std::to_string(value.cases),
            std::to_string(value.stationary_iterations),
            std::to_string(value.pcg_converged) + "/" +
                std::to_string(value.cases),
            std::to_string(value.pcg_iterations),
            experiment_support::fixed(value.operator_complexity / count, 4),
            experiment_support::fixed(
                value.interpolation_complexity / count, 4),
            experiment_support::fixed(value.setup_ms),
            experiment_support::fixed(value.stationary_ms),
            experiment_support::fixed(value.pcg_ms)});
    }

    experiment_support::Report report(
        "Multilevel interpolation comparison: V-cycle and PCG");
    report.add_summary({
        {"Version", std::string(tgi::version)},
        {"Cases", std::to_string(cases.size())},
        {"Mode", quick ? "quick" : "full"},
        {"Coarsening", "requested first H, then factor two to one unknown"},
        {"Cycle", "one forward/backward Gauss-Seidel V-cycle"},
        {"Tolerance", "1e-6"}});
    report.add_note(
        "Each interpolation policy is applied consistently at every level. "
        "The same symmetric V-cycle is measured both as a standalone "
        "stationary iteration and as an SPD preconditioner for conjugate "
        "gradients. This separates coarse-space quality from the choice of "
        "outer iteration.");
    report.add_table(
        "All multilevel cases",
        {"1/h", "1/H", "Contrast", "Topology", "Interpolation",
         "Levels", "m by level", "C_op", "C_P", "Setup ms",
         "V iterations", "V solve ms", "PCG iterations", "PCG solve ms"},
        {5, 5, 10, 20, 15, 7, 16, 8, 8, 10, 13, 11, 15, 12},
        rows, true);
    report.add_table(
        "Aggregate multilevel comparison",
        {"Interpolation", "V converged", "V iteration sum",
         "PCG converged", "PCG iteration sum", "Mean C_op", "Mean C_P",
         "Setup sum ms", "V solve sum ms", "PCG solve sum ms"},
        {15, 12, 16, 14, 18, 10, 10, 13, 14, 16}, aggregate_rows);
    report.save("multilevel_comparison");
    return 0;
}
