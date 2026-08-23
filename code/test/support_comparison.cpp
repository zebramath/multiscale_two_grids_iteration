#include "experiment/common.hpp"
#include "multigrid/adaptive_support.hpp"
#include "multigrid/reference_pruning.hpp"
#include "multigrid/residual_budget_support.hpp"

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct Candidate {
    std::string name;
    tgi::SparseMatrix prolongation;
    double build_ms = 0.0;
};

Candidate fixed_strong_candidate(
    const tgi::StructuredGrid& grid, const tgi::SparseMatrix& a,
    const tgi::InterpolationResult& local2, int threads) {
    const auto begin = experiment_support::Clock::now();
    tgi::ResidualStrongSupportOptions support_options;
    support_options.base_patch_layers = 2;
    support_options.maximum_extra_nodes_per_column = 64;
    support_options.maximum_graph_hops = 4 * grid.fine_n();
    support_options.strong_edge_fraction = 0.25;
    support_options.thread_count = threads;
    const auto support = tgi::build_residual_strong_supports(
        grid, a, local2.prolongation, support_options);
    std::vector<unsigned char> changed(
        static_cast<std::size_t>(grid.coarse_size()), 0U);
    for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
        if (support.supports[static_cast<std::size_t>(coarse)].size() >
            grid.patch_f_nodes(coarse, 2).size()) {
            changed[static_cast<std::size_t>(coarse)] = 1U;
        }
    }
    auto interpolation_options =
        experiment_support::energy_options(2, threads, 1.0e-6);
    const auto interpolation =
        tgi::refine_selected_energy_interpolation_on_supports(
            grid, a, support.supports, changed,
            local2.prolongation, interpolation_options);
    return {
        "strong-K64", interpolation.prolongation,
        local2.report.timing.total_ms + experiment_support::milliseconds(
            begin, experiment_support::Clock::now())};
}

Candidate residual_budget_candidate(
    const tgi::StructuredGrid& grid, const tgi::SparseMatrix& a,
    const tgi::InterpolationResult& local2, int threads) {
    tgi::ResidualBudgetSupportOptions options;
    options.base_patch_layers = 2;
    options.maximum_rounds = 8;
    options.maximum_extra_nodes_per_column = 128;
    options.maximum_nodes_per_round = 16;
    options.marking_fraction = 0.70;
    options.target_residual_ratio = 0.25;
    options.strength_scaling = tgi::StrengthScaling::SymmetricDiagonal;
    options.strong_edge_fraction = 0.10;
    options.thread_count = threads;
    const auto result = tgi::build_residual_budget_interpolation(
        grid, a, local2.prolongation,
        experiment_support::energy_options(2, threads, 1.0e-6),
        options);
    return {
        "residual-budget", result.prolongation,
        local2.report.timing.total_ms + result.report.total_ms};
}

} // namespace

int main(int argc, char** argv) {
    experiment_support::BasicConfig config;
    config.fine_intervals = 128;
    config.coarse_intervals = 16;
    for (int index = 1; index < argc; ++index) {
        experiment_support::parse_basic_argument(config, argv[index]);
    }
    const tgi::StructuredGrid grid =
        experiment_support::make_grid(config);
    experiment_support::Rows rows;

    for (const auto& field : experiment_support::standard_fields()) {
        const auto coefficient = experiment_support::make_field(
            grid, field, config.contrast);
        const tgi::SparseMatrix a = tgi::assemble_diffusion(
            grid, coefficient.values);
        const tgi::Vector rhs = a.multiply(
            experiment_support::manufactured_solution(grid));

        const auto geometric =
            experiment_support::geometric_interpolation(grid, a);
        const auto local2 = tgi::build_interpolation(
            grid, a,
            experiment_support::energy_options(2, config.threads));
        const auto local3 = tgi::build_interpolation(
            grid, a,
            experiment_support::energy_options(3, config.threads));
        const auto global = tgi::build_interpolation(
            grid, a,
            experiment_support::energy_options(
                0, config.threads, 1.0e-10));

        std::vector<Candidate> candidates;
        candidates.push_back({
            "geometric", geometric.prolongation,
            geometric.report.timing.total_ms});
        candidates.push_back({
            "local-2", local2.prolongation,
            local2.report.timing.total_ms});
        candidates.push_back({
            "local-3", local3.prolongation,
            local3.report.timing.total_ms});
        candidates.push_back(fixed_strong_candidate(
            grid, a, local2, config.threads));
        candidates.push_back(residual_budget_candidate(
            grid, a, local2, config.threads));

        for (const Candidate& candidate : candidates) {
            const auto matched = tgi::build_budget_matched_reference(
                grid, a, global.prolongation, candidate.prolongation,
                experiment_support::energy_options(
                    2, config.threads, 1.0e-8));
            const auto global_error =
                experiment_support::compare_prolongations_global(
                    grid, a, global.prolongation,
                    candidate.prolongation);
            const auto matched_error =
                experiment_support::compare_prolongations_global(
                    grid, a, matched.prolongation,
                    candidate.prolongation);
            const auto cycles = experiment_support::evaluate_two_grid(
                a, rhs, candidate.prolongation, config.threads,
                1.0e-6, config.max_cycles, 20);
            const auto oracle_cycles = experiment_support::evaluate_two_grid(
                a, rhs, matched.prolongation, config.threads,
                1.0e-6, config.max_cycles, 20);
            rows.push_back({
                field.name,
                candidate.name,
                std::to_string(candidate.prolongation.nnz()),
                experiment_support::scientific(
                    global_error.aggregate_relative_energy_error),
                experiment_support::scientific(
                    matched_error.aggregate_relative_energy_error),
                cycles.converged
                    ? std::to_string(cycles.cycles)
                    : "failed@" + std::to_string(cycles.cycles),
                experiment_support::fixed(
                    cycles.convergence_factor, 4),
                oracle_cycles.converged
                    ? std::to_string(oracle_cycles.cycles)
                    : "failed@" + std::to_string(oracle_cycles.cycles),
                experiment_support::fixed(candidate.build_ms),
                experiment_support::fixed(matched.build_ms)
            });
        }
        const auto global_cycles = experiment_support::evaluate_two_grid(
            a, rhs, global.prolongation, config.threads,
            1.0e-6, config.max_cycles, 20);
        rows.push_back({
            field.name, "global-reference",
            std::to_string(global.prolongation.nnz()), "0.000e+00",
            "0.000e+00",
            global_cycles.converged
                ? std::to_string(global_cycles.cycles)
                : "failed@" + std::to_string(global_cycles.cycles),
            experiment_support::fixed(
                global_cycles.convergence_factor, 4),
            "-", experiment_support::fixed(global.report.timing.total_ms),
            "-"});
    }

    const experiment_support::Row headers{
        "Field", "Method", "P nnz", "Err to global",
        "Err to oracle", "Cycles", "Rho", "Oracle cycles",
        "Build ms", "Oracle ms"};
    experiment_support::Report report(
        "Support strategies with budget-matched global pruning");
    report.add_summary({
        {"Fine grid", "h=1/" + std::to_string(config.fine_intervals)},
        {"Coarse grid", "H=1/" + std::to_string(config.coarse_intervals)},
        {"Contrast", experiment_support::scientific(config.contrast, 0)},
        {"RHS", "one manufactured solution"}
    });
    report.add_note(
        "Oracle uses the same F-entry count in every column as the candidate, "
        "selects entries from the global basis by |p_i|sqrt(A_ii), and then "
        "re-minimizes energy. It is an offline comparator, not a practical "
        "setup algorithm.");
    report.add_table(
        "Unified support comparison", headers,
        {12, 18, 10, 14, 14, 12, 9, 14, 11, 11}, rows, true);
    report.save("support_comparison");
    experiment_support::write_csv("support_comparison", headers, rows);
    return 0;
}
