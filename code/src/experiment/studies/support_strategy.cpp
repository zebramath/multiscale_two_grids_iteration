#include "experiment/study.hpp"
#include "multigrid/support_expansion.hpp"

#include <string>
#include <vector>

namespace {

experiment_support::StudyCandidate build_fixed_strong_candidate(
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
    auto options = experiment_support::energy_options(2, threads, 1.0e-6);
    options.drop_tolerance = 0.0;
    auto interpolation = tgi::refine_selected_energy_interpolation_on_supports(
        grid, a, support.supports, changed, local2.prolongation, options);
    auto candidate = experiment_support::make_candidate(
        "indicator-strong", "K=64 one-shot", std::move(interpolation));
    candidate.build_ms = local2.report.timing.total_ms +
        experiment_support::milliseconds(begin, experiment_support::Clock::now());
    return candidate;
}

experiment_support::StudyCandidate build_residual_budget_candidate(
    const tgi::StructuredGrid& grid, const tgi::SparseMatrix& a,
    const tgi::InterpolationResult& local2, int threads) {
    tgi::ResidualBudgetSupportOptions support_options;
    support_options.base_patch_layers = 2;
    support_options.maximum_rounds = 8;
    support_options.maximum_extra_nodes_per_column = 128;
    support_options.maximum_nodes_per_round = 16;
    support_options.marking_fraction = 0.70;
    support_options.target_residual_ratio = 0.25;
    support_options.refinement_tolerance = 1.0e-6;
    support_options.strength_scaling = tgi::StrengthScaling::SymmetricDiagonal;
    support_options.strong_edge_fraction = 0.10;
    support_options.thread_count = threads;

    const auto begin = experiment_support::Clock::now();
    const auto result = tgi::build_residual_budget_interpolation(
        grid, a, local2.prolongation, experiment_support::energy_options(
            2, threads, 1.0e-6), support_options);
    const double build_ms = local2.report.timing.total_ms +
        experiment_support::milliseconds(
            begin, experiment_support::Clock::now());
    return {
        "adaptive-budget", "R=8, B=128, q=16",
        result.prolongation, build_ms};
}

}

int run_support_strategy(int argc, char** argv) {
    const auto config = experiment_support::parse_config(argc, argv);
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);
    experiment_support::Rows rows;

    for (const auto& field : experiment_support::standard_fields()) {
        const auto problem = experiment_support::make_problem(
            grid, field, config);
        const auto& a = problem.matrix;
        const auto& rhs = problem.rhs;
        auto local2_options = experiment_support::energy_options(
            2, config.threads, 1.0e-6);
        local2_options.drop_tolerance = 0.0;
        auto local2 = tgi::build_interpolation(grid, a, local2_options);
        auto local2_candidate = experiment_support::make_candidate(
            "base-local", "layers=2", local2);
        rows.push_back(experiment_support::evaluate_candidate(
            field.name, a, rhs, local2_candidate, config));

        auto local3_options = experiment_support::energy_options(
            3, config.threads, 1.0e-6);
        local3_options.drop_tolerance = 0.0;
        auto layer_candidate = experiment_support::make_candidate(
            "geometric-layer", "2 -> 3",
            tgi::build_interpolation(grid, a, local3_options));
        rows.push_back(experiment_support::evaluate_candidate(
            field.name, a, rhs, layer_candidate, config));

        auto strong_candidate = build_fixed_strong_candidate(
            grid, a, local2, config.threads);
        rows.push_back(experiment_support::evaluate_candidate(
            field.name, a, rhs, strong_candidate, config));

        auto residual_budget_candidate = build_residual_budget_candidate(
            grid, a, local2, config.threads);
        rows.push_back(experiment_support::evaluate_candidate(
            field.name, a, rhs, residual_budget_candidate, config));

    }

    experiment_support::Report report(
        "Support expansion strategies: fixed and adaptive budgets");
    report.add_summary(experiment_support::fixed_study_summary(
        config, "Energy solve tolerance", "1e-6"));
    report.add_note(
        "The base, geometric-layer and indicator-strong rows use one "
        "construction pass. Adaptive-budget is the retained "
        "error-driven variant: it marks 70% of the current indicator energy, "
        "adds at most 16 strong-graph nodes per round and 128 per column, and "
        "re-solves only marked columns. The reported build time includes all "
        "support selection and refinement rounds.");
    report.add_table(
        "Fixed support strategy study", experiment_support::study_headers(),
        experiment_support::study_widths(), rows, true);
    report.save("support_strategy");
    return 0;
}
