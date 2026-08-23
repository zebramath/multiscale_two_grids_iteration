#include "experiment/study.hpp"
#include "multigrid/residual_budget_support.hpp"
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
    int changed_count = 0;
    for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
        if (support.supports[static_cast<std::size_t>(coarse)].size() >
            grid.patch_f_nodes(coarse, 2).size()) {
            changed[static_cast<std::size_t>(coarse)] = 1U;
            ++changed_count;
        }
    }
    auto options = experiment_support::energy_options(2, threads, 1.0e-6);
    options.drop_tolerance = 0.0;
    auto interpolation = tgi::refine_selected_energy_interpolation_on_supports(
        grid, a, support.supports, changed, local2.prolongation, options);
    const int total_iterations =
        interpolation.report.local_solves.total_iterations;
    auto candidate = experiment_support::make_candidate(
        "residual-strong", "K=64 one-shot", std::move(interpolation));
    candidate.mean_construction_iterations = changed_count > 0
        ? static_cast<double>(total_iterations) /
              static_cast<double>(changed_count)
        : 0.0;
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
        "residual-budget", "R=8, B=128, q=16",
        result.prolongation, build_ms, 0.0};
}

} // namespace

int main(int argc, char** argv) {
    experiment_support::BasicConfig config;
    for (int index = 1; index < argc; ++index) {
        experiment_support::parse_basic_argument(config, argv[index]);
    }
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);
    experiment_support::Rows rows;

    for (const auto& field : experiment_support::standard_fields()) {
        const auto coefficient = experiment_support::make_field(
            grid, field, config.contrast);
        const tgi::SparseMatrix a = tgi::assemble_diffusion(
            grid, coefficient.values);
        const tgi::Vector rhs = a.multiply(
            experiment_support::manufactured_solution(grid));
        const auto global = experiment_support::build_global_reference(
            grid, a, config.threads);

        auto local2_options = experiment_support::energy_options(
            2, config.threads, 1.0e-6);
        local2_options.drop_tolerance = 0.0;
        auto local2 = tgi::build_interpolation(grid, a, local2_options);
        auto local2_candidate = experiment_support::make_candidate(
            "base-local", "layers=2", local2);
        rows.push_back(experiment_support::evaluate_candidate(
            field.name, grid, a, rhs, global.prolongation,
            local2_candidate, config));

        auto local3_options = experiment_support::energy_options(
            3, config.threads, 1.0e-6);
        local3_options.drop_tolerance = 0.0;
        auto layer_candidate = experiment_support::make_candidate(
            "geometric-layer", "2 -> 3",
            tgi::build_interpolation(grid, a, local3_options));
        rows.push_back(experiment_support::evaluate_candidate(
            field.name, grid, a, rhs, global.prolongation,
            layer_candidate, config));

        auto strong_candidate = build_fixed_strong_candidate(
            grid, a, local2, config.threads);
        rows.push_back(experiment_support::evaluate_candidate(
            field.name, grid, a, rhs, global.prolongation,
            strong_candidate, config));

        auto residual_budget_candidate = build_residual_budget_candidate(
            grid, a, local2, config.threads);
        rows.push_back(experiment_support::evaluate_candidate(
            field.name, grid, a, rhs, global.prolongation,
            residual_budget_candidate, config));

        tgi::StrengthDistanceOptions distance_options;
        distance_options.coarse_candidates_per_row = 8;
        distance_options.local_tolerance = 1.0e-6;
        distance_options.local_max_iterations = 20000;
        distance_options.thread_count = config.threads;
        const auto distance = tgi::build_strength_distance_interpolation(
            grid, a, distance_options);
        experiment_support::StudyCandidate distance_candidate{
            "strength-distance", "q=8", distance.prolongation,
            distance.report.build_ms,
            distance.report.mean_construction_iterations};
        rows.push_back(experiment_support::evaluate_candidate(
            field.name, grid, a, rhs, global.prolongation,
            distance_candidate, config));
    }

    experiment_support::Report report(
        "Support expansion strategies: fixed residual and residual-budget");
    report.add_summary(experiment_support::fixed_study_summary(
        config, "Energy solve tolerance", "1e-6"));
    report.add_note(
        "The base, geometric-layer, fixed residual-strong, and strength-distance "
        "rows use one construction pass. Residual-budget is the retained "
        "error-driven variant: it marks 70% of the current residual energy, "
        "adds at most 16 strong-graph nodes per round and 128 per column, and "
        "re-solves only marked columns. The reported build time includes all "
        "support selection and refinement rounds.");
    report.add_table(
        "Fixed support strategy study", experiment_support::study_headers(),
        experiment_support::study_widths(), rows, true);
    report.save("support_strategies");
    experiment_support::write_csv(
        "support_strategies", experiment_support::study_headers(), rows);
    return 0;
}
