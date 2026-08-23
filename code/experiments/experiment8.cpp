#include "experiment/study.hpp"
#include "multigrid/frontier_gain_support.hpp"
#include "multigrid/residual_budget_support.hpp"

#include <array>
#include <string>

namespace {

experiment_support::StudyCandidate old_budget(
    const tgi::StructuredGrid& grid, const tgi::SparseMatrix& a,
    const tgi::InterpolationResult& local, int threads) {
    tgi::ResidualBudgetSupportOptions options;
    options.base_patch_layers = 2;
    options.maximum_rounds = 8;
    options.maximum_extra_nodes_per_column = 128;
    options.maximum_nodes_per_round = 16;
    options.strength_scaling = tgi::StrengthScaling::SymmetricDiagonal;
    options.strong_edge_fraction = 0.10;
    options.thread_count = threads;
    const auto result = tgi::build_residual_budget_interpolation(
        grid, a, local.prolongation,
        experiment_support::energy_options(2, threads), options);
    return {"adaptive-budget-v2", "R=8,B=128",
            result.prolongation,
            local.report.timing.total_ms + result.report.total_ms};
}

std::pair<experiment_support::StudyCandidate, tgi::FrontierGainSupportReport>
frontier_gain(
    const tgi::StructuredGrid& grid, const tgi::SparseMatrix& a,
    const tgi::InterpolationResult& local, int threads) {
    tgi::FrontierGainSupportOptions options;
    options.base_patch_layers = 2;
    options.maximum_rounds = 10;
    options.maximum_extra_nodes_per_column = 128;
    options.maximum_total_nodes_per_round =
        std::max(128, 4 * grid.coarse_size());
    options.maximum_nodes_per_column_per_round = 16;
    options.maximum_stagnant_rounds = 3;
    options.minimum_round_relative_gain = 0.002;
    options.column_marking_fraction = 0.80;
    options.candidate_marking_fraction = 0.90;
    options.strength_scaling = tgi::StrengthScaling::SymmetricDiagonal;
    options.strong_edge_fraction = 0.10;
    options.thread_count = threads;
    const auto result = tgi::build_frontier_gain_interpolation(
        grid, a, local.prolongation,
        experiment_support::energy_options(2, threads), options);
    experiment_support::StudyCandidate candidate{
        "frontier-gain-v3", "R<=10,B=128",
        result.prolongation,
        local.report.timing.total_ms + result.report.total_ms};
    return {std::move(candidate), result.report};
}

} // namespace

int main(int argc, char** argv) {
    const auto config = experiment_support::parse_config(argc, argv);
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);
    const std::array<experiment_support::FieldCase, 4> fields{
        experiment_support::standard_fields().at(2),
        experiment_support::channel_topologies().at(0),
        experiment_support::channel_topologies().at(1),
        experiment_support::channel_topologies().at(2)};
    experiment_support::Rows rows;
    const experiment_support::Row detail_headers{
        "Field", "Rounds", "Expanded cols", "Frozen cols",
        "Extra nodes", "Final/initial", "Build ms"};
    experiment_support::Rows details;

    for (const auto& field : fields) {
        const auto problem = experiment_support::make_problem(
            grid, field, config, 17);
        auto local_options = experiment_support::energy_options(
            2, config.threads, 1.0e-6);
        local_options.drop_tolerance = 0.0;
        const auto local = tgi::build_interpolation(
            grid, problem.matrix, local_options);
        const auto base = experiment_support::make_candidate(
            "base-local", "layers=2", local);
        rows.push_back(experiment_support::evaluate_candidate(
            field.name, problem.matrix, problem.rhs, base, config));
        const auto old = old_budget(
            grid, problem.matrix, local, config.threads);
        rows.push_back(experiment_support::evaluate_candidate(
            field.name, problem.matrix, problem.rhs, old, config));
        auto [frontier, support_report] = frontier_gain(
            grid, problem.matrix, local, config.threads);
        rows.push_back(experiment_support::evaluate_candidate(
            field.name, problem.matrix, problem.rhs, frontier, config));
        details.push_back({
            field.name, std::to_string(support_report.rounds),
            std::to_string(support_report.expanded_columns),
            std::to_string(support_report.frozen_columns),
            std::to_string(support_report.total_extra_nodes),
            experiment_support::fixed(
                support_report.mean_final_to_initial_ratio, 4),
            experiment_support::fixed(support_report.total_ms)});
    }

    experiment_support::Report report(
        "Sparse support expansion: global frontier gain allocation");
    report.add_summary(experiment_support::fixed_study_summary(
        config, "Coefficient seed", "17"));
    report.add_note(
        "Frontier-gain-v3 ranks connected strong-graph frontier nodes across "
        "all marked columns, allocates a global round budget, warm-starts only "
        "changed columns, and freezes columns after three low-gain rounds.");
    report.add_table(
        "Solver comparison", experiment_support::study_headers(),
        experiment_support::study_widths(), rows, true);
    report.add_table(
        "Frontier allocation diagnostics", detail_headers,
        {20, 8, 14, 12, 12, 14, 11}, details);
    report.save("experiment8");
    experiment_support::write_csv(
        "experiment8", experiment_support::study_headers(), rows);
    experiment_support::write_csv(
        "experiment8_support", detail_headers, details);
    return 0;
}
