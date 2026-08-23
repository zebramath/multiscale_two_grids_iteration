#include "experiment/study.hpp"
#include "multigrid/adaptive_global_pcg.hpp"
#include "multigrid/global_pcg_path.hpp"

#include <array>
#include <string>

int main(int argc, char** argv) {
    const auto config = experiment_support::parse_config(argc, argv);
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);
    const auto problem = experiment_support::make_problem(
        grid, experiment_support::channel_topologies().front(), config, 1);
    auto geometric = experiment_support::geometric_interpolation(
        grid, problem.matrix);
    const double geometric_ms = geometric.report.timing.total_ms;
    experiment_support::Rows rows;

    tgi::GlobalEnergyPcgPath path(
        grid, problem.matrix, geometric.prolongation, config.threads);
    constexpr std::array<int, 7> checkpoints{16, 24, 32, 40, 48, 56, 64};
    for (int steps : checkpoints) {
        path.advance_to(steps);
        const auto path_report = path.report();
        experiment_support::StudyCandidate candidate{
            "PCG-fixed", "m=" + std::to_string(steps),
            path.prolongation(0.0), geometric_ms + path_report.total_ms};
        rows.push_back(experiment_support::evaluate_candidate(
            problem.field_name, problem.matrix, problem.rhs,
            candidate, config));
    }

    tgi::AdaptiveGlobalPcgOptions options;
    options.minimum_steps = 16;
    options.maximum_steps = 64;
    options.step_increment = 4;
    options.patience = 4;
    options.probe_count = 3;
    options.power_iterations = 30;
    options.rhs_pilot_iterations = 40;
    options.rhs_tail_window = 10;
    options.thread_count = config.threads;
    options.expected_rhs = 8;
    const auto adaptive = tgi::build_adaptive_global_pcg_interpolation(
        grid, problem.matrix, geometric.prolongation, options,
        &problem.rhs);
    experiment_support::StudyCandidate adaptive_candidate{
        "PCG-adaptive",
        "selected m=" + std::to_string(adaptive.report.selected_steps),
        adaptive.prolongation,
        geometric_ms + adaptive.report.selection_wall_ms};
    rows.push_back(experiment_support::evaluate_candidate(
        problem.field_name, problem.matrix, problem.rhs,
        adaptive_candidate, config));

    const experiment_support::Row history_headers{
        "m", "rho_hat", "rho_power", "rho_rhs", "pred cycles", "P density %",
        "path ms", "probe ms", "score ms", "best"};
    experiment_support::Rows history;
    for (const auto& item : adaptive.report.history) {
        history.push_back({
            std::to_string(item.steps),
            experiment_support::fixed(item.rho_hat, 6),
            experiment_support::fixed(item.rho_power, 6),
            experiment_support::fixed(item.rho_rhs_pilot, 6),
            std::to_string(item.predicted_cycles),
            experiment_support::fixed(item.density_percent, 4),
            experiment_support::fixed(item.path_ms),
            experiment_support::fixed(item.probe_ms),
            experiment_support::fixed(item.predicted_total_ms),
            item.improved ? "yes" : "no"});
    }

    experiment_support::Report report(
        "Adaptive finite-PCG checkpoint selection");
    report.add_summary(experiment_support::fixed_study_summary(
        config, "Workload used by selector", "8 right-hand sides"));
    report.add_note(
        "The fixed rows and the adaptive search continue one PCG trajectory "
        "between checkpoints. The selector uses the maximum of three "
        "independent 30-step energy-norm power probes and a 40-cycle "
        "representative-RHS tail probe, keeps the historical "
        "best candidate (including m=0), and stops after four non-improving "
        "checkpoints.");
    report.add_table(
        "End-to-end verification", experiment_support::study_headers(),
        experiment_support::study_widths(), rows);
    report.add_table(
        "Adaptive decision trace", history_headers,
        {5, 11, 11, 11, 11, 11, 10, 10, 12, 6}, history);
    report.save("experiment7");
    experiment_support::write_csv(
        "experiment7", experiment_support::study_headers(), rows);
    experiment_support::write_csv(
        "experiment7_history", history_headers, history);
    return 0;
}
