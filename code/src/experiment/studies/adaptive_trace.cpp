#include "experiment/study.hpp"
#include "multigrid/global_pcg.hpp"

#include <array>
#include <string>

int run_adaptive_trace(int argc, char** argv) {
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
    constexpr std::array<int, 7> checkpoints{16, 24, 32, 36, 38, 40, 48};
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
    options.maximum_cycles = config.max_cycles;
    options.thread_count = config.threads;
    const auto adaptive = tgi::build_adaptive_global_pcg_interpolation(
        grid, problem.matrix, geometric.prolongation, options,
        &problem.rhs);
    rows.push_back(experiment_support::evaluate_existing_candidate(
        problem.field_name, "PCG-adaptive",
        "selected m=" + std::to_string(adaptive.report.selected_steps),
        problem.rhs, *adaptive.prolongation, *adaptive.cycle,
        geometric_ms + adaptive.report.selection_wall_ms, config));

    const experiment_support::Row history_headers{
        "m", "phase", "pilot", "rho_tail", "uncertainty",
        "pilot residual", "PCG energy residual", "pred cycles",
        "P density %", "path ms", "pilot ms", "selected"};
    experiment_support::Rows history;
    for (const auto& item : adaptive.report.history) {
        history.push_back({
            std::to_string(item.steps),
            item.phase,
            std::to_string(item.pilot_iterations),
            experiment_support::fixed(item.rho_rhs_pilot, 6),
            experiment_support::fixed(
                item.forecast_relative_uncertainty, 4),
            experiment_support::scientific(item.pilot_relative_residual, 3),
            experiment_support::scientific(
                item.preconditioned_pcg_residual, 3),
            std::to_string(item.predicted_cycles),
            experiment_support::fixed(item.density_percent, 4),
            experiment_support::fixed(item.path_ms),
            experiment_support::fixed(item.pilot_ms),
            item.selected ? "yes" : "no"});
    }

    experiment_support::Report report(
        "Adaptive finite-PCG checkpoint selection");
    report.add_summary(experiment_support::fixed_study_summary(
        config, "Selector", "midpoint adaptive screen"));
    report.add_note(
        "The selector evaluates the geometric basis, a minimum-step anchor "
        "and an interval midpoint. If the midpoint remains hard, its "
        "normalized PCG energy residual chooses one interior or forward "
        "checkpoint. No full confirmation "
        "solve is used during setup, and the selected hierarchy is reused.");
    report.add_table(
        "End-to-end verification", experiment_support::study_headers(),
        experiment_support::study_widths(), rows);
    report.add_table(
        "Adaptive decision trace", history_headers,
        {5, 10, 7, 11, 12, 15, 18, 11, 11, 10, 10, 9}, history);
    report.save("adaptive_trace");
    return 0;
}
