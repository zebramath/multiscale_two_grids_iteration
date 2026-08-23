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
    options.maximum_screening_steps = 48;
    options.screening_increment = 16;
    options.screening_pilot_iterations = 32;
    options.screening_tail_window = 8;
    options.screening_patience = 2;
    options.refinement_backtrack_steps = 10;
    options.refinement_stop_before_anchor_steps = 4;
    options.refinement_increment = 2;
    options.refinement_pilot_iterations = 64;
    options.refinement_tail_window = 16;
    options.confirmation_candidates = 2;
    options.easy_accept_cycles = 48;
    options.medium_accept_cycles = 64;
    options.maximum_confirmation_cycles = config.max_cycles;
    options.thread_count = config.threads;
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
        "m", "phase", "pilot", "rho_tail", "pilot residual",
        "pred cycles", "confirmed", "P density %", "path ms",
        "pilot ms", "confirm ms", "selected"};
    experiment_support::Rows history;
    for (const auto& item : adaptive.report.history) {
        history.push_back({
            std::to_string(item.steps),
            item.phase,
            std::to_string(item.pilot_iterations),
            experiment_support::fixed(item.rho_rhs_pilot, 6),
            experiment_support::scientific(item.pilot_relative_residual, 3),
            std::to_string(item.predicted_cycles),
            item.confirmed_cycles >= 0
                ? std::to_string(item.confirmed_cycles) : "-",
            experiment_support::fixed(item.density_percent, 4),
            experiment_support::fixed(item.path_ms),
            experiment_support::fixed(item.pilot_ms),
            experiment_support::fixed(item.confirmation_ms),
            item.selected ? "yes" : "no"});
    }

    experiment_support::Report report(
        "Adaptive finite-PCG checkpoint selection");
    report.add_summary(experiment_support::fixed_study_summary(
        config, "Selector", "screen-refine-confirm"));
    report.add_note(
        "The selector first screens coarse PCG checkpoints with 32 cycles, "
        "backtracks from the best screen anchor, refines at spacing two with "
        "64-cycle pilots, and continues only two diverse finalists to the "
        "actual solve tolerance. The final objective is measured cycles on "
        "the representative right-hand side.");
    report.add_table(
        "End-to-end verification", experiment_support::study_headers(),
        experiment_support::study_widths(), rows);
    report.add_table(
        "Adaptive decision trace", history_headers,
        {5, 9, 7, 11, 15, 11, 11, 11, 10, 10, 11, 9}, history);
    report.save("experiment7");
    experiment_support::write_csv(
        "experiment7", experiment_support::study_headers(), rows);
    experiment_support::write_csv(
        "experiment7_history", history_headers, history);
    return 0;
}
