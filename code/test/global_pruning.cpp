#include "experiment/study.hpp"
#include "multigrid/reference_pruning.hpp"

#include <array>

int main(int argc, char** argv) {
    experiment_support::BasicConfig config;
    for (int index = 1; index < argc; ++index) {
        experiment_support::parse_basic_argument(config, argv[index]);
    }
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);
    experiment_support::Rows rows;
    constexpr std::array<double, 9> thresholds{
        0.0, 1.0e-3, 3.0e-3, 1.0e-2, 2.0e-2,
        3.0e-2, 5.0e-2, 1.0e-1, 2.0e-1};

    for (const auto& field : experiment_support::standard_fields()) {
        const auto coefficient = experiment_support::make_field(
            grid, field, config.contrast);
        const tgi::SparseMatrix a = tgi::assemble_diffusion(
            grid, coefficient.values);
        const tgi::Vector rhs = a.multiply(
            experiment_support::manufactured_solution(grid));
        const auto global = experiment_support::build_global_reference(
            grid, a, config.threads);
        const int systems = global.report.local_solves.systems;
        const double mean_iterations = systems > 0
            ? static_cast<double>(global.report.local_solves.total_iterations) /
                  static_cast<double>(systems)
            : 0.0;

        for (double threshold : thresholds) {
            auto pruned = tgi::prune_global_interpolation_relative(
                grid, global.prolongation, threshold);
            experiment_support::StudyCandidate candidate{
                "global-pruned",
                "drop=" + experiment_support::scientific(threshold, 0),
                std::move(pruned.prolongation),
                global.report.timing.total_ms + pruned.pruning_ms,
                mean_iterations};
            rows.push_back(experiment_support::evaluate_candidate(
                field.name, grid, a, rhs, global.prolongation,
                candidate, config));
        }
    }

    experiment_support::Report report(
        "Sparse localization: relative magnitude pruning of the global basis");
    report.add_summary(experiment_support::fixed_study_summary(
        config, "Global solve tolerance", "1e-10"));
    report.add_note(
        "For each column, an F entry is dropped when its magnitude is at "
        "most delta times that column's maximum magnitude. C injection is "
        "always preserved. Build and total times include construction of the "
        "dense global reference, so this is a diagnostic oracle study.");
    report.add_table(
        "Global pruning study", experiment_support::study_headers(),
        experiment_support::study_widths(), rows, true);
    report.save("global_pruning");
    experiment_support::write_csv(
        "global_pruning", experiment_support::study_headers(), rows);
    return 0;
}
