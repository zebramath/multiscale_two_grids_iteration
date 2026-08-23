#include "experiment/study.hpp"
#include "multigrid/reference_pruning.hpp"

#include <array>
#include <string>

int main(int argc, char** argv) {
    experiment_support::BasicConfig config;
    for (int index = 1; index < argc; ++index) {
        experiment_support::parse_basic_argument(config, argv[index]);
    }
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);
    experiment_support::Rows rows;
    constexpr std::array<int, 5> steps{1, 2, 4, 8, 16};
    constexpr double post_prune_threshold = 1.0e-2;

    for (const auto& field : experiment_support::standard_fields()) {
        const auto coefficient = experiment_support::make_field(
            grid, field, config.contrast);
        const tgi::SparseMatrix a = tgi::assemble_diffusion(
            grid, coefficient.values);
        const tgi::Vector rhs = a.multiply(
            experiment_support::manufactured_solution(grid));
        auto global = experiment_support::build_global_reference(
            grid, a, config.threads);
        auto geometric = experiment_support::geometric_interpolation(grid, a);
        const double geometric_ms = geometric.report.timing.total_ms;

        experiment_support::StudyCandidate initial{
            "initial-geometric", "m=0", geometric.prolongation,
            geometric_ms, 0.0};
        rows.push_back(experiment_support::evaluate_candidate(
            field.name, grid, a, rhs, global.prolongation,
            initial, config));

        for (int count : steps) {
            tgi::JacobiInterpolationOptions options;
            options.steps = count;
            options.damping = 2.0 / 3.0;
            options.relative_drop_tolerance = 0.0;
            options.maximum_entries_per_row = 0;
            options.thread_count = config.threads;
            auto interpolation = tgi::build_jacobi_interpolation(
                grid, a, geometric.prolongation, options);
            const double build_ms = geometric_ms + interpolation.report.build_ms;
            experiment_support::StudyCandidate candidate{
                "Jacobi", "m=" + std::to_string(count),
                interpolation.prolongation, build_ms,
                static_cast<double>(count)};
            rows.push_back(experiment_support::evaluate_candidate(
                field.name, grid, a, rhs, global.prolongation,
                candidate, config));

            // Apply the same relative, column-wise post-pruning rule to both
            // iterative constructions.  This isolates sparsification from
            // the iteration count; the full global-pruning sweep studies the
            // threshold itself.
            const auto pruned = tgi::prune_global_interpolation_relative(
                grid, interpolation.prolongation, post_prune_threshold);
            experiment_support::StudyCandidate sparse_candidate{
                "Jacobi", "m=" + std::to_string(count) +
                    ", drop=1e-2",
                pruned.prolongation,
                build_ms + pruned.pruning_ms,
                static_cast<double>(count)};
            rows.push_back(experiment_support::evaluate_candidate(
                field.name, grid, a, rhs, global.prolongation,
                sparse_candidate, config));
        }

        for (int count : steps) {
            auto options = experiment_support::energy_options(
                0, config.threads, 1.0e-300);
            options.local_max_iterations = count;
            options.require_convergence = false;
            options.drop_tolerance = 0.0;
            auto interpolation = tgi::refine_global_energy_interpolation(
                grid, a, geometric.prolongation, options);
            auto candidate = experiment_support::make_candidate(
                "PCG", "m=" + std::to_string(count),
                std::move(interpolation));
            candidate.build_ms += geometric_ms;
            const tgi::SparseMatrix pcg_prolongation = candidate.prolongation;
            rows.push_back(experiment_support::evaluate_candidate(
                field.name, grid, a, rhs, global.prolongation,
                candidate, config));

            const auto pruned = tgi::prune_global_interpolation_relative(
                grid, pcg_prolongation, post_prune_threshold);
            experiment_support::StudyCandidate sparse_candidate{
                "PCG", "m=" + std::to_string(count) +
                    ", drop=1e-2",
                pruned.prolongation,
                candidate.build_ms + pruned.pruning_ms,
                candidate.mean_construction_iterations};
            rows.push_back(experiment_support::evaluate_candidate(
                field.name, grid, a, rhs, global.prolongation,
                sparse_candidate, config));
        }

        auto reference = experiment_support::make_candidate(
            "PCG-reference", "tol=1e-10", std::move(global));
        rows.push_back(experiment_support::evaluate_candidate(
            field.name, grid, a, rhs, reference.prolongation,
            reference, config));
    }

    experiment_support::Report report(
        "Iterative localization: Jacobi/PCG steps with and without pruning");
    report.add_summary(experiment_support::fixed_study_summary(
        config, "Initial interpolation", "geometric; iterative rows start unpruned"));
    report.add_note(
        "Jacobi and diagonal-PCG act on the same global F system and start "
        "from the same geometric interpolation. With diagonal preconditioning, "
        "one matrix action expands graph support by at most one edge. The "
        "strict PCG reference is shown only as the limiting target. For every "
        "step count, a second row applies the same post-construction relative "
        "column pruning (drop=1e-2) to expose the setup/sparsity/convergence "
        "trade-off without changing the iteration itself.");
    report.add_table(
        "Fixed-step iterative construction", experiment_support::study_headers(),
        experiment_support::study_widths(), rows, true);
    report.save("iterative_construction");
    experiment_support::write_csv(
        "iterative_construction", experiment_support::study_headers(), rows);
    return 0;
}
