#include "experiment/study.hpp"

#include <array>
#include <string>

int main(int argc, char** argv) {
    experiment_support::BasicConfig config;
    for (int index = 1; index < argc; ++index) {
        experiment_support::parse_basic_argument(config, argv[index]);
    }
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);
    experiment_support::Rows rows;
    constexpr std::array<int, 6> steps{1, 2, 4, 8, 16, 32};

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
            experiment_support::StudyCandidate candidate{
                "Jacobi", "m=" + std::to_string(count),
                std::move(interpolation.prolongation),
                geometric_ms + interpolation.report.build_ms,
                static_cast<double>(count)};
            rows.push_back(experiment_support::evaluate_candidate(
                field.name, grid, a, rhs, global.prolongation,
                candidate, config));
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
            rows.push_back(experiment_support::evaluate_candidate(
                field.name, grid, a, rhs, global.prolongation,
                candidate, config));
        }

        auto reference = experiment_support::make_candidate(
            "PCG-reference", "tol=1e-10", std::move(global));
        rows.push_back(experiment_support::evaluate_candidate(
            field.name, grid, a, rhs, reference.prolongation,
            reference, config));
    }

    experiment_support::Report report(
        "Iterative localization: fixed-step Jacobi and PCG construction");
    report.add_summary(experiment_support::fixed_study_summary(
        config, "Initial interpolation", "geometric; no pruning or row cap"));
    report.add_note(
        "Jacobi and diagonal-PCG act on the same global F system and start "
        "from the same geometric interpolation. With diagonal preconditioning, "
        "one matrix action expands graph support by at most one edge. The "
        "strict PCG reference is shown only as the limiting target.");
    report.add_table(
        "Fixed-step iterative construction", experiment_support::study_headers(),
        experiment_support::study_widths(), rows, true);
    report.save("iterative_construction");
    experiment_support::write_csv(
        "iterative_construction", experiment_support::study_headers(), rows);
    return 0;
}
