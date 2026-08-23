#include "experiment/study.hpp"

#include <array>
#include <string>
#include <utility>

int run_pcg_budget(int argc, char** argv) {
    const auto config = experiment_support::parse_config(argc, argv);
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);
    experiment_support::Rows rows;
    constexpr std::array<int, 5> steps{4, 8, 16, 32, 64};

    for (const auto& field : experiment_support::standard_fields()) {
        const auto problem = experiment_support::make_problem(
            grid, field, config);
        const auto& a = problem.matrix;
        const auto& rhs = problem.rhs;
        auto global = experiment_support::build_global_reference(
            grid, a, config.threads);
        auto geometric = experiment_support::geometric_interpolation(grid, a);
        const double geometric_ms = geometric.report.timing.total_ms;

        experiment_support::StudyCandidate initial{
            "geometric", "P_G", geometric.prolongation,
            geometric_ms};
        rows.push_back(experiment_support::evaluate_candidate(
            field.name, a, rhs, initial, config));

        for (int count : steps) {
            auto options = experiment_support::energy_options(
                0, config.threads, 0.0);
            options.local_max_iterations = count;
            options.require_convergence = false;
            options.drop_tolerance = 0.0;
            auto interpolation = tgi::refine_global_energy_interpolation(
                grid, a, geometric.prolongation, options);
            auto candidate = experiment_support::make_candidate(
                "PCG-global", "m=" + std::to_string(count),
                std::move(interpolation));
            candidate.build_ms += geometric_ms;
            rows.push_back(experiment_support::evaluate_candidate(
                field.name, a, rhs, candidate, config));
        }

        auto exact = experiment_support::make_candidate(
            "global-exact", "tol=1e-10", std::move(global));
        rows.push_back(experiment_support::evaluate_candidate(
            field.name, a, rhs, exact, config));
    }

    experiment_support::Report report(
        "Global-target finite-PCG construction");
    report.add_summary(experiment_support::fixed_study_summary(
        config, "Initial interpolation", "P_G; no pruning or row cap"));
    report.add_note(
        "Every finite candidate starts from the geometric basis and follows "
        "the global energy PCG path for 4, 8, 16, 32 or 64 steps. The exact "
        "global-energy basis is included as the limiting reference.");
    report.add_table(
        "Global-target finite-PCG construction", experiment_support::study_headers(),
        experiment_support::study_widths(), rows, true);
    report.save("pcg_budget");
    return 0;
}
