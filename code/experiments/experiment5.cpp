#include "experiment/study.hpp"
#include "multigrid/algebraic_interpolation.hpp"

#include <array>
#include <string>
#include <utility>

int main(int argc, char** argv) {
    const auto config = experiment_support::parse_config(argc, argv);
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);
    experiment_support::Rows rows;
    // These are global-F-system iteration budgets.  The larger values make
    // the approach to the global energy-minimum target visible.
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

        // P_G is the common initial guess; every finite-step row targets the
        // same global F-point equation AP=0 and is evaluated against the
        // separately computed global reference basis.
        experiment_support::StudyCandidate initial{
            "geometric", "P_G", geometric.prolongation,
            geometric_ms, 0.0};
        rows.push_back(experiment_support::evaluate_candidate(
            field.name, a, rhs, initial, config));

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
                "Jacobi-global", "m=" + std::to_string(count),
                std::move(interpolation.prolongation),
                geometric_ms + interpolation.report.build_ms,
                static_cast<double>(count)};
            rows.push_back(experiment_support::evaluate_candidate(
                field.name, a, rhs, candidate, config));
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
        "Global-target iterative construction: Jacobi and PCG");
    report.add_summary(experiment_support::fixed_study_summary(
        config, "Initial interpolation", "P_G; no pruning or row cap"));
    report.add_note(
        "Jacobi-global and PCG-global use the same geometric initial basis "
        "and target the global F-point equation AP=0. The step budgets are "
        "4, 8, 16, 32 and 64; no magnitude pruning is applied. The exact "
        "global-energy basis (PCG tolerance 1e-10) is included as the limit "
        "for comparison.");
    report.add_table(
        "Global-target fixed-step construction", experiment_support::study_headers(),
        experiment_support::study_widths(), rows, true);
    report.save("experiment5");
    experiment_support::write_csv(
        "experiment5", experiment_support::study_headers(), rows);
    return 0;
}
