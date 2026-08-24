#include "experiment/study.hpp"

#include <string>
#include <utility>

int run_pcg_dense_scan(int argc, char** argv) {
    const auto config = experiment_support::parse_config(argc, argv);
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);
    const auto& channel = experiment_support::standard_fields().at(1);
    const auto problem = experiment_support::make_problem(
        grid, channel, config);
    const auto& a = problem.matrix;
    const auto& rhs = problem.rhs;
    experiment_support::Rows rows;

    auto geometric = experiment_support::geometric_interpolation(grid, a);
    const double geometric_ms = geometric.report.timing.total_ms;
    for (int steps = 16; steps <= 64; steps += 2) {
        auto options = experiment_support::energy_options(
            0, config.threads, 0.0);
        options.local_max_iterations = steps;
        options.require_convergence = false;
        options.drop_tolerance = 0.0;
        auto interpolation = tgi::refine_global_energy_interpolation(
            grid, a, geometric.prolongation, options);
        auto candidate = experiment_support::make_candidate(
            "PCG-global", "m=" + std::to_string(steps),
            std::move(interpolation));
        candidate.build_ms += geometric_ms;
        rows.push_back(experiment_support::evaluate_candidate(
            channel.name, a, rhs, candidate, config));
    }

    auto global = experiment_support::build_global_reference(
        grid, a, config.threads);
    auto exact = experiment_support::make_candidate(
        "global-exact", "tol=1e-10", std::move(global));
    rows.push_back(experiment_support::evaluate_candidate(
        channel.name, a, rhs, exact, config));

    experiment_support::Report report(
        "Channel finite PCG-step scan");
    report.add_summary(experiment_support::fixed_study_summary(
        config, "PCG scan", "m=16,18,...,64"));
    report.add_note(
        "Each finite-step interpolation is evaluated by an independent "
        "two-grid solve. Setup ms contains interpolation construction and "
        "the production sparse Galerkin hierarchy.");
    report.add_table(
        "Channel finite-step PCG scan", experiment_support::study_headers(),
        experiment_support::study_widths(), rows);
    report.save("pcg_dense_scan");
    return 0;
}
