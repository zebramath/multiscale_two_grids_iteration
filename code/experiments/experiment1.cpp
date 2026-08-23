#include "experiment/study.hpp"

#include <string>

int main(int argc, char** argv) {
    const auto config = experiment_support::parse_config(argc, argv);
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);
    experiment_support::Rows rows;

    for (const auto& field : experiment_support::standard_fields()) {
        const auto problem = experiment_support::make_problem(
            grid, field, config);
        const auto& a = problem.matrix;
        const auto& rhs = problem.rhs;
        auto global = experiment_support::build_global_reference(
            grid, a, config.threads);

        auto geometric = experiment_support::make_candidate(
            "geometric", "P_G",
            experiment_support::geometric_interpolation(grid, a));
        rows.push_back(experiment_support::evaluate_candidate(
            field.name, a, rhs, geometric, config));

        for (int layers : {2, 3, 4}) {
            auto options = experiment_support::energy_options(
                layers, config.threads, 1.0e-10);
            options.drop_tolerance = 0.0;
            auto candidate = experiment_support::make_candidate(
                "local-energy", "layers=" + std::to_string(layers),
                tgi::build_interpolation(grid, a, options));
            rows.push_back(experiment_support::evaluate_candidate(
                field.name, a, rhs, candidate, config));
        }
        auto global_candidate = experiment_support::make_candidate(
            "global-energy", "layers=inf", std::move(global));
        rows.push_back(experiment_support::evaluate_candidate(
            field.name, a, rhs, global_candidate, config));
    }

    experiment_support::Report report(
        "Spatial localization: P_G, local support radius and global reference");
    report.add_summary(experiment_support::fixed_study_summary(
        config, "Energy solve tolerance", "1e-10 for every support"));
    report.add_note(
        "P_G is the pure geometric bilinear baseline. The energy rows only "
        "change support radius; their common PCG tolerance makes algebraic "
        "stopping error negligible relative to localization.");
    report.add_table(
        "Support radius study", experiment_support::study_headers(),
        experiment_support::study_widths(), rows, true);
    report.save("experiment1");
    experiment_support::write_csv(
        "experiment1", experiment_support::study_headers(), rows);
    return 0;
}

