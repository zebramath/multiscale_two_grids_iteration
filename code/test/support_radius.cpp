#include "experiment/study.hpp"

#include <string>

int main(int argc, char** argv) {
    experiment_support::BasicConfig config;
    for (int index = 1; index < argc; ++index) {
        experiment_support::parse_basic_argument(config, argv[index]);
    }
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);
    experiment_support::Rows rows;

    for (const auto& field : experiment_support::standard_fields()) {
        const auto coefficient = experiment_support::make_field(
            grid, field, config.contrast);
        const tgi::SparseMatrix a = tgi::assemble_diffusion(
            grid, coefficient.values);
        const tgi::Vector rhs = a.multiply(
            experiment_support::manufactured_solution(grid));
        auto global = experiment_support::build_global_reference(
            grid, a, config.threads);

        auto geometric = experiment_support::make_candidate(
            "geometric", "P_G",
            experiment_support::geometric_interpolation(grid, a));
        rows.push_back(experiment_support::evaluate_candidate(
            field.name, grid, a, rhs, geometric, config));

        for (int layers : {2, 3, 4}) {
            auto options = experiment_support::energy_options(
                layers, config.threads, 1.0e-10);
            options.drop_tolerance = 0.0;
            auto candidate = experiment_support::make_candidate(
                "local-energy", "layers=" + std::to_string(layers),
                tgi::build_interpolation(grid, a, options));
            rows.push_back(experiment_support::evaluate_candidate(
                field.name, grid, a, rhs, candidate, config));
        }
        auto global_candidate = experiment_support::make_candidate(
            "global-energy", "layers=inf", std::move(global));
        rows.push_back(experiment_support::evaluate_candidate(
            field.name, grid, a, rhs, global_candidate, config));
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
    report.save("support_radius");
    experiment_support::write_csv(
        "support_radius", experiment_support::study_headers(), rows);
    return 0;
}
