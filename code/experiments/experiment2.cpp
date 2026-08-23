#include "experiment/study.hpp"

#include <array>

int main(int argc, char** argv) {
    const auto config = experiment_support::parse_config(argc, argv);
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);
    experiment_support::Rows rows;
    constexpr std::array<double, 6> tolerances{
        1.0e-2, 3.0e-3, 1.0e-3, 3.0e-4, 1.0e-4, 1.0e-10};

    for (const auto& field : experiment_support::standard_fields()) {
        const auto problem = experiment_support::make_problem(
            grid, field, config);
        const auto& a = problem.matrix;
        const auto& rhs = problem.rhs;
        for (double tolerance : tolerances) {
            auto options = experiment_support::energy_options(
                4, config.threads, tolerance);
            options.drop_tolerance = 0.0;
            auto candidate = experiment_support::make_candidate(
                "local-energy-4",
                "tol=" + experiment_support::scientific(tolerance, 0),
                tgi::build_interpolation(grid, a, options));
            rows.push_back(experiment_support::evaluate_candidate(
                field.name, a, rhs, candidate, config));
        }
    }

    experiment_support::Report report(
        "Algebraic localization: PCG tolerance on a fixed local support");
    report.add_summary(experiment_support::fixed_study_summary(
        config, "Fixed support", "four coarse patch layers"));
    report.add_note(
        "The support is identical in every row. The 1e-10 case is the "
        "reference-accurate local minimizer, not a direct solve.");
    report.add_table(
        "PCG tolerance study", experiment_support::study_headers(),
        experiment_support::study_widths(), rows, true);
    report.save("experiment2");
    experiment_support::write_csv(
        "experiment2", experiment_support::study_headers(), rows);
    return 0;
}

