#include "experiment/study.hpp"

#include <array>

int main(int argc, char** argv) {
    experiment_support::BasicConfig config;
    for (int index = 1; index < argc; ++index) {
        experiment_support::parse_basic_argument(config, argv[index]);
    }
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);
    experiment_support::Rows rows;
    constexpr std::array<double, 6> tolerances{
        1.0e-2, 3.0e-3, 1.0e-3, 3.0e-4, 1.0e-4, 1.0e-10};

    for (const auto& field : experiment_support::standard_fields()) {
        const auto coefficient = experiment_support::make_field(
            grid, field, config.contrast);
        const tgi::SparseMatrix a = tgi::assemble_diffusion(
            grid, coefficient.values);
        const tgi::Vector rhs = a.multiply(
            experiment_support::manufactured_solution(grid));
        for (double tolerance : tolerances) {
            auto options = experiment_support::energy_options(
                4, config.threads, tolerance);
            options.drop_tolerance = 0.0;
            auto candidate = experiment_support::make_candidate(
                "local-energy-4",
                "tol=" + experiment_support::scientific(tolerance, 0),
                tgi::build_interpolation(grid, a, options));
            rows.push_back(experiment_support::evaluate_candidate(
                field.name, grid, a, rhs, candidate, config));
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
    report.save("solver_tolerance");
    experiment_support::write_csv(
        "solver_tolerance", experiment_support::study_headers(), rows);
    return 0;
}
