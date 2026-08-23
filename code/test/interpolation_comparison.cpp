#include "experiment/common.hpp"
#include "multigrid/algebraic_interpolation.hpp"

#include <string>
#include <utility>
#include <vector>

namespace {

struct Candidate {
    std::string name;
    tgi::SparseMatrix prolongation;
    double build_ms = 0.0;
    double f_residual = 0.0;
};

double f_residual(
    const tgi::StructuredGrid& grid, const tgi::SparseMatrix& a,
    const tgi::SparseMatrix& p, int threads) {
    return tgi::algebraic_interpolation_detail::scaled_f_residual(
        grid, a, p, threads);
}

} // namespace

int main(int argc, char** argv) {
    experiment_support::BasicConfig config;
    for (int index = 1; index < argc; ++index) {
        experiment_support::parse_basic_argument(config, argv[index]);
    }
    const tgi::StructuredGrid grid =
        experiment_support::make_grid(config);
    experiment_support::Rows rows;

    for (const auto& field : experiment_support::standard_fields()) {
        const auto coefficient = experiment_support::make_field(
            grid, field, config.contrast);
        const tgi::SparseMatrix a = tgi::assemble_diffusion(
            grid, coefficient.values);
        const tgi::Vector rhs = a.multiply(
            experiment_support::manufactured_solution(grid));
        const auto geometric =
            experiment_support::geometric_interpolation(grid, a);
        const auto local2 = tgi::build_interpolation(
            grid, a,
            experiment_support::energy_options(2, config.threads));
        const auto global = tgi::build_interpolation(
            grid, a,
            experiment_support::energy_options(
                0, config.threads, 1.0e-10));

        std::vector<Candidate> candidates;
        candidates.push_back({
            "geometric", geometric.prolongation,
            geometric.report.timing.total_ms,
            f_residual(grid, a, geometric.prolongation, config.threads)});
        for (int steps : {1, 2, 4, 8}) {
            tgi::JacobiInterpolationOptions options;
            options.steps = steps;
            options.damping = 2.0 / 3.0;
            options.maximum_entries_per_row = 8;
            options.relative_drop_tolerance = 1.0e-10;
            options.thread_count = config.threads;
            const auto interpolation = tgi::build_jacobi_interpolation(
                grid, a, geometric.prolongation, options);
            candidates.push_back({
                "Jacobi-" + std::to_string(steps),
                interpolation.prolongation,
                geometric.report.timing.total_ms +
                    interpolation.report.build_ms,
                interpolation.report.final_f_residual});
        }
        tgi::StrengthDistanceOptions distance_options;
        distance_options.coarse_candidates_per_row = 8;
        distance_options.thread_count = config.threads;
        const auto distance =
            tgi::build_strength_distance_interpolation(
                grid, a, distance_options);
        candidates.push_back({
            "strength-distance", distance.prolongation,
            distance.report.build_ms, distance.report.final_f_residual});
        candidates.push_back({
            "local-2", local2.prolongation,
            local2.report.timing.total_ms,
            f_residual(grid, a, local2.prolongation, config.threads)});
        candidates.push_back({
            "global-reference", global.prolongation,
            global.report.timing.total_ms,
            f_residual(grid, a, global.prolongation, config.threads)});

        for (const Candidate& candidate : candidates) {
            const auto error =
                experiment_support::compare_prolongations_global(
                    grid, a, global.prolongation,
                    candidate.prolongation);
            const auto cycles = experiment_support::evaluate_two_grid(
                a, rhs, candidate.prolongation, config.threads,
                1.0e-6, config.max_cycles, 20);
            rows.push_back({
                field.name,
                candidate.name,
                std::to_string(candidate.prolongation.nnz()),
                experiment_support::scientific(candidate.f_residual),
                experiment_support::scientific(
                    error.aggregate_relative_energy_error),
                cycles.converged
                    ? std::to_string(cycles.cycles)
                    : "failed@" + std::to_string(cycles.cycles),
                experiment_support::fixed(
                    cycles.convergence_factor, 4),
                experiment_support::fixed(candidate.build_ms),
                experiment_support::fixed(cycles.total_ms)
            });
        }
    }

    const experiment_support::Row headers{
        "Field", "Interpolation", "P nnz", "F residual",
        "Energy error", "Cycles", "Rho", "Build ms", "Solve ms"};
    experiment_support::Report report(
        "Geometric, Jacobi-relaxed, algebraic-distance and energy interpolation");
    report.add_summary({
        {"Fine grid", "h=1/" + std::to_string(config.fine_intervals)},
        {"Coarse grid", "H=1/" + std::to_string(config.coarse_intervals)},
        {"Contrast", experiment_support::scientific(config.contrast, 0)},
        {"Jacobi damping", "2/3"},
        {"Jacobi row cap", "8"},
        {"Strength candidates", "8 per F row"}
    });
    report.add_note(
        "Jacobi relaxes AP=0 only on F rows and restores exact C-point "
        "injection after every step. Strength-distance uses only the matrix "
        "graph and the prescribed C set; it is included as a deliberately "
        "algebraic support selector followed by local energy minimization.");
    report.add_table(
        "Interpolation comparison", headers,
        {12, 20, 10, 13, 13, 12, 9, 11, 11}, rows, true);
    report.save("interpolation_comparison");
    experiment_support::write_csv(
        "interpolation_comparison", headers, rows);
    return 0;
}
