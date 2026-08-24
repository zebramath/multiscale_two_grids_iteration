#include "experiment/study.hpp"
#include "multigrid/global_pcg.hpp"

#include <array>
#include <string>

namespace {

struct PathCase {
    int fine;
    int coarse;
    double contrast;
    std::uint64_t seed;
    experiment_support::FieldCase field;
};

experiment_support::Row path_row(
    const PathCase& item, const std::string& method, int steps,
    const tgi::SparseMatrix& matrix, const tgi::Vector& rhs,
    const tgi::SparseMatrix& prolongation, double exact_energy,
    double initial_excess, int threads, int maximum_cycles) {
    const tgi::TwoGridCycle cycle(matrix, prolongation, 1, threads);
    const auto solved = tgi::solve_two_grid(
        rhs, cycle, 1.0e-6, maximum_cycles);
    const double excess = initial_excess > 0.0
        ? (cycle.setup_report().interpolation_energy - exact_energy) /
            initial_excess
        : 0.0;
    return {
        std::to_string(item.fine),
        std::to_string(item.coarse),
        item.field.name,
        experiment_support::scientific(item.contrast, 0),
        method,
        steps < 0 ? "exact" : std::to_string(steps),
        experiment_support::scientific(excess, 3),
        experiment_support::fixed(
            experiment_support::interpolation_density_percent(
                prolongation), 4),
        solved.converged ? std::to_string(solved.cycles)
            : "failed@" + std::to_string(solved.cycles)};
}

}

int run_finite_path_evidence(int argc, char** argv) {
    const auto base = experiment_support::parse_config(argc, argv);
    const auto& topologies = experiment_support::channel_topologies();
    const std::array<PathCase, 3> cases{
        PathCase{64, 16, 1.0e4, 1, topologies[0]},
        PathCase{128, 16, 1.0e4, 1, topologies[0]},
        PathCase{128, 16, 1.0e6, 1, topologies[0]}};
    constexpr std::array<int, 7> checkpoints{
        12, 20, 28, 36, 44, 52, 60};
    constexpr int maximum_cycles = 6000;
    experiment_support::Rows rows;

    for (const auto& item : cases) {
        experiment_support::BasicConfig config = base;
        config.fine_intervals = item.fine;
        config.coarse_intervals = item.coarse;
        config.contrast = item.contrast;
        config.max_cycles = maximum_cycles;
        const tgi::StructuredGrid grid = experiment_support::make_grid(config);
        const auto problem = experiment_support::make_problem(
            grid, item.field, config, item.seed);
        const auto geometric = experiment_support::geometric_interpolation(
            grid, problem.matrix);
        const auto exact = experiment_support::build_global_reference(
            grid, problem.matrix, config.threads);
        const tgi::TwoGridCycle exact_cycle(
            problem.matrix, exact.prolongation, 1, config.threads);
        const tgi::TwoGridCycle geometric_cycle(
            problem.matrix, geometric.prolongation, 1, config.threads);
        const double exact_energy =
            exact_cycle.setup_report().interpolation_energy;
        const double initial_excess =
            geometric_cycle.setup_report().interpolation_energy -
            exact_energy;

        rows.push_back(path_row(
            item, "geometric", 0, problem.matrix, problem.rhs,
            geometric.prolongation, exact_energy, initial_excess,
            config.threads, maximum_cycles));
        tgi::GlobalEnergyPcgPath path(
            grid, problem.matrix, geometric.prolongation, config.threads);
        for (int steps : checkpoints) {
            path.advance_to(steps);
            rows.push_back(path_row(
                item, "finite-PCG", steps, problem.matrix, problem.rhs,
                path.prolongation(0.0), exact_energy, initial_excess,
                config.threads, maximum_cycles));
        }
        rows.push_back(path_row(
            item, "global-exact", -1, problem.matrix, problem.rhs,
            exact.prolongation, exact_energy, initial_excess,
            config.threads, maximum_cycles));
    }

    experiment_support::Report report(
        "Finite-PCG path evidence for nonmonotone two-grid behavior");
    report.add_summary({
        {"Version", std::string(tgi::version)},
        {"Cases", std::to_string(cases.size())},
        {"Checkpoints", "m=0,12,20,...,60 and exact"},
        {"Solve tolerance", "1e-6"}});
    report.add_note(
        "The three cross-channel cases isolate scale and contrast. Only the "
        "energy excess, interpolation density and independently measured "
        "two-grid cycles are retained: energy decreases along the PCG path, "
        "whereas the iteration count can attain a much better finite minimum.");
    report.add_table(
        "Finite path comparison",
        {"1/h", "1/H", "Topology", "Contrast", "Method", "m", "Energy excess",
         "P density %", "Cycles"},
        {5, 5, 20, 10, 14, 7, 13, 11, 12}, rows, true);
    report.save("finite_path_evidence");
    return 0;
}
