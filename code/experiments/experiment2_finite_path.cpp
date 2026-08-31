#include "experiment/problem.hpp"
#include "experiment/reporting.hpp"
#include "multigrid/global_pcg.hpp"
#include "version.hpp"

#include <array>
#include <string>

namespace {

struct PathCase {
    int fine;
    int coarse;
    double contrast;
    experiment_support::FieldCase field;
};

experiment_support::Row path_row(
    const PathCase& item, const std::string& method, int steps,
    const tgi::SparseMatrix& matrix, const tgi::Vector& rhs,
    const tgi::SparseMatrix& prolongation, double reference_energy,
    double initial_excess, int threads, int maximum_cycles) {
    const tgi::TwoGridCycle cycle(matrix, prolongation, 1, threads);
    const auto solved = tgi::solve_two_grid(
        rhs, cycle, 1.0e-6, maximum_cycles);
    const double excess = initial_excess > 0.0
        ? (cycle.setup_report().interpolation_energy - reference_energy) /
            initial_excess
        : 0.0;
    return {
        std::to_string(item.fine),
        std::to_string(item.coarse),
        item.field.name,
        experiment_support::scientific(item.contrast, 0),
        method,
        steps < 0 ? "tol=1e-10" : std::to_string(steps),
        experiment_support::scientific(excess, 3),
        experiment_support::fixed(
            experiment_support::interpolation_density_percent(
                prolongation), 4),
        std::to_string(solved.cycles),
        tgi::two_grid_status_name(solved.status),
        experiment_support::scientific(solved.relative_residual, 2),
        experiment_support::fixed(solved.effective_factor, 6),
        experiment_support::fixed(solved.tail_factor, 6)};
}

}

int main(int argc, char** argv) {
    experiment_support::BasicConfig base;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument.rfind("--threads=", 0) == 0) {
            base.threads = std::stoi(argument.substr(10));
        }
    }
    const auto& topologies = experiment_support::channel_topologies();
    const std::array<PathCase, 3> cases{
        PathCase{64, 16, 1.0e4, topologies[0]},
        PathCase{128, 16, 1.0e4, topologies[0]},
        PathCase{128, 16, 1.0e6, topologies[0]}};
    constexpr int maximum_cycles =
        experiment_support::maximum_two_grid_cycles;
    experiment_support::Rows rows;

    for (const auto& item : cases) {
        experiment_support::BasicConfig config = base;
        config.fine_intervals = item.fine;
        config.coarse_intervals = item.coarse;
        config.contrast = item.contrast;
        const tgi::StructuredGrid grid = experiment_support::make_grid(config);
        const auto problem = experiment_support::make_problem(
            grid, item.field, config);
        const auto geometric = tgi::build_geometric_interpolation(grid);
        const auto reference = experiment_support::build_global_reference(
            grid, problem.matrix, config.threads);
        const tgi::TwoGridCycle reference_cycle(
            problem.matrix, reference.prolongation, 1, config.threads);
        const tgi::TwoGridCycle geometric_cycle(
            problem.matrix, geometric.prolongation, 1, config.threads);
        const double reference_energy =
            reference_cycle.setup_report().interpolation_energy;
        const double initial_excess =
            geometric_cycle.setup_report().interpolation_energy -
            reference_energy;

        rows.push_back(path_row(
            item, "geometric", 0, problem.matrix, problem.rhs,
            geometric.prolongation, reference_energy, initial_excess,
            config.threads, experiment_support::maximum_geometric_cycles));
        tgi::GlobalEnergyPcgPath path(
            grid, problem.matrix, geometric.prolongation, config.threads);
        for (int numerator = 1; numerator <= 8; ++numerator) {
            const int steps =
                (numerator * item.fine + 8) / 16;
            path.advance_to(steps);
            rows.push_back(path_row(
                item, "finite-PCG", steps, problem.matrix, problem.rhs,
                path.prolongation(), reference_energy, initial_excess,
                config.threads, maximum_cycles));
        }
        rows.push_back(path_row(
            item, "global-reference", -1, problem.matrix, problem.rhs,
            reference.prolongation, reference_energy, initial_excess,
            config.threads, maximum_cycles));
    }

    experiment_support::Report report(
        "Finite-PCG path evidence for nonmonotone two-grid behavior");
    report.add_summary({
        {"Version", std::string(tgi::version)},
        {"Threads", std::to_string(base.threads)},
        {"Cases", std::to_string(cases.size())},
        {"Checkpoints", "m/(1/h)=0,1/16,...,1/2 and reference"},
        {"Solve tolerance", "1e-6"},
        {"Maximum cycles", "finite/reference 20000; geometric 30000"}});
    report.add_note(
        "The three cross-channel cases isolate scale and contrast. Only the "
        "energy excess, interpolation density and independently measured "
        "two-grid cycles are retained: energy decreases along the PCG path, "
        "whereas the iteration count can attain a much better finite minimum.");
    report.add_table(
        "Finite path comparison",
        {"1/h", "1/H", "Topology", "Contrast", "Method", "m",
         "Energy excess", "P density %", "Cycles", "Status",
         "Final relres", "Effective factor", "Tail factor"},
        {5, 5, 20, 10, 14, 7, 13, 11, 8, 10, 12, 16, 11}, rows, true);
    report.save("experiment2_finite_path");
    return 0;
}
