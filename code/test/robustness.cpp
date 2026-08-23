#include "experiment/common.hpp"
#include "multigrid/adaptive_support.hpp"
#include "multigrid/algebraic_interpolation.hpp"

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct Scale {
    int fine_intervals;
    int coarse_intervals;
};

struct Candidate {
    std::string name;
    tgi::SparseMatrix prolongation;
    double build_ms;
};

Candidate strong_k64(
    const tgi::StructuredGrid& grid, const tgi::SparseMatrix& a,
    const tgi::InterpolationResult& local2, int threads) {
    const auto begin = experiment_support::Clock::now();
    tgi::ResidualStrongSupportOptions options;
    options.base_patch_layers = 2;
    options.maximum_extra_nodes_per_column = 64;
    options.maximum_graph_hops = 4 * grid.fine_n();
    options.strong_edge_fraction = 0.25;
    options.thread_count = threads;
    const auto support = tgi::build_residual_strong_supports(
        grid, a, local2.prolongation, options);
    std::vector<unsigned char> changed(
        static_cast<std::size_t>(grid.coarse_size()), 0U);
    for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
        changed[static_cast<std::size_t>(coarse)] =
            support.supports[static_cast<std::size_t>(coarse)].size() >
                    grid.patch_f_nodes(coarse, 2).size()
                ? 1U
                : 0U;
    }
    const auto interpolation =
        tgi::refine_selected_energy_interpolation_on_supports(
            grid, a, support.supports, changed, local2.prolongation,
            experiment_support::energy_options(2, threads));
    return {
        "strong-K64", interpolation.prolongation,
        local2.report.timing.total_ms + experiment_support::milliseconds(
            begin, experiment_support::Clock::now())};
}

} // namespace

int main(int argc, char** argv) {
    int threads = 4;
    int max_cycles = 4000;
    bool quick = false;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument.rfind("--threads=", 0) == 0) {
            threads = std::stoi(argument.substr(10));
        } else if (argument.rfind("--max-cycles=", 0) == 0) {
            max_cycles = std::stoi(argument.substr(13));
        } else if (argument == "--quick") {
            quick = true;
        } else {
            throw std::invalid_argument(
                "usage: robustness [--quick] [--threads=N] "
                "[--max-cycles=N]");
        }
    }
    const std::vector<Scale> scales = quick
        ? std::vector<Scale>{{64, 8}}
        : std::vector<Scale>{{64, 8}, {128, 16}};
    const std::vector<double> contrasts{1.0e2, 1.0e4, 1.0e6};
    experiment_support::Rows rows;

    for (const Scale& scale : scales) {
        const tgi::StructuredGrid grid(
            scale.fine_intervals - 1,
            scale.fine_intervals / scale.coarse_intervals);
        for (double contrast : contrasts) {
            for (const auto& field : experiment_support::standard_fields()) {
                const auto coefficient = experiment_support::make_field(
                    grid, field, contrast);
                const tgi::SparseMatrix a = tgi::assemble_diffusion(
                    grid, coefficient.values);
                const tgi::Vector rhs = a.multiply(
                    experiment_support::manufactured_solution(grid));
                const auto geometric =
                    experiment_support::geometric_interpolation(grid, a);
                const auto local2 = tgi::build_interpolation(
                    grid, a,
                    experiment_support::energy_options(2, threads));
                const auto local3 = tgi::build_interpolation(
                    grid, a,
                    experiment_support::energy_options(3, threads));
                tgi::JacobiInterpolationOptions jacobi_options;
                jacobi_options.steps = 4;
                jacobi_options.maximum_entries_per_row = 8;
                jacobi_options.thread_count = threads;
                const auto jacobi = tgi::build_jacobi_interpolation(
                    grid, a, geometric.prolongation, jacobi_options);
                const auto strong = strong_k64(
                    grid, a, local2, threads);
                const auto global = tgi::build_interpolation(
                    grid, a,
                    experiment_support::energy_options(
                        0, threads, 1.0e-10));

                const std::vector<Candidate> candidates{
                    {"geometric", geometric.prolongation,
                     geometric.report.timing.total_ms},
                    {"Jacobi-4", jacobi.prolongation,
                     geometric.report.timing.total_ms +
                         jacobi.report.build_ms},
                    {"local-3", local3.prolongation,
                     local3.report.timing.total_ms},
                    strong,
                    {"global-reference", global.prolongation,
                     global.report.timing.total_ms}
                };
                for (const Candidate& candidate : candidates) {
                    const auto cycles = experiment_support::evaluate_two_grid(
                        a, rhs, candidate.prolongation, threads,
                        1.0e-6, max_cycles);
                    rows.push_back({
                        std::to_string(scale.fine_intervals) + "/" +
                            std::to_string(scale.coarse_intervals),
                        experiment_support::scientific(contrast, 0),
                        field.name,
                        candidate.name,
                        std::to_string(candidate.prolongation.nnz()),
                        cycles.converged
                            ? std::to_string(cycles.cycles)
                            : "failed@" + std::to_string(cycles.cycles),
                        experiment_support::fixed(candidate.build_ms),
                        experiment_support::fixed(cycles.total_ms)
                    });
                }
            }
        }
    }

    const experiment_support::Row headers{
        "Fine/Coarse", "Contrast", "Field", "Method", "P nnz",
        "Cycles", "Build ms", "Solve ms"};
    experiment_support::Report report(
        "Contrast and scale robustness with one deterministic RHS");
    report.add_summary({
        {"Scales", std::to_string(scales.size())},
        {"Contrasts", "1e2, 1e4, 1e6"},
        {"Seed", "1"},
        {"RHS", "one manufactured solution"}
    });
    report.add_note(
        "This script intentionally removes repeated seeds and multiple RHS. "
        "It keeps the two dimensions that changed the qualitative conclusion: "
        "coefficient contrast and physical grid scale.");
    report.add_table(
        "Robustness comparison", headers,
        {11, 10, 12, 18, 10, 12, 11, 11}, rows, true);
    report.save("robustness");
    experiment_support::write_csv("robustness", headers, rows);
    return 0;
}
