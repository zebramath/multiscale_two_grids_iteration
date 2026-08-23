#include "experiment/common.hpp"
#include "experiment/test_field_dataset.hpp"
#include "multigrid/energy_interpolation.hpp"
#include "multigrid/two_grid_solver.hpp"

#include <cmath>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Config {
    int fine_intervals = 128;
    int coarse_intervals = 16;
    double exact_tolerance = 1.0e-10;
    double local_tolerance = 1.0e-3;
    int local_max_iterations = 20000;
    int threads = 4;
    double outer_tolerance = 1.0e-6;
    int outer_max_cycles = 20000;
    std::filesystem::path fields = "models/test_fields.tgi";
};

struct Method {
    std::string name;
    tgi::InterpolationStrategy strategy;
    int patch_layers;
    double tolerance;
    bool use_reference = false;
};

double half_trace(const tgi::SparseMatrix& matrix) {
    double trace = 0.0;
    for (int row = 0; row < matrix.rows(); ++row) {
        for (int position =
                 matrix.row_ptr()[static_cast<std::size_t>(row)];
             position <
                 matrix.row_ptr()[static_cast<std::size_t>(row) + 1U];
             ++position) {
            if (matrix.col_idx()[static_cast<std::size_t>(position)] == row) {
                trace += matrix.values()[static_cast<std::size_t>(position)];
                break;
            }
        }
    }
    return 0.5 * trace;
}

} // namespace

int main(int argc, char** argv) {
    Config cfg;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument.rfind("--fields=", 0) == 0) {
            cfg.fields = argument.substr(9);
        } else {
            throw std::invalid_argument(
                "usage: experiment1 [--fields=PATH]");
        }
    }

    const auto dataset = experiment_support::load_test_fields(
        cfg.fields, cfg.fine_intervals, cfg.coarse_intervals);
    const int ratio = cfg.fine_intervals / cfg.coarse_intervals;
    const tgi::StructuredGrid grid(cfg.fine_intervals - 1, ratio);
    const std::vector<Method> methods{
        {"P_G", tgi::InterpolationStrategy::GeometricBilinear, 0, 0.0},
        {"P_E,2 exact", tgi::InterpolationStrategy::LocalEnergyMinimum,
         2, cfg.exact_tolerance},
        {"P_E,2 tol=1e-3", tgi::InterpolationStrategy::LocalEnergyMinimum,
         2, cfg.local_tolerance},
        {"P_E,3 exact", tgi::InterpolationStrategy::LocalEnergyMinimum,
         3, cfg.exact_tolerance},
        {"P_E,3 tol=1e-3", tgi::InterpolationStrategy::LocalEnergyMinimum,
         3, cfg.local_tolerance},
        {"P_E,4 exact", tgi::InterpolationStrategy::LocalEnergyMinimum,
         4, cfg.exact_tolerance},
        {"P_E,4 tol=1e-3", tgi::InterpolationStrategy::LocalEnergyMinimum,
         4, cfg.local_tolerance},
        {"P_E,inf exact", tgi::InterpolationStrategy::GlobalEnergyMinimum,
         0, cfg.exact_tolerance, true},
        {"P_E,inf tol=1e-3",
         tgi::InterpolationStrategy::GlobalEnergyMinimum,
         0, cfg.local_tolerance}
    };

    experiment_support::Rows rows;
    for (const auto& field : dataset.fields) {
        const tgi::SparseMatrix a =
            tgi::assemble_diffusion(grid, field.coefficient.values);
        const tgi::Vector rhs =
            a.multiply(experiment_support::manufactured_solution(grid));

        tgi::InterpolationOptions reference_options;
        reference_options.strategy =
            tgi::InterpolationStrategy::GlobalEnergyMinimum;
        reference_options.local_tolerance = cfg.exact_tolerance;
        reference_options.local_max_iterations =
            cfg.local_max_iterations;
        reference_options.thread_count = cfg.threads;
        const auto reference_begin = experiment_support::Clock::now();
        const tgi::InterpolationResult reference =
            tgi::build_interpolation(grid, a, reference_options);
        const double reference_build_ms =
            experiment_support::milliseconds(
                reference_begin, experiment_support::Clock::now());
        const tgi::TwoGridCycle reference_cycle(
            a, reference.prolongation, 1, cfg.threads);
        const double reference_energy =
            half_trace(reference_cycle.coarse_matrix());

        for (const Method& method : methods) {
            tgi::InterpolationResult computed;
            const tgi::InterpolationResult* interpolation = &reference;
            double build_ms = reference_build_ms;
            if (!method.use_reference) {
                tgi::InterpolationOptions options;
                options.strategy = method.strategy;
                options.patch_layers = method.patch_layers;
                options.local_tolerance = method.tolerance;
                options.local_max_iterations =
                    cfg.local_max_iterations;
                options.thread_count = cfg.threads;
                const auto build_begin =
                    experiment_support::Clock::now();
                computed = tgi::build_interpolation(grid, a, options);
                build_ms = experiment_support::milliseconds(
                    build_begin, experiment_support::Clock::now());
                interpolation = &computed;
            }

            const auto evaluate = [&](const tgi::TwoGridCycle& cycle) {
                const auto solve_begin =
                    experiment_support::Clock::now();
                const tgi::TwoGridIterationResult solve =
                    tgi::solve_two_grid(
                        a, rhs, cycle, cfg.outer_tolerance,
                        cfg.outer_max_cycles);
                const double solve_ms =
                    experiment_support::milliseconds(
                        solve_begin,
                        experiment_support::Clock::now());
                rows.push_back({
                    field.name,
                    method.name,
                    experiment_support::fixed(
                        half_trace(cycle.coarse_matrix()) /
                            reference_energy,
                        6),
                    solve.converged
                        ? std::to_string(solve.cycles)
                        : "failed@" + std::to_string(solve.cycles),
                    experiment_support::scientific(
                        solve.relative_residual),
                    std::to_string(
                        interpolation->prolongation.nnz()),
                    std::to_string(
                        cycle.setup_report().coarse_nnz),
                    std::to_string(
                        cycle.setup_report().factor_nnz),
                    experiment_support::fixed(
                        build_ms +
                        cycle.setup_report().total_ms +
                        solve_ms)
                });
            };
            if (method.use_reference) {
                evaluate(reference_cycle);
            } else {
                const tgi::TwoGridCycle cycle(
                    a, interpolation->prolongation, 1, cfg.threads);
                evaluate(cycle);
            }
        }
    }

    experiment_support::Report report(
        "Experiment 1 - Localization across shared high-contrast fields");
    report.add_summary({
        {"Dataset", cfg.fields.string()},
        {"Fine grid", "h = 1/" + std::to_string(cfg.fine_intervals)},
        {"Coarse grid", "H = 1/" + std::to_string(cfg.coarse_intervals)},
        {"H/h", std::to_string(ratio)},
        {"Fields", std::to_string(dataset.fields.size())},
        {"Local tolerance",
         experiment_support::scientific(cfg.local_tolerance, 1)},
        {"Outer tolerance",
         experiment_support::scientific(cfg.outer_tolerance, 1)}
    });
    report.add_note(
        "For each support, exact uses tolerance 1e-10 and tol uses 1e-3. "
        "All methods share the same manufactured solution, Galerkin "
        "construction, two-grid cycle, and stopping rule.");
    report.add_table(
        "Localization, convergence, sparsity, and timing",
        {"Field", "Method", "J/Jglobal", "Cycles", "Final residual",
         "P nnz", "Ac nnz", "L nnz", "Total ms"},
        {20, 20, 11, 12, 15, 10, 10, 10, 10},
        rows, true);
    report.save("experiment1");
    return 0;
}
