#include "experiment/common.hpp"
#include "experiment/metrics.hpp"
#include "experiment/test_field_dataset.hpp"

#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Config {
    int fine_intervals = 128;
    int coarse_intervals = 16;
    int patch_layers = 3;
    std::vector<double> tolerances{
        1.0e-2, 3.0e-3, 1.0e-3,
        3.0e-4, 1.0e-4, 1.0e-10
    };
    double reference_tolerance = 1.0e-10;
    int local_max_iterations = 20000;
    int threads = 4;
    double outer_tolerance = 1.0e-6;
    int outer_max_cycles = 20000;
    std::filesystem::path fields = "models/test_fields.tgi";
};

} // namespace

int main(int argc, char** argv) {
    Config cfg;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument.rfind("--fields=", 0) == 0) {
            cfg.fields = argument.substr(9);
        } else {
            throw std::invalid_argument(
                "usage: experiment2 [--fields=PATH]");
        }
    }

    const auto dataset = experiment_support::load_test_fields(
        cfg.fields, cfg.fine_intervals, cfg.coarse_intervals);
    const int ratio = cfg.fine_intervals / cfg.coarse_intervals;
    const tgi::StructuredGrid grid(cfg.fine_intervals - 1, ratio);
    experiment_support::Rows rows;

    for (const auto& field : dataset.fields) {
        const tgi::SparseMatrix a =
            tgi::assemble_diffusion(grid, field.coefficient.values);
        const tgi::Vector rhs =
            a.multiply(experiment_support::manufactured_solution(grid));

        tgi::InterpolationOptions reference_options;
        reference_options.strategy =
            tgi::InterpolationStrategy::LocalEnergyMinimum;
        reference_options.patch_layers = cfg.patch_layers;
        reference_options.local_tolerance = cfg.reference_tolerance;
        reference_options.local_max_iterations = cfg.local_max_iterations;
        reference_options.thread_count = cfg.threads;
        const tgi::InterpolationResult reference =
            tgi::build_interpolation(grid, a, reference_options);

        for (double tolerance : cfg.tolerances) {
            const auto build_begin = experiment_support::Clock::now();
            tgi::InterpolationResult candidate;
            double build_ms = 0.0;
            if (tolerance == cfg.reference_tolerance) {
                candidate = reference;
                build_ms = reference.report.timing.total_ms;
            } else {
                tgi::InterpolationOptions options = reference_options;
                options.local_tolerance = tolerance;
                candidate = tgi::build_interpolation(grid, a, options);
                build_ms = experiment_support::milliseconds(
                    build_begin, experiment_support::Clock::now());
            }
            const experiment_support::BasisMetrics basis =
                experiment_support::compare_prolongations(
                    grid, a, reference.prolongation,
                    candidate.prolongation, cfg.patch_layers);
            const auto residual = experiment_support::summarize(
                basis.relative_residual);
            const auto l2 = experiment_support::summarize(
                basis.relative_l2_error);
            const auto energy = experiment_support::summarize(
                basis.relative_energy_error);
            const experiment_support::CycleMetrics cycle =
                experiment_support::evaluate_two_grid(
                    a, rhs, candidate.prolongation, cfg.threads,
                    cfg.outer_tolerance, cfg.outer_max_cycles);
            rows.push_back({
                field.name,
                experiment_support::scientific(tolerance, 1),
                experiment_support::scientific(residual.mean),
                experiment_support::scientific(l2.mean),
                experiment_support::scientific(energy.mean),
                std::to_string(candidate.prolongation.nnz()),
                cycle.converged
                    ? std::to_string(cycle.cycles)
                    : "failed@" + std::to_string(cycle.cycles),
                experiment_support::scientific(cycle.relative_residual),
                experiment_support::fixed(build_ms + cycle.total_ms)
            });
        }
    }

    experiment_support::Report report(
        "Experiment 2 - Local tolerance across shared high-contrast fields");
    report.add_summary({
        {"Dataset", cfg.fields.string()},
        {"Fine grid", "h = 1/" + std::to_string(cfg.fine_intervals)},
        {"Coarse grid", "H = 1/" + std::to_string(cfg.coarse_intervals)},
        {"H/h", std::to_string(ratio)},
        {"Patch layers", std::to_string(cfg.patch_layers)},
        {"Fields", std::to_string(dataset.fields.size())},
        {"Reference tolerance",
         experiment_support::scientific(cfg.reference_tolerance, 1)}
    });
    report.add_table(
        "Basis accuracy, convergence, and timing",
        {"Field", "Local tol", "Residual mean", "L2 mean",
         "Energy mean", "P nnz", "Cycles", "Final residual", "Total ms"},
        {20, 10, 13, 11, 13, 10, 12, 15, 10},
        rows, true);
    report.save("experiment2");
    return 0;
}
