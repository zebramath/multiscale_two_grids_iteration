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
    double global_tolerance = 1.0e-10;
    std::vector<double> drop_tolerances{
        0.0, 1.0e-3, 3.0e-3, 1.0e-2, 2.0e-2,
        3.0e-2, 5.0e-2, 1.0e-1, 2.0e-1
    };
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
                "usage: experiment3 [--fields=PATH]");
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

        tgi::InterpolationOptions options;
        options.strategy =
            tgi::InterpolationStrategy::GlobalEnergyMinimum;
        options.patch_layers = 0;
        options.local_tolerance = cfg.global_tolerance;
        options.local_max_iterations = cfg.local_max_iterations;
        options.thread_count = cfg.threads;
        const tgi::InterpolationResult base =
            tgi::build_interpolation(grid, a, options);
        const double base_nnz =
            static_cast<double>(base.prolongation.nnz());

        for (double threshold : cfg.drop_tolerances) {
            const auto drop_begin = experiment_support::Clock::now();
            const tgi::SparseMatrix candidate =
                experiment_support::drop_prolongation(
                    base.prolongation, threshold);
            const double drop_ms = experiment_support::milliseconds(
                drop_begin, experiment_support::Clock::now());
            const experiment_support::GlobalProlongationMetrics quality =
                experiment_support::compare_prolongations_global(
                    grid, a, base.prolongation,
                    candidate);
            const auto l2 = experiment_support::summarize(
                quality.basis.relative_l2_error);
            const auto energy = experiment_support::summarize(
                quality.basis.relative_energy_error);
            const experiment_support::CycleMetrics cycle =
                experiment_support::evaluate_two_grid(
                    a, rhs, candidate, cfg.threads,
                    cfg.outer_tolerance, cfg.outer_max_cycles);
            rows.push_back({
                field.name,
                experiment_support::scientific(threshold, 1),
                std::to_string(candidate.nnz()),
                experiment_support::fixed(
                    100.0 * static_cast<double>(candidate.nnz()) /
                        base_nnz,
                    2) + "%",
                std::to_string(cycle.coarse_nnz),
                std::to_string(cycle.factor_nnz),
                experiment_support::scientific(l2.mean),
                experiment_support::scientific(energy.mean),
                cycle.converged
                    ? std::to_string(cycle.cycles)
                    : "failed@" + std::to_string(cycle.cycles),
                experiment_support::scientific(cycle.relative_residual),
                experiment_support::fixed(
                    drop_ms + cycle.total_ms)
            });
        }
    }

    experiment_support::Report report(
        "Experiment 3 - Aggressive global prolongation truncation");
    report.add_summary({
        {"Dataset", cfg.fields.string()},
        {"Fine grid", "h = 1/" + std::to_string(cfg.fine_intervals)},
        {"Coarse grid", "H = 1/" + std::to_string(cfg.coarse_intervals)},
        {"H/h", std::to_string(ratio)},
        {"Base interpolation", "P_E,inf exact"},
        {"Global tolerance",
         experiment_support::scientific(cfg.global_tolerance, 1)},
        {"Fields", std::to_string(dataset.fields.size())}
    });
    report.add_note(
        "Thresholds span mild compression through loss of convergence. "
        "Timing includes truncation and two-grid setup/solve, excluding the "
        "shared global interpolation construction.");
    report.add_table(
        "Sparsity, accuracy, convergence, and timing",
        {"Field", "Drop", "P nnz", "Retained", "Ac nnz", "L nnz",
         "Global L2", "Global energy", "Cycles", "Final residual",
         "Truncate+2G ms"},
        {20, 9, 10, 10, 10, 10, 11, 13, 12, 15, 15},
        rows, true);
    report.save("experiment3");
    return 0;
}
