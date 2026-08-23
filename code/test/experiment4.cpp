#include "experiment/common.hpp"
#include "experiment/metrics.hpp"
#include "experiment/test_field_dataset.hpp"
#include "multigrid/adaptive_support.hpp"
#include "multigrid/energy_interpolation.hpp"

#include <filesystem>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct Config {
    int fine_intervals = 128;
    int coarse_intervals = 16;
    int base_patch_layers = 2;
    double exact_tolerance = 1.0e-10;
    double local_tolerance = 1.0e-3;
    double strong_edge_fraction = 0.25;
    int local_max_iterations = 20000;
    int threads = 4;
    double outer_tolerance = 1.0e-6;
    int outer_max_cycles = 20000;
    std::filesystem::path fields = "models/test_fields.tgi";
};

struct Candidate {
    std::string name;
    tgi::SparseMatrix prolongation;
    double construction_ms = 0.0;
    int expanded_columns = 0;
    int extra_nodes = 0;
};

tgi::InterpolationOptions interpolation_options(
    const Config& cfg, int patch_layers, double tolerance,
    bool global = false) {
    tgi::InterpolationOptions options;
    options.strategy = global
        ? tgi::InterpolationStrategy::GlobalEnergyMinimum
        : tgi::InterpolationStrategy::LocalEnergyMinimum;
    options.patch_layers = patch_layers;
    options.local_tolerance = tolerance;
    options.local_max_iterations = cfg.local_max_iterations;
    options.thread_count = cfg.threads;
    return options;
}

Candidate build_standard(
    const std::string& name,
    const tgi::StructuredGrid& grid,
    const tgi::SparseMatrix& a,
    const tgi::InterpolationOptions& options) {
    const auto begin = experiment_support::Clock::now();
    const tgi::InterpolationResult interpolation =
        tgi::build_interpolation(grid, a, options);
    return {
        name,
        interpolation.prolongation,
        experiment_support::milliseconds(
            begin, experiment_support::Clock::now()),
        0,
        0
    };
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
                "usage: experiment4 [--fields=PATH]");
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
        std::vector<Candidate> candidates;

        candidates.push_back(build_standard(
            "P_E,inf exact", grid, a,
            interpolation_options(
                cfg, 0, cfg.exact_tolerance, true)));
        const tgi::SparseMatrix reference =
            candidates.front().prolongation;
        candidates.push_back(build_standard(
            "P_E,2 tol", grid, a,
            interpolation_options(
                cfg, cfg.base_patch_layers,
                cfg.local_tolerance)));
        const tgi::SparseMatrix local2 =
            candidates.back().prolongation;
        const double local2_ms =
            candidates.back().construction_ms;
        for (int budget : std::vector<int>{32, 64, 128}) {
            tgi::ResidualStrongSupportOptions support_options;
            support_options.base_patch_layers =
                cfg.base_patch_layers;
            support_options.maximum_extra_nodes_per_column = budget;
            support_options.maximum_graph_hops = 4 * grid.fine_n();
            support_options.strong_edge_fraction =
                cfg.strong_edge_fraction;
            support_options.thread_count = cfg.threads;
            const tgi::AdaptiveSupportResult support =
                tgi::build_residual_strong_supports(
                    grid, a, local2, support_options);
            const auto begin = experiment_support::Clock::now();
            const tgi::InterpolationResult interpolation =
                tgi::build_energy_interpolation_on_supports(
                    grid, a, support.supports,
                    interpolation_options(
                        cfg, cfg.base_patch_layers,
                        cfg.local_tolerance));
            candidates.push_back({
                "strong K=" + std::to_string(budget),
                interpolation.prolongation,
                local2_ms + support.report.selection_ms +
                    experiment_support::milliseconds(
                        begin, experiment_support::Clock::now()),
                support.report.expanded_columns,
                support.report.total_extra_nodes
            });
        }

        for (const Candidate& candidate : candidates) {
            const auto quality =
                experiment_support::compare_prolongations_global(
                    grid, a, reference, candidate.prolongation);
            const auto cycle =
                experiment_support::evaluate_two_grid(
                    a, rhs, candidate.prolongation, cfg.threads,
                    cfg.outer_tolerance, cfg.outer_max_cycles);
            rows.push_back({
                field.name,
                candidate.name,
                experiment_support::scientific(
                    quality.aggregate_relative_energy_error),
                std::to_string(candidate.expanded_columns),
                std::to_string(candidate.extra_nodes),
                std::to_string(candidate.prolongation.nnz()),
                std::to_string(cycle.coarse_nnz),
                std::to_string(cycle.factor_nnz),
                cycle.converged
                    ? std::to_string(cycle.cycles)
                    : "failed@" + std::to_string(cycle.cycles),
                experiment_support::scientific(cycle.relative_residual),
                experiment_support::fixed(
                    candidate.construction_ms + cycle.total_ms)
            });
        }
    }

    experiment_support::Report report(
        "Experiment 4 - Adaptive nonlocal supports across shared fields");
    report.add_summary({
        {"Dataset", cfg.fields.string()},
        {"Fine grid", "h = 1/" + std::to_string(cfg.fine_intervals)},
        {"Coarse grid", "H = 1/" + std::to_string(cfg.coarse_intervals)},
        {"H/h", std::to_string(ratio)},
        {"Base patch layers", std::to_string(cfg.base_patch_layers)},
        {"Local tolerance",
         experiment_support::scientific(cfg.local_tolerance, 1)},
        {"Fields", std::to_string(dataset.fields.size())}
    });
    report.add_note(
        "Strong supports use AP residuals and coefficient-weighted graph "
        "paths. Only K=32, 64, and 128 residual-driven expansions are kept.");
    report.add_table(
        "Support quality, sparsity, convergence, and timing",
        {"Field", "Method", "Energy error", "Expanded", "Extra nodes",
         "P nnz", "Ac nnz", "L nnz", "Cycles", "Final residual",
         "Total ms"},
        {20, 18, 13, 10, 12, 10, 10, 10, 12, 15, 10},
        rows, true);
    report.save("experiment4");
    return 0;
}
