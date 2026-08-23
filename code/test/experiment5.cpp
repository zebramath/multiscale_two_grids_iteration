#include "experiment/common.hpp"
#include "experiment/metrics.hpp"
#include "experiment/test_field_dataset.hpp"
#include "multigrid/adaptive_support.hpp"
#include "multigrid/energy_interpolation.hpp"
#include "multigrid/residual_budget_support.hpp"
#include "multigrid/support_pruning.hpp"

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
    int rounds = 0;
    int expanded_columns = 0;
    int extra_nodes = 0;
    double residual_ratio = 1.0;
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

std::vector<unsigned char> expanded_mask(
    const tgi::StructuredGrid& grid,
    const std::vector<std::vector<int>>& supports,
    int base_layers) {
    std::vector<unsigned char> mask(
        static_cast<std::size_t>(grid.coarse_size()), 0U);
    for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
        const std::size_t base_size =
            grid.patch_f_nodes(coarse, base_layers).size();
        if (supports[static_cast<std::size_t>(coarse)].size() > base_size) {
            mask[static_cast<std::size_t>(coarse)] = 1U;
        }
    }
    return mask;
}

Candidate build_budget_candidate(
    const std::string& name, const Config& cfg,
    const tgi::StructuredGrid& grid, const tgi::SparseMatrix& a,
    const tgi::SparseMatrix& local,
    double local_build_ms, tgi::StrengthScaling scaling) {
    tgi::ResidualBudgetSupportOptions support;
    support.base_patch_layers = cfg.base_patch_layers;
    support.maximum_rounds = 8;
    support.maximum_extra_nodes_per_column = 128;
    support.maximum_nodes_per_round = 16;
    support.marking_fraction = 0.70;
    support.target_residual_ratio = 0.25;
    support.strength_scaling = scaling;
    support.strong_edge_fraction = scaling == tgi::StrengthScaling::SymmetricDiagonal
        ? 0.10
        : 0.25;
    support.thread_count = cfg.threads;
    const auto adaptive = tgi::build_residual_budget_interpolation(
        grid, a, local,
        interpolation_options(
            cfg, cfg.base_patch_layers, cfg.local_tolerance),
        support);
    return {
        name,
        adaptive.prolongation,
        local_build_ms + adaptive.report.total_ms,
        adaptive.report.rounds,
        adaptive.report.expanded_columns,
        adaptive.report.total_extra_nodes,
        adaptive.report.mean_final_to_initial_ratio
    };
}

} // namespace

int main(int argc, char** argv) {
    Config cfg;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument.rfind("--fields=", 0) == 0) {
            cfg.fields = argument.substr(9);
        } else if (argument.rfind("--threads=", 0) == 0) {
            cfg.threads = std::stoi(argument.substr(10));
        } else if (argument.rfind("--max-cycles=", 0) == 0) {
            cfg.outer_max_cycles = std::stoi(argument.substr(13));
        } else {
            throw std::invalid_argument(
                "usage: experiment5 [--fields=PATH] [--threads=N] "
                "[--max-cycles=N]");
        }
    }

    const auto dataset = experiment_support::load_test_fields(
        cfg.fields, cfg.fine_intervals, cfg.coarse_intervals);
    const int ratio = cfg.fine_intervals / cfg.coarse_intervals;
    const tgi::StructuredGrid grid(cfg.fine_intervals - 1, ratio);
    experiment_support::Rows rows;

    for (const auto& field : dataset.fields) {
        const tgi::SparseMatrix a = tgi::assemble_diffusion(
            grid, field.coefficient.values);
        const tgi::Vector rhs = a.multiply(
            experiment_support::manufactured_solution(grid));

        const auto reference = tgi::build_interpolation(
            grid, a, interpolation_options(
                cfg, 0, cfg.exact_tolerance, true));
        const auto local_begin = experiment_support::Clock::now();
        const auto local = tgi::build_interpolation(
            grid, a, interpolation_options(
                cfg, cfg.base_patch_layers, cfg.local_tolerance));
        const double local_ms = experiment_support::milliseconds(
            local_begin, experiment_support::Clock::now());

        std::vector<Candidate> candidates;
        candidates.push_back({
            "P_E,2 tol", local.prolongation, local_ms, 0, 0, 0, 1.0});

        tgi::ResidualStrongSupportOptions fixed_options;
        fixed_options.base_patch_layers = cfg.base_patch_layers;
        fixed_options.maximum_extra_nodes_per_column = 64;
        fixed_options.maximum_graph_hops = 4 * grid.fine_n();
        fixed_options.strong_edge_fraction = 0.25;
        fixed_options.thread_count = cfg.threads;
        const auto fixed = tgi::build_residual_strong_supports(
            grid, a, local.prolongation, fixed_options);
        const auto fixed_mask = expanded_mask(
            grid, fixed.supports, cfg.base_patch_layers);
        const auto fixed_refine_begin = experiment_support::Clock::now();
        auto fixed_refine_options = interpolation_options(
            cfg, cfg.base_patch_layers, cfg.local_tolerance);
        fixed_refine_options.local_tolerance = 1.0e-6;
        const auto fixed_interpolation =
            tgi::refine_selected_energy_interpolation_on_supports(
                grid, a, fixed.supports, fixed_mask,
                local.prolongation,
                fixed_refine_options);
        const double fixed_refine_ms = experiment_support::milliseconds(
            fixed_refine_begin, experiment_support::Clock::now());
        candidates.push_back({
            "fixed-global K=64",
            fixed_interpolation.prolongation,
            local_ms + fixed.report.selection_ms + fixed_refine_ms,
            1,
            fixed.report.expanded_columns,
            fixed.report.total_extra_nodes,
            1.0
        });

        fixed_options.maximum_extra_nodes_per_column = 128;
        const auto expanded128 = tgi::build_residual_strong_supports(
            grid, a, local.prolongation, fixed_options);
        const auto expanded128_mask = expanded_mask(
            grid, expanded128.supports, cfg.base_patch_layers);
        const auto expanded128_begin = experiment_support::Clock::now();
        const auto expanded128_interpolation =
            tgi::refine_selected_energy_interpolation_on_supports(
                grid, a, expanded128.supports, expanded128_mask,
                local.prolongation, fixed_refine_options);
        const double expanded128_refine_ms =
            experiment_support::milliseconds(
                expanded128_begin, experiment_support::Clock::now());
        const double expanded128_build_ms =
            local_ms + expanded128.report.selection_ms +
            expanded128_refine_ms;
        candidates.push_back({
            "fixed-global K=128",
            expanded128_interpolation.prolongation,
            expanded128_build_ms,
            1,
            expanded128.report.expanded_columns,
            expanded128.report.total_extra_nodes,
            1.0
        });

        tgi::SupportPruningOptions pruning_options;
        pruning_options.base_patch_layers = cfg.base_patch_layers;
        pruning_options.maximum_extra_nodes_per_column = 64;
        pruning_options.relative_magnitude_threshold = 1.0e-4;
        pruning_options.refinement_tolerance = 1.0e-6;
        pruning_options.thread_count = cfg.threads;
        const auto pruned = tgi::prune_energy_supports_by_magnitude(
            grid, a, expanded128.supports,
            expanded128_interpolation.prolongation,
            fixed_refine_options, pruning_options);
        candidates.push_back({
            "K128-prune64",
            pruned.prolongation,
            expanded128_build_ms + pruned.report.total_ms,
            2,
            pruned.report.pruned_columns,
            pruned.report.retained_extra_nodes,
            1.0
        });

        candidates.push_back(build_budget_candidate(
            "budget-global", cfg, grid, a, local.prolongation,
            local_ms, tgi::StrengthScaling::GlobalMaximum));
        candidates.push_back(build_budget_candidate(
            "budget-row", cfg, grid, a, local.prolongation,
            local_ms, tgi::StrengthScaling::RowMaximum));
        candidates.push_back(build_budget_candidate(
            "budget-symmetric", cfg, grid, a, local.prolongation,
            local_ms, tgi::StrengthScaling::SymmetricDiagonal));

        for (const Candidate& candidate : candidates) {
            const auto quality =
                experiment_support::compare_prolongations_global(
                    grid, a, reference.prolongation,
                    candidate.prolongation);
            const auto cycle = experiment_support::evaluate_two_grid(
                a, rhs, candidate.prolongation, cfg.threads,
                cfg.outer_tolerance, cfg.outer_max_cycles);
            rows.push_back({
                field.name,
                candidate.name,
                std::to_string(candidate.rounds),
                std::to_string(candidate.expanded_columns),
                std::to_string(candidate.extra_nodes),
                experiment_support::fixed(candidate.residual_ratio, 3),
                experiment_support::scientific(
                    quality.aggregate_relative_energy_error),
                std::to_string(candidate.prolongation.nnz()),
                cycle.converged
                    ? std::to_string(cycle.cycles)
                    : "failed@" + std::to_string(cycle.cycles),
                experiment_support::fixed(candidate.construction_ms),
                experiment_support::fixed(cycle.coarse_setup_ms),
                experiment_support::fixed(cycle.solve_ms),
                experiment_support::fixed(
                    candidate.construction_ms + cycle.total_ms)
            });
        }
    }

    experiment_support::Report report(
        "Experiment 5 - Residual-budget support selection and ablations");
    report.add_summary({
        {"Dataset", cfg.fields.string()},
        {"Fine grid", "h = 1/" + std::to_string(cfg.fine_intervals)},
        {"Coarse grid", "H = 1/" + std::to_string(cfg.coarse_intervals)},
        {"H/h", std::to_string(ratio)},
        {"Threads", std::to_string(cfg.threads)},
        {"Target residual ratio", "0.25"},
        {"Maximum extra nodes", "128"}
    });
    report.add_note(
        "Budget variants recompute AP after each support update, assign each "
        "round to columns covering 70% of the unresolved residual energy, "
        "follow residual-seeded strong paths, and warm-start only changed "
        "columns. Timings are split into "
        "interpolation/support construction, Galerkin/factorization, and solve.");
    const experiment_support::Row headers{
        "Field", "Method", "Rounds", "Expanded", "Extra", "Res ratio",
        "Energy err", "P nnz", "Cycles", "Build ms", "Coarse ms",
        "Solve ms", "Total ms"};
    report.add_table(
        "Support efficiency, quality, and end-to-end cost",
        headers,
        {20, 20, 7, 9, 8, 10, 11, 10, 12, 10, 10, 10, 10},
        rows, true);
    report.save("experiment5");
    experiment_support::write_csv("experiment5", headers, rows);
    return 0;
}
