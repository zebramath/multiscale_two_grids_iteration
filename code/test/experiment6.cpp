#include "experiment/common.hpp"
#include "experiment/test_field_dataset.hpp"
#include "multigrid/adaptive_support.hpp"
#include "multigrid/energy_interpolation.hpp"
#include "multigrid/residual_budget_support.hpp"
#include "multigrid/support_pruning.hpp"
#include "multigrid/two_grid_solver.hpp"

#include <algorithm>
#include <filesystem>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct Config {
    int fine_intervals = 128;
    int coarse_intervals = 16;
    int threads = 4;
    int spectral_iterations = 80;
    int outer_max_cycles = 20000;
    double outer_tolerance = 1.0e-6;
    std::filesystem::path fields = "models/test_fields.tgi";
};

struct Candidate {
    std::string name;
    tgi::SparseMatrix prolongation;
    double build_ms = 0.0;
};

tgi::InterpolationOptions energy_options(
    int layers, double tolerance, int threads, bool global = false) {
    tgi::InterpolationOptions options;
    options.strategy = global
        ? tgi::InterpolationStrategy::GlobalEnergyMinimum
        : tgi::InterpolationStrategy::LocalEnergyMinimum;
    options.patch_layers = layers;
    options.local_tolerance = tolerance;
    options.local_max_iterations = 20000;
    options.thread_count = threads;
    return options;
}

Candidate build_standard(
    const std::string& name, const tgi::StructuredGrid& grid,
    const tgi::SparseMatrix& a,
    const tgi::InterpolationOptions& options) {
    const auto begin = experiment_support::Clock::now();
    const auto interpolation = tgi::build_interpolation(grid, a, options);
    return {
        name,
        interpolation.prolongation,
        experiment_support::milliseconds(
            begin, experiment_support::Clock::now())
    };
}

std::vector<unsigned char> expanded_mask(
    const tgi::StructuredGrid& grid,
    const std::vector<std::vector<int>>& supports,
    int base_layers) {
    std::vector<unsigned char> mask(
        static_cast<std::size_t>(grid.coarse_size()), 0U);
    for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
        if (supports[static_cast<std::size_t>(coarse)].size() >
            grid.patch_f_nodes(coarse, base_layers).size()) {
            mask[static_cast<std::size_t>(coarse)] = 1U;
        }
    }
    return mask;
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
        } else if (argument.rfind("--spectral-iters=", 0) == 0) {
            cfg.spectral_iterations = std::stoi(argument.substr(17));
        } else if (argument.rfind("--max-cycles=", 0) == 0) {
            cfg.outer_max_cycles = std::stoi(argument.substr(13));
        } else {
            throw std::invalid_argument(
                "usage: experiment6 [--fields=PATH] [--threads=N] "
                "[--spectral-iters=N] [--max-cycles=N]");
        }
    }

    const auto dataset = experiment_support::load_test_fields(
        cfg.fields, cfg.fine_intervals, cfg.coarse_intervals);
    const int ratio = cfg.fine_intervals / cfg.coarse_intervals;
    const tgi::StructuredGrid grid(cfg.fine_intervals - 1, ratio);
    experiment_support::Rows rows;
    experiment_support::Rows rhs_rows;

    for (const auto& field : dataset.fields) {
        const auto a = tgi::assemble_diffusion(
            grid, field.coefficient.values);
        std::vector<tgi::Vector> right_hand_sides;
        std::vector<std::string> rhs_names;
        right_hand_sides.push_back(a.multiply(
            experiment_support::manufactured_solution(grid)));
        rhs_names.push_back("manufactured");
        for (std::uint64_t seed : {17ULL, 31ULL, 47ULL}) {
            right_hand_sides.push_back(a.multiply(
                experiment_support::random_spectral_solution(
                    grid, seed)));
            rhs_names.push_back("spectral-" + std::to_string(seed));
        }
        right_hand_sides.push_back(experiment_support::point_source_rhs(
            grid, grid.fine_n() / 2, grid.fine_n() / 2));
        rhs_names.push_back("point-center");
        right_hand_sides.push_back(experiment_support::point_source_rhs(
            grid, grid.fine_n() / 4, 3 * grid.fine_n() / 4));
        rhs_names.push_back("point-quarter");

        std::vector<Candidate> candidates;
        tgi::InterpolationOptions geometric;
        geometric.strategy = tgi::InterpolationStrategy::GeometricBilinear;
        geometric.thread_count = cfg.threads;
        candidates.push_back(build_standard(
            "P_G", grid, a, geometric));
        candidates.push_back(build_standard(
            "P_E,2 tol", grid, a,
            energy_options(2, 1.0e-3, cfg.threads)));
        const Candidate local2 = candidates.back();
        candidates.push_back(build_standard(
            "P_E,3 tol", grid, a,
            energy_options(3, 1.0e-3, cfg.threads)));

        tgi::ResidualStrongSupportOptions fixed_options;
        fixed_options.base_patch_layers = 2;
        fixed_options.maximum_extra_nodes_per_column = 64;
        fixed_options.maximum_graph_hops = 4 * grid.fine_n();
        fixed_options.strong_edge_fraction = 0.25;
        fixed_options.thread_count = cfg.threads;
        const auto fixed64_begin = experiment_support::Clock::now();
        const auto fixed64 = tgi::build_residual_strong_supports(
            grid, a, local2.prolongation, fixed_options);
        auto refine_options = energy_options(2, 1.0e-6, cfg.threads);
        const auto fixed64_interpolation =
            tgi::refine_selected_energy_interpolation_on_supports(
                grid, a, fixed64.supports,
                expanded_mask(grid, fixed64.supports, 2),
                local2.prolongation, refine_options);
        candidates.push_back({
            "fixed K=64",
            fixed64_interpolation.prolongation,
            local2.build_ms + experiment_support::milliseconds(
                fixed64_begin, experiment_support::Clock::now())
        });

        fixed_options.maximum_extra_nodes_per_column = 128;
        const auto prune_begin = experiment_support::Clock::now();
        const auto fixed128 = tgi::build_residual_strong_supports(
            grid, a, local2.prolongation, fixed_options);
        const auto fixed128_interpolation =
            tgi::refine_selected_energy_interpolation_on_supports(
                grid, a, fixed128.supports,
                expanded_mask(grid, fixed128.supports, 2),
                local2.prolongation, refine_options);
        tgi::SupportPruningOptions pruning;
        pruning.base_patch_layers = 2;
        pruning.maximum_extra_nodes_per_column = 64;
        pruning.relative_magnitude_threshold = 1.0e-4;
        pruning.refinement_tolerance = 1.0e-6;
        pruning.thread_count = cfg.threads;
        const auto pruned = tgi::prune_energy_supports_by_magnitude(
            grid, a, fixed128.supports,
            fixed128_interpolation.prolongation,
            refine_options, pruning);
        candidates.push_back({
            "K128-prune64",
            pruned.prolongation,
            local2.build_ms + experiment_support::milliseconds(
                prune_begin, experiment_support::Clock::now())
        });

        tgi::ResidualBudgetSupportOptions support;
        support.base_patch_layers = 2;
        support.maximum_rounds = 8;
        support.maximum_extra_nodes_per_column = 128;
        support.maximum_nodes_per_round = 16;
        support.marking_fraction = 0.70;
        support.target_residual_ratio = 0.25;
        support.strength_scaling = tgi::StrengthScaling::RowMaximum;
        support.strong_edge_fraction = 0.25;
        support.thread_count = cfg.threads;
        const auto adaptive = tgi::build_residual_budget_interpolation(
            grid, a, local2.prolongation,
            energy_options(2, 1.0e-3, cfg.threads), support);
        candidates.push_back({
            "budget-row",
            adaptive.prolongation,
            local2.build_ms + adaptive.report.total_ms
        });
        candidates.push_back(build_standard(
            "P_E,inf ref", grid, a,
            energy_options(0, 1.0e-10, cfg.threads, true)));

        for (const Candidate& candidate : candidates) {
            const auto setup_begin = experiment_support::Clock::now();
            const tgi::TwoGridCycle cycle(
                a, candidate.prolongation, 1, cfg.threads);
            const double setup_ms = experiment_support::milliseconds(
                setup_begin, experiment_support::Clock::now());

            double rho = 0.0;
            for (std::uint64_t seed : {20260803ULL, 20260817ULL}) {
                rho = std::max(
                    rho,
                    cycle.estimate_convergence_factor(
                        cfg.spectral_iterations, seed));
            }

            const auto solves_begin = experiment_support::Clock::now();
            std::vector<int> cycles;
            int failed = 0;
            for (std::size_t rhs_index = 0;
                 rhs_index < right_hand_sides.size(); ++rhs_index) {
                const auto rhs_begin = experiment_support::Clock::now();
                const auto solved = tgi::solve_two_grid(
                    a, right_hand_sides[rhs_index], cycle,
                    cfg.outer_tolerance,
                    cfg.outer_max_cycles);
                const double rhs_ms = experiment_support::milliseconds(
                    rhs_begin, experiment_support::Clock::now());
                cycles.push_back(solved.cycles);
                if (!solved.converged) ++failed;
                rhs_rows.push_back({
                    field.name,
                    candidate.name,
                    rhs_names[rhs_index],
                    std::to_string(solved.cycles),
                    solved.converged ? "1" : "0",
                    experiment_support::scientific(
                        solved.relative_residual),
                    experiment_support::fixed(rhs_ms)
                });
            }
            const double solve_ms = experiment_support::milliseconds(
                solves_begin, experiment_support::Clock::now());
            const double mean_cycles = std::accumulate(
                cycles.begin(), cycles.end(), 0.0) /
                static_cast<double>(cycles.size());
            const int maximum_cycles = *std::max_element(
                cycles.begin(), cycles.end());

            rows.push_back({
                field.name,
                candidate.name,
                experiment_support::scientific(rho),
                experiment_support::fixed(mean_cycles, 1),
                std::to_string(maximum_cycles),
                std::to_string(failed),
                std::to_string(candidate.prolongation.nnz()),
                experiment_support::fixed(candidate.build_ms),
                experiment_support::fixed(setup_ms),
                experiment_support::fixed(solve_ms),
                experiment_support::fixed(
                    candidate.build_ms + setup_ms + solve_ms)
            });
        }
    }

    experiment_support::Report report(
        "Experiment 6 - Spectral proxy and multiple right-hand sides");
    report.add_summary({
        {"Dataset", cfg.fields.string()},
        {"Fine grid", "h = 1/" + std::to_string(cfg.fine_intervals)},
        {"Coarse grid", "H = 1/" + std::to_string(cfg.coarse_intervals)},
        {"Right-hand sides", "6"},
        {"Spectral starts", "2"},
        {"Power iterations", std::to_string(cfg.spectral_iterations)},
        {"Maximum cycles", std::to_string(cfg.outer_max_cycles)},
        {"Threads", std::to_string(cfg.threads)}
    });
    report.add_note(
        "The RHS set contains the original manufactured solution, three "
        "random sine mixtures, and two point sources. Rho is the largest "
        "energy-norm power-iteration estimate over two deterministic starts.");
    const experiment_support::Row headers{
        "Field", "Method", "Rho", "Mean cycles", "Max cycles", "Failed",
        "P nnz", "Build ms", "Setup ms", "6 solves ms", "Total ms"};
    report.add_table(
        "RHS-independent and multi-query evaluation",
        headers,
        {20, 18, 11, 10, 9, 8, 10, 10, 10, 10, 10},
        rows, true);
    report.save("experiment6");
    experiment_support::write_csv("experiment6", headers, rows);
    experiment_support::write_csv(
        "experiment6_rhs",
        {"Field", "Method", "RHS", "Cycles", "Converged",
         "Final residual", "Solve ms"},
        rhs_rows);
    return 0;
}
