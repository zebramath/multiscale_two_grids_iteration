#include "experiment/common.hpp"
#include "multigrid/adaptive_support.hpp"
#include "multigrid/energy_interpolation.hpp"
#include "multigrid/residual_budget_support.hpp"
#include "multigrid/two_grid_solver.hpp"
#include "pde/diffusion_problem.hpp"

#include <cstdint>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Config {
    int fine_intervals = 64;
    int coarse_intervals = 8;
    int seed_count = 3;
    int threads = 4;
    int max_cycles = 10000;
    std::vector<double> contrasts{1.0e2, 1.0e4, 1.0e6};
};

struct FieldSpec {
    std::string name;
    tgi::CoefficientDistribution distribution;
};

struct Aggregate {
    std::string field;
    std::string contrast;
    std::string method;
    int samples = 0;
    int failures = 0;
    double cycle_sum = 0.0;
    double cycle_squared_sum = 0.0;
    int maximum_cycles = 0;
    double extra_sum = 0.0;
    double total_ms_sum = 0.0;
};

tgi::InterpolationOptions energy_options(
    int layers, int threads) {
    tgi::InterpolationOptions options;
    options.strategy = tgi::InterpolationStrategy::LocalEnergyMinimum;
    options.patch_layers = layers;
    options.local_tolerance = 1.0e-3;
    options.local_max_iterations = 20000;
    options.thread_count = threads;
    return options;
}

} // namespace

int main(int argc, char** argv) {
    Config cfg;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument.rfind("--fine=", 0) == 0) {
            cfg.fine_intervals = std::stoi(argument.substr(7));
        } else if (argument.rfind("--coarse=", 0) == 0) {
            cfg.coarse_intervals = std::stoi(argument.substr(9));
        } else if (argument.rfind("--seeds=", 0) == 0) {
            cfg.seed_count = std::stoi(argument.substr(8));
        } else if (argument.rfind("--threads=", 0) == 0) {
            cfg.threads = std::stoi(argument.substr(10));
        } else if (argument.rfind("--contrast=", 0) == 0) {
            cfg.contrasts = {std::stod(argument.substr(11))};
        } else if (argument.rfind("--max-cycles=", 0) == 0) {
            cfg.max_cycles = std::stoi(argument.substr(13));
        } else {
            throw std::invalid_argument(
                "usage: experiment7 [--fine=N] [--coarse=N] [--seeds=N] "
                "[--contrast=X] [--threads=N] [--max-cycles=N]");
        }
    }
    if (cfg.fine_intervals <= 0 || cfg.coarse_intervals <= 0 ||
        cfg.fine_intervals % cfg.coarse_intervals != 0 ||
        cfg.seed_count <= 0) {
        throw std::invalid_argument("invalid sweep dimensions");
    }

    const int ratio = cfg.fine_intervals / cfg.coarse_intervals;
    const tgi::StructuredGrid grid(cfg.fine_intervals - 1, ratio);
    const std::vector<FieldSpec> fields{
        {"continuous", tgi::CoefficientDistribution::RandomContinuous},
        {"channel", tgi::CoefficientDistribution::ChannelizedBinary},
        {"block_random", tgi::CoefficientDistribution::RandomBinaryCheckerboard}
    };
    experiment_support::Rows rows;
    std::map<std::string, Aggregate> aggregates;

    for (double contrast : cfg.contrasts) {
        for (int seed = 1; seed <= cfg.seed_count; ++seed) {
            for (const FieldSpec& specification : fields) {
                tgi::CoefficientOptions coefficient_options;
                coefficient_options.distribution = specification.distribution;
                coefficient_options.contrast = contrast;
                coefficient_options.seed = static_cast<std::uint64_t>(seed);
                coefficient_options.checkerboard_block_size = ratio;
                coefficient_options.channel_background_block_size = ratio;
                coefficient_options.channel_width_fine_cells = 2;
                const auto coefficient = tgi::make_coefficient(
                    grid, coefficient_options);
                const auto a = tgi::assemble_diffusion(
                    grid, coefficient.values);
                const auto rhs = a.multiply(
                    experiment_support::manufactured_solution(grid));

                const auto local2_begin = experiment_support::Clock::now();
                const auto local2 = tgi::build_interpolation(
                    grid, a, energy_options(2, cfg.threads));
                const double local2_ms = experiment_support::milliseconds(
                    local2_begin, experiment_support::Clock::now());

                const auto local3_begin = experiment_support::Clock::now();
                const auto local3 = tgi::build_interpolation(
                    grid, a, energy_options(3, cfg.threads));
                const double local3_ms = experiment_support::milliseconds(
                    local3_begin, experiment_support::Clock::now());

                tgi::ResidualStrongSupportOptions fixed_support;
                fixed_support.base_patch_layers = 2;
                fixed_support.maximum_extra_nodes_per_column = 64;
                fixed_support.maximum_graph_hops = 4 * grid.fine_n();
                fixed_support.strong_edge_fraction = 0.25;
                fixed_support.thread_count = cfg.threads;
                const auto fixed_begin = experiment_support::Clock::now();
                const auto fixed = tgi::build_residual_strong_supports(
                    grid, a, local2.prolongation, fixed_support);
                std::vector<unsigned char> fixed_mask(
                    static_cast<std::size_t>(grid.coarse_size()), 0U);
                for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
                    if (fixed.supports[static_cast<std::size_t>(coarse)].size() >
                        grid.patch_f_nodes(coarse, 2).size()) {
                        fixed_mask[static_cast<std::size_t>(coarse)] = 1U;
                    }
                }
                auto fixed_options = energy_options(2, cfg.threads);
                fixed_options.local_tolerance = 1.0e-6;
                const auto fixed_interpolation =
                    tgi::refine_selected_energy_interpolation_on_supports(
                        grid, a, fixed.supports, fixed_mask,
                        local2.prolongation, fixed_options);
                const double fixed_ms = local2_ms +
                    experiment_support::milliseconds(
                        fixed_begin, experiment_support::Clock::now());

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
                const auto adaptive =
                    tgi::build_residual_budget_interpolation(
                        grid, a, local2.prolongation,
                        energy_options(2, cfg.threads), support);

                struct Method {
                    std::string name;
                    const tgi::SparseMatrix* p;
                    double build_ms;
                    int extras;
                };
                const std::vector<Method> methods{
                    {"local-2", &local2.prolongation, local2_ms, 0},
                    {"local-3", &local3.prolongation, local3_ms, 0},
                    {"fixed-K64", &fixed_interpolation.prolongation,
                     fixed_ms, fixed.report.total_extra_nodes},
                    {"budget-row", &adaptive.prolongation,
                     local2_ms + adaptive.report.total_ms,
                     adaptive.report.total_extra_nodes}
                };
                for (const Method& method : methods) {
                    const auto begin = experiment_support::Clock::now();
                    const tgi::TwoGridCycle cycle(
                        a, *method.p, 1, cfg.threads);
                    const auto solved = tgi::solve_two_grid(
                        a, rhs, cycle, 1.0e-6, cfg.max_cycles);
                    const double total_ms = method.build_ms +
                        experiment_support::milliseconds(
                            begin, experiment_support::Clock::now());
                    const std::string contrast_text =
                        experiment_support::scientific(contrast, 0);
                    const std::string aggregate_key = specification.name +
                        "|" + contrast_text + "|" + method.name;
                    auto& aggregate = aggregates[aggregate_key];
                    aggregate.field = specification.name;
                    aggregate.contrast = contrast_text;
                    aggregate.method = method.name;
                    ++aggregate.samples;
                    if (!solved.converged) ++aggregate.failures;
                    aggregate.cycle_sum += static_cast<double>(solved.cycles);
                    aggregate.cycle_squared_sum +=
                        static_cast<double>(solved.cycles) *
                        static_cast<double>(solved.cycles);
                    aggregate.maximum_cycles = std::max(
                        aggregate.maximum_cycles, solved.cycles);
                    aggregate.extra_sum += static_cast<double>(method.extras);
                    aggregate.total_ms_sum += total_ms;
                    rows.push_back({
                        specification.name,
                        contrast_text,
                        std::to_string(seed),
                        method.name,
                        std::to_string(method.extras),
                        std::to_string(method.p->nnz()),
                        solved.converged
                            ? std::to_string(solved.cycles)
                            : "failed@" + std::to_string(solved.cycles),
                        experiment_support::fixed(method.build_ms),
                        experiment_support::fixed(total_ms)
                    });
                }
            }
        }
    }

    experiment_support::Report report(
        "Experiment 7 - Seed, contrast, and scale-ready robustness sweep");
    report.add_summary({
        {"Fine grid", "h = 1/" + std::to_string(cfg.fine_intervals)},
        {"Coarse grid", "H = 1/" + std::to_string(cfg.coarse_intervals)},
        {"H/h", std::to_string(ratio)},
        {"Seeds", std::to_string(cfg.seed_count)},
        {"Contrasts", std::to_string(cfg.contrasts.size())},
        {"Maximum cycles", std::to_string(cfg.max_cycles)},
        {"Threads", std::to_string(cfg.threads)}
    });
    report.add_note(
        "Use --fine, --coarse, --seeds, and --contrast to create targeted "
        "sweeps. The default 64/8 configuration is intentionally inexpensive; "
        "repeat at 128/16 and 256/16 for scale studies.");
    const experiment_support::Row headers{
        "Field", "Contrast", "Seed", "Method", "Extra", "P nnz",
        "Cycles", "Build ms", "Total ms"};
    report.add_table(
        "Per-instance robustness results",
        headers,
        {16, 10, 6, 14, 8, 10, 12, 10, 10},
        rows, true);
    experiment_support::Rows summary_rows;
    for (const auto& [key, aggregate] : aggregates) {
        (void)key;
        const double count = static_cast<double>(aggregate.samples);
        const double mean = aggregate.cycle_sum / count;
        const double variance = std::max(
            0.0,
            aggregate.cycle_squared_sum / count - mean * mean);
        summary_rows.push_back({
            aggregate.field,
            aggregate.contrast,
            aggregate.method,
            std::to_string(aggregate.samples),
            std::to_string(aggregate.failures),
            experiment_support::fixed(mean, 1),
            experiment_support::fixed(std::sqrt(variance), 1),
            std::to_string(aggregate.maximum_cycles),
            experiment_support::fixed(aggregate.extra_sum / count, 1),
            experiment_support::fixed(aggregate.total_ms_sum / count)
        });
    }
    const experiment_support::Row summary_headers{
        "Field", "Contrast", "Method", "Samples", "Failed",
        "Mean cycles", "Std cycles", "Worst cycles", "Mean extra",
        "Mean total ms"};
    report.add_table(
        "Across-seed summary",
        summary_headers,
        {16, 10, 14, 8, 8, 12, 11, 12, 11, 13},
        summary_rows, true);
    report.save("experiment7");
    experiment_support::write_csv("experiment7", headers, rows);
    experiment_support::write_csv(
        "experiment7_summary", summary_headers, summary_rows);
    return 0;
}
