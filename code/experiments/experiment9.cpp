#include "experiment/study.hpp"
#include "multigrid/adaptive_global_pcg.hpp"

#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct SweepCase {
    int fine = 32;
    int coarse = 8;
    double contrast = 1.0e4;
    std::uint64_t seed = 1;
    experiment_support::FieldCase field;
};

experiment_support::Row measure(
    const SweepCase& item, const std::string& method,
    const std::string& parameter, const tgi::SparseMatrix& a,
    const tgi::Vector& rhs, const tgi::SparseMatrix& p,
    double interpolation_ms, int threads, int max_cycles,
    const std::string& selector_cycles = "-") {
    const auto metric = experiment_support::evaluate_two_grid(
        a, rhs, p, threads, 1.0e-6, max_cycles);
    return {
        std::to_string(item.fine), std::to_string(item.coarse),
        std::to_string(item.fine / item.coarse),
        experiment_support::scientific(item.contrast, 0),
        std::to_string(item.seed), item.field.name, method, parameter,
        experiment_support::fixed(
            experiment_support::interpolation_density_percent(p), 4),
        selector_cycles,
        metric.converged ? std::to_string(metric.cycles)
            : "failed@" + std::to_string(metric.cycles),
        experiment_support::fixed(
            interpolation_ms + metric.coarse_setup_ms),
        experiment_support::fixed(interpolation_ms + metric.total_ms)};
}

}

int main(int argc, char** argv) {
    bool quick = false;
    int threads = 4;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--quick") {
            quick = true;
        } else if (argument.rfind("--threads=", 0) == 0) {
            threads = std::stoi(argument.substr(10));
        }
    }

    const std::array<std::pair<int, int>, 3> grids{
        std::pair<int, int>{32, 8}, {64, 8}, {64, 16}};
    const std::array<double, 3> contrasts{1.0e2, 1.0e4, 1.0e6};
    const std::array<std::uint64_t, 2> seeds{1, 17};
    const auto& topologies = experiment_support::channel_topologies();
    std::vector<SweepCase> cases;
    int topology = 0;
    for (const auto& [fine, coarse] : grids) {
        for (double contrast : contrasts) {
            for (std::uint64_t seed : seeds) {
                cases.push_back(
                    {fine, coarse, contrast, seed,
                     topologies[static_cast<std::size_t>(topology) %
                                topologies.size()]});
                ++topology;
            }
        }
    }
    if (quick) {
        std::vector<SweepCase> selected;
        for (std::size_t index : {0U, 4U, 6U, 10U, 12U, 16U}) {
            selected.push_back(cases[index]);
        }
        cases = std::move(selected);
    }

    const experiment_support::Row headers{
        "1/h", "1/H", "H/h", "Contrast", "Seed", "Topology",
        "Method", "Parameter", "P density %", "selector estimate", "Cycles",
        "Setup ms", "Total ms"};
    experiment_support::Rows rows;
    const experiment_support::Row history_headers{
        "1/h", "1/H", "Contrast", "Seed", "Topology", "m",
        "Pilot residual", "Tail rho", "Forecast", "Selected"};
    experiment_support::Rows history_rows;
    constexpr int max_cycles = 12000;
    int case_index = 0;
    for (const SweepCase& item : cases) {
        ++case_index;
        experiment_support::progress(
            "experiment9 case " + std::to_string(case_index) + "/" +
            std::to_string(cases.size()) + ": " + item.field.name);
        experiment_support::BasicConfig config;
        config.fine_intervals = item.fine;
        config.coarse_intervals = item.coarse;
        config.contrast = item.contrast;
        config.threads = threads;
        config.max_cycles = max_cycles;
        const tgi::StructuredGrid grid = experiment_support::make_grid(config);
        const auto problem = experiment_support::make_problem(
            grid, item.field, config, item.seed);
        const auto geometric = experiment_support::geometric_interpolation(
            grid, problem.matrix);
        const double geometric_ms = geometric.report.timing.total_ms;
        rows.push_back(measure(
            item, "geometric", "P_G", problem.matrix, problem.rhs,
            geometric.prolongation, geometric_ms, threads, max_cycles));

        auto fixed_options = experiment_support::energy_options(
            0, threads, 0.0);
        fixed_options.local_max_iterations = 40;
        fixed_options.require_convergence = false;
        fixed_options.drop_tolerance = 0.0;
        const auto fixed = tgi::refine_global_energy_interpolation(
            grid, problem.matrix, geometric.prolongation, fixed_options);
        rows.push_back(measure(
            item, "PCG-fixed", "m=40", problem.matrix, problem.rhs,
            fixed.prolongation,
            geometric_ms + fixed.report.timing.total_ms,
            threads, max_cycles));

        tgi::AdaptiveGlobalPcgOptions adaptive_options;
        adaptive_options.minimum_steps = 12;
        adaptive_options.maximum_steps = 56;
        adaptive_options.maximum_cycles = max_cycles;
        adaptive_options.thread_count = threads;
        const auto adaptive = tgi::build_adaptive_global_pcg_interpolation(
            grid, problem.matrix, geometric.prolongation, adaptive_options,
            &problem.rhs);
        for (const auto& checkpoint : adaptive.report.history) {
            history_rows.push_back({
                std::to_string(item.fine), std::to_string(item.coarse),
                experiment_support::scientific(item.contrast, 0),
                std::to_string(item.seed), item.field.name,
                std::to_string(checkpoint.steps),
                experiment_support::scientific(
                    checkpoint.pilot_relative_residual, 3),
                experiment_support::fixed(checkpoint.rho_rhs_pilot, 6),
                std::to_string(checkpoint.predicted_cycles),
                checkpoint.selected ? "yes" : "no"});
        }
        rows.push_back(measure(
            item, "PCG-adaptive",
            "m=" + std::to_string(adaptive.report.selected_steps),
            problem.matrix, problem.rhs, adaptive.prolongation,
            geometric_ms + adaptive.report.selection_wall_ms,
            threads, max_cycles,
            std::to_string(adaptive.report.estimated_selected_cycles)));

    }

    experiment_support::Report report(
        "Robustness sweep over seeds, channel topology, contrast and H/h");
    report.add_summary({
        {"Version", std::string(tgi::version)},
        {"Cases", std::to_string(cases.size())},
        {"Mode", quick ? "quick" : "full"},
        {"Solve tolerance", "1e-6"},
        {"Maximum cycles", std::to_string(max_cycles)}});
    report.add_note(
        "The full matrix contains 18 coefficient problems: three grids, "
        "three contrasts, two seeds, and four channel topologies assigned in "
        "balanced rotation. Each problem compares the geometric basis, fixed "
        "m=40 PCG, and the midpoint adaptive selector without case-specific "
        "tuning.");
    report.add_table(
        "Robustness matrix", headers,
        {5, 5, 5, 10, 6, 20, 14, 10, 11, 10, 12, 10, 10}, rows, true);
    report.save("experiment9");
    experiment_support::write_csv("experiment9", headers, rows);
    experiment_support::write_csv(
        "experiment9_history", history_headers, history_rows);
    return 0;
}
