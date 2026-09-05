#include "experiment/problem.hpp"
#include "experiment/reporting.hpp"
#include "multigrid/global_pcg.hpp"
#include "version.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace {

struct PathPoint {
    int steps = 0;
    double energy = 0.0;
    double effective_factor = 1.0;
};

struct ScanSummary {
    int minimum_factor_steps = 0;
    double minimum_factor = 1.0;
    int adaptive_steps = 0;
    double adaptive_factor = 1.0;
    double reference_factor = 1.0;
    double reference_energy = 0.0;
    double maximum_energy_increment = 0.0;
    bool energy_monotone = false;
};

PathPoint measure(
    int steps, const tgi::SparseMatrix& matrix, const tgi::Vector& rhs,
    const tgi::SparseMatrix& prolongation, int threads,
    int observation_limit) {
    const tgi::TwoGridCycle cycle(matrix, prolongation, 1, threads);
    const auto solved = tgi::solve_two_grid(
        rhs, cycle, 1.0e-6, observation_limit);
    return {
        steps, cycle.setup_report().interpolation_energy,
        solved.effective_factor};
}

ScanSummary scan_case(
    const experiment_support::BasicConfig& config,
    const experiment_support::FieldCase& field,
    const std::string& output_name) {
    constexpr int observation_limit = 12000;
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);
    const auto problem = experiment_support::make_problem(grid, field, config);
    const auto initial = tgi::build_geometric_interpolation(grid);
    tgi::GlobalEnergyPcgPath path(
        grid, problem.matrix, initial.prolongation, config.threads);

    std::vector<PathPoint> points;
    points.reserve(static_cast<std::size_t>(config.fine_intervals));
    for (int steps = 1; steps <= config.fine_intervals; ++steps) {
        experiment_support::progress(
            field.name + " path scan " + std::to_string(steps) + "/" +
            std::to_string(config.fine_intervals));
        path.advance_to(steps);
        points.push_back(measure(
            steps, problem.matrix, problem.rhs, path.prolongation(),
            config.threads, observation_limit));
    }

    path.advance_until_relative_residual(1.0e-10);
    const PathPoint reference = measure(
        0, problem.matrix, problem.rhs, path.prolongation(), config.threads,
        experiment_support::maximum_two_grid_cycles);
    const double energy_scale = points.front().energy - reference.energy;
    experiment_support::Rows rows;
    rows.reserve(points.size());
    for (const PathPoint& point : points) {
        rows.push_back({
            std::to_string(point.steps),
            experiment_support::scientific(point.energy, 12),
            experiment_support::scientific(
                (point.energy - reference.energy) / energy_scale, 12),
            experiment_support::fixed(point.effective_factor, 9)});
    }
    experiment_support::save_csv(
        output_name,
        {"m", "Interpolation energy", "Normalized energy excess",
         "Effective factor"},
        rows);

    const auto minimum = std::min_element(
        points.begin(), points.end(),
        [](const PathPoint& left, const PathPoint& right) {
            return left.effective_factor < right.effective_factor;
        });
    const int adaptive_steps =
        tgi::adaptive_global_pcg_detail::select_steps(grid, problem.matrix);
    const auto adaptive = std::find_if(
        points.begin(), points.end(),
        [adaptive_steps](const PathPoint& point) {
            return point.steps == adaptive_steps;
        });
    double maximum_increment = -std::numeric_limits<double>::infinity();
    for (std::size_t index = 1; index < points.size(); ++index) {
        maximum_increment = std::max(
            maximum_increment,
            points[index].energy - points[index - 1U].energy);
    }
    const double monotonicity_tolerance =
        100.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::abs(points.front().energy));
    return {
        minimum->steps, minimum->effective_factor,
        adaptive_steps, adaptive->effective_factor,
        reference.effective_factor, reference.energy,
        maximum_increment, maximum_increment <= monotonicity_tolerance};
}

}

int main(int argc, char** argv) {
    int threads = 4;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument.rfind("--threads=", 0) == 0) {
            threads = std::stoi(argument.substr(10));
        }
    }

    experiment_support::BasicConfig config;
    config.fine_intervals = 128;
    config.coarse_intervals = 16;
    config.contrast = 1.0e4;
    config.threads = threads;
    const auto& topologies = experiment_support::channel_topologies();
    const ScanSummary cross = scan_case(
        config, topologies[0], "experiment2_cross_channel_path");
    const ScanSummary ring = scan_case(
        config, topologies[5], "experiment2_winding_ring_path");

    experiment_support::Rows summary_rows;
    for (const auto& item : {
             std::pair<std::string, ScanSummary>{"cross-channel", cross},
             {"winding-ring", ring}}) {
        const ScanSummary& value = item.second;
        summary_rows.push_back({
            item.first,
            std::to_string(value.minimum_factor_steps),
            experiment_support::fixed(value.minimum_factor, 7),
            std::to_string(value.adaptive_steps),
            experiment_support::fixed(value.adaptive_factor, 7),
            experiment_support::fixed(value.reference_factor, 7),
            experiment_support::scientific(value.reference_energy, 7),
            experiment_support::scientific(
                value.maximum_energy_increment, 7),
            value.energy_monotone ? "yes" : "no"});
    }

    experiment_support::Report report(
        "Finite-PCG energy and effective-factor path scans");
    report.add_summary({
        {"Version", std::string(tgi::version)},
        {"Problems", "128/16, contrast 1e4"},
        {"Scanned interval", "m=1,...,128"},
        {"Threads", std::to_string(threads)},
        {"Solve tolerance", "1e-6"},
        {"Reference column tolerance", "1e-10"}});
    report.add_table(
        "Path-scan summary",
        {"Topology", "Min m", "Min rho_eff", "Adaptive m",
         "Adaptive rho_eff", "Reference rho_eff", "Reference J",
         "Max delta J", "J monotone"},
        {15, 6, 11, 10, 16, 17, 14, 14, 10}, summary_rows);
    report.add_note(
        "For each topology, every integer m in the scanned interval is "
        "evaluated on one Jacobi-PCG path. J(W)=0.5 tr(P^T A P), and the "
        "normalized energy excess is (J(W_m)-J(W_ref))/(J(W_1)-J(W_ref)). "
        "rho_eff is computed from the observed residual history and is "
        "distinct from the energy-norm two-grid factor rho_TG. Minima are "
        "minima only within m=1,...,128. Max delta J is max_m "
        "[J(W_{m+1})-J(W_m)].");
    report.save("experiment2_step_scan");
    return 0;
}
