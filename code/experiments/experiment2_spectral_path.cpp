#include "experiment/problem.hpp"
#include "experiment/reporting.hpp"
#include "multigrid/global_pcg.hpp"
#include "multigrid/spectral_diagnostics.hpp"

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr int scan_cycle_limit = 12000;

struct PathPoint {
    int steps = 0;
    double rho_tg = 1.0;
    double spectral_stage_difference = 0.0;
    double effective_factor = 1.0;
    int cycles = 0;
    std::string status;
};

struct ScanSummary {
    int minimum_spectral_steps = 0;
    double minimum_spectral_factor = 1.0;
    int minimum_effective_steps = 0;
    double minimum_effective_factor = 1.0;
    double endpoint_spectral_factor = 1.0;
    double endpoint_stage_difference = 0.0;
    double endpoint_effective_factor = 1.0;
    int endpoint_cycles = 0;
    std::string endpoint_status;
    double maximum_stage_difference = 0.0;
};

PathPoint measure(
    int steps, const tgi::SparseMatrix& matrix, const tgi::Vector& rhs,
    const tgi::SparseMatrix& prolongation, int threads,
    int spectral_iterations) {
    const tgi::TwoGridCycle cycle(matrix, prolongation, 1, threads);
    const auto spectral = tgi::estimate_two_grid_spectral_radius(
        matrix, cycle, 0x243f6a8885a308d3ULL,
        spectral_iterations);
    const auto solved = tgi::solve_two_grid(
        rhs, cycle, 1.0e-6, scan_cycle_limit);
    return {
        steps, spectral.rho, spectral.stage_difference,
        solved.effective_factor, solved.cycles,
        tgi::stationary_status_name(solved.status)};
}

ScanSummary scan_case(
    const experiment_support::BasicConfig& config,
    const experiment_support::FieldCase& field,
    const std::string& output_name, int spectral_iterations,
    int maximum_steps) {
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);
    const auto problem = experiment_support::make_problem(grid, field, config);
    const auto initial = tgi::build_geometric_interpolation(grid);
    tgi::GlobalEnergyPcgPath path(
        grid, problem.matrix, initial, config.threads);

    std::vector<PathPoint> points;
    points.reserve(static_cast<std::size_t>(maximum_steps));
    for (int steps = 1; steps <= maximum_steps; ++steps) {
        experiment_support::progress(
            field.name + " spectral path " + std::to_string(steps) + "/" +
            std::to_string(maximum_steps));
        path.advance_to(steps);
        const tgi::SparseMatrix prolongation = path.prolongation();
        points.push_back(measure(
            steps, problem.matrix, problem.rhs, prolongation,
            config.threads, spectral_iterations));
    }

    tgi::GlobalEnergyPcgPath endpoint_path(
        grid, problem.matrix, initial, config.threads);
    endpoint_path.advance_until_relative_residual(1.0e-10);
    const tgi::SparseMatrix endpoint = endpoint_path.prolongation();
    const PathPoint endpoint_point = measure(
        0, problem.matrix, problem.rhs, endpoint, config.threads,
        spectral_iterations);

    experiment_support::Rows rows;
    rows.reserve(points.size());
    for (const PathPoint& point : points) {
        rows.push_back({
            std::to_string(point.steps),
            experiment_support::fixed(point.rho_tg, 10),
            experiment_support::scientific(
                point.spectral_stage_difference, 5),
            experiment_support::fixed(point.effective_factor, 10),
            std::to_string(point.cycles), point.status});
    }
    experiment_support::save_csv(
        output_name,
        {"m", "rho_TG", "spectral_stage_difference", "rho_eff",
         "cycles", "status"},
        rows);

    const auto spectral_minimum = std::min_element(
        points.begin(), points.end(),
        [](const PathPoint& left, const PathPoint& right) {
            return left.rho_tg < right.rho_tg;
        });
    const auto effective_minimum = std::min_element(
        points.begin(), points.end(),
        [](const PathPoint& left, const PathPoint& right) {
            return left.effective_factor < right.effective_factor;
        });
    double maximum_stage_difference = 0.0;
    for (const PathPoint& point : points) {
        maximum_stage_difference = std::max(
            maximum_stage_difference, point.spectral_stage_difference);
    }
    return {
        spectral_minimum->steps, spectral_minimum->rho_tg,
        effective_minimum->steps, effective_minimum->effective_factor,
        endpoint_point.rho_tg, endpoint_point.spectral_stage_difference,
        endpoint_point.effective_factor, endpoint_point.cycles,
        endpoint_point.status,
        maximum_stage_difference};
}

}

int main(int argc, char** argv) {
    int threads = 4;
    int spectral_iterations = 200;
    bool quick = false;
    std::string topology = "both";
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--quick") quick = true;
        if (argument.rfind("--threads=", 0) == 0) {
            threads = std::stoi(argument.substr(10));
        }
        if (argument.rfind("--spectral-iterations=", 0) == 0) {
            spectral_iterations = std::stoi(argument.substr(22));
        }
        if (argument.rfind("--topology=", 0) == 0) {
            topology = argument.substr(11);
        }
    }

    experiment_support::BasicConfig config;
    config.fine_intervals = quick ? 32 : 128;
    config.coarse_intervals = quick ? 8 : 16;
    config.contrast = 1.0e4;
    config.threads = threads;
    const int maximum_steps = config.fine_intervals;
    const auto& topologies = experiment_support::channel_topologies();
    if (topology != "both" && topology != "cross" && topology != "ring") {
        throw std::invalid_argument(
            "--topology must be cross, ring, or both");
    }
    std::vector<std::pair<std::string, ScanSummary>> summaries;
    if (topology == "both" || topology == "cross") {
        summaries.push_back({
            "cross-channel",
            scan_case(
                config, topologies[0],
                "experiment2_cross_channel_spectral_path",
                spectral_iterations, maximum_steps)});
    }
    if (topology == "both" || topology == "ring") {
        summaries.push_back({
            "winding-ring",
            scan_case(
                config, topologies[5],
                "experiment2_winding_ring_spectral_path",
                spectral_iterations, maximum_steps)});
    }
    experiment_support::Rows summary_rows;
    experiment_support::Rows endpoint_rows;
    for (const auto& item : summaries) {
        const ScanSummary& value = item.second;
        summary_rows.push_back({
            item.first,
            std::to_string(value.minimum_spectral_steps),
            experiment_support::fixed(value.minimum_spectral_factor, 9),
            std::to_string(value.minimum_effective_steps),
            experiment_support::fixed(value.minimum_effective_factor, 9),
            experiment_support::fixed(value.endpoint_spectral_factor, 9),
            experiment_support::fixed(value.endpoint_effective_factor, 9),
            experiment_support::scientific(
                value.maximum_stage_difference, 3)});
        endpoint_rows.push_back({
            item.first,
            experiment_support::fixed(value.endpoint_spectral_factor, 10),
            experiment_support::scientific(
                value.endpoint_stage_difference, 5),
            experiment_support::fixed(value.endpoint_effective_factor, 10),
            std::to_string(value.endpoint_cycles), value.endpoint_status});
    }
    experiment_support::save_csv(
        "experiment2_spectral_endpoints",
        {"topology", "rho_TG", "spectral_stage_difference", "rho_eff",
         "cycles", "status"}, endpoint_rows);

    experiment_support::Report report(
        "True two-grid spectral-radius paths along finite PCG");
    report.add_summary({
        {"Mode", quick ? "quick" : "full"},
        {"Problems", quick ? "32/8, contrast 1e4"
                            : "128/16, contrast 1e4"},
        {"Scanned interval", "m=1,...," + std::to_string(maximum_steps)},
        {"Lanczos iterations", std::to_string(spectral_iterations)},
        {"Threads", std::to_string(threads)},
        {"RHS solve tolerance", "1e-6"},
        {"Endpoint column tolerance", "1e-10"}});
    report.add_table(
        "Path-scan summary",
        {"Topology", "Min rho m", "Min rho_TG", "Min eff m",
         "Min rho_eff", "Endpoint rho_TG", "Endpoint rho_eff",
         "Max stage diff"},
        {15, 10, 13, 10, 13, 16, 17, 15}, summary_rows);
    report.add_note(
        "rho_TG is obtained by cold-start matrix-free Lanczos iteration in "
        "the A-energy norm for the symmetric forward-GS/coarse/backward-GS "
        "error propagator.  The stage difference is the absolute change "
        "between the halfway and final Rayleigh-factor estimates.  rho_eff "
        "is a separate constant-RHS residual-history statistic.");
    report.save(
        topology == "both" ? "experiment2_spectral_path"
                           : "experiment2_spectral_path_" + topology);
    return 0;
}
