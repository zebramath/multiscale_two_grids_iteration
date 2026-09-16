#include "experiment/reporting.hpp"
#include "multigrid/global_pcg.hpp"
#include "multigrid/spectral_diagnostics.hpp"
#include "pde/diffusion_problem.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <string>
namespace {
using Clock = std::chrono::steady_clock;
struct ScaleCase {
    std::string regime;
    int fine_intervals = 0;
    int coarse_intervals = 0;
};
double elapsed_ms(const Clock::time_point& begin) {
    return std::chrono::duration<double, std::milli>(
               Clock::now() - begin)
        .count();
}
int center_coarse_column(const tgi::StructuredGrid& grid) {
    const int center = grid.coarse_n() / 2;
    return grid.coarse_id(center, center);
}
void run_checkpoint(
    const ScaleCase& item, const std::string& rule, int steps,
    const tgi::StructuredGrid& grid, const tgi::SparseMatrix& matrix,
    const tgi::Vector& rhs, const tgi::SparseMatrix& initial,
    int threads, int spectral_iterations, experiment_support::Rows& rows) {
    const auto begin = Clock::now();
    tgi::GlobalEnergyPcgPath path(grid, matrix, initial, threads);
    path.advance_to(steps);
    const auto propagation = path.column_propagation_report(
        center_coarse_column(grid), 1.0e-14);
    const tgi::SparseMatrix prolongation = path.prolongation();
    const tgi::TwoGridCycle cycle(matrix, prolongation, 1, threads);
    const double setup_ms = elapsed_ms(begin);
    const auto spectral = tgi::estimate_two_grid_spectral_radius(
        matrix, cycle, 0x13198a2e03707344ULL, spectral_iterations);
    const auto solved = tgi::solve_two_grid(
        rhs, cycle, 1.0e-6, 20000);
    const double correction_percent = 100.0 *
        static_cast<double>(propagation.correction_support) /
        static_cast<double>(propagation.f_unknowns);
    const double reachable_percent = 100.0 *
        static_cast<double>(propagation.reachable_support) /
        static_cast<double>(propagation.f_unknowns);
    rows.push_back({
        item.regime, std::to_string(item.fine_intervals),
        std::to_string(item.coarse_intervals),
        std::to_string(grid.ratio()), rule, std::to_string(steps),
        experiment_support::fixed(
            static_cast<double>(steps) * grid.h(), 6),
        experiment_support::fixed(correction_percent, 4),
        experiment_support::fixed(reachable_percent, 4),
        std::to_string(propagation.maximum_correction_distance),
        std::to_string(propagation.theoretical_distance_bound),
        std::to_string(propagation.bound_violations),
        experiment_support::fixed(
            experiment_support::interpolation_density_percent(prolongation),
            5),
        std::to_string(cycle.coarse_matrix().nnz()),
        experiment_support::fixed(spectral.rho, 9),
        experiment_support::scientific(spectral.stage_difference, 3),
        std::to_string(solved.cycles),
        experiment_support::fixed(solved.effective_factor, 9),
        tgi::stationary_status_name(solved.status),
        experiment_support::fixed(setup_ms, 2)});
}
}
int main(int argc, char** argv) {
    int threads = 4;
    int spectral_iterations = 160;
    bool quick = false;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--quick") quick = true;
        if (argument.rfind("--threads=", 0) == 0) {
            threads = std::stoi(argument.substr(10));
        }
        if (argument.rfind("--spectral-iterations=", 0) == 0) {
            spectral_iterations = std::stoi(argument.substr(22));
        }
    }
    const std::array<int, 4> resolutions{{32, 64, 128, 256}};
    const std::size_t count = quick ? 2U : resolutions.size();
    experiment_support::Rows rows;
    for (const std::string& regime : {
             std::string("fixed-H"), std::string("fixed-q")}) {
        for (std::size_t index = 0; index < count; ++index) {
            const int fine = resolutions[index];
            const int coarse = regime == "fixed-H" ? 8 : fine / 8;
            const ScaleCase item{regime, fine, coarse};
            experiment_support::progress(
                regime + " scale " + std::to_string(fine) + "/" +
                std::to_string(coarse));
            const tgi::StructuredGrid grid(fine - 1, fine / coarse);
            tgi::FixedPhysicalCoefficientOptions options;
            options.distribution =
                tgi::CoefficientDistribution::ChannelizedBinary;
            options.contrast = 1.0e4;
            options.seed = 1;
            options.background_blocks_per_direction = 8;
            options.channel_width = 1.0 / 16.0;
            const auto coefficient =
                tgi::make_fixed_physical_coefficient(grid, options);
            const auto matrix =
                tgi::assemble_diffusion(grid, coefficient.values);
            const tgi::Vector rhs(
                static_cast<std::size_t>(grid.fine_size()), 1.0);
            const auto geometric = tgi::build_geometric_interpolation(grid);
            const int q_steps = grid.ratio();
            const int h_steps = std::max(1, fine / 4);
            run_checkpoint(
                item, "m=q", q_steps, grid, matrix, rhs,
                geometric, threads, spectral_iterations, rows);
            run_checkpoint(
                item, "m=(1/h)/4", h_steps, grid, matrix, rhs,
                geometric, threads, spectral_iterations, rows);
        }
    }
    experiment_support::save_csv(
        "experiment3_scaling_propagation",
        {"regime", "1/h", "1/H", "q", "rule", "m", "m*h",
         "correction_support_percent", "reachable_support_percent",
         "max_correction_distance", "theoretical_bound",
         "bound_violations", "P_density_percent", "Ac_nnz", "rho_TG",
         "spectral_stage_difference", "cycles", "rho_eff", "status",
         "setup_ms"},
        rows);
    experiment_support::Report report(
        "Fixed-H/fixed-q scaling and finite propagation");
    report.add_summary({
        {"Mode", quick ? "quick" : "full"},
        {"Resolutions", quick ? "32, 64" : "32, 64, 128, 256"},
        {"Fixed-H family", "1/H=8"},
        {"Fixed-q family", "q=8"},
        {"Coefficient field", "fixed physical cross-channel, contrast 1e4"},
        {"Lanczos iterations", std::to_string(spectral_iterations)},
        {"Threads", std::to_string(threads)}});
    report.add_note(
        "The two families separate growth of the coarsening ratio q from "
        "growth of the global resolution 1/h.  Each family compares m=q "
        "with m=(1/h)/4 on the same fixed-physical coefficient rule.  The "
        "support diagnostic uses the central coarse column and graph "
        "distance in B=A_FF.  In exact arithmetic the PCG correction after "
        "m steps stays inside the radius-(m-1) neighborhood of the initial "
        "residual; bound_violations should therefore be zero.  A 1e-14 "
        "threshold identifies floating-point support.");
    report.add_table(
        "Scaling measurements",
        {"Regime", "1/h", "1/H", "q", "Rule", "m", "m h",
         "Corr %", "Reach %", "dmax", "bound", "viol", "P %",
         "Ac nnz", "rho_TG", "stage d", "cycles", "rho_eff",
         "status", "setup ms"},
        {9, 5, 5, 4, 11, 5, 8, 8, 8, 6, 6, 5, 8, 9, 11, 10, 8,
         11, 10, 10},
        rows, true);
    report.save("experiment3_scaling_propagation");
    return 0;
}
