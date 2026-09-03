#include "experiment/problem.hpp"
#include "experiment/reporting.hpp"
#include "multigrid/global_pcg.hpp"
#include "version.hpp"

#include <array>
#include <cstddef>
#include <string>

namespace {

struct RefinementLevel {
    int fine;
    int coarse;
};

experiment_support::Row measurement_row(
    int fine, int coarse, const std::string& method,
    const std::string& parameter, const tgi::SparseMatrix& prolongation,
    const tgi::TwoGridCycle& cycle, const tgi::Vector& rhs,
    int maximum_cycles) {
    const auto solved = tgi::solve_two_grid(
        rhs, cycle, 1.0e-6, maximum_cycles);
    return {
        std::to_string(fine),
        std::to_string(coarse),
        std::to_string(fine / coarse),
        method,
        parameter,
        experiment_support::scientific(
            cycle.setup_report().interpolation_energy, 4),
        experiment_support::fixed(
            experiment_support::interpolation_density_percent(
                prolongation), 4),
        std::to_string(solved.cycles),
        tgi::stationary_status_name(solved.status),
        experiment_support::scientific(solved.relative_residual, 2),
        experiment_support::fixed(solved.effective_factor, 6),
        experiment_support::fixed(solved.tail_factor, 6)};
}

}

int main(int argc, char** argv) {
    int threads = 4;
    bool quick = false;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--quick") quick = true;
        if (argument.rfind("--threads=", 0) == 0) {
            threads = std::stoi(argument.substr(10));
        }
    }

    const std::array<RefinementLevel, 3> levels{{
        {32, 4}, {64, 8}, {128, 16}
    }};
    const std::size_t level_count = quick ? 2U : levels.size();
    constexpr double contrast = 1.0e4;
    constexpr double physical_channel_width = 1.0 / 16.0;
    constexpr int physical_background_blocks = 8;
    constexpr std::uint64_t seed = 1;
    const auto field_case =
        experiment_support::channel_topologies().front();
    experiment_support::Rows rows;
    tgi::Vector previous_coefficient;
    int previous_intervals = 0;
    int shared_node_mismatches = 0;

    for (std::size_t level_index = 0;
         level_index < level_count; ++level_index) {
        const RefinementLevel level = levels[level_index];
        experiment_support::progress(
            "fixed-physical refinement " +
            std::to_string(level_index + 1U) + "/" +
            std::to_string(level_count) + ": n=" +
            std::to_string(level.fine));
        experiment_support::BasicConfig config;
        config.fine_intervals = level.fine;
        config.coarse_intervals = level.coarse;
        config.contrast = contrast;
        config.threads = threads;
        const tgi::StructuredGrid grid =
            experiment_support::make_grid(config);
        tgi::FixedPhysicalCoefficientOptions options;
        options.distribution = field_case.distribution;
        options.contrast = contrast;
        options.seed = seed;
        options.channel_width = physical_channel_width;
        options.background_blocks_per_direction =
            physical_background_blocks;
        const auto coefficient =
            tgi::make_fixed_physical_coefficient(grid, options);

        if (!previous_coefficient.empty()) {
            const int refinement_ratio = level.fine / previous_intervals;
            const int previous_n = previous_intervals - 1;
            for (int iy = 0; iy < previous_n; ++iy) {
                for (int ix = 0; ix < previous_n; ++ix) {
                    const int previous_id = iy * previous_n + ix;
                    const int fine_ix = refinement_ratio * (ix + 1) - 1;
                    const int fine_iy = refinement_ratio * (iy + 1) - 1;
                    const int fine_id = grid.fine_id(fine_ix, fine_iy);
                    if (previous_coefficient[
                            static_cast<std::size_t>(previous_id)] !=
                        coefficient.values[
                            static_cast<std::size_t>(fine_id)]) {
                        ++shared_node_mismatches;
                    }
                }
            }
        }
        previous_coefficient = coefficient.values;
        previous_intervals = level.fine;

        const tgi::SparseMatrix matrix =
            tgi::assemble_diffusion(grid, coefficient.values);
        const tgi::Vector rhs(
            static_cast<std::size_t>(grid.fine_size()), 1.0);
        const auto initial = tgi::build_geometric_interpolation(grid);

        const auto adaptive = tgi::build_adaptive_global_pcg_interpolation(
            grid, matrix, initial.prolongation, threads);
        rows.push_back(measurement_row(
            level.fine, level.coarse, "adaptive",
            "m=" + std::to_string(adaptive.report.selected_steps) +
                "; m/n=" + experiment_support::fixed(
                    static_cast<double>(adaptive.report.selected_steps) /
                        static_cast<double>(level.fine), 3),
            *adaptive.prolongation, *adaptive.cycle, rhs,
            experiment_support::maximum_two_grid_cycles));

        const auto reference = experiment_support::build_global_reference(
            grid, matrix, threads);
        const tgi::TwoGridCycle reference_cycle(
            matrix, reference.prolongation, 1, threads);
        rows.push_back(measurement_row(
            level.fine, level.coarse, "global-reference", "tol=1e-10",
            reference.prolongation, reference_cycle, rhs,
            experiment_support::maximum_two_grid_cycles));

    }

    experiment_support::Report report(
        "Fixed-physical coefficient under nested grid refinement");
    report.add_summary({
        {"Version", std::string(tgi::version)},
        {"Mode", quick ? "quick" : "full"},
        {"Levels", std::to_string(level_count)},
        {"Fine/coarse pairs", quick ? "32/4, 64/8" :
             "32/4, 64/8, 128/16"},
        {"Coarsening ratio", "8 on every level"},
        {"Contrast", "1e4"},
        {"Topology", field_case.name},
        {"Physical channel width", "1/16"},
        {"Physical background partition", "8 x 8"},
        {"Shared-node mismatches",
         std::to_string(shared_node_mismatches)},
        {"Threads", std::to_string(threads)},
        {"Solve tolerance", "1e-6"}});
    report.add_note(
        "The channel width, channel centerlines, random-block partition and "
        "seed are defined in physical coordinates and held fixed while both "
        "fine and coarse meshes are dyadically refined. Equality on all "
        "shared nested-grid nodes is checked before solving. This three-level "
        "sequence complements the cell-count-defined coefficient families in "
        "experiments 1--5 and quantifies finite-grid refinement behavior.");
    report.add_table(
        "Nested-refinement comparison",
        {"1/h", "1/H", "H/h", "Method", "Parameter", "Energy",
         "P density %", "Cycles", "Status", "Final relres",
         "Eff factor", "Tail factor"},
        {5, 5, 5, 18, 22, 13, 11, 8, 10, 12, 11, 11}, rows, true);
    report.save("experiment6_fixed_physical_refinement");
    return shared_node_mismatches == 0 ? 0 : 1;
}
