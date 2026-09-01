#include "experiment/problem.hpp"
#include "experiment/reporting.hpp"
#include "multigrid/global_pcg.hpp"
#include "multigrid/multilevel_solver.hpp"
#include "version.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace {

enum class TransferMethod { Adaptive, GlobalReference, Geometric };

struct MultilevelCase {
    std::string name;
    std::array<int, 3> intervals;
    double contrast;
    tgi::CoefficientDistribution distribution;
};

struct Hierarchy {
    std::vector<tgi::SparseMatrix> matrices;
    std::vector<tgi::SparseMatrix> prolongations;
    std::vector<std::string> level_parameters;
};

const char* method_name(TransferMethod method) {
    switch (method) {
        case TransferMethod::Adaptive:
            return "adaptive";
        case TransferMethod::GlobalReference:
            return "global-reference";
        case TransferMethod::Geometric:
            return "geometric";
    }
    return "unknown";
}

Hierarchy build_hierarchy(
    const MultilevelCase& item, const tgi::SparseMatrix& fine_matrix,
    TransferMethod method, int threads) {
    Hierarchy hierarchy;
    hierarchy.matrices.push_back(fine_matrix);
    for (std::size_t level = 0; level < 2U; ++level) {
        const int fine_intervals = item.intervals[level];
        const int coarse_intervals = item.intervals[level + 1U];
        const tgi::StructuredGrid grid(
            fine_intervals - 1, fine_intervals / coarse_intervals);
        const auto geometric = tgi::build_geometric_interpolation(grid);
        tgi::SparseMatrix prolongation;
        std::string parameter;
        if (method == TransferMethod::Adaptive) {
            const int steps = tgi::adaptive_global_pcg_detail::select_steps(
                grid, hierarchy.matrices.back());
            tgi::GlobalEnergyPcgPath path(
                grid, hierarchy.matrices.back(), geometric.prolongation,
                threads);
            path.advance_to(steps);
            prolongation = path.prolongation();
            parameter = "m=" + std::to_string(steps);
        } else if (method == TransferMethod::GlobalReference) {
            const auto reference = experiment_support::build_global_reference(
                grid, hierarchy.matrices.back(), threads);
            prolongation = reference.prolongation;
            parameter = "tol=1e-10";
        } else {
            prolongation = geometric.prolongation;
            parameter = "P_G";
        }
        hierarchy.level_parameters.push_back(parameter);
        hierarchy.prolongations.push_back(std::move(prolongation));
        hierarchy.matrices.push_back(tgi::galerkin_coarse_operator(
            hierarchy.matrices.back(), hierarchy.prolongations.back(),
            threads));
    }
    return hierarchy;
}

std::string level_pair(const std::vector<std::string>& values) {
    return values[0] + "/" + values[1];
}

std::string density_pair(
    const tgi::SparseMatrix& fine_prolongation,
    const tgi::SparseMatrix& coarse_prolongation) {
    return experiment_support::fixed(
               experiment_support::interpolation_density_percent(
                   fine_prolongation), 3) +
        "/" + experiment_support::fixed(
            experiment_support::interpolation_density_percent(
                coarse_prolongation), 3);
}

experiment_support::Row measurement_row(
    const MultilevelCase& item, TransferMethod method,
    const tgi::Vector& rhs, Hierarchy hierarchy, int threads) {
    const tgi::TwoGridCycle exact_two_grid(
        hierarchy.matrices[0], hierarchy.prolongations[0], 1, threads);
    const auto two_grid = tgi::solve_two_grid(
        rhs, exact_two_grid, 1.0e-6,
        method == TransferMethod::Geometric
            ? experiment_support::maximum_geometric_cycles
            : experiment_support::maximum_two_grid_cycles);

    const tgi::MultilevelVCycle multilevel(
        std::move(hierarchy.matrices),
        std::move(hierarchy.prolongations), 1, threads);
    const auto v_cycle = tgi::solve_multilevel(
        rhs, multilevel, 1.0e-6,
        method == TransferMethod::Geometric
            ? experiment_support::maximum_geometric_cycles
            : experiment_support::maximum_two_grid_cycles);
    const double cycle_ratio = two_grid.cycles > 0
        ? static_cast<double>(v_cycle.cycles) /
              static_cast<double>(two_grid.cycles)
        : 0.0;
    return {
        item.name,
        std::to_string(item.intervals[0]) + "/" +
            std::to_string(item.intervals[1]) + "/" +
            std::to_string(item.intervals[2]),
        method_name(method),
        level_pair(hierarchy.level_parameters),
        density_pair(
            multilevel.prolongation(0), multilevel.prolongation(1)),
        experiment_support::fixed(multilevel.operator_complexity(), 4),
        experiment_support::fixed(
            multilevel.interpolation_complexity(), 4),
        std::to_string(two_grid.cycles),
        std::to_string(v_cycle.cycles),
        experiment_support::fixed(cycle_ratio, 3),
        tgi::stationary_status_name(v_cycle.status),
        experiment_support::scientific(v_cycle.relative_residual, 2),
        experiment_support::fixed(v_cycle.effective_factor, 6)};
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

    const auto& topologies = experiment_support::channel_topologies();
    const std::array<MultilevelCase, 2> cases{{
        {"cross-medium", {64, 16, 8}, 1.0e4,
         tgi::CoefficientDistribution::ChannelizedBinary},
        {"ring-center", {128, 16, 8}, 1.0e4,
         tgi::CoefficientDistribution::WindingRingBinary}}};
    const std::size_t case_count = quick ? 1U : cases.size();
    experiment_support::Rows rows;
    for (std::size_t case_index = 0; case_index < case_count; ++case_index) {
        const MultilevelCase& item = cases[case_index];
        experiment_support::BasicConfig config;
        config.fine_intervals = item.intervals[0];
        config.coarse_intervals = item.intervals[1];
        config.contrast = item.contrast;
        config.threads = threads;
        const tgi::StructuredGrid fine_grid =
            experiment_support::make_grid(config);
        experiment_support::FieldCase field;
        field.name = item.distribution ==
                tgi::CoefficientDistribution::ChannelizedBinary
            ? topologies[0].name : topologies[5].name;
        field.distribution = item.distribution;
        const auto problem = experiment_support::make_problem(
            fine_grid, field, config);
        for (const TransferMethod method : {
                 TransferMethod::Adaptive,
                 TransferMethod::GlobalReference,
                 TransferMethod::Geometric}) {
            experiment_support::progress(
                "multilevel " + std::to_string(case_index + 1U) + "/" +
                std::to_string(case_count) + ": " + item.name + ", " +
                method_name(method));
            rows.push_back(measurement_row(
                item, method, problem.rhs,
                build_hierarchy(item, problem.matrix, method, threads),
                threads));
        }
    }

    experiment_support::Report report(
        "Three-level V-cycle feasibility study");
    report.add_summary({
        {"Version", std::string(tgi::version)},
        {"Mode", quick ? "quick" : "full"},
        {"Cases", std::to_string(case_count)},
        {"Hierarchy", quick ? "64/16/8" : "64/16/8, 128/16/8"},
        {"Contrast", "1e4"},
        {"Smoother", "1 forward + 1 backward Gauss-Seidel"},
        {"Coarsest solve", "exact sparse Cholesky"},
        {"Setup threads", std::to_string(threads)},
        {"Solve tolerance", "1e-6"}});
    report.add_note(
        "The experiment applies the same interpolation construction "
        "independently on both Galerkin transitions. Exact two-grid cycles "
        "with the first interpolation provide a controlled reference for "
        "replacing the exact coarse solve by one recursive V-cycle. The two "
        "hierarchies quantify finite-level feasibility and recursion cost.");
    report.add_table(
        "Exact two-grid versus recursive V-cycle",
        {"Case", "Levels", "Method", "Level params", "P density %",
         "C_A", "C_P", "TG cyc", "V cyc", "V/TG", "Status",
         "Final relres", "Eff factor"},
        {13, 10, 18, 21, 15, 8, 8, 7, 7, 7, 10, 12, 11}, rows, true);
    report.save("experiment7_multilevel_pilot");
    return 0;
}
