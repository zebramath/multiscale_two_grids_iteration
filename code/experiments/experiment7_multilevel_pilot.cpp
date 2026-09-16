#include "experiment/problem.hpp"
#include "experiment/reporting.hpp"
#include "multigrid/global_pcg.hpp"
#include "multigrid/multilevel_solver.hpp"
#include <array>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>
namespace {
enum class TransferMethod { FiniteOneThird, EnergyEndpoint };
struct MultilevelCase {
    std::string name;
    std::array<int, 3> intervals;
    tgi::CoefficientDistribution distribution;
};
struct Hierarchy {
    std::vector<tgi::SparseMatrix> matrices;
    std::vector<tgi::SparseMatrix> prolongations;
    std::vector<std::string> parameters;
};
const char* method_name(TransferMethod method) {
    return method == TransferMethod::FiniteOneThird
        ? "finite-1/3" : "energy-endpoint";
}
Hierarchy build_hierarchy(
    const MultilevelCase& item, const tgi::SparseMatrix& fine_matrix,
    TransferMethod method, int threads) {
    Hierarchy hierarchy;
    hierarchy.matrices.push_back(fine_matrix);
    for (std::size_t level = 0; level < 2U; ++level) {
        const int fine = item.intervals[level];
        const int coarse = item.intervals[level + 1U];
        const tgi::StructuredGrid grid(fine - 1, fine / coarse);
        const auto initial = tgi::build_geometric_interpolation(grid);
        tgi::GlobalEnergyPcgPath path(
            grid, hierarchy.matrices.back(), initial, threads);
        if (method == TransferMethod::FiniteOneThird) {
            const int steps = experiment_support::fixed_path_steps(
                fine, 1, 3);
            path.advance_to(steps);
            hierarchy.parameters.push_back("m=" + std::to_string(steps));
        } else {
            path.advance_until_relative_residual(1.0e-10);
            hierarchy.parameters.push_back("relres<=1e-10");
        }
        hierarchy.prolongations.push_back(path.prolongation());
        hierarchy.matrices.push_back(tgi::galerkin_coarse_operator(
            hierarchy.matrices.back(), hierarchy.prolongations.back(),
            threads));
    }
    return hierarchy;
}
std::string pair_text(const std::vector<std::string>& values) {
    return values.at(0) + "/" + values.at(1);
}
std::string density_pair(
    const tgi::SparseMatrix& first, const tgi::SparseMatrix& second) {
    return experiment_support::fixed(
               experiment_support::interpolation_density_percent(first), 3) +
        "/" + experiment_support::fixed(
            experiment_support::interpolation_density_percent(second), 3);
}
experiment_support::Row measure(
    const MultilevelCase& item, TransferMethod method,
    const tgi::Vector& rhs, Hierarchy hierarchy, int threads) {
    const tgi::TwoGridCycle exact_two_grid(
        hierarchy.matrices[0], hierarchy.prolongations[0], 1, threads);
    const auto two_grid = tgi::solve_two_grid(
        rhs, exact_two_grid, 1.0e-6,
        experiment_support::maximum_two_grid_cycles);
    const std::string parameters = pair_text(hierarchy.parameters);
    const std::string densities = density_pair(
        hierarchy.prolongations[0], hierarchy.prolongations[1]);
    const tgi::MultilevelVCycle multilevel(
        std::move(hierarchy.matrices),
        std::move(hierarchy.prolongations), 1, threads);
    const auto v_cycle = tgi::solve_multilevel(
        rhs, multilevel, 1.0e-6,
        experiment_support::maximum_two_grid_cycles);
    const double ratio = two_grid.cycles > 0
        ? static_cast<double>(v_cycle.cycles) /
              static_cast<double>(two_grid.cycles)
        : 0.0;
    return {
        item.name,
        std::to_string(item.intervals[0]) + "/" +
            std::to_string(item.intervals[1]) + "/" +
            std::to_string(item.intervals[2]),
        method_name(method), parameters, densities,
        experiment_support::fixed(multilevel.operator_complexity(), 4),
        experiment_support::fixed(
            multilevel.interpolation_complexity(), 4),
        std::to_string(two_grid.cycles),
        std::to_string(v_cycle.cycles),
        experiment_support::fixed(ratio, 3),
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
    const std::array<MultilevelCase, 2> cases{{
        {"cross-medium", {64, 16, 8},
         tgi::CoefficientDistribution::ChannelizedBinary},
        {"ring-center", {128, 16, 8},
         tgi::CoefficientDistribution::WindingRingBinary}}};
    const std::size_t count = quick ? 1U : cases.size();
    experiment_support::Rows rows;
    for (std::size_t index = 0; index < count; ++index) {
        const MultilevelCase& item = cases[index];
        experiment_support::BasicConfig config;
        config.fine_intervals = item.intervals[0];
        config.coarse_intervals = item.intervals[1];
        config.contrast = 1.0e4;
        config.threads = threads;
        const tgi::StructuredGrid grid =
            experiment_support::make_grid(config);
        experiment_support::FieldCase field;
        field.name = item.distribution ==
                tgi::CoefficientDistribution::ChannelizedBinary
            ? "cross-channel" : "winding-ring";
        field.distribution = item.distribution;
        const auto problem = experiment_support::make_problem(
            grid, field, config);
        for (const TransferMethod method : {
                 TransferMethod::FiniteOneThird,
                 TransferMethod::EnergyEndpoint}) {
            experiment_support::progress(
                "multilevel " + std::to_string(index + 1U) + "/" +
                std::to_string(count) + ": " + item.name + ", " +
                method_name(method));
            rows.push_back(measure(
                item, method, problem.rhs,
                build_hierarchy(item, problem.matrix, method, threads),
                threads));
        }
    }
    experiment_support::Report report(
        "Three-level V-cycle feasibility study");
    report.add_summary({
        {"Mode", quick ? "quick" : "full"},
        {"Cases", std::to_string(count)},
        {"Hierarchy", quick ? "64/16/8" : "64/16/8, 128/16/8"},
        {"Contrast", "1e4"},
        {"Smoother", "1 forward + 1 backward Gauss-Seidel"},
        {"Coarsest solve", "exact sparse Cholesky"},
        {"Threads", std::to_string(threads)},
        {"Solve tolerance", "1e-6"}});
    report.add_note(
        "The same interpolation construction is applied independently on "
        "both Galerkin transitions. Exact two-grid cycles using the first "
        "transfer provide the controlled reference for replacing the exact "
        "coarse solve by one recursive V-cycle. This experiment is a "
        "feasibility check; the proved results remain two-grid statements.");
    report.add_table(
        "Exact two-grid versus recursive V-cycle",
        {"Case", "Levels", "Method", "Level parameters", "P density %",
         "C_A", "C_P", "TG cycles", "V cycles", "V/TG", "Status",
         "Final relres", "Eff factor"},
        {13, 10, 18, 29, 15, 8, 8, 9, 9, 7, 10, 12, 11}, rows, true);
    report.save("experiment7_multilevel_pilot");
    return 0;
}
