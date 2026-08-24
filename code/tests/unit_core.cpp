#include "experiment/study.hpp"
#include "multigrid/support_expansion.hpp"
#include "multigrid/energy_interpolation.hpp"
#include "multigrid/multilevel_hierarchy.hpp"
#include "multigrid/two_grid_solver.hpp"
#include "pde/diffusion_problem.hpp"
#include "version.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

template <class Action>
void require_throws(Action&& action, const std::string& message) {
    try {
        action();
    } catch (const std::exception&) {
        return;
    }
    throw std::runtime_error(message);
}

double row_sum(const tgi::SparseMatrix& matrix, int row) {
    double sum = 0.0;
    for (int position = matrix.row_ptr()[static_cast<std::size_t>(row)];
         position < matrix.row_ptr()[static_cast<std::size_t>(row) + 1U];
         ++position) {
        sum += matrix.values()[static_cast<std::size_t>(position)];
    }
    return sum;
}

}

int main() {
    require(tgi::version == "4.0.0", "wrong package version");
    require(tgi::version_major == 4 && tgi::version_minor == 0 &&
                tgi::version_patch == 0,
            "inconsistent numeric package version");
    const experiment_support::Row expected_headers{
        "Field", "Method", "Parameter", "P density %",
        "Setup ms", "Total ms", "Cycles"};
    require(experiment_support::study_headers() == expected_headers,
            "public output schema is not the compact schema");
    for (const auto& header : experiment_support::study_headers()) {
        std::string lower = header;
        std::transform(
            lower.begin(), lower.end(), lower.begin(),
            [](unsigned char value) {
                return static_cast<char>(std::tolower(value));
            });
        require(lower.find("residual") == std::string::npos,
                "forbidden residual metric in public output");
    }

    const tgi::StructuredGrid grid(15, 4);
    tgi::CoefficientOptions coefficient_options;
    coefficient_options.distribution =
        tgi::CoefficientDistribution::ChannelizedBinary;
    coefficient_options.contrast = 1.0e4;
    coefficient_options.channel_background_block_size = 4;
    const auto coefficient = tgi::make_coefficient(grid, coefficient_options);
    std::vector<tgi::Vector> complex_topologies;
    for (const auto distribution : {
             tgi::CoefficientDistribution::BranchingChannelsBinary,
             tgi::CoefficientDistribution::WindingRingBinary}) {
        coefficient_options.distribution = distribution;
        const auto field = tgi::make_coefficient(grid, coefficient_options);
        require(field.actual_contrast == coefficient_options.contrast,
                "complex topology lost the requested contrast");
        require(std::find(field.values.begin(), field.values.end(), 1.0) !=
                    field.values.end(),
                "complex topology has no background region");
        require(std::find(
                    field.values.begin(), field.values.end(),
                    coefficient_options.contrast) != field.values.end(),
                "complex topology has no high-conductivity region");
        complex_topologies.push_back(field.values);
    }
    require(tgi::norm2(tgi::subtract(
                complex_topologies[0], complex_topologies[1])) > 0.0,
            "complex coefficient topologies are indistinguishable");
    const tgi::SparseMatrix a = tgi::assemble_diffusion(
        grid, coefficient.values);
    require(a.rows() == grid.fine_size(), "wrong matrix dimension");
    experiment_support::BasicConfig problem_config;
    problem_config.fine_intervals = 16;
    problem_config.coarse_intervals = 4;
    const auto problem = experiment_support::make_problem(
        grid, experiment_support::standard_fields().front(), problem_config);
    require(std::all_of(
                problem.rhs.begin(), problem.rhs.end(),
                [](double value) { return value == 1.0; }),
            "experiment right-hand side is not constant one");
    for (double diagonal : a.diagonal()) {
        require(diagonal > 0.0, "diffusion diagonal is not positive");
    }

    tgi::InterpolationOptions geometric_options;
    geometric_options.strategy =
        tgi::InterpolationStrategy::GeometricBilinear;
    const auto geometric = tgi::build_interpolation(
        grid, a, geometric_options);
    tgi::InterpolationOptions local_options;
    local_options.strategy =
        tgi::InterpolationStrategy::LocalEnergyMinimum;
    local_options.patch_layers = 1;
    local_options.local_tolerance = 1.0e-8;
    local_options.local_max_iterations = 40000;
    local_options.thread_count = 2;
    const auto local = tgi::build_interpolation(grid, a, local_options);
    require(local.prolongation.cols() == grid.coarse_size(),
            "wrong local interpolation dimension");
    tgi::InterpolationOptions unfinished_local = local_options;
    unfinished_local.local_tolerance = 0.0;
    unfinished_local.local_max_iterations = 0;
    unfinished_local.require_convergence = false;
    const auto unfinished = tgi::build_interpolation(
        grid, a, unfinished_local);
    require(unfinished.report.local_solves.failed_systems > 0,
            "nonconverged local solves were not reported");
    unfinished_local.require_convergence = true;
    require_throws(
        [&]() { (void)tgi::build_interpolation(grid, a, unfinished_local); },
        "required local convergence was not enforced");

    tgi::StrengthDistanceOptions distance_options;
    distance_options.coarse_candidates_per_row = 3;
    distance_options.thread_count = 2;
    const auto distance = tgi::build_strength_distance_interpolation(
        grid, a, distance_options);
    for (int fine = 0; fine < grid.fine_size(); ++fine) {
        if (!grid.is_coarse_node(fine)) continue;
        require(std::abs(row_sum(distance.prolongation, fine) - 1.0) <
                    1.0e-12,
                "strength-distance interpolation lost C-point injection");
    }

    tgi::InterpolationOptions global_options = local_options;
    global_options.strategy =
        tgi::InterpolationStrategy::GlobalEnergyMinimum;
    global_options.patch_layers = 0;
    global_options.local_tolerance = 1.0e-10;
    const auto global = tgi::build_interpolation(grid, a, global_options);
    const auto pruned = tgi::prune_global_interpolation_relative(
        grid, global.prolongation, 1.0e-2);
    require(pruned.prolongation.nnz() <= global.prolongation.nnz(),
            "relative pruning increased interpolation density");

    tgi::ResidualStrongSupportOptions strong_options;
    strong_options.base_patch_layers = 1;
    strong_options.maximum_extra_nodes_per_column = 8;
    strong_options.thread_count = 2;
    const auto strong = tgi::build_residual_strong_supports(
        grid, a, local.prolongation, strong_options);
    require(strong.supports.size() ==
                static_cast<std::size_t>(grid.coarse_size()),
            "strong support count is wrong");

    tgi::InterpolationOptions fixed_step_options = global_options;
    fixed_step_options.local_tolerance = 0.0;
    fixed_step_options.local_max_iterations = 2;
    fixed_step_options.require_convergence = false;
    fixed_step_options.drop_tolerance = 0.0;
    const auto pcg2 = tgi::refine_global_energy_interpolation(
        grid, a, geometric.prolongation, fixed_step_options);
    require(pcg2.report.local_solves.max_iterations == 2,
            "fixed-step PCG did not honor the iteration budget");
    require(pcg2.report.local_solves.total_iterations ==
                2 * pcg2.report.local_solves.systems,
            "fixed-step PCG stopped before its explicit budget");

    const tgi::Vector rhs(
        static_cast<std::size_t>(grid.fine_size()), 1.0);
    const tgi::TwoGridCycle cycle(a, global.prolongation, 1, 2);
    const auto& setup = cycle.setup_report();
    require(setup.interpolation_energy > 0.0,
            "interpolation energy diagnostic is not positive");
    double coarse_trace = 0.0;
    for (double value : cycle.coarse_matrix().diagonal()) {
        coarse_trace += value;
    }
    require(std::abs(coarse_trace - setup.interpolation_energy) <=
                1.0e-12 * coarse_trace,
            "interpolation energy does not match trace(P^T A P)");
    tgi::Vector coarse_rhs(
        static_cast<std::size_t>(cycle.coarse_size()), 1.0);
    tgi::Vector coarse_solution;
    tgi::Vector coarse_work;
    cycle.solve_coarse_system(
        coarse_rhs, coarse_solution, coarse_work);
    const tgi::Vector coarse_product =
        cycle.coarse_matrix().multiply(coarse_solution);
    const double coarse_relative_residual =
        tgi::norm2(tgi::subtract(coarse_rhs, coarse_product)) /
        tgi::norm2(coarse_rhs);
    require(coarse_relative_residual < 1.0e-10,
            "coarse direct solve is not numerically accurate");
    const auto solved = tgi::solve_two_grid(rhs, cycle, 1.0e-6, 1000);
    require(solved.converged, "two-grid solve did not converge");
    const auto zero_budget = tgi::solve_two_grid(
        rhs, cycle, 1.0e-6, 0);
    require(!zero_budget.converged && zero_budget.cycles == 0,
            "two-grid solve ignored its zero-cycle hard stop");

    tgi::MultilevelHierarchyOptions hierarchy_options;
    hierarchy_options.first_coarse_intervals = 4;
    hierarchy_options.thread_count = 2;
    const auto hierarchy = tgi::build_multilevel_hierarchy(
        grid.intervals(), a, rhs,
        tgi::MultilevelInterpolationPolicy::Geometric,
        hierarchy_options);
    require(hierarchy.cycle->setup_report().levels == 3,
            "multilevel hierarchy has the wrong number of levels");
    require(hierarchy.cycle->setup_report().operator_complexity > 1.0 &&
                hierarchy.cycle->setup_report()
                    .interpolation_complexity > 0.0,
            "multilevel complexity report is invalid");
    const tgi::Vector preconditioned = hierarchy.cycle->apply(rhs);
    require(tgi::dot(rhs, preconditioned) > 0.0,
            "multilevel V-cycle is not a positive preconditioner");
    const auto multilevel = tgi::solve_multilevel(
        rhs, *hierarchy.cycle, 1.0e-6, 1000);
    require(multilevel.converged,
            "standalone multilevel iteration did not converge");
    const auto multilevel_pcg = tgi::solve_preconditioned_cg(
        a, rhs, *hierarchy.cycle, 1.0e-6, 200);
    require(multilevel_pcg.converged,
            "multilevel-preconditioned CG did not converge");
    return 0;
}
