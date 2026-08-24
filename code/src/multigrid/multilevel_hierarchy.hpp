#pragma once

#include "multigrid/global_pcg.hpp"
#include "multigrid/multilevel_solver.hpp"

#include <chrono>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace tgi {

enum class MultilevelInterpolationPolicy {
    Geometric,
    AdaptiveGlobalPcg,
    ExactGlobalEnergy
};

struct MultilevelHierarchyOptions {
    int first_coarse_intervals = 16;
    int smoothing_steps = 1;
    int thread_count = 1;
    int maximum_pilot_cycles = 6000;
    AdaptiveGlobalPcgOptions first_level_adaptive;
};

struct MultilevelHierarchyReport {
    double total_setup_ms = 0.0;
    std::vector<int> intervals;
    std::vector<int> selected_steps;
};

struct MultilevelHierarchyResult {
    std::unique_ptr<MultilevelCycle> cycle;
    MultilevelHierarchyReport hierarchy;
};

inline InterpolationResult geometric_level_interpolation(
    const StructuredGrid& grid, const SparseMatrix& matrix) {
    InterpolationOptions options;
    options.strategy = InterpolationStrategy::GeometricBilinear;
    return build_interpolation(grid, matrix, options);
}

inline InterpolationResult exact_level_interpolation(
    const StructuredGrid& grid, const SparseMatrix& matrix,
    int thread_count) {
    InterpolationOptions options;
    options.strategy = InterpolationStrategy::GlobalEnergyMinimum;
    options.local_tolerance = 1.0e-10;
    options.local_max_iterations = 40000;
    options.thread_count = thread_count;
    options.drop_tolerance = 0.0;
    options.require_convergence = true;
    return build_interpolation(grid, matrix, options);
}

inline MultilevelHierarchyResult build_multilevel_hierarchy(
    int fine_intervals, const SparseMatrix& fine_matrix,
    const Vector& representative_rhs,
    MultilevelInterpolationPolicy policy,
    MultilevelHierarchyOptions options = {}) {
    const auto begin = std::chrono::steady_clock::now();
    if (fine_intervals < 4 || options.first_coarse_intervals < 2 ||
        fine_intervals % options.first_coarse_intervals != 0 ||
        representative_rhs.size() !=
            static_cast<std::size_t>(fine_matrix.rows())) {
        throw std::invalid_argument(
            "invalid multilevel hierarchy grid or representative RHS");
    }
    options.thread_count = std::max(1, options.thread_count);

    std::vector<SparseMatrix> matrices;
    std::vector<SparseMatrix> prolongations;
    matrices.push_back(fine_matrix);
    Vector level_rhs = representative_rhs;
    int intervals = fine_intervals;
    bool first_level = true;
    MultilevelHierarchyReport hierarchy_report;
    hierarchy_report.intervals.push_back(intervals);

    while (intervals > 2) {
        const int ratio = first_level
            ? intervals / options.first_coarse_intervals
            : 2;
        if (ratio < 2 || intervals % ratio != 0) {
            throw std::invalid_argument(
                "multilevel hierarchy requires nested integer coarsening");
        }
        const StructuredGrid grid(intervals - 1, ratio);
        const auto geometric = geometric_level_interpolation(
            grid, matrices.back());
        SparseMatrix prolongation;
        int selected_steps = 0;

        if (policy == MultilevelInterpolationPolicy::Geometric) {
            prolongation = geometric.prolongation;
        } else if (policy ==
                   MultilevelInterpolationPolicy::ExactGlobalEnergy) {
            prolongation = exact_level_interpolation(
                grid, matrices.back(), options.thread_count).prolongation;
            selected_steps = -1;
        } else {
            AdaptiveGlobalPcgOptions adaptive =
                options.first_level_adaptive;
            adaptive.thread_count = options.thread_count;
            adaptive.maximum_cycles = options.maximum_pilot_cycles;
            if (!first_level) {
                adaptive.minimum_steps = 4;
                adaptive.maximum_steps = 24;
                adaptive.step_quantum = 2;
                adaptive.pilot_iterations = 8;
                adaptive.tail_window = 4;
            }
            const auto selected =
                build_adaptive_global_pcg_interpolation(
                    grid, matrices.back(), geometric.prolongation,
                    adaptive, &level_rhs);
            prolongation = *selected.prolongation;
            selected_steps = selected.report.selected_steps;
        }

        SparseMatrix next_matrix = galerkin_sparse(
            matrices.back(), prolongation, 0.0, options.thread_count);
        Vector next_rhs;
        prolongation.transpose_multiply(level_rhs, next_rhs);
        prolongations.push_back(std::move(prolongation));
        matrices.push_back(std::move(next_matrix));
        level_rhs = std::move(next_rhs);
        hierarchy_report.selected_steps.push_back(selected_steps);
        intervals /= ratio;
        hierarchy_report.intervals.push_back(intervals);
        first_level = false;
    }

    auto cycle = std::make_unique<MultilevelCycle>(
        std::move(matrices), std::move(prolongations),
        options.smoothing_steps, options.thread_count);
    hierarchy_report.total_setup_ms =
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - begin).count();
    return {std::move(cycle), std::move(hierarchy_report)};
}

}
