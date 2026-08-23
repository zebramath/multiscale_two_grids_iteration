#pragma once

#include "multigrid/energy_interpolation.hpp"
#include "multigrid/two_grid_solver.hpp"
#include "pde/diffusion_problem.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

namespace tgi {

struct AlgebraicInterpolationReport {
    double build_ms = 0.0;
    double final_f_residual = 0.0;
    double mean_construction_iterations = 0.0;
};

struct AlgebraicInterpolationResult {
    SparseMatrix prolongation;
    AlgebraicInterpolationReport report;
};

struct JacobiInterpolationOptions {
    int steps = 4;
    double damping = 2.0 / 3.0;
    double relative_drop_tolerance = 1.0e-10;
    int maximum_entries_per_row = 8;
    int thread_count = 1;
};

struct StrengthDistanceOptions {
    int coarse_candidates_per_row = 4;
    double minimum_strength = 1.0e-12;
    double local_tolerance = 1.0e-6;
    int local_max_iterations = 40000;
    int thread_count = 1;
};

namespace algebraic_interpolation_detail {

using Clock = std::chrono::steady_clock;

inline double milliseconds(Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

inline int coarse_id_from_fine(
    const StructuredGrid& grid, int fine) {
    const auto [ix, iy] = grid.fine_coords(fine);
    return grid.coarse_id(
        (ix + 1) / grid.ratio() - 1,
        (iy + 1) / grid.ratio() - 1);
}

inline void sparsify_row(
    std::vector<std::pair<int, double>>& entries,
    double relative_drop, int maximum_entries) {
    if (entries.empty()) return;
    double maximum = 0.0;
    for (const auto& entry : entries) {
        maximum = std::max(maximum, std::abs(entry.second));
    }
    const double threshold = relative_drop * maximum;
    entries.erase(
        std::remove_if(
            entries.begin(), entries.end(),
            [&](const auto& entry) {
                return std::abs(entry.second) <= threshold;
            }),
        entries.end());
    if (maximum_entries > 0 &&
        entries.size() > static_cast<std::size_t>(maximum_entries)) {
        std::nth_element(
            entries.begin(), entries.begin() + maximum_entries,
            entries.end(), [](const auto& lhs, const auto& rhs) {
                return std::abs(lhs.second) > std::abs(rhs.second);
            });
        entries.resize(static_cast<std::size_t>(maximum_entries));
    }
    std::sort(entries.begin(), entries.end(), [](const auto& lhs,
                                                  const auto& rhs) {
        return lhs.first < rhs.first;
    });
}

inline double scaled_f_residual(
    const StructuredGrid& grid, const SparseMatrix& a,
    const SparseMatrix& p, int thread_count) {
    const SparseMatrix ap = sparse_multiply(a, p, 0.0, thread_count);
    const Vector diagonal = a.diagonal();
    double squared = 0.0;
    for (int row = 0; row < grid.fine_size(); ++row) {
        if (grid.is_coarse_node(row)) continue;
        for (int position = ap.row_ptr()[static_cast<std::size_t>(row)];
             position < ap.row_ptr()[static_cast<std::size_t>(row) + 1U];
             ++position) {
            const double value =
                ap.values()[static_cast<std::size_t>(position)];
            squared += value * value /
                diagonal[static_cast<std::size_t>(row)];
        }
    }
    return std::sqrt(squared);
}

inline SparseMatrix jacobi_step(
    const StructuredGrid& grid, const SparseMatrix& a,
    const SparseMatrix& p, const Vector& diagonal,
    const JacobiInterpolationOptions& options) {
    const SparseMatrix ap = sparse_multiply(
        a, p, 0.0, options.thread_count);
    std::vector<Triplet> triplets;
    triplets.reserve(
        p.nnz() + static_cast<std::size_t>(2 * grid.fine_size()));
    std::vector<std::pair<int, double>> row;
    for (int fine = 0; fine < grid.fine_size(); ++fine) {
        if (grid.is_coarse_node(fine)) {
            triplets.push_back(
                {fine, coarse_id_from_fine(grid, fine), 1.0});
            continue;
        }
        row.clear();
        int p_position = p.row_ptr()[static_cast<std::size_t>(fine)];
        int ap_position = ap.row_ptr()[static_cast<std::size_t>(fine)];
        const int p_end = p.row_ptr()[static_cast<std::size_t>(fine) + 1U];
        const int ap_end = ap.row_ptr()[static_cast<std::size_t>(fine) + 1U];
        while (p_position < p_end || ap_position < ap_end) {
            const int p_col = p_position < p_end
                ? p.col_idx()[static_cast<std::size_t>(p_position)]
                : std::numeric_limits<int>::max();
            const int ap_col = ap_position < ap_end
                ? ap.col_idx()[static_cast<std::size_t>(ap_position)]
                : std::numeric_limits<int>::max();
            const int column = std::min(p_col, ap_col);
            double value = 0.0;
            if (p_col == column) {
                value += p.values()[static_cast<std::size_t>(p_position++)];
            }
            if (ap_col == column) {
                value -= options.damping *
                    ap.values()[static_cast<std::size_t>(ap_position++)] /
                    diagonal[static_cast<std::size_t>(fine)];
            }
            row.emplace_back(column, value);
        }
        sparsify_row(
            row, options.relative_drop_tolerance,
            options.maximum_entries_per_row);
        for (const auto& [column, value] : row) {
            triplets.push_back({fine, column, value});
        }
    }
    return SparseMatrix(
        grid.fine_size(), grid.coarse_size(), triplets, 0.0);
}

inline void insert_candidate(
    std::vector<std::pair<double, int>>& candidates,
    std::pair<double, int> candidate, int limit) {
    if (candidates.size() < static_cast<std::size_t>(limit)) {
        candidates.push_back(candidate);
        std::sort(candidates.begin(), candidates.end());
        return;
    }
    if (candidate < candidates.back()) {
        candidates.back() = candidate;
        std::sort(candidates.begin(), candidates.end());
    }
}

} // namespace algebraic_interpolation_detail

// F-point weighted Jacobi applied to AP=0 while C rows remain exact injection.
// Without dropping, the limit is the ideal interpolation
// P_F=-A_FF^{-1}A_FC for the prescribed coarse variables.
inline AlgebraicInterpolationResult build_jacobi_interpolation(
    const StructuredGrid& grid, const SparseMatrix& a,
    const SparseMatrix& initial,
    const JacobiInterpolationOptions& options = {}) {
    if (a.rows() != grid.fine_size() || a.cols() != grid.fine_size() ||
        initial.rows() != grid.fine_size() ||
        initial.cols() != grid.coarse_size() || options.steps < 0 ||
        !(options.damping > 0.0 && options.damping < 2.0) ||
        options.relative_drop_tolerance < 0.0 ||
        options.maximum_entries_per_row < 0) {
        throw std::invalid_argument(
            "build_jacobi_interpolation: invalid input or options");
    }
    const auto begin = algebraic_interpolation_detail::Clock::now();
    const Vector diagonal = a.diagonal();
    SparseMatrix current = initial;
    for (int step = 0; step < options.steps; ++step) {
        current = algebraic_interpolation_detail::jacobi_step(
            grid, a, current, diagonal, options);
    }
    AlgebraicInterpolationReport report;
    report.final_f_residual =
        algebraic_interpolation_detail::scaled_f_residual(
            grid, a, current, options.thread_count);
    report.mean_construction_iterations =
        static_cast<double>(options.steps);
    report.build_ms = algebraic_interpolation_detail::milliseconds(
        begin, algebraic_interpolation_detail::Clock::now());
    return {std::move(current), report};
}

// Coarse variables are prescribed. Each F row selects the q nearest coarse
// variables in the graph metric sqrt(A_ii A_jj)/|A_ij|. These row selections
// define column supports; weights are then obtained by constrained local energy
// minimization, rather than by an ad-hoc inverse-distance formula.
inline AlgebraicInterpolationResult build_strength_distance_interpolation(
    const StructuredGrid& grid, const SparseMatrix& a,
    const StrengthDistanceOptions& options = {}) {
    if (a.rows() != grid.fine_size() || a.cols() != grid.fine_size() ||
        options.coarse_candidates_per_row <= 0 ||
        !(options.minimum_strength > 0.0) ||
        !(options.local_tolerance > 0.0) ||
        options.local_max_iterations <= 0 || options.thread_count <= 0) {
        throw std::invalid_argument(
            "build_strength_distance_interpolation: invalid input or options");
    }
    const auto begin = algebraic_interpolation_detail::Clock::now();
    const int candidate_count = std::min(
        options.coarse_candidates_per_row, grid.coarse_size());
    const Vector diagonal = a.diagonal();
    std::vector<std::vector<std::pair<double, int>>> nearest(
        static_cast<std::size_t>(grid.fine_size()));
    using QueueEntry = std::pair<double, int>;

    for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
        std::vector<double> distance(
            static_cast<std::size_t>(grid.fine_size()),
            std::numeric_limits<double>::infinity());
        std::priority_queue<
            QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>>
            queue;
        const int source = grid.coarse_fine_id(coarse);
        distance[static_cast<std::size_t>(source)] = 0.0;
        queue.push({0.0, source});
        while (!queue.empty()) {
            const auto [current_distance, node] = queue.top();
            queue.pop();
            if (current_distance !=
                distance[static_cast<std::size_t>(node)]) {
                continue;
            }
            for (int position =
                     a.row_ptr()[static_cast<std::size_t>(node)];
                 position <
                     a.row_ptr()[static_cast<std::size_t>(node) + 1U];
                 ++position) {
                const int neighbor =
                    a.col_idx()[static_cast<std::size_t>(position)];
                if (neighbor == node) continue;
                const double edge = std::abs(
                    a.values()[static_cast<std::size_t>(position)]);
                const double strength = edge / std::sqrt(
                    diagonal[static_cast<std::size_t>(node)] *
                    diagonal[static_cast<std::size_t>(neighbor)]);
                const double cost = 1.0 /
                    std::max(strength, options.minimum_strength);
                const double proposed = current_distance + cost;
                if (proposed < distance[static_cast<std::size_t>(neighbor)]) {
                    distance[static_cast<std::size_t>(neighbor)] = proposed;
                    queue.push({proposed, neighbor});
                }
            }
        }
        for (int fine = 0; fine < grid.fine_size(); ++fine) {
            algebraic_interpolation_detail::insert_candidate(
                nearest[static_cast<std::size_t>(fine)],
                {distance[static_cast<std::size_t>(fine)], coarse},
                candidate_count);
        }
    }

    std::vector<std::vector<int>> supports(
        static_cast<std::size_t>(grid.coarse_size()));
    for (int fine = 0; fine < grid.fine_size(); ++fine) {
        if (grid.is_coarse_node(fine)) continue;
        const auto& candidates = nearest[static_cast<std::size_t>(fine)];
        for (const auto& [distance, coarse] : candidates) {
            (void)distance;
            supports[static_cast<std::size_t>(coarse)].push_back(fine);
        }
    }
    for (auto& support : supports) {
        std::sort(support.begin(), support.end());
        support.erase(
            std::unique(support.begin(), support.end()), support.end());
    }
    InterpolationOptions interpolation_options;
    interpolation_options.strategy =
        InterpolationStrategy::LocalEnergyMinimum;
    interpolation_options.local_tolerance = options.local_tolerance;
    interpolation_options.local_max_iterations =
        options.local_max_iterations;
    interpolation_options.thread_count = options.thread_count;
    const InterpolationResult interpolation =
        build_energy_interpolation_on_supports(
            grid, a, supports, interpolation_options);
    SparseMatrix prolongation = interpolation.prolongation;
    AlgebraicInterpolationReport report;
    report.final_f_residual =
        algebraic_interpolation_detail::scaled_f_residual(
            grid, a, prolongation, 1);
    if (interpolation.report.local_solves.systems > 0) {
        report.mean_construction_iterations =
            static_cast<double>(
                interpolation.report.local_solves.total_iterations) /
            static_cast<double>(
                interpolation.report.local_solves.systems);
    }
    report.build_ms = algebraic_interpolation_detail::milliseconds(
        begin, algebraic_interpolation_detail::Clock::now());
    return {std::move(prolongation), report};
}

} // namespace tgi
