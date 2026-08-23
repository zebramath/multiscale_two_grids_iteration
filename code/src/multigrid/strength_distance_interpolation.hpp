#pragma once

#include "multigrid/energy_interpolation.hpp"
#include "pde/diffusion_problem.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

namespace tgi {

struct StrengthDistanceOptions {
    int coarse_candidates_per_row = 4;
    double minimum_strength = 1.0e-12;
    double local_tolerance = 1.0e-6;
    int local_max_iterations = 40000;
    int thread_count = 1;
};

struct StrengthDistanceResult {
    SparseMatrix prolongation;
    double build_ms = 0.0;
};

namespace strength_distance_detail {

using Clock = std::chrono::steady_clock;

inline void insert_candidate(
    std::vector<std::pair<double, int>>& candidates,
    std::pair<double, int> candidate, int limit) {
    if (candidates.size() < static_cast<std::size_t>(limit)) {
        candidates.push_back(candidate);
        std::sort(candidates.begin(), candidates.end());
    } else if (candidate < candidates.back()) {
        candidates.back() = candidate;
        std::sort(candidates.begin(), candidates.end());
    }
}

}

inline StrengthDistanceResult build_strength_distance_interpolation(
    const StructuredGrid& grid, const SparseMatrix& a,
    const StrengthDistanceOptions& options = {}) {
    const auto begin = strength_distance_detail::Clock::now();
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
            if (current_distance != distance[static_cast<std::size_t>(node)])
                continue;
            for (int position = a.row_ptr()[static_cast<std::size_t>(node)];
                 position < a.row_ptr()[static_cast<std::size_t>(node) + 1U];
                 ++position) {
                const int neighbor =
                    a.col_idx()[static_cast<std::size_t>(position)];
                if (neighbor == node) continue;
                const double edge = std::abs(
                    a.values()[static_cast<std::size_t>(position)]);
                const double strength = edge / std::sqrt(
                    diagonal[static_cast<std::size_t>(node)] *
                    diagonal[static_cast<std::size_t>(neighbor)]);
                const double proposed = current_distance +
                    1.0 / std::max(strength, options.minimum_strength);
                if (proposed < distance[static_cast<std::size_t>(neighbor)]) {
                    distance[static_cast<std::size_t>(neighbor)] = proposed;
                    queue.push({proposed, neighbor});
                }
            }
        }
        for (int fine = 0; fine < grid.fine_size(); ++fine) {
            strength_distance_detail::insert_candidate(
                nearest[static_cast<std::size_t>(fine)],
                {distance[static_cast<std::size_t>(fine)], coarse},
                candidate_count);
        }
    }

    std::vector<std::vector<int>> supports(
        static_cast<std::size_t>(grid.coarse_size()));
    for (int fine = 0; fine < grid.fine_size(); ++fine) {
        if (grid.is_coarse_node(fine)) continue;
        for (const auto& [distance, coarse] :
             nearest[static_cast<std::size_t>(fine)]) {
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
    InterpolationResult interpolation = build_energy_interpolation_on_supports(
        grid, a, supports, interpolation_options);
    const double build_ms = std::chrono::duration<double, std::milli>(
        strength_distance_detail::Clock::now() - begin).count();
    return {std::move(interpolation.prolongation), build_ms};
}

}
