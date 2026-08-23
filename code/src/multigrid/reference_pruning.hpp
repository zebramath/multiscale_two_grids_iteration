#pragma once

#include "multigrid/energy_interpolation.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace tgi {

struct MatchedReferenceResult {
    SparseMatrix prolongation;
    std::vector<std::vector<int>> supports;
    double build_ms = 0.0;
};

inline std::vector<int> interpolation_f_entries_per_column(
    const StructuredGrid& grid, const SparseMatrix& prolongation,
    int thread_count = 1) {
    if (prolongation.rows() != grid.fine_size() ||
        prolongation.cols() != grid.coarse_size()) {
        throw std::invalid_argument(
            "interpolation_f_entries_per_column: incompatible dimensions");
    }
    const SparseMatrix transpose = prolongation.transpose(thread_count);
    std::vector<int> counts(
        static_cast<std::size_t>(grid.coarse_size()), 0);
    for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
        for (int position =
                 transpose.row_ptr()[static_cast<std::size_t>(coarse)];
             position <
                 transpose.row_ptr()[static_cast<std::size_t>(coarse) + 1U];
             ++position) {
            const int fine =
                transpose.col_idx()[static_cast<std::size_t>(position)];
            if (!grid.is_coarse_node(fine)) {
                ++counts[static_cast<std::size_t>(coarse)];
            }
        }
    }
    return counts;
}

inline std::vector<std::vector<int>> select_reference_supports(
    const StructuredGrid& grid, const SparseMatrix& a,
    const SparseMatrix& global_reference,
    const std::vector<int>& f_entries_per_column,
    int thread_count = 1) {
    if (a.rows() != grid.fine_size() || a.cols() != grid.fine_size() ||
        global_reference.rows() != grid.fine_size() ||
        global_reference.cols() != grid.coarse_size() ||
        f_entries_per_column.size() !=
            static_cast<std::size_t>(grid.coarse_size())) {
        throw std::invalid_argument(
            "select_reference_supports: incompatible dimensions");
    }
    const SparseMatrix transpose = global_reference.transpose(thread_count);
    const Vector diagonal = a.diagonal();
    std::vector<std::vector<int>> supports(
        static_cast<std::size_t>(grid.coarse_size()));
    std::vector<std::pair<double, int>> ranking;
    for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
        ranking.clear();
        for (int position =
                 transpose.row_ptr()[static_cast<std::size_t>(coarse)];
             position <
                 transpose.row_ptr()[static_cast<std::size_t>(coarse) + 1U];
             ++position) {
            const int fine =
                transpose.col_idx()[static_cast<std::size_t>(position)];
            if (grid.is_coarse_node(fine)) continue;
            const double value =
                transpose.values()[static_cast<std::size_t>(position)];
            const double score = std::abs(value) *
                std::sqrt(diagonal[static_cast<std::size_t>(fine)]);
            ranking.push_back({score, fine});
        }
        const int requested = std::max(
            0, f_entries_per_column[static_cast<std::size_t>(coarse)]);
        const int keep = std::min(
            requested, static_cast<int>(ranking.size()));
        if (keep < static_cast<int>(ranking.size())) {
            std::nth_element(
                ranking.begin(), ranking.begin() + keep, ranking.end(),
                [](const auto& lhs, const auto& rhs) {
                    if (lhs.first != rhs.first) return lhs.first > rhs.first;
                    return lhs.second < rhs.second;
                });
            ranking.resize(static_cast<std::size_t>(keep));
        }
        auto& support = supports[static_cast<std::size_t>(coarse)];
        support.reserve(static_cast<std::size_t>(keep));
        for (const auto& entry : ranking) support.push_back(entry.second);
        std::sort(support.begin(), support.end());
    }
    return supports;
}

// Sparse oracle used only for comparison: select the same number of F entries
// per column as a candidate from the global basis, then re-minimize energy on
// those supports. It is not available in a practical setup because it uses the
// global reference.
inline MatchedReferenceResult build_budget_matched_reference(
    const StructuredGrid& grid, const SparseMatrix& a,
    const SparseMatrix& global_reference,
    const SparseMatrix& budget_source,
    InterpolationOptions interpolation_options) {
    const auto begin = std::chrono::steady_clock::now();
    const std::vector<int> counts = interpolation_f_entries_per_column(
        grid, budget_source, interpolation_options.thread_count);
    auto supports = select_reference_supports(
        grid, a, global_reference, counts,
        interpolation_options.thread_count);
    interpolation_options.strategy =
        InterpolationStrategy::LocalEnergyMinimum;
    interpolation_options.local_tolerance = std::min(
        interpolation_options.local_tolerance, 1.0e-8);
    const InterpolationResult interpolation =
        build_energy_interpolation_on_supports(
            grid, a, supports, interpolation_options);
    const double build_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - begin).count();
    return {
        interpolation.prolongation, std::move(supports), build_ms};
}

} // namespace tgi
