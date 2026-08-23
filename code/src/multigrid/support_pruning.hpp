#pragma once

#include "multigrid/energy_interpolation.hpp"
#include "pde/diffusion_problem.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace tgi {

struct SupportPruningOptions {
    int base_patch_layers = 2;
    int maximum_extra_nodes_per_column = 64;
    double relative_magnitude_threshold = 0.0;
    double refinement_tolerance = 1.0e-6;
    int thread_count = 1;
};

struct SupportPruningReport {
    int pruned_columns = 0;
    int input_extra_nodes = 0;
    int retained_extra_nodes = 0;
    int removed_extra_nodes = 0;
    double selection_ms = 0.0;
    double refinement_ms = 0.0;
    double total_ms = 0.0;
};

struct SupportPruningResult {
    SparseMatrix prolongation;
    std::vector<std::vector<int>> supports;
    SupportPruningReport report;
};

inline SupportPruningResult prune_energy_supports_by_magnitude(
    const StructuredGrid& grid, const SparseMatrix& a,
    const std::vector<std::vector<int>>& expanded_supports,
    const SparseMatrix& expanded_prolongation,
    const InterpolationOptions& interpolation_options,
    const SupportPruningOptions& options = {}) {
    using Clock = std::chrono::steady_clock;
    const auto milliseconds = [](Clock::time_point begin,
                                 Clock::time_point end) {
        return std::chrono::duration<double, std::milli>(end - begin).count();
    };
    if (expanded_supports.size() !=
            static_cast<std::size_t>(grid.coarse_size()) ||
        expanded_prolongation.rows() != grid.fine_size() ||
        expanded_prolongation.cols() != grid.coarse_size() ||
        options.base_patch_layers <= 0 ||
        options.maximum_extra_nodes_per_column < 0 ||
        !(options.relative_magnitude_threshold >= 0.0) ||
        !(options.refinement_tolerance > 0.0)) {
        throw std::invalid_argument(
            "prune_energy_supports_by_magnitude: invalid input");
    }

    const auto total_begin = Clock::now();
    const auto selection_begin = Clock::now();
    SupportPruningResult result;
    result.supports.resize(static_cast<std::size_t>(grid.coarse_size()));
    const SparseMatrix transpose = expanded_prolongation.transpose(
        options.thread_count);
    std::vector<unsigned char> is_base(
        static_cast<std::size_t>(grid.fine_size()), 0U);
    Vector value(static_cast<std::size_t>(grid.fine_size()), 0.0);
    std::vector<int> active_values;
    std::vector<unsigned char> refine_column(
        static_cast<std::size_t>(grid.coarse_size()), 0U);

    for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
        const std::size_t offset = static_cast<std::size_t>(coarse);
        std::vector<int> base = grid.patch_f_nodes(
            coarse, options.base_patch_layers);
        for (int node : base) {
            is_base[static_cast<std::size_t>(node)] = 1U;
        }
        active_values.clear();
        for (int position = transpose.row_ptr()[offset];
             position < transpose.row_ptr()[offset + 1U]; ++position) {
            const int node = transpose.col_idx()[
                static_cast<std::size_t>(position)];
            if (grid.is_coarse_node(node)) continue;
            value[static_cast<std::size_t>(node)] = std::abs(
                transpose.values()[static_cast<std::size_t>(position)]);
            active_values.push_back(node);
        }

        std::vector<int> extras;
        extras.reserve(expanded_supports[offset].size());
        for (int node : expanded_supports[offset]) {
            if (is_base[static_cast<std::size_t>(node)] == 0U) {
                extras.push_back(node);
            }
        }
        result.report.input_extra_nodes +=
            static_cast<int>(extras.size());
        std::sort(extras.begin(), extras.end(), [&](int lhs, int rhs) {
            const double lhs_value = value[static_cast<std::size_t>(lhs)];
            const double rhs_value = value[static_cast<std::size_t>(rhs)];
            if (lhs_value != rhs_value) return lhs_value > rhs_value;
            return lhs < rhs;
        });
        const double maximum = extras.empty()
            ? 0.0
            : value[static_cast<std::size_t>(extras.front())];
        std::vector<int> retained;
        retained.reserve(static_cast<std::size_t>(
            options.maximum_extra_nodes_per_column));
        for (int node : extras) {
            if (static_cast<int>(retained.size()) >=
                options.maximum_extra_nodes_per_column) {
                break;
            }
            if (maximum > 0.0 &&
                value[static_cast<std::size_t>(node)] <
                    options.relative_magnitude_threshold * maximum) {
                break;
            }
            retained.push_back(node);
        }
        result.report.retained_extra_nodes +=
            static_cast<int>(retained.size());
        result.report.removed_extra_nodes +=
            static_cast<int>(extras.size() - retained.size());
        if (retained.size() != extras.size()) {
            refine_column[offset] = 1U;
            ++result.report.pruned_columns;
        }
        auto& support = result.supports[offset];
        support = std::move(base);
        support.insert(support.end(), retained.begin(), retained.end());
        std::sort(support.begin(), support.end());
        support.erase(std::unique(support.begin(), support.end()), support.end());

        for (int node : support) {
            is_base[static_cast<std::size_t>(node)] = 0U;
        }
        for (int node : active_values) {
            value[static_cast<std::size_t>(node)] = 0.0;
        }
    }
    result.report.selection_ms = milliseconds(
        selection_begin, Clock::now());

    const auto refinement_begin = Clock::now();
    InterpolationOptions refine_options = interpolation_options;
    refine_options.local_tolerance = std::min(
        refine_options.local_tolerance, options.refinement_tolerance);
    refine_options.thread_count = options.thread_count;
    result.prolongation = refine_selected_energy_interpolation_on_supports(
        grid, a, result.supports, refine_column,
        expanded_prolongation, refine_options).prolongation;
    result.report.refinement_ms = milliseconds(
        refinement_begin, Clock::now());
    result.report.total_ms = milliseconds(total_begin, Clock::now());
    return result;
}

} // namespace tgi
