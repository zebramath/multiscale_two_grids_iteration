#pragma once

#include "pde/diffusion_problem.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace tgi {

struct RelativePruningResult {
    SparseMatrix prolongation;
    double pruning_ms = 0.0;
};

inline RelativePruningResult prune_global_interpolation_relative(
    const StructuredGrid& grid, const SparseMatrix& global_reference,
    double relative_threshold) {
    const auto begin = std::chrono::steady_clock::now();
    Vector column_max(
        static_cast<std::size_t>(grid.coarse_size()), 0.0);
    for (int row = 0; row < global_reference.rows(); ++row) {
        for (int position =
                 global_reference.row_ptr()[static_cast<std::size_t>(row)];
             position < global_reference.row_ptr()[
                 static_cast<std::size_t>(row) + 1U];
             ++position) {
            const int column = global_reference.col_idx()[
                static_cast<std::size_t>(position)];
            column_max[static_cast<std::size_t>(column)] = std::max(
                column_max[static_cast<std::size_t>(column)],
                std::abs(global_reference.values()[
                    static_cast<std::size_t>(position)]));
        }
    }

    std::vector<int> row_ptr(
        static_cast<std::size_t>(grid.fine_size()) + 1U, 0);
    std::vector<int> col_idx;
    Vector values;
    col_idx.reserve(global_reference.nnz());
    values.reserve(global_reference.nnz());
    for (int row = 0; row < global_reference.rows(); ++row) {
        const bool coarse_row = grid.is_coarse_node(row);
        for (int position =
                 global_reference.row_ptr()[static_cast<std::size_t>(row)];
             position < global_reference.row_ptr()[
                 static_cast<std::size_t>(row) + 1U];
             ++position) {
            const int column = global_reference.col_idx()[
                static_cast<std::size_t>(position)];
            const double value = global_reference.values()[
                static_cast<std::size_t>(position)];
            const double threshold = relative_threshold *
                column_max[static_cast<std::size_t>(column)];
            if (coarse_row || std::abs(value) > threshold) {
                col_idx.push_back(column);
                values.push_back(value);
            }
        }
        row_ptr[static_cast<std::size_t>(row) + 1U] =
            static_cast<int>(values.size());
    }
    SparseMatrix prolongation(
        grid.fine_size(), grid.coarse_size(), std::move(row_ptr),
        std::move(col_idx), std::move(values));
    const double pruning_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - begin).count();
    return {std::move(prolongation), pruning_ms};
}

}
