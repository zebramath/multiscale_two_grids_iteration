#pragma once

#include "multigrid/two_grid_solver.hpp"
#include "pde/diffusion_problem.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <exception>
#include <functional>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

namespace tgi {

struct ResidualStrongSupportOptions {
    int base_patch_layers = 2;
    int maximum_extra_nodes_per_column = 64;
    int maximum_graph_hops = 256;
    double strong_edge_fraction = 0.25;
    int thread_count = 1;
};

struct FixedSupportExpansionReport {
    double selection_ms = 0.0;
    double maximum_off_diagonal = 0.0;
    double strong_edge_threshold = 0.0;
    int expanded_columns = 0;
    int total_extra_nodes = 0;
    int maximum_extra_nodes = 0;
    double mean_scaled_column_residual = 0.0;
    double maximum_scaled_column_residual = 0.0;
};

struct FixedSupportExpansionResult {
    std::vector<std::vector<int>> supports;
    FixedSupportExpansionReport report;
};

inline FixedSupportExpansionResult build_residual_strong_supports(
    const StructuredGrid& grid, const SparseMatrix& a,
    const SparseMatrix& base_prolongation,
    const ResidualStrongSupportOptions& options = {});

namespace support_expansion_detail {

using Clock = std::chrono::steady_clock;

inline double milliseconds(Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

struct QueueEntry {
    double score = 0.0;
    int hops = 0;
    int node = -1;
};

struct QueueLess {
    bool operator()(const QueueEntry& lhs,
                    const QueueEntry& rhs) const {
        if (lhs.score != rhs.score) return lhs.score < rhs.score;
        if (lhs.hops != rhs.hops) return lhs.hops > rhs.hops;
        return lhs.node > rhs.node;
    }
};

inline double maximum_off_diagonal(const SparseMatrix& a) {
    double maximum = 0.0;
    for (int row = 0; row < a.rows(); ++row) {
        for (int position =
                 a.row_ptr()[static_cast<std::size_t>(row)];
             position <
                 a.row_ptr()[static_cast<std::size_t>(row) + 1U];
             ++position) {
            if (a.col_idx()[static_cast<std::size_t>(position)] != row) {
                maximum = std::max(
                    maximum,
                    std::abs(a.values()[static_cast<std::size_t>(position)]));
            }
        }
    }
    return maximum;
}

inline SparseMatrix build_strong_graph(
    const SparseMatrix& a, double threshold) {
    std::vector<int> row_ptr(
        static_cast<std::size_t>(a.rows()) + 1U, 0);
    std::vector<int> col_idx;
    Vector weights;
    col_idx.reserve(a.nnz());
    weights.reserve(a.nnz());
    for (int row = 0; row < a.rows(); ++row) {
        for (int position =
                 a.row_ptr()[static_cast<std::size_t>(row)];
             position <
                 a.row_ptr()[static_cast<std::size_t>(row) + 1U];
             ++position) {
            const int col =
                a.col_idx()[static_cast<std::size_t>(position)];
            const double edge = std::abs(
                a.values()[static_cast<std::size_t>(position)]);
            if (col != row && edge >= threshold) {
                col_idx.push_back(col);
                weights.push_back(edge);
            }
        }
        row_ptr[static_cast<std::size_t>(row) + 1U] =
            static_cast<int>(weights.size());
    }
    return SparseMatrix(
        a.rows(), a.cols(), std::move(row_ptr),
        std::move(col_idx), std::move(weights));
}

}


inline FixedSupportExpansionResult build_residual_strong_supports(
    const StructuredGrid& grid, const SparseMatrix& a,
    const SparseMatrix& base_prolongation,
    const ResidualStrongSupportOptions& options) {
    if (a.rows() != grid.fine_size() ||
        a.cols() != grid.fine_size() ||
        base_prolongation.rows() != grid.fine_size() ||
        base_prolongation.cols() != grid.coarse_size()) {
        throw std::invalid_argument(
            "build_residual_strong_supports: incompatible dimensions");
    }
    if (options.base_patch_layers <= 0 ||
        options.maximum_extra_nodes_per_column < 0 ||
        options.maximum_graph_hops < 0 ||
        !(options.strong_edge_fraction >= 0.0)) {
        throw std::invalid_argument(
            "build_residual_strong_supports: invalid options");
    }

    FixedSupportExpansionResult result;
    result.supports.resize(
        static_cast<std::size_t>(grid.coarse_size()));
    const auto begin = support_expansion_detail::Clock::now();
    result.report.maximum_off_diagonal =
        support_expansion_detail::maximum_off_diagonal(a);
    result.report.strong_edge_threshold =
        options.strong_edge_fraction *
        result.report.maximum_off_diagonal;
    if (!(result.report.maximum_off_diagonal > 0.0)) {
        throw std::runtime_error(
            "build_residual_strong_supports: matrix has no graph edges");
    }
    const SparseMatrix strong_graph =
        support_expansion_detail::build_strong_graph(
            a, result.report.strong_edge_threshold);

    const SparseMatrix ap = sparse_multiply(
        a, base_prolongation, 0.0, options.thread_count);
    const SparseMatrix ap_transpose =
        ap.transpose(options.thread_count);
    const Vector diagonal = a.diagonal();
    std::vector<double> column_residuals(
        static_cast<std::size_t>(grid.coarse_size()), 0.0);
    for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
        double squared = 0.0;
        for (int position =
                 ap_transpose.row_ptr()[static_cast<std::size_t>(coarse)];
             position <
                 ap_transpose.row_ptr()[
                     static_cast<std::size_t>(coarse) + 1U];
             ++position) {
            const int node =
                ap_transpose.col_idx()[static_cast<std::size_t>(position)];
            if (grid.is_coarse_node(node)) continue;
            const double scaled =
                ap_transpose.values()[static_cast<std::size_t>(position)] /
                diagonal[static_cast<std::size_t>(node)];
            squared += scaled * scaled;
        }
        const double value = std::sqrt(squared);
        column_residuals[static_cast<std::size_t>(coarse)] = value;
        result.report.mean_scaled_column_residual += value;
        result.report.maximum_scaled_column_residual = std::max(
            result.report.maximum_scaled_column_residual, value);
    }
    result.report.mean_scaled_column_residual /=
        static_cast<double>(grid.coarse_size());
    std::vector<unsigned char> coarse_mask(
        static_cast<std::size_t>(grid.fine_size()), 0U);
    for (int coarse = 0; coarse < grid.coarse_size(); ++coarse) {
        coarse_mask[static_cast<std::size_t>(
            grid.coarse_fine_id(coarse))] = 1U;
    }
    std::vector<int> extras_per_column(
        static_cast<std::size_t>(grid.coarse_size()), 0);
    const unsigned int requested = options.thread_count > 0
        ? static_cast<unsigned int>(options.thread_count)
        : std::thread::hardware_concurrency();
    const int worker_count = std::max(
        1, std::min(
               grid.coarse_size(),
               static_cast<int>(requested == 0U ? 1U : requested)));
    std::atomic<int> next_coarse{0};
    std::exception_ptr worker_error;
    std::mutex error_mutex;

    auto worker = [&]() {
        Vector residual(
            static_cast<std::size_t>(grid.fine_size()), 0.0);
        std::vector<int> residual_active;
        std::vector<unsigned char> in_support(
            static_cast<std::size_t>(grid.fine_size()), 0U);
        std::vector<unsigned char> visited(
            static_cast<std::size_t>(grid.fine_size()), 0U);
        std::vector<int> visited_nodes;
        try {
            while (true) {
                const int coarse = next_coarse.fetch_add(
                    1, std::memory_order_relaxed);
                if (coarse >= grid.coarse_size()) break;
                auto& support =
                    result.supports[static_cast<std::size_t>(coarse)];
                grid.patch_f_nodes(
                    coarse, options.base_patch_layers, support);
                for (int node : support) {
                    in_support[static_cast<std::size_t>(node)] = 1U;
                }

                residual_active.clear();
                for (int position =
                         ap_transpose.row_ptr()[
                             static_cast<std::size_t>(coarse)];
                     position <
                         ap_transpose.row_ptr()[
                             static_cast<std::size_t>(coarse) + 1U];
                     ++position) {
                    const int node =
                        ap_transpose.col_idx()[
                            static_cast<std::size_t>(position)];
                    if (coarse_mask[
                            static_cast<std::size_t>(node)] != 0U) {
                        continue;
                    }
                    const double scaled =
                        ap_transpose.values()[
                            static_cast<std::size_t>(position)] /
                        diagonal[static_cast<std::size_t>(node)];
                    residual[static_cast<std::size_t>(node)] = scaled;
                    residual_active.push_back(node);
                }

                int extras = 0;
                if (column_residuals[
                        static_cast<std::size_t>(coarse)] >
                        0.0 &&
                    options.maximum_extra_nodes_per_column > 0) {
                    std::priority_queue<
                        support_expansion_detail::QueueEntry,
                        std::vector<
                            support_expansion_detail::QueueEntry>,
                        support_expansion_detail::QueueLess> queue;

                    for (int node : support) {
                        const double source_score = std::abs(
                            residual[static_cast<std::size_t>(node)]);
                        for (int position =
                                 strong_graph.row_ptr()[
                                     static_cast<std::size_t>(node)];
                             position <
                                 strong_graph.row_ptr()[
                                     static_cast<std::size_t>(node)
                                     + 1U];
                             ++position) {
                            const int neighbor =
                                strong_graph.col_idx()[
                                    static_cast<std::size_t>(
                                        position)];
                            if (in_support[
                                    static_cast<std::size_t>(
                                        neighbor)] != 0U) {
                                continue;
                            }
                            const double edge =
                                strong_graph.values()[
                                    static_cast<std::size_t>(
                                        position)];
                            const double neighbor_score = std::abs(
                                residual[static_cast<std::size_t>(
                                    neighbor)]);
                            queue.push({
                                std::max(
                                    source_score, neighbor_score) +
                                    edge /
                                        result.report
                                            .maximum_off_diagonal,
                                1,
                                neighbor
                            });
                        }
                    }

                    while (!queue.empty() &&
                           extras <
                               options
                                   .maximum_extra_nodes_per_column) {
                        const auto entry = queue.top();
                        queue.pop();
                        if (entry.hops >
                            options.maximum_graph_hops) {
                            continue;
                        }
                        const std::size_t offset =
                            static_cast<std::size_t>(entry.node);
                        if (visited[offset] != 0U ||
                            in_support[offset] != 0U) {
                            continue;
                        }
                        visited[offset] = 1U;
                        visited_nodes.push_back(entry.node);
                        if (coarse_mask[offset] == 0U) {
                            support.push_back(entry.node);
                            in_support[offset] = 1U;
                            ++extras;
                        }
                        for (int position =
                                 strong_graph.row_ptr()[offset];
                             position <
                                 strong_graph.row_ptr()[offset + 1U];
                             ++position) {
                            const int neighbor =
                                strong_graph.col_idx()[
                                    static_cast<std::size_t>(
                                        position)];
                            const std::size_t neighbor_offset =
                                static_cast<std::size_t>(neighbor);
                            if (visited[neighbor_offset] != 0U ||
                                in_support[neighbor_offset] != 0U) {
                                continue;
                            }
                            const double edge =
                                strong_graph.values()[
                                    static_cast<std::size_t>(
                                        position)];
                            queue.push({
                                entry.score * 0.999 +
                                    1.0e-3 * edge /
                                        result.report
                                            .maximum_off_diagonal,
                                entry.hops + 1,
                                neighbor
                            });
                        }
                    }
                }

                std::sort(support.begin(), support.end());
                support.erase(
                    std::unique(support.begin(), support.end()),
                    support.end());
                extras_per_column[
                    static_cast<std::size_t>(coarse)] = extras;

                for (int node : support) {
                    in_support[static_cast<std::size_t>(node)] = 0U;
                }
                for (int node : residual_active) {
                    residual[static_cast<std::size_t>(node)] = 0.0;
                }
                for (int node : visited_nodes) {
                    visited[static_cast<std::size_t>(node)] = 0U;
                }
                visited_nodes.clear();
            }
        } catch (...) {
            std::lock_guard<std::mutex> lock(error_mutex);
            if (!worker_error) {
                worker_error = std::current_exception();
            }
            next_coarse.store(grid.coarse_size());
        }
    };

    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(worker_count));
    for (int index = 0; index < worker_count; ++index) {
        workers.emplace_back(worker);
    }
    for (auto& thread : workers) thread.join();
    if (worker_error) std::rethrow_exception(worker_error);

    for (int extras : extras_per_column) {
        if (extras > 0) ++result.report.expanded_columns;
        result.report.total_extra_nodes += extras;
        result.report.maximum_extra_nodes = std::max(
            result.report.maximum_extra_nodes, extras);
    }
    result.report.selection_ms =
        support_expansion_detail::milliseconds(
            begin, support_expansion_detail::Clock::now());
    return result;
}

}
