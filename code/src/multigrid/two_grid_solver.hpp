#pragma once

#include "core/linear_algebra.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace tgi {

struct CoarseSetupReport {
    double interpolation_energy = 0.0;
};

class TwoGridCycle {
public:
    struct Workspace {
        Vector coarse_rhs;
        Vector coarse_solution;
        Vector coarse_work;
    };

    TwoGridCycle(const SparseMatrix& a, const SparseMatrix& p,
                 int smoothing_steps = 1, int setup_threads = 1);

    double iterate(const Vector& rhs, Vector& solution, Vector& residual,
                   Workspace& workspace) const;
    void solve_coarse_system(const Vector& rhs, Vector& solution,
                             Vector& work) const {
        coarse_solver_.solve(rhs, solution, work);
    }
    const SparseMatrix& coarse_matrix() const { return coarse_matrix_; }
    const CoarseSetupReport& setup_report() const { return setup_report_; }

private:
    void solve_gauss_seidel_sweep(const Vector& rhs, bool forward,
                                  Vector& solution) const;
    void restrict_fine_residual(const Vector& rhs, const Vector& solution,
                                Vector& coarse_rhs) const;

    const SparseMatrix& a_;
    const SparseMatrix& p_;
    Vector inverse_diagonal_;
    std::vector<int> diagonal_position_;
    int smoothing_steps_;
    int application_threads_ = 1;
    SparseMatrix p_transpose_;
    SparseMatrix coarse_matrix_;
    SparseCholesky coarse_solver_;
    CoarseSetupReport setup_report_;
};

enum class TwoGridIterationStatus {
    Converged,
    SlowAtLimit,
    Diverged
};

inline const char* two_grid_status_name(TwoGridIterationStatus status) {
    switch (status) {
        case TwoGridIterationStatus::Converged:
            return "converged";
        case TwoGridIterationStatus::SlowAtLimit:
            return "slow-limit";
        case TwoGridIterationStatus::Diverged:
            return "diverged";
    }
    return "unknown";
}

struct TwoGridIterationResult {
    Vector solution;
    int cycles = 0;
    double relative_residual = 0.0;
    double best_relative_residual = 1.0;
    double effective_factor = 1.0;
    double tail_factor = 1.0;
    TwoGridIterationStatus status = TwoGridIterationStatus::SlowAtLimit;
    bool converged = false;
};

inline TwoGridIterationResult solve_two_grid(
    const Vector& rhs, const TwoGridCycle& cycle,
    double relative_tolerance = 1e-8, int max_cycles = 40000);

template <class Cycle>
inline TwoGridIterationResult solve_stationary_cycles(
    const Vector& rhs, const Cycle& cycle, double relative_tolerance,
    int max_cycles, const char* method_name);

namespace two_grid_solver_detail {

inline void nested_dissection_rectangle(int side, int separator_width,
                                        int xmin, int xmax,
                                        int ymin, int ymax,
                                        std::vector<int>& permutation) {
    const int width = xmax - xmin;
    const int height = ymax - ymin;
    if (width <= 0 || height <= 0) return;
    if (width <= 2 * separator_width + 4 ||
        height <= 2 * separator_width + 4) {
        for (int y = ymin; y < ymax; ++y) {
            for (int x = xmin; x < xmax; ++x) {
                permutation.push_back(y * side + x);
            }
        }
        return;
    }
    if (width >= height) {
        const int separator_begin =
            xmin + (width - separator_width) / 2;
        const int separator_end =
            separator_begin + separator_width;
        nested_dissection_rectangle(
            side, separator_width, xmin, separator_begin,
            ymin, ymax, permutation);
        nested_dissection_rectangle(
            side, separator_width, separator_end, xmax,
            ymin, ymax, permutation);
        for (int y = ymin; y < ymax; ++y) {
            for (int x = separator_begin; x < separator_end; ++x) {
                permutation.push_back(y * side + x);
            }
        }
    } else {
        const int separator_begin =
            ymin + (height - separator_width) / 2;
        const int separator_end =
            separator_begin + separator_width;
        nested_dissection_rectangle(
            side, separator_width, xmin, xmax,
            ymin, separator_begin, permutation);
        nested_dissection_rectangle(
            side, separator_width, xmin, xmax,
            separator_end, ymax, permutation);
        for (int y = separator_begin; y < separator_end; ++y) {
            for (int x = xmin; x < xmax; ++x) {
                permutation.push_back(y * side + x);
            }
        }
    }
}

inline std::vector<int> coarse_ordering(const SparseMatrix& matrix) {
    const int size = matrix.rows();
    const int side = static_cast<int>(
        std::llround(std::sqrt(static_cast<double>(size))));
    if (side * side != size) return {};
    int separator_width = 1;
    for (int row = 0; row < size; ++row) {
        const int row_x = row % side;
        const int row_y = row / side;
        for (int position =
                 matrix.row_ptr()[static_cast<std::size_t>(row)];
             position <
                 matrix.row_ptr()[static_cast<std::size_t>(row) + 1U];
             ++position) {
            const int col =
                matrix.col_idx()[static_cast<std::size_t>(position)];
            separator_width = std::max(
                separator_width,
                std::max(std::abs(row_x - col % side),
                         std::abs(row_y - col / side)));
        }
    }
    std::vector<int> permutation;
    permutation.reserve(static_cast<std::size_t>(size));
    nested_dissection_rectangle(
        side, separator_width, 0, side, 0, side, permutation);
    return permutation;
}

inline SparseMatrix symmetrize_upper(const SparseMatrix& upper) {
    std::vector<int> lower_counts(
        static_cast<std::size_t>(upper.rows()), 0);
    std::vector<int> row_ptr(
        static_cast<std::size_t>(upper.rows()) + 1U, 0);
    for (int row = 0; row < upper.rows(); ++row) {
        int row_entries = 0;
        for (int position =
                 upper.row_ptr()[static_cast<std::size_t>(row)];
             position <
                 upper.row_ptr()[static_cast<std::size_t>(row) + 1U];
             ++position) {
            const int col =
                upper.col_idx()[static_cast<std::size_t>(position)];
            ++row_entries;
            if (col != row) {
                ++lower_counts[static_cast<std::size_t>(col)];
            }
        }
        row_ptr[static_cast<std::size_t>(row) + 1U] = row_entries;
    }
    for (int row = 0; row < upper.rows(); ++row) {
        row_ptr[static_cast<std::size_t>(row) + 1U] +=
            row_ptr[static_cast<std::size_t>(row)] +
            lower_counts[static_cast<std::size_t>(row)];
    }

    std::vector<int> lower_next = row_ptr;
    std::vector<int> upper_next(static_cast<std::size_t>(upper.rows()));
    for (int row = 0; row < upper.rows(); ++row) {
        upper_next[static_cast<std::size_t>(row)] =
            row_ptr[static_cast<std::size_t>(row)] +
            lower_counts[static_cast<std::size_t>(row)];
    }
    std::vector<int> col_idx(static_cast<std::size_t>(row_ptr.back()));
    Vector values(static_cast<std::size_t>(row_ptr.back()));
    for (int row = 0; row < upper.rows(); ++row) {
        for (int position =
                 upper.row_ptr()[static_cast<std::size_t>(row)];
             position <
                 upper.row_ptr()[static_cast<std::size_t>(row) + 1U];
             ++position) {
            const int col =
                upper.col_idx()[static_cast<std::size_t>(position)];
            const double value =
                upper.values()[static_cast<std::size_t>(position)];
            const int target =
                upper_next[static_cast<std::size_t>(row)]++;
            col_idx[static_cast<std::size_t>(target)] = col;
            values[static_cast<std::size_t>(target)] = value;
            if (col != row) {
                const int mirror =
                    lower_next[static_cast<std::size_t>(col)]++;
                col_idx[static_cast<std::size_t>(mirror)] = row;
                values[static_cast<std::size_t>(mirror)] = value;
            }
        }
    }
    return SparseMatrix(
        upper.rows(), upper.cols(), std::move(row_ptr),
        std::move(col_idx), std::move(values));
}

inline SparseMatrix multiply_sparse_matrices(
    const SparseMatrix& lhs, const SparseMatrix& rhs,
    double drop_tolerance, int thread_count, bool upper_triangle_only);

inline void multiply_parallel(const SparseMatrix& matrix, const Vector& x,
                              Vector& result, int thread_count) {
#if defined(_OPENMP)
    matrix.multiply(x, result, thread_count);
#else
    if (thread_count <= 1) {
        matrix.multiply(x, result);
        return;
    }
    result.resize(static_cast<std::size_t>(matrix.rows()));
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(thread_count));
    for (int worker = 0; worker < thread_count; ++worker) {
        workers.emplace_back([&, worker]() {
            const int first = matrix.rows() * worker / thread_count;
            const int last = matrix.rows() * (worker + 1) / thread_count;
            for (int row = first; row < last; ++row) {
                double sum = 0.0;
                for (int position =
                         matrix.row_ptr()[static_cast<std::size_t>(row)];
                     position <
                         matrix.row_ptr()[static_cast<std::size_t>(row) + 1U];
                     ++position) {
                    sum +=
                        matrix.values()[static_cast<std::size_t>(position)] *
                        x[static_cast<std::size_t>(
                            matrix.col_idx()[
                                static_cast<std::size_t>(position)])];
                }
                result[static_cast<std::size_t>(row)] = sum;
            }
        });
    }
    for (auto& thread : workers) thread.join();
#endif
}

inline void multiply_add_parallel(const SparseMatrix& matrix, double alpha,
                                  const Vector& x, Vector& result,
                                  int thread_count) {
#if defined(_OPENMP)
    matrix.multiply_add(alpha, x, result, thread_count);
#else
    if (thread_count <= 1) {
        matrix.multiply_add(alpha, x, result);
        return;
    }
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(thread_count));
    for (int worker = 0; worker < thread_count; ++worker) {
        workers.emplace_back([&, worker]() {
            const int first = matrix.rows() * worker / thread_count;
            const int last = matrix.rows() * (worker + 1) / thread_count;
            for (int row = first; row < last; ++row) {
                double sum = 0.0;
                for (int position =
                         matrix.row_ptr()[static_cast<std::size_t>(row)];
                     position <
                         matrix.row_ptr()[static_cast<std::size_t>(row) + 1U];
                     ++position) {
                    sum +=
                        matrix.values()[static_cast<std::size_t>(position)] *
                        x[static_cast<std::size_t>(
                            matrix.col_idx()[
                                static_cast<std::size_t>(position)])];
                }
                result[static_cast<std::size_t>(row)] += alpha * sum;
            }
        });
    }
    for (auto& thread : workers) thread.join();
#endif
}

inline void gauss_seidel_sweep(
    const SparseMatrix& matrix, const Vector& inverse_diagonal,
    const std::vector<int>& diagonal_position, const Vector& rhs,
    bool forward, Vector& solution) {
    const int begin = forward ? 0 : matrix.rows() - 1;
    const int end = forward ? matrix.rows() : -1;
    const int increment = forward ? 1 : -1;
    for (int row = begin; row != end; row += increment) {
        double value = rhs[static_cast<std::size_t>(row)];
        const int diagonal =
            diagonal_position[static_cast<std::size_t>(row)];
        for (int position = matrix.row_ptr()[static_cast<std::size_t>(row)];
             position < diagonal; ++position) {
            value -= matrix.values()[static_cast<std::size_t>(position)] *
                solution[static_cast<std::size_t>(
                    matrix.col_idx()[static_cast<std::size_t>(position)])];
        }
        for (int position = diagonal + 1;
             position < matrix.row_ptr()[static_cast<std::size_t>(row) + 1U];
             ++position) {
            value -= matrix.values()[static_cast<std::size_t>(position)] *
                solution[static_cast<std::size_t>(
                    matrix.col_idx()[static_cast<std::size_t>(position)])];
        }
        solution[static_cast<std::size_t>(row)] = value *
            inverse_diagonal[static_cast<std::size_t>(row)];
    }
}

}

inline SparseMatrix two_grid_solver_detail::multiply_sparse_matrices(
    const SparseMatrix& lhs, const SparseMatrix& rhs,
    double drop_tolerance, int thread_count,
    bool upper_triangle_only) {
    unsigned int requested = thread_count > 0
        ? static_cast<unsigned int>(thread_count)
        : std::thread::hardware_concurrency();
    if (requested == 0U) requested = 1U;
    const int worker_count = std::max(
        1, std::min(lhs.rows(), static_cast<int>(requested)));
    const bool use_dense_upper_path = upper_triangle_only;

    struct RowBlockProduct {
        int first_row = 0;
        int last_row = 0;
        std::vector<int> row_offsets;
        std::vector<int> columns;
        std::vector<double> entries;
    };
    std::vector<RowBlockProduct> products(
        static_cast<std::size_t>(worker_count));
    std::exception_ptr worker_error;
    std::mutex error_mutex;

    auto worker = [&](int worker_id) {
        try {
            RowBlockProduct& product =
                products[static_cast<std::size_t>(worker_id)];
            product.first_row =
                lhs.rows() * worker_id / worker_count;
            product.last_row =
                lhs.rows() * (worker_id + 1) / worker_count;
            const int local_rows =
                product.last_row - product.first_row;
            product.row_offsets.assign(
                static_cast<std::size_t>(local_rows) + 1U, 0);
            const std::size_t reserve_entries =
                static_cast<std::size_t>(local_rows) *
                static_cast<std::size_t>(
                    std::min(rhs.cols(), 64));
            product.columns.reserve(reserve_entries);
            product.entries.reserve(reserve_entries);

            std::vector<double> accumulator(
                static_cast<std::size_t>(rhs.cols()), 0.0);
            std::vector<int> marker(
                static_cast<std::size_t>(rhs.cols()), -1);
            std::vector<int> touched;
            touched.reserve(static_cast<std::size_t>(
                std::min(rhs.cols(), 256)));
            for (int row = product.first_row;
                 row < product.last_row; ++row) {
                if (use_dense_upper_path) {
                    std::fill(
                        accumulator.begin() + row,
                        accumulator.end(), 0.0);
                    for (int lhs_position =
                             lhs.row_ptr()[static_cast<std::size_t>(row)];
                         lhs_position <
                             lhs.row_ptr()[
                                 static_cast<std::size_t>(row) + 1U];
                         ++lhs_position) {
                        const int inner =
                            lhs.col_idx()[
                                static_cast<std::size_t>(lhs_position)];
                        const double lhs_value =
                            lhs.values()[
                                static_cast<std::size_t>(lhs_position)];
                        const int rhs_row_begin =
                            rhs.row_ptr()[static_cast<std::size_t>(inner)];
                        const int rhs_end =
                            rhs.row_ptr()[
                                static_cast<std::size_t>(inner) + 1U];
                        const int rhs_begin = static_cast<int>(
                            std::lower_bound(
                                rhs.col_idx().begin() + rhs_row_begin,
                                rhs.col_idx().begin() + rhs_end,
                                row) -
                            rhs.col_idx().begin());
                        for (int rhs_position = rhs_begin;
                             rhs_position < rhs_end;
                             ++rhs_position) {
                            const int col =
                                rhs.col_idx()[
                                    static_cast<std::size_t>(rhs_position)];
                            accumulator[static_cast<std::size_t>(col)] +=
                                lhs_value *
                                rhs.values()[
                                    static_cast<std::size_t>(rhs_position)];
                        }
                    }
                    for (int col = row; col < rhs.cols(); ++col) {
                        const double value =
                            accumulator[static_cast<std::size_t>(col)];
                        if (std::abs(value) > drop_tolerance) {
                            product.columns.push_back(col);
                            product.entries.push_back(value);
                        }
                    }
                    product.row_offsets[
                        static_cast<std::size_t>(
                            row - product.first_row) + 1U] =
                        static_cast<int>(product.entries.size());
                    continue;
                }

                touched.clear();
                for (int lhs_position =
                         lhs.row_ptr()[static_cast<std::size_t>(row)];
                     lhs_position <
                         lhs.row_ptr()[static_cast<std::size_t>(row) + 1U];
                     ++lhs_position) {
                    const int inner =
                        lhs.col_idx()[static_cast<std::size_t>(lhs_position)];
                    const double lhs_value =
                        lhs.values()[static_cast<std::size_t>(lhs_position)];
                    const int rhs_begin =
                        rhs.row_ptr()[static_cast<std::size_t>(inner)];
                    const int rhs_end =
                        rhs.row_ptr()[static_cast<std::size_t>(inner) + 1U];
                    const int first_rhs = upper_triangle_only
                        ? static_cast<int>(
                              std::lower_bound(
                                  rhs.col_idx().begin() + rhs_begin,
                                  rhs.col_idx().begin() + rhs_end,
                                  row) -
                              rhs.col_idx().begin())
                        : rhs_begin;
                    for (int rhs_position = first_rhs;
                         rhs_position < rhs_end;
                         ++rhs_position) {
                        const int col =
                            rhs.col_idx()[static_cast<std::size_t>(rhs_position)];
                        if (marker[static_cast<std::size_t>(col)] != row) {
                            marker[static_cast<std::size_t>(col)] = row;
                            accumulator[static_cast<std::size_t>(col)] = 0.0;
                            touched.push_back(col);
                        }
                        accumulator[static_cast<std::size_t>(col)] +=
                            lhs_value *
                            rhs.values()[static_cast<std::size_t>(rhs_position)];
                    }
                }
                const int first_candidate =
                    upper_triangle_only ? row : 0;
                const int candidate_count =
                    rhs.cols() - first_candidate;
                const bool dense_scan =
                    candidate_count > 0 &&
                    touched.size() * 2U >=
                        static_cast<std::size_t>(candidate_count);
                if (dense_scan) {
                    for (int col = first_candidate;
                         col < rhs.cols(); ++col) {
                        if (marker[static_cast<std::size_t>(col)] != row) {
                            continue;
                        }
                        const double value =
                            accumulator[static_cast<std::size_t>(col)];
                        if (std::abs(value) > drop_tolerance) {
                            product.columns.push_back(col);
                            product.entries.push_back(value);
                        }
                    }
                } else {
                    std::sort(touched.begin(), touched.end());
                    for (int col : touched) {
                        const double value =
                            accumulator[static_cast<std::size_t>(col)];
                        if (std::abs(value) > drop_tolerance) {
                            product.columns.push_back(col);
                            product.entries.push_back(value);
                        }
                    }
                }
                product.row_offsets[
                    static_cast<std::size_t>(
                        row - product.first_row) + 1U] =
                    static_cast<int>(product.entries.size());
            }
        } catch (...) {
            std::lock_guard<std::mutex> lock(error_mutex);
            if (!worker_error) worker_error = std::current_exception();
        }
    };

    if (worker_count == 1) {
        worker(0);
    } else {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) num_threads(worker_count)
        for (int worker_id = 0;
             worker_id < worker_count; ++worker_id) {
            worker(worker_id);
        }
#else
        std::vector<std::thread> workers;
        workers.reserve(static_cast<std::size_t>(worker_count));
        for (int worker_id = 0; worker_id < worker_count; ++worker_id) {
            workers.emplace_back(worker, worker_id);
        }
        for (auto& thread : workers) thread.join();
#endif
    }
    if (worker_error) std::rethrow_exception(worker_error);

    if (worker_count == 1) {
        RowBlockProduct& product = products.front();
        SparseMatrix result(
            lhs.rows(), rhs.cols(), std::move(product.row_offsets),
            std::move(product.columns), std::move(product.entries));
        return upper_triangle_only
            ? symmetrize_upper(result)
            : std::move(result);
    }

    std::vector<int> row_ptr(
        static_cast<std::size_t>(lhs.rows()) + 1U, 0);
    for (const RowBlockProduct& product : products) {
        for (int row = product.first_row; row < product.last_row; ++row) {
            const int local_row = row - product.first_row;
            const int count =
                product.row_offsets[
                    static_cast<std::size_t>(local_row) + 1U] -
                product.row_offsets[static_cast<std::size_t>(local_row)];
            row_ptr[static_cast<std::size_t>(row) + 1U] =
                row_ptr[static_cast<std::size_t>(row)] + count;
        }
    }
    std::vector<int> col_idx(
        static_cast<std::size_t>(row_ptr.back()));
    std::vector<double> values(
        static_cast<std::size_t>(row_ptr.back()));
    for (RowBlockProduct& product : products) {
        const std::size_t begin = static_cast<std::size_t>(
            row_ptr[static_cast<std::size_t>(product.first_row)]);
        std::move(product.columns.begin(), product.columns.end(),
                  col_idx.begin() + static_cast<std::ptrdiff_t>(begin));
        std::move(product.entries.begin(), product.entries.end(),
                  values.begin() + static_cast<std::ptrdiff_t>(begin));
    }
    SparseMatrix result(lhs.rows(), rhs.cols(), std::move(row_ptr),
                        std::move(col_idx), std::move(values));
    return upper_triangle_only
        ? symmetrize_upper(result)
        : std::move(result);
}

inline SparseMatrix galerkin_coarse_operator(
    const SparseMatrix& fine_matrix, const SparseMatrix& prolongation,
    int setup_threads = 1) {
    if (fine_matrix.rows() != fine_matrix.cols() ||
        fine_matrix.rows() != prolongation.rows() ||
        prolongation.cols() <= 0) {
        throw std::invalid_argument(
            "incompatible dimensions in Galerkin coarse operator");
    }
    const SparseMatrix restriction =
        prolongation.transpose(setup_threads);
    const SparseMatrix action =
        two_grid_solver_detail::multiply_sparse_matrices(
            fine_matrix, prolongation, 0.0, setup_threads, false);
    return two_grid_solver_detail::multiply_sparse_matrices(
        restriction, action, 0.0, setup_threads, true);
}

inline TwoGridCycle::TwoGridCycle(const SparseMatrix& a, const SparseMatrix& p,
                                  int smoothing_steps, int setup_threads)
    : a_(a), p_(p), smoothing_steps_(smoothing_steps) {
    inverse_diagonal_.resize(static_cast<std::size_t>(a_.rows()));
    diagonal_position_.resize(static_cast<std::size_t>(a_.rows()));
    for (int row = 0; row < a_.rows(); ++row) {
        const std::size_t i = static_cast<std::size_t>(row);
        int position = a_.row_ptr()[i];
        const int end = a_.row_ptr()[i + 1U];
        while (position < end &&
               a_.col_idx()[static_cast<std::size_t>(position)] < row) {
            ++position;
        }
        if (position == end ||
            a_.col_idx()[static_cast<std::size_t>(position)] != row) {
            throw std::runtime_error("fine matrix is missing a diagonal entry");
        }
        const double diagonal =
            a_.values()[static_cast<std::size_t>(position)];
        if (!(diagonal > 0.0)) {
            throw std::runtime_error(
                "fine matrix has nonpositive diagonal");
        }
        inverse_diagonal_[i] = 1.0 / diagonal;
        diagonal_position_[i] = position;
    }

    p_transpose_ = p_.transpose(setup_threads);
    const SparseMatrix ap = two_grid_solver_detail::multiply_sparse_matrices(
        a_, p_, 0.0, setup_threads, false);
    coarse_matrix_ = two_grid_solver_detail::multiply_sparse_matrices(
        p_transpose_, ap, 0.0, setup_threads, true);
    for (double value : coarse_matrix_.diagonal()) {
        setup_report_.interpolation_energy += value;
    }

    const std::vector<int> ordering =
        two_grid_solver_detail::coarse_ordering(coarse_matrix_);
    coarse_solver_.factorize(coarse_matrix_, ordering);
    unsigned int requested = setup_threads > 0
        ? static_cast<unsigned int>(setup_threads)
        : std::thread::hardware_concurrency();
    if (requested == 0U) requested = 1U;
#if defined(_OPENMP)
    constexpr std::size_t parallel_threshold = 40000U;
#else
    constexpr std::size_t parallel_threshold = 200000U;
#endif
    application_threads_ = p_.nnz() >= parallel_threshold
        ? std::max(
              1, std::min(
                     {8, static_cast<int>(requested),
                      p_.rows(), p_.cols()}))
        : 1;
}

inline void TwoGridCycle::solve_gauss_seidel_sweep(
    const Vector& rhs, bool forward, Vector& solution) const {
    two_grid_solver_detail::gauss_seidel_sweep(
        a_, inverse_diagonal_, diagonal_position_, rhs, forward, solution);
}

inline void TwoGridCycle::restrict_fine_residual(
    const Vector& rhs, const Vector& solution, Vector& coarse_rhs) const {
    coarse_rhs.assign(static_cast<std::size_t>(p_.cols()), 0.0);
    for (int row = 0; row < a_.rows(); ++row) {
        double product = 0.0;
        for (int position =
                 a_.row_ptr()[static_cast<std::size_t>(row)];
             position <
                 a_.row_ptr()[static_cast<std::size_t>(row) + 1U];
             ++position) {
            product += a_.values()[static_cast<std::size_t>(position)] *
                       solution[static_cast<std::size_t>(
                           a_.col_idx()[static_cast<std::size_t>(position)])];
        }
        const double residual =
            rhs[static_cast<std::size_t>(row)] - product;
        for (int position =
                 p_.row_ptr()[static_cast<std::size_t>(row)];
             position <
                 p_.row_ptr()[static_cast<std::size_t>(row) + 1U];
             ++position) {
            coarse_rhs[static_cast<std::size_t>(
                p_.col_idx()[static_cast<std::size_t>(position)])] +=
                p_.values()[static_cast<std::size_t>(position)] * residual;
        }
    }
}

inline double TwoGridCycle::iterate(
    const Vector& rhs, Vector& solution, Vector& residual,
    Workspace& workspace) const {
    for (int step = 0; step < smoothing_steps_; ++step) {
        solve_gauss_seidel_sweep(rhs, true, solution);
    }

    if (application_threads_ == 1) {
        restrict_fine_residual(rhs, solution, workspace.coarse_rhs);
    } else {
        a_.residual_squared(
            solution, rhs, residual, application_threads_);
        two_grid_solver_detail::multiply_parallel(
            p_transpose_, residual, workspace.coarse_rhs,
            application_threads_);
    }
    coarse_solver_.solve(
        workspace.coarse_rhs, workspace.coarse_solution,
        workspace.coarse_work);
    two_grid_solver_detail::multiply_add_parallel(
        p_, 1.0, workspace.coarse_solution, solution,
        application_threads_);

    for (int step = 0; step < smoothing_steps_; ++step) {
        solve_gauss_seidel_sweep(rhs, false, solution);
    }
    return a_.residual_squared(
        solution, rhs, residual, application_threads_);
}

template <class Cycle>
inline TwoGridIterationResult solve_stationary_cycles(
    const Vector& rhs, const Cycle& cycle, double relative_tolerance,
    int max_cycles, const char* method_name) {
    if (!(relative_tolerance > 0.0) ||
        !(relative_tolerance < 1.0) || max_cycles <= 0) {
        throw std::invalid_argument(
            std::string("invalid ") + method_name + " stopping criteria");
    }
    constexpr int tail_window = 32;
    TwoGridIterationResult result;
    result.solution.assign(rhs.size(), 0.0);
    Vector residual = rhs;
    typename Cycle::Workspace workspace;
    std::array<double, static_cast<std::size_t>(tail_window + 1)>
        recent_residuals{};
    const double initial_norm = norm2(residual);
    if (initial_norm == 0.0) {
        result.converged = true;
        result.relative_residual = 0.0;
        result.best_relative_residual = 0.0;
        result.effective_factor = 0.0;
        result.tail_factor = 0.0;
        result.status = TwoGridIterationStatus::Converged;
        return result;
    }
    recent_residuals[0] = 1.0;
    for (int iteration = 0; iteration < max_cycles; ++iteration) {
        const double residual_squared =
            cycle.iterate(rhs, result.solution, residual, workspace);
        result.cycles = iteration + 1;
        if (!(residual_squared >= 0.0) ||
            !std::isfinite(residual_squared)) {
            result.relative_residual =
                std::numeric_limits<double>::infinity();
            result.status = TwoGridIterationStatus::Diverged;
            break;
        }
        result.relative_residual =
            std::sqrt(residual_squared) / initial_norm;
        result.best_relative_residual = std::min(
            result.best_relative_residual, result.relative_residual);
        recent_residuals[static_cast<std::size_t>(
            result.cycles % (tail_window + 1))] =
            result.relative_residual;
        if (result.relative_residual <= relative_tolerance) {
            result.converged = true;
            result.status = TwoGridIterationStatus::Converged;
            break;
        }
        if (result.cycles >= tail_window &&
            result.relative_residual > 1.0e8 &&
            result.relative_residual >
                10.0 * result.best_relative_residual) {
            result.status = TwoGridIterationStatus::Diverged;
            break;
        }
    }
    if (result.cycles > 0 && std::isfinite(result.relative_residual)) {
        if (result.relative_residual > 0.0) {
            result.effective_factor = std::exp(
                std::log(result.relative_residual) /
                static_cast<double>(result.cycles));
        } else {
            result.effective_factor = 0.0;
        }
        const int span = std::min(result.cycles, tail_window);
        const double anchor = recent_residuals[static_cast<std::size_t>(
            (result.cycles - span) % (tail_window + 1))];
        if (anchor > 0.0 && result.relative_residual > 0.0) {
            result.tail_factor = std::exp(
                std::log(result.relative_residual / anchor) /
                static_cast<double>(span));
        }
    }
    if (!result.converged &&
        result.status != TwoGridIterationStatus::Diverged) {
        if (result.tail_factor > 1.001 &&
            result.relative_residual >= 1.0 &&
            result.relative_residual >
                10.0 * result.best_relative_residual) {
            result.status = TwoGridIterationStatus::Diverged;
        } else {
            result.status = TwoGridIterationStatus::SlowAtLimit;
        }
    }
    return result;
}

inline TwoGridIterationResult solve_two_grid(
    const Vector& rhs, const TwoGridCycle& cycle,
    double relative_tolerance, int max_cycles) {
    return solve_stationary_cycles(
        rhs, cycle, relative_tolerance, max_cycles, "two-grid");
}

}
