#include "experiment/study.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using LongVector = std::vector<long double>;

LongVector dense_long_double_solve(
    const tgi::SparseMatrix& matrix, const tgi::Vector& rhs) {
    const int n = matrix.rows();
    std::vector<long double> lower(
        static_cast<std::size_t>(n) * static_cast<std::size_t>(n), 0.0L);
    const auto at = [n](int row, int col) {
        return static_cast<std::size_t>(row) *
                   static_cast<std::size_t>(n) +
               static_cast<std::size_t>(col);
    };
    for (int row = 0; row < n; ++row) {
        for (int position =
                 matrix.row_ptr()[static_cast<std::size_t>(row)];
             position <
                 matrix.row_ptr()[static_cast<std::size_t>(row) + 1U];
             ++position) {
            const int col =
                matrix.col_idx()[static_cast<std::size_t>(position)];
            if (col <= row) {
                lower[at(row, col)] = static_cast<long double>(
                    matrix.values()[static_cast<std::size_t>(position)]);
            }
        }
    }

    for (int row = 0; row < n; ++row) {
        for (int col = 0; col <= row; ++col) {
            long double value = lower[at(row, col)];
            for (int inner = 0; inner < col; ++inner) {
                value -= lower[at(row, inner)] * lower[at(col, inner)];
            }
            if (row == col) {
                if (!(value > 0.0L) || !std::isfinite(value)) {
                    throw std::runtime_error(
                        "independent dense Cholesky found a non-SPD coarse matrix");
                }
                lower[at(row, col)] = std::sqrt(value);
            } else {
                lower[at(row, col)] = value / lower[at(col, col)];
            }
        }
    }

    LongVector solution(static_cast<std::size_t>(n));
    for (int row = 0; row < n; ++row) {
        long double value = static_cast<long double>(
            rhs[static_cast<std::size_t>(row)]);
        for (int col = 0; col < row; ++col) {
            value -= lower[at(row, col)] *
                     solution[static_cast<std::size_t>(col)];
        }
        solution[static_cast<std::size_t>(row)] =
            value / lower[at(row, row)];
    }
    for (int row = n - 1; row >= 0; --row) {
        long double value = solution[static_cast<std::size_t>(row)];
        for (int col = row + 1; col < n; ++col) {
            value -= lower[at(col, row)] *
                     solution[static_cast<std::size_t>(col)];
        }
        solution[static_cast<std::size_t>(row)] =
            value / lower[at(row, row)];
    }
    return solution;
}

long double norm2_long(const LongVector& values) {
    long double squared = 0.0L;
    for (long double value : values) squared += value * value;
    return std::sqrt(squared);
}

void verify_coarse_solver(const tgi::TwoGridCycle& cycle) {
    const int n = cycle.coarse_size();
    tgi::Vector rhs(static_cast<std::size_t>(n));
    for (int index = 0; index < n; ++index) {
        const long double x = static_cast<long double>(index + 1);
        rhs[static_cast<std::size_t>(index)] = static_cast<double>(
            std::sin(0.37L * x) + 0.31L * std::cos(0.11L * x));
    }

    tgi::Vector sparse_solution;
    tgi::Vector work;
    cycle.solve_coarse_system(rhs, sparse_solution, work);
    const LongVector dense_solution =
        dense_long_double_solve(cycle.coarse_matrix(), rhs);

    const tgi::Vector product =
        cycle.coarse_matrix().multiply(sparse_solution);
    const double relative_residual =
        tgi::norm2(tgi::subtract(rhs, product)) / tgi::norm2(rhs);

    LongVector difference(static_cast<std::size_t>(n));
    for (int index = 0; index < n; ++index) {
        difference[static_cast<std::size_t>(index)] =
            static_cast<long double>(
                sparse_solution[static_cast<std::size_t>(index)]) -
            dense_solution[static_cast<std::size_t>(index)];
    }
    const long double dense_norm = std::max(
        norm2_long(dense_solution),
        std::numeric_limits<long double>::min());
    const long double relative_difference =
        norm2_long(difference) / dense_norm;

    if (!(relative_residual <= 1.0e-10) ||
        !(relative_difference <= 1.0e-8L)) {
        throw std::runtime_error(
            "coarse solve cross-check failed: sparse Cholesky and "
            "independent long-double dense Cholesky disagree");
    }
}

experiment_support::Row evaluate_verified_candidate(
    const std::string& field, const tgi::SparseMatrix& matrix,
    const tgi::Vector& rhs,
    const experiment_support::StudyCandidate& candidate,
    const experiment_support::BasicConfig& config) {
    const tgi::TwoGridCycle validation_cycle(
        matrix, candidate.prolongation, 1, config.threads);
    verify_coarse_solver(validation_cycle);
    return experiment_support::evaluate_candidate(
        field, matrix, rhs, candidate, config);
}

}

int run_pcg_dense_scan(int argc, char** argv) {
    const auto config = experiment_support::parse_config(argc, argv);
    const tgi::StructuredGrid grid = experiment_support::make_grid(config);
    const auto& channel = experiment_support::standard_fields().at(1);
    const auto problem = experiment_support::make_problem(
        grid, channel, config);
    const auto& a = problem.matrix;
    const auto& rhs = problem.rhs;
    experiment_support::Rows rows;

    auto geometric = experiment_support::geometric_interpolation(grid, a);
    const double geometric_ms = geometric.report.timing.total_ms;
    for (int steps = 16; steps <= 64; steps += 2) {
        auto options = experiment_support::energy_options(
            0, config.threads, 0.0);
        options.local_max_iterations = steps;
        options.require_convergence = false;
        options.drop_tolerance = 0.0;
        auto interpolation = tgi::refine_global_energy_interpolation(
            grid, a, geometric.prolongation, options);
        auto candidate = experiment_support::make_candidate(
            "PCG-global", "m=" + std::to_string(steps),
            std::move(interpolation));
        candidate.build_ms += geometric_ms;
        rows.push_back(evaluate_verified_candidate(
            channel.name, a, rhs, candidate, config));
    }

    auto global = experiment_support::build_global_reference(
        grid, a, config.threads);
    auto exact = experiment_support::make_candidate(
        "global-exact", "tol=1e-10", std::move(global));
    rows.push_back(evaluate_verified_candidate(
        channel.name, a, rhs, exact, config));

    experiment_support::Report report(
        "Channel dense PCG-step scan and coarse-solve cross-check");
    report.add_summary(experiment_support::fixed_study_summary(
        config, "PCG scan", "m=16,18,...,64"));
    report.add_note(
        "Every row passed an independent long-double dense-Cholesky "
        "cross-check of the Galerkin coarse solve. Validation work is "
        "excluded from Setup ms and Total ms, and the production sparse "
        "coarse solver is unchanged.");
    report.add_table(
        "Channel finite-step PCG scan", experiment_support::study_headers(),
        experiment_support::study_widths(), rows);
    report.save("pcg_dense_scan");
    return 0;
}
