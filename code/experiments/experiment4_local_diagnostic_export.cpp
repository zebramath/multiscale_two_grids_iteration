#include "experiment/problem.hpp"
#include "experiment/reporting.hpp"
#include "multigrid/global_pcg.hpp"

#include <cstddef>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

namespace {

void write_sparse(
    std::ofstream& stream, const std::string& kind, int step,
    const tgi::SparseMatrix& matrix) {
    for (int row = 0; row < matrix.rows(); ++row) {
        for (int position =
                 matrix.row_ptr()[static_cast<std::size_t>(row)];
             position <
                 matrix.row_ptr()[static_cast<std::size_t>(row) + 1U];
             ++position) {
            stream << kind << ',' << step << ',' << row << ','
                   << matrix.col_idx()[static_cast<std::size_t>(position)]
                   << ','
                   << matrix.values()[static_cast<std::size_t>(position)]
                   << '\n';
        }
    }
}

}

int main(int argc, char** argv) {
    int threads = 4;
    int maximum_steps = 80;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument.rfind("--threads=", 0) == 0) {
            threads = std::stoi(argument.substr(10));
        }
        if (argument.rfind("--maximum-steps=", 0) == 0) {
            maximum_steps = std::stoi(argument.substr(16));
        }
    }

    experiment_support::BasicConfig config;
    config.fine_intervals = 16;
    config.coarse_intervals = 4;
    config.contrast = 1.0e4;
    config.threads = threads;
    const auto& field = experiment_support::channel_topologies().front();
    const auto grid = experiment_support::make_grid(config);
    const auto problem = experiment_support::make_problem(
        grid, field, config, 1);
    const auto geometric = tgi::build_geometric_interpolation(grid);

    const auto output_path = experiment_support::results_directory() /
        "experiment4_local_diagnostic_matrices.csv";
    std::ofstream stream(output_path);
    if (!stream) {
        throw std::runtime_error(
            "local diagnostic export open failure: " + output_path.string());
    }
    stream << std::setprecision(17);
    stream << "kind,m,row,col,value\n";
    write_sparse(stream, "A", -2, problem.matrix);

    tgi::GlobalEnergyPcgPath path(
        grid, problem.matrix, geometric.prolongation, threads);
    write_sparse(stream, "P", 0, geometric.prolongation);
    for (int step = 1; step <= maximum_steps; ++step) {
        experiment_support::progress(
            "local diagnostic export " + std::to_string(step) + "/" +
            std::to_string(maximum_steps));
        path.advance_to(step);
        write_sparse(stream, "P", step, path.prolongation());
    }
    path.advance_until_relative_residual(1.0e-12, 40000);
    write_sparse(stream, "P_endpoint", -1, path.prolongation());
    if (!stream) {
        throw std::runtime_error(
            "failed while writing local diagnostic matrices");
    }

    return 0;
}
