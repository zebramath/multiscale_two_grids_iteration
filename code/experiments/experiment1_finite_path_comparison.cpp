#include "experiment/comparison_cases.hpp"
#include "experiment/reporting.hpp"
#include "multigrid/global_pcg.hpp"
#include <array>
#include <chrono>
#include <cmath>
#include <map>
#include <memory>
#include <string>
namespace {
struct Aggregate {
    int cases = 0;
    int converged = 0;
    long long cycles = 0;
    double setup_ms = 0.0;
    double density = 0.0;
};
double elapsed_ms(const std::chrono::steady_clock::time_point& start) {
    return std::chrono::duration<double, std::milli>(
               std::chrono::steady_clock::now() - start)
        .count();
}
void append_measurement(
    const experiment_support::ComparisonCase& item,
    const std::string& method, const std::string& parameter,
    const tgi::SparseMatrix& matrix, const tgi::Vector& rhs,
    const tgi::SparseMatrix& prolongation, int threads,
    double setup_ms, experiment_support::Rows& rows,
    std::map<std::string, Aggregate>& aggregates) {
    const tgi::TwoGridCycle cycle(matrix, prolongation, 1, threads);
    const auto solved = tgi::solve_two_grid(
        rhs, cycle, 1.0e-6,
        experiment_support::maximum_two_grid_cycles);
    const double density =
        experiment_support::interpolation_density_percent(prolongation);
    rows.push_back({
        item.axis,
        std::to_string(item.fine),
        std::to_string(item.coarse),
        experiment_support::scientific(item.contrast, 0),
        item.field.name,
        method,
        parameter,
        experiment_support::fixed(density, 4),
        std::to_string(cycle.coarse_matrix().nnz()),
        experiment_support::fixed(setup_ms, 2),
        std::to_string(solved.cycles),
        tgi::stationary_status_name(solved.status),
        experiment_support::fixed(solved.effective_factor, 6),
        experiment_support::fixed(solved.tail_factor, 6)});
    Aggregate& aggregate = aggregates[method];
    ++aggregate.cases;
    aggregate.setup_ms += setup_ms;
    aggregate.density += density;
    if (solved.converged) {
        ++aggregate.converged;
        aggregate.cycles += solved.cycles;
    }
}
}
int main(int argc, char** argv) {
    int threads = 4;
    bool quick = false;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--quick") quick = true;
        if (argument.rfind("--threads=", 0) == 0) {
            threads = std::stoi(argument.substr(10));
        }
    }
    const auto cases = experiment_support::comparison_cases(quick);
    experiment_support::Rows rows;
    std::map<std::string, Aggregate> aggregates;
    const std::array<std::array<int, 2>, 3> fractions{{
        {{1, 4}}, {{1, 3}}, {{1, 2}}
    }};
    for (std::size_t case_index = 0; case_index < cases.size(); ++case_index) {
        const auto& item = cases[case_index];
        experiment_support::progress(
            "finite-path comparison " + std::to_string(case_index + 1U) +
            "/" + std::to_string(cases.size()) + ": " + item.field.name);
        const auto config = experiment_support::comparison_config(item, threads);
        const auto grid = experiment_support::make_grid(config);
        const auto problem = experiment_support::make_problem(
            grid, item.field, config);
        const auto initial = tgi::build_geometric_interpolation(grid);
        auto start = std::chrono::steady_clock::now();
        tgi::GlobalEnergyPcgPath path(
            grid, problem.matrix, initial, threads);
        double cumulative_pcg_ms = elapsed_ms(start);
        for (const auto& fraction : fractions) {
            const int steps = experiment_support::fixed_path_steps(
                item.fine, fraction[0], fraction[1]);
            start = std::chrono::steady_clock::now();
            path.advance_to(steps);
            cumulative_pcg_ms += elapsed_ms(start);
            start = std::chrono::steady_clock::now();
            const auto prolongation = path.prolongation();
            const double snapshot_ms = elapsed_ms(start);
            const std::string method = "finite-" +
                std::to_string(fraction[0]) + "/" +
                std::to_string(fraction[1]);
            append_measurement(
                item, method,
                "m=" + std::to_string(steps),
                problem.matrix, problem.rhs, prolongation, threads,
                cumulative_pcg_ms + snapshot_ms, rows, aggregates);
        }
        start = std::chrono::steady_clock::now();
        tgi::GlobalEnergyPcgPath endpoint_path(
            grid, problem.matrix, initial, threads);
        endpoint_path.advance_until_relative_residual(1.0e-10);
        const auto endpoint = endpoint_path.prolongation();
        const double endpoint_setup_ms = elapsed_ms(start);
        append_measurement(
            item, "energy-endpoint", "column relres<=1e-10",
            problem.matrix, problem.rhs, endpoint, threads,
            endpoint_setup_ms, rows, aggregates);
    }
    experiment_support::Rows aggregate_rows;
    for (const std::string& method : {
             std::string("finite-1/4"), std::string("finite-1/3"),
             std::string("finite-1/2"), std::string("energy-endpoint")}) {
        const Aggregate& value = aggregates[method];
        aggregate_rows.push_back({
            method,
            std::to_string(value.converged) + "/" +
                std::to_string(value.cases),
            std::to_string(value.cycles),
            experiment_support::fixed(
                value.density / static_cast<double>(value.cases), 4),
            experiment_support::fixed(value.setup_ms, 1)});
    }
    experiment_support::Report report(
        "Fixed O(1/h) checkpoints along the energy-minimization path");
    report.add_summary({
        {"Cases", std::to_string(cases.size())},
        {"Mode", quick ? "quick" : "full"},
        {"Threads", std::to_string(threads)},
        {"Solve tolerance", "1e-6"},
        {"Endpoint column tolerance", "1e-10"}});
    report.add_note(
        "The three finite checkpoints m=(1/h)/4, (1/h)/3 and (1/h)/2 "
        "are reported separately. The endpoint follows the same Jacobi-PCG "
        "column paths from geometric "
        "interpolation until every column reaches relative residual 1e-10.");
    report.add_table(
        "All two-grid cases",
        {"Axis", "1/h", "1/H", "Contrast", "Topology", "Path point",
         "Parameter", "P density %", "Ac nnz", "Setup ms", "Cycles",
         "Status", "Effective factor", "Tail factor"},
        {12, 5, 5, 10, 20, 17, 22, 11, 9, 10, 8, 10, 16, 11},
        rows, true);
    report.add_table(
        "Aggregate summary",
        {"Path point", "Converged", "Cycle sum", "Mean density %",
         "Setup sum ms"},
        {17, 10, 12, 14, 13}, aggregate_rows);
    report.save("experiment1_finite_path_comparison");
    return 0;
}
