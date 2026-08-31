#include "experiment/problem.hpp"
#include "experiment/reporting.hpp"
#include "multigrid/global_pcg.hpp"
#include "version.hpp"

#include <array>
#include <string>

namespace {

struct OracleCase {
    const char* split;
    int fine;
    int coarse;
    double contrast;
    experiment_support::FieldCase field;
};

struct CycleMeasurement {
    int cycles = 0;
    double effective_factor = 1.0;
    tgi::TwoGridIterationStatus status =
        tgi::TwoGridIterationStatus::SlowAtLimit;
    bool converged = false;
};

CycleMeasurement cycles_for(
    const tgi::SparseMatrix& a, const tgi::Vector& rhs,
    const tgi::SparseMatrix& p, int threads, int maximum) {
    const tgi::TwoGridCycle cycle(a, p, 1, threads);
    const auto solved = tgi::solve_two_grid(
        rhs, cycle, 1.0e-6, maximum);
    return {solved.cycles, solved.effective_factor,
            solved.status, solved.converged};
}

CycleMeasurement cycles_for(
    const tgi::Vector& rhs, const tgi::TwoGridCycle& cycle, int maximum) {
    const auto solved = tgi::solve_two_grid(
        rhs, cycle, 1.0e-6, maximum);
    return {solved.cycles, solved.effective_factor,
            solved.status, solved.converged};
}

}

int main(int argc, char** argv) {
    int threads = 4;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument.rfind("--threads=", 0) == 0)
            threads = std::stoi(argument.substr(10));
    }
    const auto& topology = experiment_support::channel_topologies();
    const std::array<OracleCase, 10> cases{{
        {"design", 32, 8, 1.0e4, topology[0]},
        {"design", 64, 8, 1.0e4, topology[0]},
        {"design", 64, 16, 1.0e4, topology[0]},
        {"design", 128, 16, 1.0e4, topology[0]},
        {"design", 128, 16, 1.0e2, topology[0]},
        {"design", 128, 16, 1.0e6, topology[0]},
        {"design", 128, 16, 1.0e4, topology[5]},
        {"validation", 88, 11, 1.0e4, topology[2]},
        {"validation", 120, 15, 1.0e2, topology[4]},
        {"validation", 152, 19, 1.0e6, topology[3]}
    }};
    constexpr int maximum_cycles =
        experiment_support::maximum_two_grid_cycles;
    constexpr int oracle_candidate_limit = 6000;
    const experiment_support::Row headers{
        "Split", "1/h", "1/H", "Contrast", "Topology",
        "Selected m", "m/(1/h)", "Cycles", "Eff factor", "Status",
        "Oracle m", "Oracle m/(1/h)", "Oracle cycles", "Oracle factor",
        "Oracle status", "Gap %"};
    experiment_support::Rows rows;

    int case_index = 0;
    for (const OracleCase& item : cases) {
        ++case_index;
        experiment_support::progress(
            "oracle case " + std::to_string(case_index) + "/" +
            std::to_string(cases.size()) + ": " + item.field.name);
        experiment_support::BasicConfig config;
        config.fine_intervals = item.fine;
        config.coarse_intervals = item.coarse;
        config.contrast = item.contrast;
        config.threads = threads;
        const tgi::StructuredGrid grid = experiment_support::make_grid(config);
        const auto problem = experiment_support::make_problem(
            grid, item.field, config);
        const auto geometric = tgi::build_geometric_interpolation(grid);

        int oracle_steps = 0;
        CycleMeasurement oracle = cycles_for(
            problem.matrix, problem.rhs, geometric.prolongation,
            threads, maximum_cycles);
        tgi::GlobalEnergyPcgPath path(
            grid, problem.matrix, geometric.prolongation, threads);
        int first_oracle_step = std::max(2, item.fine / 8);
        if (first_oracle_step % 2 != 0) ++first_oracle_step;
        for (int steps = first_oracle_step;
             steps <= item.fine / 2; steps += 2) {
            path.advance_to(steps);
            const tgi::SparseMatrix candidate = path.prolongation();
            const CycleMeasurement measured = cycles_for(
                problem.matrix, problem.rhs, candidate,
                threads, oracle_candidate_limit);
            if (measured.converged &&
                (!oracle.converged || measured.cycles < oracle.cycles)) {
                oracle = measured;
                oracle_steps = steps;
            }
        }

        const auto adaptive = tgi::build_adaptive_global_pcg_interpolation(
            grid, problem.matrix, geometric.prolongation, threads);
        const CycleMeasurement adaptive_cycles = cycles_for(
            problem.rhs, *adaptive.cycle, maximum_cycles);
        const double gap = adaptive_cycles.converged && oracle.converged
            ? 100.0 * static_cast<double>(
                  adaptive_cycles.cycles - oracle.cycles) /
                  static_cast<double>(oracle.cycles)
            : 0.0;
        rows.push_back({
            item.split, std::to_string(item.fine),
            std::to_string(item.coarse),
            experiment_support::scientific(item.contrast, 0),
            item.field.name,
            std::to_string(adaptive.report.selected_steps),
            experiment_support::fixed(
                static_cast<double>(adaptive.report.selected_steps) /
                    static_cast<double>(item.fine), 3),
            std::to_string(adaptive_cycles.cycles),
            experiment_support::fixed(
                adaptive_cycles.effective_factor, 6),
            tgi::two_grid_status_name(adaptive_cycles.status),
            std::to_string(oracle_steps),
            experiment_support::fixed(
                static_cast<double>(oracle_steps) /
                    static_cast<double>(item.fine), 3),
            std::to_string(oracle.cycles),
            experiment_support::fixed(oracle.effective_factor, 6),
            tgi::two_grid_status_name(oracle.status),
            adaptive_cycles.converged && oracle.converged
                ? experiment_support::fixed(gap, 2)
                : "--"});
    }

    experiment_support::Report report(
        "Adaptive PCG versus an offline step-two oracle");
    report.add_summary({
        {"Version", std::string(tgi::version)},
        {"Threads", std::to_string(threads)},
        {"Oracle candidates", "m=0 and (1/h)/8,...,(1/h)/2 by 2"},
        {"Solve tolerance", "1e-6"},
        {"Maximum cycles", std::to_string(maximum_cycles)}});
    report.add_note(
        "The step-two oracle is evaluation-only. Adaptive uses (1/h)/8 when "
        "1/H<=8; otherwise it maps the matrix diagonal ratio to (1/h)/4, "
        "(1/h)/3 or (1/h)/2. It uses matrix and grid information but no "
        "contrast or topology label. Oracle candidates are screened to 6000 cycles; a "
        "candidate still above tolerance at that point cannot beat any "
        "reported oracle minimum. The final three cases are a post-freeze "
        "validation split with unseen resolutions and topology/contrast "
        "combinations.");
    report.add_table(
        "Representative oracle gaps", headers,
        {10, 5, 5, 10, 20, 10, 9, 8, 11, 10, 9, 16, 14, 13, 13, 8}, rows,
        true);
    report.save("experiment3_oracle_validation");
    return 0;
}
