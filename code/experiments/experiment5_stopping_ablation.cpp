#include "experiment/problem.hpp"
#include "experiment/reporting.hpp"
#include "multigrid/global_pcg.hpp"
#include "version.hpp"

#include <array>
#include <map>
#include <memory>
#include <string>

namespace {

struct AblationCase {
    int fine;
    int coarse;
    double contrast;
    experiment_support::FieldCase field;
};

struct Measurement {
    std::string policy;
    std::string parameter;
    tgi::GlobalPcgPathReport path;
    double energy = 0.0;
    double density = 0.0;
    tgi::StationaryIterationResult solved;
};

Measurement measure(
    const std::string& policy, const std::string& parameter,
    const tgi::GlobalPcgPathReport& path,
    const tgi::SparseMatrix& prolongation, const tgi::TwoGridCycle& cycle,
    const tgi::Vector& rhs) {
    return {
        policy, parameter, path,
        cycle.setup_report().interpolation_energy,
        experiment_support::interpolation_density_percent(prolongation),
        tgi::solve_two_grid(
            rhs, cycle, 1.0e-6,
            experiment_support::maximum_two_grid_cycles)};
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

    const auto& topology = experiment_support::channel_topologies();
    const std::array<AblationCase, 6> all_cases{{
        {32, 8, 1.0e4, topology[0]},
        {64, 16, 1.0e4, topology[0]},
        {128, 16, 1.0e2, topology[0]},
        {128, 16, 1.0e4, topology[0]},
        {128, 16, 1.0e6, topology[0]},
        {128, 16, 1.0e4, topology[5]}
    }};
    const std::size_t case_count = quick ? 2U : all_cases.size();
    constexpr double residual_tolerance = 1.0e-2;
    experiment_support::Rows rows;
    std::map<std::string, int> converged_counts;
    std::map<std::string, long long> recorded_cycle_sums;
    std::map<std::string, long long> column_iteration_sums;
    int adaptive_cycle_wins_fixed = 0;
    int adaptive_cycle_wins_residual = 0;

    for (std::size_t case_index = 0; case_index < case_count; ++case_index) {
        const AblationCase& item = all_cases[case_index];
        experiment_support::progress(
            "stopping ablation " + std::to_string(case_index + 1U) + "/" +
            std::to_string(case_count) + ": " + item.field.name);
        experiment_support::BasicConfig config;
        config.fine_intervals = item.fine;
        config.coarse_intervals = item.coarse;
        config.contrast = item.contrast;
        config.threads = threads;
        const tgi::StructuredGrid grid = experiment_support::make_grid(config);
        const auto problem = experiment_support::make_problem(
            grid, item.field, config);
        const auto geometric = tgi::build_geometric_interpolation(grid);

        std::array<Measurement, 3> measurements;

        const auto adaptive = tgi::build_adaptive_global_pcg_interpolation(
            grid, problem.matrix, geometric.prolongation, threads);
        measurements[0] = measure(
            "adaptive",
            "m=" + std::to_string(adaptive.report.selected_steps),
            adaptive.report.path, *adaptive.prolongation, *adaptive.cycle,
            problem.rhs);

        const int fixed_steps =
            tgi::adaptive_global_pcg_detail::scaled_checkpoint(
                item.fine, 1, 4);
        tgi::GlobalEnergyPcgPath fixed_path(
            grid, problem.matrix, geometric.prolongation, threads);
        fixed_path.advance_to(fixed_steps);
        auto fixed_prolongation = std::make_shared<tgi::SparseMatrix>(
            fixed_path.prolongation());
        const tgi::TwoGridCycle fixed_cycle(
            problem.matrix, *fixed_prolongation, 1, threads);
        measurements[1] = measure(
            "fixed-step", "m=round(n/4)", fixed_path.report(),
            *fixed_prolongation, fixed_cycle, problem.rhs);

        tgi::GlobalEnergyPcgPath residual_path(
            grid, problem.matrix, geometric.prolongation, threads);
        const auto residual_report =
            residual_path.advance_until_relative_residual(
                residual_tolerance);
        auto residual_prolongation = std::make_shared<tgi::SparseMatrix>(
            residual_path.prolongation());
        const tgi::TwoGridCycle residual_cycle(
            problem.matrix, *residual_prolongation, 1, threads);
        measurements[2] = measure(
            "fixed-residual", "column relres<=1e-2", residual_report,
            *residual_prolongation, residual_cycle, problem.rhs);

        if (measurements[0].solved.cycles <=
            measurements[1].solved.cycles) {
            ++adaptive_cycle_wins_fixed;
        }
        if (measurements[0].solved.cycles <=
            measurements[2].solved.cycles) {
            ++adaptive_cycle_wins_residual;
        }

        for (const Measurement& value : measurements) {
            if (value.solved.converged) {
                ++converged_counts[value.policy];
            }
            recorded_cycle_sums[value.policy] += value.solved.cycles;
            column_iteration_sums[value.policy] +=
                value.path.total_iterations;
            rows.push_back({
                std::to_string(item.fine),
                std::to_string(item.coarse),
                experiment_support::scientific(item.contrast, 0),
                item.field.name,
                value.policy,
                value.parameter,
                std::to_string(value.path.minimum_iterations) + "/" +
                    std::to_string(value.path.maximum_iterations),
                std::to_string(value.path.total_iterations),
                experiment_support::scientific(
                    value.path.maximum_relative_residual, 2),
                experiment_support::scientific(value.energy, 4),
                experiment_support::fixed(value.density, 4),
                std::to_string(value.solved.cycles),
                tgi::stationary_status_name(value.solved.status),
                experiment_support::fixed(
                    value.solved.effective_factor, 6)});
        }
    }

    experiment_support::Rows summary_rows;
    for (const std::string& policy :
         {std::string("adaptive"), std::string("fixed-step"),
          std::string("fixed-residual")}) {
        summary_rows.push_back({
            policy,
            std::to_string(converged_counts[policy]) + "/" +
                std::to_string(case_count),
            std::to_string(recorded_cycle_sums[policy]),
            std::to_string(column_iteration_sums[policy])});
    }

    experiment_support::Report report(
        "Stopping-policy ablation for the adaptive selector");
    report.add_summary({
        {"Version", std::string(tgi::version)},
        {"Mode", quick ? "quick" : "full"},
        {"Cases", std::to_string(case_count)},
        {"Threads", std::to_string(threads)},
        {"Fixed-step baseline", "m=round((1/h)/4)"},
        {"Residual baseline", "per-column relative residual 1e-2"},
        {"Solve tolerance", "1e-6"},
        {"Adaptive <= fixed-step cycles",
         std::to_string(adaptive_cycle_wins_fixed) + "/" +
             std::to_string(case_count)},
        {"Adaptive <= residual cycles",
         std::to_string(adaptive_cycle_wins_residual) + "/" +
             std::to_string(case_count)}});
    report.add_note(
        "All three policies start from the same geometric interpolation and "
        "use the same global Jacobi-PCG column equations. Fixed-step applies "
        "one normalized checkpoint to every problem. Fixed-residual stops "
        "each column independently at the same relative residual. The "
        "six-case comparison tests scale, contrast and topology without "
        "retuning either baseline. Recorded cycle sums include cycles actually "
        "executed by slow-limit cases; no unconverged time is extrapolated. "
        "Column iteration totals are deterministic setup-work proxies, and "
        "wall time is not used in this ablation.");
    report.add_table(
        "Per-case stopping-policy comparison",
        {"1/h", "1/H", "Contrast", "Topology", "Policy", "Parameter",
         "Column m min/max", "Column iter sum", "Max column relres",
         "Energy", "P density %", "Cycles", "Status", "Eff factor"},
        {5, 5, 10, 20, 16, 22, 16, 15, 17, 13, 11, 8, 10, 11},
        rows, true);
    report.add_table(
        "Aggregate ablation summary",
        {"Policy", "Converged", "Recorded cycle sum",
         "Column iteration sum"},
        {16, 11, 19, 20}, summary_rows);
    report.save("experiment5_stopping_ablation");
    return 0;
}
