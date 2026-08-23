#include "experiment/reporting.hpp"
#include "experiment/studies.hpp"

int main(int argc, char** argv) {
    experiment_support::begin_result_group("experiment2_pcg_path");
    if (run_pcg_budget(argc, argv) != 0) return 1;
    if (run_pcg_dense_scan(argc, argv) != 0) return 1;
    return run_pilot_diagnostic(argc, argv);
}
