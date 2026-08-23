#include "experiment/reporting.hpp"
#include "experiment/studies.hpp"

int main(int argc, char** argv) {
    experiment_support::begin_result_group("experiment5_diagnostics");
    if (run_diagnostic_cases(argc, argv) != 0) return 1;
    return run_stress_cases(argc, argv);
}
