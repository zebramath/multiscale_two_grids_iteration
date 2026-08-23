#include "experiment/reporting.hpp"
#include "experiment/studies.hpp"

int main(int argc, char** argv) {
    experiment_support::begin_result_group("experiment3_adaptive_oracle");
    if (run_adaptive_trace(argc, argv) != 0) return 1;
    return run_oracle_quality(argc, argv);
}
