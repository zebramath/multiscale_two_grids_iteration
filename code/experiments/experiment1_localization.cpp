#include "experiment/reporting.hpp"
#include "experiment/studies.hpp"

int main(int argc, char** argv) {
    experiment_support::begin_result_group("experiment1_localization");
    if (run_localization_radius(argc, argv) != 0) return 1;
    if (run_local_tolerance(argc, argv) != 0) return 1;
    if (run_global_pruning(argc, argv) != 0) return 1;
    return run_support_strategy(argc, argv);
}
