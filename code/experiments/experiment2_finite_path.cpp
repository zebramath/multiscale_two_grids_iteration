#include "experiment/reporting.hpp"
#include "experiment/studies.hpp"

int main(int argc, char** argv) {
    experiment_support::begin_result_group("experiment2_finite_path");
    return run_finite_path_evidence(argc, argv);
}
