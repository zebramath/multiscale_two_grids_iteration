#include "experiment/reporting.hpp"
#include "experiment/studies.hpp"

int main(int argc, char** argv) {
    experiment_support::begin_result_group("experiment6_workload");
    return run_workload(argc, argv);
}
