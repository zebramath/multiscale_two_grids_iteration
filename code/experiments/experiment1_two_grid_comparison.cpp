#include "experiment/reporting.hpp"
#include "experiment/studies.hpp"

int main(int argc, char** argv) {
    experiment_support::begin_result_group(
        "experiment1_two_grid_comparison");
    return run_two_grid_comparison(argc, argv);
}
