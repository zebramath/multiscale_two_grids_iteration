#include "experiment/reporting.hpp"
#include "experiment/studies.hpp"

int main(int argc, char** argv) {
    experiment_support::begin_result_group(
        "experiment4_multilevel_comparison");
    return run_multilevel_comparison(argc, argv);
}
