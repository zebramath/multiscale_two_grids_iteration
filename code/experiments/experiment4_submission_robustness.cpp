#include "experiment/reporting.hpp"
#include "experiment/studies.hpp"

int main(int argc, char** argv) {
    experiment_support::begin_result_group(
        "experiment4_submission_robustness");
    return run_submission_robustness(argc, argv);
}
