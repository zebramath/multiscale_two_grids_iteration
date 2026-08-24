#include "experiment/reporting.hpp"
#include "experiment/studies.hpp"

int main(int argc, char** argv) {
    experiment_support::begin_result_group(
        "experiment3_oracle_validation");
    return run_oracle_quality(argc, argv);
}
