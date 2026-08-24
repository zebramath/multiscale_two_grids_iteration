#pragma once

#include "experiment/problem.hpp"
#include "experiment/reporting.hpp"
#include "version.hpp"

#include <string>

namespace experiment_support {

inline double interpolation_density_percent(
    const tgi::SparseMatrix& prolongation) {
    const double entries = static_cast<double>(prolongation.rows()) *
        static_cast<double>(prolongation.cols());
    return entries > 0.0
        ? 100.0 * static_cast<double>(prolongation.nnz()) / entries
        : 0.0;
}

}
