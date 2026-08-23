#pragma once

#include "pde/diffusion_problem.hpp"

#include <cmath>
namespace experiment_support {

inline tgi::Vector manufactured_solution(
    const tgi::StructuredGrid& grid) {
    constexpr double pi = 3.141592653589793238462643383279502884;
    tgi::Vector exact(static_cast<std::size_t>(grid.fine_size()));
    for (int id = 0; id < grid.fine_size(); ++id) {
        const auto [ix, iy] = grid.fine_coords(id);
        const double x = static_cast<double>(ix + 1) * grid.h();
        const double y = static_cast<double>(iy + 1) * grid.h();
        exact[static_cast<std::size_t>(id)] =
            std::sin(pi * x) * std::sin(pi * y) +
            0.23 * std::sin(3.0 * pi * x) *
                std::sin(2.0 * pi * y) +
            0.11 * std::sin(7.0 * pi * x) *
                std::sin(5.0 * pi * y);
    }
    return exact;
}

} // namespace experiment_support
