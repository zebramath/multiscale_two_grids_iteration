#pragma once

#include "pde/diffusion_problem.hpp"

#include <cmath>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>
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

inline tgi::Vector random_spectral_solution(
    const tgi::StructuredGrid& grid, std::uint64_t seed,
    int mode_count = 12, int maximum_frequency = 12) {
    constexpr double pi = 3.141592653589793238462643383279502884;
    std::mt19937_64 generator(seed);
    std::uniform_int_distribution<int> frequency(1, maximum_frequency);
    std::normal_distribution<double> amplitude(0.0, 1.0);
    struct Mode {
        int kx;
        int ky;
        double value;
    };
    std::vector<Mode> modes;
    modes.reserve(static_cast<std::size_t>(mode_count));
    for (int index = 0; index < mode_count; ++index) {
        const int kx = frequency(generator);
        const int ky = frequency(generator);
        const double decay = std::sqrt(
            static_cast<double>(kx * kx + ky * ky));
        modes.push_back({kx, ky, amplitude(generator) / decay});
    }
    tgi::Vector exact(
        static_cast<std::size_t>(grid.fine_size()), 0.0);
    for (int id = 0; id < grid.fine_size(); ++id) {
        const auto [ix, iy] = grid.fine_coords(id);
        const double x = static_cast<double>(ix + 1) * grid.h();
        const double y = static_cast<double>(iy + 1) * grid.h();
        double value = 0.0;
        for (const Mode& mode : modes) {
            value += mode.value *
                std::sin(static_cast<double>(mode.kx) * pi * x) *
                std::sin(static_cast<double>(mode.ky) * pi * y);
        }
        exact[static_cast<std::size_t>(id)] = value;
    }
    const double norm = tgi::norm2(exact);
    if (norm > 0.0) tgi::scale(1.0 / norm, exact);
    return exact;
}

inline tgi::Vector point_source_rhs(
    const tgi::StructuredGrid& grid, int ix, int iy) {
    if (ix < 0 || ix >= grid.fine_n() ||
        iy < 0 || iy >= grid.fine_n()) {
        throw std::out_of_range("point source is outside the fine grid");
    }
    tgi::Vector rhs(
        static_cast<std::size_t>(grid.fine_size()), 0.0);
    rhs[static_cast<std::size_t>(grid.fine_id(ix, iy))] = 1.0;
    return rhs;
}

} // namespace experiment_support
