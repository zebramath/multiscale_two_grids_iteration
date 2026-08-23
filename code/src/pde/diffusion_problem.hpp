#pragma once

#include "core/linear_algebra.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

namespace tgi {

class StructuredGrid {
public:
    StructuredGrid(int fine_interior_points, int coarsening_ratio);

    int fine_n() const { return fine_n_; }
    int ratio() const { return ratio_; }
    int intervals() const { return fine_n_ + 1; }
    int fine_size() const { return fine_n_ * fine_n_; }
    int coarse_n() const { return coarse_n_; }
    int coarse_size() const { return coarse_n_ * coarse_n_; }
    double h() const { return 1.0 / static_cast<double>(intervals()); }

    int fine_id(int ix, int iy) const { return iy * fine_n_ + ix; }
    std::pair<int, int> fine_coords(int id) const;
    int coarse_id(int cx, int cy) const { return cy * coarse_n_ + cx; }
    std::pair<int, int> coarse_coords(int id) const;
    int coarse_fine_id(int coarse_id) const;
    bool is_coarse_node(int fine_id) const;
    std::vector<int> all_f_nodes() const;
    std::vector<int> patch_f_nodes(int coarse_id, int patch_layers) const;
    void patch_f_nodes(int coarse_id, int patch_layers,
                       std::vector<int>& nodes) const;

private:
    int fine_n_;
    int ratio_;
    int coarse_n_;
};

enum class CoefficientDistribution {
    Analytic,
    RandomContinuous,
    RandomBinaryCheckerboard,
    ChannelizedBinary
};

struct CoefficientOptions {
    CoefficientDistribution distribution = CoefficientDistribution::Analytic;
    double contrast = 1e4;
    std::uint64_t seed = 1;
    int random_modes = 32;
    int minimum_frequency = 2;
    int maximum_frequency = 12;
    double spectral_decay = 1.1;
    int checkerboard_block_size = 8;
    int channel_background_block_size = 8;
    int channel_width_fine_cells = 2;
};

struct CoefficientField {
    Vector values;
    double minimum = 0.0;
    double maximum = 0.0;
    double actual_contrast = 0.0;
};

inline CoefficientField make_coefficient(const StructuredGrid& grid,
                                  const CoefficientOptions& options);
inline SparseMatrix assemble_diffusion(const StructuredGrid& grid,
                                const Vector& coefficient);


inline StructuredGrid::StructuredGrid(int fine_interior_points, int coarsening_ratio)
    : fine_n_(fine_interior_points), ratio_(coarsening_ratio),
      coarse_n_(0) {
    if (fine_n_ <= 0 || ratio_ <= 0 ||
        (fine_n_ + 1) % ratio_ != 0 ||
        (fine_n_ + 1) / ratio_ <= 1) {
        throw std::invalid_argument(
            "StructuredGrid: incompatible fine size and coarsening ratio");
    }
    coarse_n_ = (fine_n_ + 1) / ratio_ - 1;
}

inline std::pair<int, int> StructuredGrid::fine_coords(int id) const {
    return {id % fine_n_, id / fine_n_};
}

inline std::pair<int, int> StructuredGrid::coarse_coords(int id) const {
    return {id % coarse_n_, id / coarse_n_};
}

inline int StructuredGrid::coarse_fine_id(int id) const {
    const auto [cx, cy] = coarse_coords(id);
    return fine_id((cx + 1) * ratio_ - 1, (cy + 1) * ratio_ - 1);
}

inline bool StructuredGrid::is_coarse_node(int id) const {
    const auto [ix, iy] = fine_coords(id);
    return ((ix + 1) % ratio_ == 0) && ((iy + 1) % ratio_ == 0);
}

inline std::vector<int> StructuredGrid::all_f_nodes() const {
    std::vector<int> nodes;
    nodes.reserve(static_cast<std::size_t>(fine_size() - coarse_size()));
    for (int fine = 0; fine < fine_size(); ++fine) {
        if (!is_coarse_node(fine)) nodes.push_back(fine);
    }
    return nodes;
}

inline std::vector<int> StructuredGrid::patch_f_nodes(int id, int patch_layers) const {
    std::vector<int> nodes;
    patch_f_nodes(id, patch_layers, nodes);
    return nodes;
}

inline void StructuredGrid::patch_f_nodes(int id, int patch_layers,
                                   std::vector<int>& nodes) const {
    const auto [center_x, center_y] = fine_coords(coarse_fine_id(id));
    const int radius = patch_layers * ratio_;
    const int xmin = std::max(0, center_x - radius);
    const int xmax = std::min(fine_n_ - 1, center_x + radius);
    const int ymin = std::max(0, center_y - radius);
    const int ymax = std::min(fine_n_ - 1, center_y + radius);
    nodes.clear();
    nodes.reserve(static_cast<std::size_t>((xmax - xmin + 1) * (ymax - ymin + 1)));
    for (int iy = ymin; iy <= ymax; ++iy) {
        for (int ix = xmin; ix <= xmax; ++ix) {
            const int fine = fine_id(ix, iy);
            if (!is_coarse_node(fine)) nodes.push_back(fine);
        }
    }
}

namespace diffusion_problem_detail {

constexpr double pi = 3.141592653589793238462643383279502884;

struct FourierMode {
    int kx;
    int ky;
    double amplitude;
    double phase;
};

inline std::vector<FourierMode> make_random_modes(
    int count, int minimum_frequency, int maximum_frequency,
    double spectral_decay, std::uint64_t seed) {
    std::mt19937_64 generator(seed);
    std::uniform_int_distribution<int> wave_number(
        minimum_frequency, maximum_frequency);
    std::uniform_real_distribution<double> phase(0.0, 2.0 * pi);
    std::vector<FourierMode> modes;
    modes.reserve(static_cast<std::size_t>(count));
    for (int index = 0; index < count; ++index) {
        const int kx = wave_number(generator);
        const int ky = wave_number(generator);
        const double wave_norm = std::sqrt(static_cast<double>(kx * kx + ky * ky));
        modes.push_back(
            {kx, ky, 1.0 / std::pow(wave_norm, spectral_decay),
             phase(generator)});
    }
    return modes;
}

inline double harmonic(double lhs, double rhs) {
    return 2.0 * lhs * rhs / (lhs + rhs);
}

inline std::uint64_t mix_bits(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
}

inline bool is_high_conductivity_channel(double x, double y,
                                  double width) {
    const double horizontal_center =
        9.0 / 32.0 + 0.007 * std::sin(4.0 * pi * x + 0.30);
    const double vertical_center =
        23.0 / 32.0 + 0.007 * std::sin(3.0 * pi * y + 0.80);

    const bool horizontal =
        std::abs(y - horizontal_center) <= 0.5 * width;
    const bool vertical =
        std::abs(x - vertical_center) <= 0.5 * width;
    return horizontal || vertical;
}

}

inline CoefficientField make_coefficient(const StructuredGrid& grid,
                                  const CoefficientOptions& options) {
    if (options.distribution ==
        CoefficientDistribution::ChannelizedBinary) {
        if (options.channel_width_fine_cells <= 0) {
            throw std::invalid_argument(
                "channel width must be positive");
        }
        if (options.channel_background_block_size <= 0) {
            throw std::invalid_argument(
                "channelized inclusion block size must be positive");
        }
        CoefficientField field;
        field.values.resize(static_cast<std::size_t>(grid.fine_size()));
        const double width =
            static_cast<double>(options.channel_width_fine_cells) *
            grid.h();
        for (int id = 0; id < grid.fine_size(); ++id) {
            const auto [ix, iy] = grid.fine_coords(id);
            const double x = static_cast<double>(ix + 1) * grid.h();
            const double y = static_cast<double>(iy + 1) * grid.h();
            const int block_x =
                (ix + 1) / options.channel_background_block_size;
            const int block_y =
                (iy + 1) / options.channel_background_block_size;
            const std::uint64_t hash = diffusion_problem_detail::mix_bits(
                options.seed ^
                (static_cast<std::uint64_t>(
                     static_cast<std::uint32_t>(block_x)) << 32U) ^
                static_cast<std::uint64_t>(
                    static_cast<std::uint32_t>(block_y)));
            const bool high_inclusion = (hash & 1ULL) != 0ULL;
            field.values[static_cast<std::size_t>(id)] =
                (diffusion_problem_detail::is_high_conductivity_channel(x, y, width) ||
                 high_inclusion)
                    ? options.contrast
                    : 1.0;
        }
        const auto [field_min, field_max] =
            std::minmax_element(field.values.begin(), field.values.end());
        field.minimum = *field_min;
        field.maximum = *field_max;
        field.actual_contrast = field.maximum / field.minimum;
        return field;
    }

    if (options.distribution ==
        CoefficientDistribution::RandomBinaryCheckerboard) {
        if (options.checkerboard_block_size <= 0) {
            throw std::invalid_argument(
                "checkerboard block size must be positive");
        }
        CoefficientField field;
        field.values.resize(static_cast<std::size_t>(grid.fine_size()));
        for (int id = 0; id < grid.fine_size(); ++id) {
            const auto [ix, iy] = grid.fine_coords(id);
            const int block_x =
                (ix + 1) / options.checkerboard_block_size;
            const int block_y =
                (iy + 1) / options.checkerboard_block_size;
            const std::uint64_t hash = diffusion_problem_detail::mix_bits(
                options.seed ^
                (static_cast<std::uint64_t>(
                     static_cast<std::uint32_t>(block_x)) << 32U) ^
                static_cast<std::uint64_t>(
                    static_cast<std::uint32_t>(block_y)));
            field.values[static_cast<std::size_t>(id)] =
                (hash & 1ULL) == 0ULL ? 1.0 : options.contrast;
        }
        const auto [field_min, field_max] =
            std::minmax_element(field.values.begin(), field.values.end());
        field.minimum = *field_min;
        field.maximum = *field_max;
        field.actual_contrast = field.maximum / field.minimum;
        return field;
    }

    Vector raw(static_cast<std::size_t>(grid.fine_size()), 0.0);
    const std::vector<diffusion_problem_detail::FourierMode> modes =
        options.distribution == CoefficientDistribution::RandomContinuous
            ? diffusion_problem_detail::make_random_modes(
                  options.random_modes, options.minimum_frequency,
                  options.maximum_frequency, options.spectral_decay,
                  options.seed)
            : std::vector<diffusion_problem_detail::FourierMode>{};

    for (int id = 0; id < grid.fine_size(); ++id) {
        const auto [ix, iy] = grid.fine_coords(id);
        const double x = static_cast<double>(ix + 1) * grid.h();
        const double y = static_cast<double>(iy + 1) * grid.h();
        double value = 0.0;
        if (options.distribution == CoefficientDistribution::Analytic) {
            value = std::sin(2.0 * diffusion_problem_detail::pi * x) *
                    std::sin(2.0 * diffusion_problem_detail::pi * y);
        } else {
            for (const auto& mode : modes) {
                value += mode.amplitude *
                    std::cos(2.0 * diffusion_problem_detail::pi *
                             (static_cast<double>(mode.kx) * x +
                              static_cast<double>(mode.ky) * y) +
                             mode.phase);
            }
        }
        raw[static_cast<std::size_t>(id)] = value;
    }

    const auto [minimum_it, maximum_it] =
        std::minmax_element(raw.begin(), raw.end());
    const double raw_range = *maximum_it - *minimum_it;

    CoefficientField field;
    field.values.resize(raw.size());
    const double log_contrast = std::log(options.contrast);
    for (std::size_t i = 0; i < raw.size(); ++i) {
        const double normalized = (raw[i] - *minimum_it) / raw_range;
        field.values[i] = std::exp(log_contrast * normalized);
    }
    const auto [field_min, field_max] =
        std::minmax_element(field.values.begin(), field.values.end());
    field.minimum = *field_min;
    field.maximum = *field_max;
    field.actual_contrast = field.maximum / field.minimum;
    return field;
}

inline SparseMatrix assemble_diffusion(const StructuredGrid& grid,
                                const Vector& coefficient) {
    if (coefficient.size() !=
        static_cast<std::size_t>(grid.fine_size())) {
        throw std::invalid_argument(
            "assemble_diffusion: coefficient size mismatch");
    }
    const double inverse_h2 = 1.0 / (grid.h() * grid.h());
    std::vector<int> row_ptr(
        static_cast<std::size_t>(grid.fine_size()) + 1U, 0);
    std::vector<int> col_idx;
    Vector values;
    col_idx.reserve(static_cast<std::size_t>(5 * grid.fine_size()));
    values.reserve(static_cast<std::size_t>(5 * grid.fine_size()));

    for (int id = 0; id < grid.fine_size(); ++id) {
        const auto [ix, iy] = grid.fine_coords(id);
        const double center =
            coefficient[static_cast<std::size_t>(id)];
        double diagonal = 0.0;

        const auto face_value = [&](int neighbor) {
            return diffusion_problem_detail::harmonic(
                center, coefficient[static_cast<std::size_t>(neighbor)]);
        };
        const auto add_neighbor = [&](int neighbor) {
            const double face = face_value(neighbor);
            col_idx.push_back(neighbor);
            values.push_back(-face * inverse_h2);
            diagonal += face * inverse_h2;
        };

        // Row-major numbering gives the sorted CSR order
        // up, left, diagonal, right, down.
        if (iy > 0) add_neighbor(id - grid.fine_n());
        else diagonal += center * inverse_h2;
        if (ix > 0) add_neighbor(id - 1);
        else diagonal += center * inverse_h2;

        col_idx.push_back(id);
        values.push_back(0.0);
        const std::size_t diagonal_position = values.size() - 1U;

        if (ix + 1 < grid.fine_n()) add_neighbor(id + 1);
        else diagonal += center * inverse_h2;
        if (iy + 1 < grid.fine_n()) add_neighbor(id + grid.fine_n());
        else diagonal += center * inverse_h2;

        values[diagonal_position] = diagonal;
        row_ptr[static_cast<std::size_t>(id) + 1U] =
            static_cast<int>(values.size());
    }
    return SparseMatrix(
        grid.fine_size(), grid.fine_size(), std::move(row_ptr),
        std::move(col_idx), std::move(values));
}

}
