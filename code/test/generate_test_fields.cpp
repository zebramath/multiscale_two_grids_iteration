#include "experiment/reporting.hpp"
#include "experiment/test_field_dataset.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>

namespace {

struct Config {
    int fine_intervals = 128;
    int coarse_intervals = 16;
    double contrast = 1.0e4;
    std::uint64_t seed = 2;
    std::filesystem::path output = "models/test_fields.tgi";
};

std::size_t high_phase_nodes(const tgi::CoefficientField& field) {
    const double threshold =
        std::sqrt(field.minimum * field.maximum);
    return static_cast<std::size_t>(std::count_if(
        field.values.begin(), field.values.end(),
        [threshold](double value) { return value > threshold; }));
}

} // namespace

int main(int argc, char** argv) {
    Config cfg;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument.rfind("--output=", 0) == 0) {
            cfg.output = argument.substr(9);
        } else if (argument.rfind("--seed=", 0) == 0) {
            cfg.seed = static_cast<std::uint64_t>(
                std::stoull(argument.substr(7)));
        } else {
            throw std::invalid_argument(
                "usage: generate_test_fields [--output=PATH] [--seed=N]");
        }
    }

    const int ratio = cfg.fine_intervals / cfg.coarse_intervals;
    const tgi::StructuredGrid grid(cfg.fine_intervals - 1, ratio);
    experiment_support::TestFieldDataset dataset;
    dataset.fine_intervals = cfg.fine_intervals;
    dataset.coarse_intervals = cfg.coarse_intervals;

    tgi::CoefficientOptions random;
    random.distribution =
        tgi::CoefficientDistribution::RandomContinuous;
    random.contrast = cfg.contrast;
    random.seed = cfg.seed;
    random.random_modes = 32;
    random.minimum_frequency = 2;
    random.maximum_frequency = 12;
    random.spectral_decay = 1.1;
    dataset.fields.push_back({
        "random_continuous", tgi::make_coefficient(grid, random)});

    tgi::CoefficientOptions channel;
    channel.distribution =
        tgi::CoefficientDistribution::ChannelizedBinary;
    channel.contrast = cfg.contrast;
    channel.seed = cfg.seed;
    channel.channel_background_block_size = ratio;
    channel.channel_width_fine_cells = 2;
    dataset.fields.push_back({
        "channelized_binary", tgi::make_coefficient(grid, channel)});

    tgi::CoefficientOptions checkerboard;
    checkerboard.distribution =
        tgi::CoefficientDistribution::RandomBinaryCheckerboard;
    checkerboard.contrast = cfg.contrast;
    checkerboard.seed = cfg.seed;
    checkerboard.checkerboard_block_size = ratio;
    dataset.fields.push_back({
        "coarse_checkerboard", tgi::make_coefficient(grid, checkerboard)});

    experiment_support::save_test_fields(cfg.output, dataset);

    experiment_support::Rows rows;
    for (const auto& field : dataset.fields) {
        rows.push_back({
            field.name,
            experiment_support::integer(field.coefficient.actual_contrast),
            experiment_support::fixed(
                100.0 * static_cast<double>(high_phase_nodes(
                    field.coefficient)) /
                    static_cast<double>(field.coefficient.values.size()),
                2) + "%"
        });
    }
    experiment_support::Report report(
        "Shared high-contrast test-field dataset");
    report.add_summary({
        {"Dataset", cfg.output.string()},
        {"Fine grid", "h = 1/" + std::to_string(cfg.fine_intervals)},
        {"Coarse grid", "H = 1/" + std::to_string(cfg.coarse_intervals)},
        {"H/h", std::to_string(ratio)},
        {"Fields", std::to_string(dataset.fields.size())},
        {"Seed", std::to_string(cfg.seed)}
    });
    report.add_table(
        "Coefficient fields",
        {"Name", "Contrast", "High phase"},
        {22, 10, 12}, rows);
    report.save("generate_test_fields");
    return 0;
}
