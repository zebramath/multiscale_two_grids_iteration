#pragma once

#include "pde/diffusion_problem.hpp"

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace experiment_support {

struct NamedCoefficientField {
    std::string name;
    tgi::CoefficientField coefficient;
};

struct TestFieldDataset {
    int fine_intervals = 0;
    int coarse_intervals = 0;
    std::vector<NamedCoefficientField> fields;
};

namespace test_field_dataset_detail {

template <class T>
inline void write_scalar(std::ofstream& stream, const T& value) {
    stream.write(
        reinterpret_cast<const char*>(&value),
        static_cast<std::streamsize>(sizeof(T)));
}

template <class T>
inline T read_scalar(std::ifstream& stream) {
    T value{};
    stream.read(
        reinterpret_cast<char*>(&value),
        static_cast<std::streamsize>(sizeof(T)));
    return value;
}

inline void require_stream(
    const std::ios& stream, const std::filesystem::path& path) {
    if (!stream) {
        throw std::runtime_error(
            "invalid test-field dataset: " + path.string());
    }
}

inline void validate(const TestFieldDataset& dataset) {
    if (dataset.fine_intervals <= 1 ||
        dataset.coarse_intervals <= 0 ||
        dataset.fine_intervals % dataset.coarse_intervals != 0 ||
        dataset.fields.empty()) {
        throw std::invalid_argument("invalid test-field dataset metadata");
    }
    const std::size_t expected_values =
        static_cast<std::size_t>(dataset.fine_intervals - 1) *
        static_cast<std::size_t>(dataset.fine_intervals - 1);
    std::unordered_set<std::string> names;
    for (const NamedCoefficientField& field : dataset.fields) {
        if (field.name.empty() || !names.insert(field.name).second ||
            field.coefficient.values.size() != expected_values ||
            !(field.coefficient.minimum > 0.0) ||
            !(field.coefficient.maximum >= field.coefficient.minimum) ||
            !(field.coefficient.actual_contrast >= 1.0)) {
            throw std::invalid_argument("invalid named coefficient field");
        }
    }
}

} // namespace test_field_dataset_detail

inline void save_test_fields(
    const std::filesystem::path& path,
    const TestFieldDataset& dataset) {
    test_field_dataset_detail::validate(dataset);
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream stream(path, std::ios::binary);
    if (!stream) {
        throw std::runtime_error(
            "cannot create test-field dataset: " + path.string());
    }

    constexpr char magic[8] = {'T', 'G', 'I', 'F', 'L', 'D', 'S', '1'};
    stream.write(magic, 8);
    test_field_dataset_detail::write_scalar(stream, std::uint32_t{1});
    test_field_dataset_detail::write_scalar(
        stream, static_cast<std::int32_t>(dataset.fine_intervals));
    test_field_dataset_detail::write_scalar(
        stream, static_cast<std::int32_t>(dataset.coarse_intervals));
    test_field_dataset_detail::write_scalar(
        stream, static_cast<std::uint32_t>(dataset.fields.size()));

    for (const NamedCoefficientField& field : dataset.fields) {
        test_field_dataset_detail::write_scalar(
            stream, static_cast<std::uint32_t>(field.name.size()));
        stream.write(
            field.name.data(),
            static_cast<std::streamsize>(field.name.size()));
        test_field_dataset_detail::write_scalar(
            stream, field.coefficient.minimum);
        test_field_dataset_detail::write_scalar(
            stream, field.coefficient.maximum);
        test_field_dataset_detail::write_scalar(
            stream, field.coefficient.actual_contrast);
        test_field_dataset_detail::write_scalar(
            stream,
            static_cast<std::uint64_t>(field.coefficient.values.size()));
        stream.write(
            reinterpret_cast<const char*>(field.coefficient.values.data()),
            static_cast<std::streamsize>(
                field.coefficient.values.size() * sizeof(double)));
    }
    if (!stream) {
        throw std::runtime_error(
            "failed to write test-field dataset: " + path.string());
    }
}

inline TestFieldDataset load_test_fields(
    const std::filesystem::path& path,
    int expected_fine_intervals,
    int expected_coarse_intervals) {
    std::ifstream stream(path, std::ios::binary);
    if (!stream) {
        throw std::runtime_error(
            "cannot open test-field dataset " + path.string() +
            "; run generate_test_fields first");
    }

    char magic[8]{};
    stream.read(magic, 8);
    constexpr char expected_magic[8] =
        {'T', 'G', 'I', 'F', 'L', 'D', 'S', '1'};
    if (!std::equal(std::begin(magic), std::end(magic),
                    std::begin(expected_magic))) {
        throw std::runtime_error(
            "unsupported test-field dataset: " + path.string());
    }
    const std::uint32_t version =
        test_field_dataset_detail::read_scalar<std::uint32_t>(stream);
    if (version != 1U) {
        throw std::runtime_error("unsupported test-field dataset version");
    }

    TestFieldDataset dataset;
    dataset.fine_intervals =
        test_field_dataset_detail::read_scalar<std::int32_t>(stream);
    dataset.coarse_intervals =
        test_field_dataset_detail::read_scalar<std::int32_t>(stream);
    const std::uint32_t count =
        test_field_dataset_detail::read_scalar<std::uint32_t>(stream);
    if (dataset.fine_intervals != expected_fine_intervals ||
        dataset.coarse_intervals != expected_coarse_intervals ||
        count == 0U || count > 1024U) {
        throw std::runtime_error(
            "test-field dataset grid or field count does not match");
    }

    const std::uint64_t expected_values =
        static_cast<std::uint64_t>(dataset.fine_intervals - 1) *
        static_cast<std::uint64_t>(dataset.fine_intervals - 1);
    dataset.fields.reserve(count);
    for (std::uint32_t index = 0; index < count; ++index) {
        const std::uint32_t name_size =
            test_field_dataset_detail::read_scalar<std::uint32_t>(stream);
        if (name_size == 0U || name_size > 1024U) {
            throw std::runtime_error("invalid test-field name");
        }
        NamedCoefficientField field;
        field.name.resize(name_size);
        stream.read(
            field.name.data(),
            static_cast<std::streamsize>(name_size));
        field.coefficient.minimum =
            test_field_dataset_detail::read_scalar<double>(stream);
        field.coefficient.maximum =
            test_field_dataset_detail::read_scalar<double>(stream);
        field.coefficient.actual_contrast =
            test_field_dataset_detail::read_scalar<double>(stream);
        const std::uint64_t value_count =
            test_field_dataset_detail::read_scalar<std::uint64_t>(stream);
        if (value_count != expected_values ||
            value_count >
                std::numeric_limits<std::size_t>::max() / sizeof(double)) {
            throw std::runtime_error("invalid test-field value count");
        }
        field.coefficient.values.resize(
            static_cast<std::size_t>(value_count));
        stream.read(
            reinterpret_cast<char*>(field.coefficient.values.data()),
            static_cast<std::streamsize>(value_count * sizeof(double)));
        dataset.fields.push_back(std::move(field));
    }
    test_field_dataset_detail::require_stream(stream, path);
    if (stream.peek() != std::char_traits<char>::eof()) {
        throw std::runtime_error(
            "test-field dataset has trailing data: " + path.string());
    }
    test_field_dataset_detail::validate(dataset);
    return dataset;
}

} // namespace experiment_support
