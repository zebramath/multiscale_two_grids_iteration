#pragma once

#include "core/linear_algebra.hpp"

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifndef TGI_RESULTS_DIR
#define TGI_RESULTS_DIR "results"
#endif

namespace experiment_support {

using Row = std::vector<std::string>;
using Rows = std::vector<Row>;
using Summary = std::vector<std::pair<std::string, std::string>>;

inline std::string fixed(double value, int precision = 3) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(precision) << value;
    return stream.str();
}

inline std::string scientific(double value, int precision = 3) {
    std::ostringstream stream;
    stream << std::scientific << std::setprecision(precision) << value;
    return stream.str();
}

inline void progress(const std::string& message) {
    std::cerr << "[progress] " << message << std::endl;
}

inline double interpolation_density_percent(
    const tgi::SparseMatrix& prolongation) {
    const double entries = static_cast<double>(prolongation.rows()) *
        static_cast<double>(prolongation.cols());
    return 100.0 * static_cast<double>(prolongation.nnz()) / entries;
}

inline std::filesystem::path results_directory() {
    const char* runtime_directory = std::getenv("TGI_RESULTS_DIR");
    const std::filesystem::path directory =
        runtime_directory != nullptr && runtime_directory[0] != '\0'
            ? runtime_directory
            : TGI_RESULTS_DIR;
    std::filesystem::create_directories(directory);
    return directory;
}

inline void save_csv(
    const std::string& name, const Row& headers, const Rows& rows) {
    const std::filesystem::path path =
        results_directory() / (name + ".csv");
    std::ofstream stream(path);
    auto write_row = [&stream](const Row& row) {
        for (std::size_t index = 0; index < row.size(); ++index) {
            if (index != 0U) stream << ',';
            stream << row[index];
        }
        stream << '\n';
    };
    write_row(headers);
    for (const Row& row : rows) write_row(row);
    if (!stream) {
        throw std::runtime_error(
            "cannot write result file " + path.string());
    }
}

class Report {
public:
    explicit Report(const std::string& title) {
        text_ << title << '\n' << std::string(title.size(), '=') << "\n\n";
    }

    void add_summary(const Summary& values) {
        std::size_t label_width = 0;
        for (const auto& item : values) {
            label_width = std::max(label_width, item.first.size());
        }
        for (const auto& [label, value] : values) {
            text_ << std::left << std::setw(static_cast<int>(label_width))
                  << label << " : " << value << '\n';
        }
        text_ << '\n';
    }

    void add_note(const std::string& note) {
        text_ << "Note: " << note << "\n\n";
    }

    void add_table(const std::string& section,
                   const Row& headers,
                   const std::vector<int>& widths,
                   const Rows& rows,
                   bool separate_first_column_groups = false) {
        text_ << section << '\n' << std::string(section.size(), '-') << '\n';
        write_text_row(headers, widths);
        for (std::size_t i = 0; i < widths.size(); ++i) {
            if (i != 0) text_ << "  ";
            text_ << std::string(static_cast<std::size_t>(widths[i]), '-');
        }
        text_ << '\n';
        std::string previous_group;
        for (const Row& row : rows) {
            if (separate_first_column_groups && !previous_group.empty() &&
                row.front() != previous_group) {
                text_ << '\n';
            }
            write_text_row(row, widths);
            previous_group = row.front();
        }
        text_ << '\n';
    }

    void save(const std::string& name) const {
        const std::filesystem::path path =
            results_directory() / (name + ".txt");
        std::ofstream stream(path);
        if (!stream || !(stream << text_.str())) {
            throw std::runtime_error(
                "cannot write result file " + path.string());
        }
        std::cout << text_.str() << "Saved: " << path.string() << '\n';
    }

private:
    void write_text_row(const Row& row, const std::vector<int>& widths) {
        for (std::size_t i = 0; i < row.size(); ++i) {
            if (i != 0) text_ << "  ";
            text_ << (i == 0 ? std::left : std::right)
                  << std::setw(widths[i]) << row[i];
        }
        text_ << '\n';
    }

    std::ostringstream text_;
};

}
