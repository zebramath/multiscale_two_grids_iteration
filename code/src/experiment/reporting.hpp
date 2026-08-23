#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
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

using Clock = std::chrono::steady_clock;
using Row = std::vector<std::string>;
using Rows = std::vector<Row>;
using Summary = std::vector<std::pair<std::string, std::string>>;

inline double milliseconds(Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

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

inline std::string integer(double value) {
    return std::to_string(static_cast<long long>(std::llround(value)));
}

inline void progress(const std::string& message) {
    std::cerr << "[progress] " << message << std::endl;
}

inline std::filesystem::path write_result(
    const std::string& name, const std::string& contents) {
    const std::filesystem::path directory = TGI_RESULTS_DIR;
    std::filesystem::create_directories(directory);
    const std::filesystem::path path = directory / (name + ".txt");
    std::ofstream stream(path);
    if (!stream || !(stream << contents)) {
        throw std::runtime_error(
            "cannot write result file " + path.string());
    }
    return path;
}

inline std::filesystem::path write_csv(
    const std::string& name, const Row& headers, const Rows& rows) {
    const std::filesystem::path directory = TGI_RESULTS_DIR;
    std::filesystem::create_directories(directory);
    const std::filesystem::path path = directory / (name + ".csv");
    std::ofstream stream(path);
    const auto write_cell = [&](const std::string& cell) {
        const bool quote = cell.find_first_of(",\"\n\r") !=
            std::string::npos;
        if (!quote) {
            stream << cell;
            return;
        }
        stream << '"';
        for (char character : cell) {
            if (character == '"') stream << '"';
            stream << character;
        }
        stream << '"';
    };
    const auto write_row = [&](const Row& row) {
        for (std::size_t index = 0; index < row.size(); ++index) {
            if (index != 0U) stream << ',';
            write_cell(row[index]);
        }
        stream << '\n';
    };
    write_row(headers);
    for (const Row& row : rows) write_row(row);
    if (!stream) {
        throw std::runtime_error(
            "cannot write CSV result file " + path.string());
    }
    return path;
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
            if (separate_first_column_groups && !row.empty() &&
                !previous_group.empty() &&
                row.front() != previous_group) {
                text_ << '\n';
            }
            write_text_row(row, widths);
            if (!row.empty()) previous_group = row.front();
        }
        text_ << '\n';
    }

    void save(const std::string& name) const {
        const std::filesystem::path path = write_result(name, text_.str());
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
