#pragma once

#include <fstream>
#include <string>
#include <vector>

namespace cardillo { namespace io {

class CsvWriter {
public:
    CsvWriter() = default;
    CsvWriter(const std::string& path, const std::vector<std::string>& header) {
        open(path, header);
    }

    bool isOpen() const { return static_cast<bool>(m_stream); }

    void open(const std::string& path, const std::vector<std::string>& header) {
        m_stream.open(path);
        if (!m_stream) return;
        for (size_t i = 0; i < header.size(); ++i) {
            m_stream << header[i];
            if (i + 1 < header.size()) m_stream << ",";
        }
        m_stream << "\n";
    }

    template <typename... Args>
    void writeRow(const Args&... args) {
        writeRowImpl(args...);
    }

private:
    void writeRowImpl() {
        m_stream << "\n";
    }

    template <typename T, typename... Rest>
    void writeRowImpl(const T& value, const Rest&... rest) {
        m_stream << value;
        if constexpr (sizeof...(Rest) > 0) {
            m_stream << ",";
            writeRowImpl(rest...);
        } else {
            m_stream << "\n";
        }
    }

    std::ofstream m_stream;
};

}} // namespace cardillo::io
