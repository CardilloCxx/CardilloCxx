#pragma once
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace cardillo { namespace misc {

class ProgressBar {
public:
    using clock = std::chrono::steady_clock;

    explicit ProgressBar(std::size_t total,
                         std::ostream& os = std::cout,
                         std::size_t barBlocks = 50, // number of block chars inside the bar after percent marker
                         bool enabled = true)
        : m_total(total), m_os(os), m_barBlocks(barBlocks), m_enabled(enabled) {
            
        m_start = clock::now();
        m_lastDraw = m_start;
    }

    // Backward-compatible overload: third argument as 'enabled' (bool)
    explicit ProgressBar(std::size_t total,
                         std::ostream& os,
                         bool enabled)
        : ProgressBar(total, os, (std::size_t)24, enabled) {}

    void set_description(const std::string& desc) { m_desc = desc; }
    void set_postfix(const std::string& postfix) { m_postfix = postfix; }

    void update(std::size_t n = 1) { m_current += n; maybeDraw(); }
    void advance_to(std::size_t n) { m_current = n; maybeDraw(); }

    void refresh() { draw(false); }

    void close() {
        if (!m_enabled || m_closed) return;
        m_current = m_total; // snap to total
        draw(true);
        m_os << std::endl;
        m_closed = true;
    }

    ~ProgressBar() { if (!m_closed) close(); }

private:
    std::size_t m_total{0};
    std::size_t m_current{0};
    std::ostream& m_os;
    std::size_t m_barBlocks{24};
    double m_minRefresh{0.01}; // kept internal constant refresh cadence
    bool m_enabled{true};
    bool m_closed{false};
    std::string m_desc;
    std::string m_postfix;
    clock::time_point m_start;
    clock::time_point m_lastDraw;

    static std::string fmt_mmss(double sec) {
        if (sec < 0) sec = 0;
        int m = static_cast<int>(sec / 60.0); sec -= m * 60.0;
        int s = static_cast<int>(sec + 0.5);
        std::ostringstream oss; oss << std::setfill('0') << std::setw(2) << m << ":" << std::setw(2) << s; return oss.str();
    }

    void maybeDraw() {
        if (!m_enabled) return;
        auto now = clock::now();
        double sinceLast = std::chrono::duration<double>(now - m_lastDraw).count();
        if (sinceLast >= m_minRefresh || m_current == m_total) draw(false);
    }

    void draw(bool force) {
        if (!m_enabled) return;
        auto now = clock::now();
        m_lastDraw = now;
        double elapsed = std::chrono::duration<double>(now - m_start).count();
        double frac = (m_total > 0) ? static_cast<double>(m_current) / static_cast<double>(m_total) : 0.0;
        if (frac < 0) frac = 0; if (frac > 1) frac = 1;
        // Determine filled cells (integer percent of barBlocks, clamp to [0, barBlocks])
        std::size_t filled = (m_barBlocks > 0)
            ? static_cast<std::size_t>(frac * static_cast<double>(m_barBlocks))
            : 0;
        if (filled > m_barBlocks) filled = m_barBlocks;
        double it_per_sec = (elapsed > 1e-9) ? static_cast<double>(m_current) / elapsed : 0.0;
        double eta = (it_per_sec > 1e-12 && m_total > 0) ? (static_cast<double>(m_total - m_current) / it_per_sec) : 0.0;

        int pctInt = static_cast<int>(frac * 100.0 + 0.5);
        std::ostringstream line;
        if (!m_desc.empty()) line << m_desc << ' ';
        line << pctInt << "%|";

        std::string bar; bar.reserve(m_barBlocks);
        for (std::size_t i = 0; i < m_barBlocks; ++i) {
            if (i < filled) {
                bar += std::string("█");
            } else {
                bar += ' ';
            }
        }
        line << bar;
        line << " | " << m_current << '/' << m_total << ' ';
        // timing segment: [elapsed<ETA, rate it/s]
        line << '[' << fmt_mmss(elapsed) << '<' << fmt_mmss(eta) << ", ";
        line << std::fixed << std::setprecision(2) << it_per_sec << " it/s]";
        if (!m_postfix.empty()) line << ' ' << m_postfix;

        m_os << '\r' << line.str();
        m_os.flush();
    }
};

}} // namespace cardillo::misc
