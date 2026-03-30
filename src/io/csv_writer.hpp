#pragma once

#include <fstream>
#include <string>
#include <vector>
#include <entt/entt.hpp>
#include "../physics/ecs_types.hpp"
#include "../config/config.hpp"
#include "../misc/types.hpp"

namespace cardillo { namespace io {

class CsvWriter {
public:
    CsvWriter() = default;
    CsvWriter(const std::string& path, const std::vector<std::string>& header) {
        open(path, header);
    }

    // Construct with config: store pending path for lazy opening using config values
    CsvWriter(const config::Config& cfg, const std::string& suffix = "_tracked.csv") {
        m_pending_path = cfg.output_folder + "/" + cfg.output_filename_prefix + suffix;
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

    // Public API for pipeline-level writer: write tracked state (opens lazily)
    void writeTrackedState(real_t t, const entt::registry& reg) {
        if (!hasTracked(reg)) return;
        if (!isOpen() && !m_pending_path.empty()) {
            open(m_pending_path, trackedHeader());
        }
        if (!isOpen()) return;
        writeEntities(t, reg);
    }

private:
    static bool hasTracked(const entt::registry& reg) {
        for (auto e : reg.view<cardillo::C_TrackTag>()) { (void)e; return true; }
        return false;
    }

    static std::vector<std::string> trackedHeader() {
        return std::vector<std::string>{
            "t",
            "name",
            "px","py","pz",
            "vx","vy","vz",
            "wx","wy","wz",
            "euler_x","euler_y","euler_z"
        };
    }

    
    void writeEntities(real_t t, const entt::registry& reg) {
        auto view = reg.view<cardillo::C_TrackTag, cardillo::C_Position3, cardillo::C_LinearVelocity3, cardillo::C_AngularVelocity3, cardillo::C_Orientation>();
        for (auto e : view) {
            const auto& tag   = view.get<cardillo::C_TrackTag>(e);
            const auto& pos   = view.get<cardillo::C_Position3>(e).value;
            const auto& v     = view.get<cardillo::C_LinearVelocity3>(e).value;
            const auto& w     = view.get<cardillo::C_AngularVelocity3>(e).value;
            const auto& euler = view.get<cardillo::C_Orientation>(e).value.toRotationMatrix().eulerAngles(0, 1, 2);
            writeRow(t,
                     tag.name,
                     pos.x(), pos.y(), pos.z(),
                     v.x(), v.y(), v.z(),
                     w.x(), w.y(), w.z(),
                     euler.x(), euler.y(), euler.z());
        }
    }

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
    std::string m_pending_path;
};

}} // namespace cardillo::io
