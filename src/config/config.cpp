#include "config.hpp"
#include <fstream>
#include <sstream>
#include <cctype>
#include <algorithm>

namespace cardillo::config {

static inline std::string trim(const std::string& s) {
    size_t b = 0, e = s.size();
    while (b < e && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(s[e-1]))) --e;
    return s.substr(b, e - b);
}

static inline bool iequals(const std::string& a, const std::string& b) {
    return std::equal(a.begin(), a.end(), b.begin(), b.end(), [](char c1, char c2){
        return std::tolower(static_cast<unsigned char>(c1)) == std::tolower(static_cast<unsigned char>(c2));
    });
}

Config ConfigReader::fromFile(const std::string& path) {
    Config cfg; // start with defaults from header
    std::ifstream f(path);
    if (!f) return cfg; // return defaults if file can't be opened
    std::string line;
    while (std::getline(f, line)) {
        std::string s = trim(line);
        if (s.empty() || s[0] == '#') continue;
        auto pos = s.find('=');
        if (pos == std::string::npos) continue;
        std::string key = trim(s.substr(0, pos));
        std::string val = trim(s.substr(pos + 1));
        
        if (key == "pj.max_iterations") {
            try { cfg.pj_max_iterations = std::max(1, std::stoi(val)); cfg.has_pj_max_iterations = true; } catch (...) {}
        }
        else if (key == "pj.tol_abs") {
            try { cfg.pj_tol_abs = static_cast<real_t>(std::stod(val)); cfg.has_pj_tol_abs = true; } catch (...) {}
        }
        else if (key == "pj.relaxation") {
            try { cfg.pj_relaxation = static_cast<real_t>(std::stod(val)); cfg.has_pj_relaxation = true; } catch (...) {}
        }
        else if (key == "pj.alpha") {
            try { cfg.pj_alpha = static_cast<real_t>(std::stod(val)); cfg.has_pj_alpha = true; } catch (...) {}
        }
        else if (key == "pj.compliance") {
            try { cfg.pj_compliance = static_cast<real_t>(std::stod(val)); cfg.has_pj_compliance = true; } catch (...) {}
        }
        else if (key == "sim.T") {
            try { cfg.sim_T = static_cast<real_t>(std::stod(val)); } catch (...) {}
        }
        else if (key == "sim.dt") {
            try { cfg.sim_dt = static_cast<real_t>(std::stod(val)); } catch (...) {}
        }
        else if (key == "output.interval_steps") {
            try { cfg.output_interval_steps = std::max(1, std::stoi(val)); } catch (...) {}
        }
    }
    return cfg;
}

}
