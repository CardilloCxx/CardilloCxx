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
        if (s.empty()) continue;
        // Strip inline comments starting with '#'
        auto hashPos = s.find('#');
        if (hashPos != std::string::npos) {
            s = trim(s.substr(0, hashPos));
        }
        if (s.empty()) continue;
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
        
        else if (key == "pj.nesterov") {
            std::string v = val; std::transform(v.begin(), v.end(), v.begin(), [](unsigned char c){ return (char)std::tolower(c); });
            cfg.pj_nesterov = (v == "1" || v == "true" || v == "yes" || v == "on");
        }
        else if (key == "pj.nesterov_beta_threshold") {
            try { cfg.pj_nesterov_beta_threshold = static_cast<real_t>(std::stod(val)); } catch (...) {}
        }
        else if (key == "pj.nesterov_restart_limit") {
            try { cfg.pj_nesterov_restart_limit = std::max(0, std::stoi(val)); } catch (...) {}
        }
        else if (key == "pj.warmstart") {
            std::string v = val; std::transform(v.begin(), v.end(), v.begin(), [](unsigned char c){ return (char)std::tolower(c); });
            cfg.pj_warmstart = (v == "1" || v == "true" || v == "yes" || v == "on");
        }
        else if (key == "debug.rb") {
            std::string v = val; std::transform(v.begin(), v.end(), v.begin(), [](unsigned char c){ return (char)std::tolower(c); });
            cfg.debug_rb = (v == "1" || v == "true" || v == "yes" || v == "on");
        }
        else if (key == "debug.pj") {
            std::string v = val; std::transform(v.begin(), v.end(), v.begin(), [](unsigned char c){ return (char)std::tolower(c); });
            cfg.debug_pj = (v == "1" || v == "true" || v == "yes" || v == "on");
        }
        else if (key == "debug.mesh") {
            std::string v = val; std::transform(v.begin(), v.end(), v.begin(), [](unsigned char c){ return (char)std::tolower(c); });
            cfg.debug_mesh = (v == "1" || v == "true" || v == "yes" || v == "on");
        }
        else if (key == "sim.T") {
            try { cfg.sim_T = static_cast<real_t>(std::stod(val)); } catch (...) {}
        }
        else if (key == "sim.dt") {
            try { cfg.sim_dt = static_cast<real_t>(std::stod(val)); } catch (...) {}
        }
        else if (key == "sim.gravity") {
            std::string v = val;
            // Allow commas or spaces as separators
            std::replace(v.begin(), v.end(), ',', ' ');
            std::istringstream ss(v);
            double gx=0, gy=0, gz=-9.81;
            if (ss >> gx >> gy >> gz) {
                cfg.sim_gravity = Vector3r((real_t)gx, (real_t)gy, (real_t)gz);
            }
        }
        else if (key == "output.write_contacts") {
            cfg.output_write_contacts = (iequals(val, "1") || iequals(val, "true") || iequals(val, "yes") || iequals(val, "on"));
        }
        else if (key == "output.interval_steps") {
            try { cfg.output_interval_steps = std::max(1, std::stoi(val)); } catch (...) {}
        }
        else if (key == "output.heightfield_stride") {
            try { cfg.output_heightfield_stride = std::max(1, std::stoi(val)); } catch (...) {}
        }
        else if (key == "output.contacts_body_vectors") {
            std::string v = val; std::transform(v.begin(), v.end(), v.begin(), [](unsigned char c){ return (char)std::tolower(c); });
            cfg.output_contacts_body_vectors = (v == "1" || v == "true" || v == "yes" || v == "on");
        }
        else if (key == "output.folder") {
            cfg.output_folder = val;
        }
        else if (key == "output.filename_prefix") {
            cfg.output_filename_prefix = val;
        }
        else if (key == "collision.broadphase") {
            cfg.collision_broadphase = val;
        }
        else if (key == "collision.max_raw_contacts") {
            try { cfg.collision_max_raw_contacts = std::max(1, std::stoi(val)); } catch (...) {}
        }
        else if (key == "collision.max_patches") {
            try {
                int v = std::max(0, std::stoi(val));
                cfg.collision_max_patches = static_cast<std::size_t>(v);
            } catch (...) {}
        }
        else if (key == "collision.use_patch_vertices") {
            std::string v = val; std::transform(v.begin(), v.end(), v.begin(), [](unsigned char c){ return (char)std::tolower(c); });
            cfg.collision_use_patch_vertices = (v == "1" || v == "true" || v == "yes" || v == "on");
        }
        else if (key == "collision.match_max_dist") {
            try {
                cfg.collision_match_max_dist = static_cast<real_t>(std::stod(val));
                if (cfg.collision_match_max_dist < (real_t)0) cfg.collision_match_max_dist = (real_t)0;
            } catch (...) {}
        }
        else if (key == "collision.min_pair_contact_distance") {
            try {
                cfg.collision_min_pair_contact_distance = static_cast<real_t>(std::stod(val));
                if (cfg.collision_min_pair_contact_distance < (real_t)0) cfg.collision_min_pair_contact_distance = (real_t)0;
            } catch (...) {}
        }
        else if (key == "collision.security_margin") {
            try {
                cfg.collision_security_margin = static_cast<real_t>(std::stod(val));
            } catch (...) {}
        }
        else if (key == "friction.enable") {
            std::string v = val; std::transform(v.begin(), v.end(), v.begin(), [](unsigned char c){ return (char)std::tolower(c); });
            cfg.friction_enable = (v == "1" || v == "true" || v == "yes" || v == "on");
        }
        else if (key == "friction.combine") {
            std::string v = val; std::transform(v.begin(), v.end(), v.begin(), [](unsigned char c){ return (char)std::tolower(c); });
            if (v == "min" || v == "arith" || v == "geom") cfg.friction_combine = v;
        }
        else if (key == "friction.default_mu") {
            try { cfg.friction_default_mu = static_cast<real_t>(std::stod(val)); } catch (...) {}
        }
        else if (key == "scene.name") {
            cfg.scene_name = val;
        }
    }
    return cfg;
}

}
