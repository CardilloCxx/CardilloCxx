#include "config.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace cardillo::config {

namespace {

std::string trim(const std::string& s) {
    size_t b = 0, e = s.size();
    while (b < e && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
    return s.substr(b, e - b);
}

std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return (char)std::tolower(c); });
    return s;
}

// Thrown by the value parsers; ConfigReader::fromFile() adds file:line and the key.
struct BadValue : std::runtime_error {
    using std::runtime_error::runtime_error;
};

real_t parseReal(const std::string& val) {
    std::size_t pos = 0;
    double x = 0;
    try {
        x = std::stod(val, &pos);
    } catch (const std::exception&) {
        throw BadValue("expected a number, got '" + val + "'");
    }
    if (pos != val.size()) throw BadValue("expected a number, got '" + val + "'");
    return static_cast<real_t>(x);
}

int parseInt(const std::string& val) {
    std::size_t pos = 0;
    int x = 0;
    try {
        x = std::stoi(val, &pos);
    } catch (const std::exception&) {
        throw BadValue("expected an integer, got '" + val + "'");
    }
    if (pos != val.size()) throw BadValue("expected an integer, got '" + val + "'");
    return x;
}

bool parseBool(const std::string& val) {
    const std::string v = lower(val);
    if (v == "1" || v == "true" || v == "yes" || v == "on") return true;
    if (v == "0" || v == "false" || v == "no" || v == "off") return false;
    throw BadValue("expected a boolean (true/false, yes/no, on/off, 1/0), got '" + val + "'");
}

// Returns the lower-cased value if it is one of `allowed`.
std::string parseChoice(const std::string& val, std::initializer_list<const char*> allowed) {
    const std::string v = lower(val);
    std::string list;
    for (const char* a : allowed) {
        if (v == a) return v;
        list += (list.empty() ? "" : ", ") + std::string(a);
    }
    throw BadValue("expected one of [" + list + "], got '" + val + "'");
}

void requireAtLeast(real_t x, real_t min, const std::string& val) {
    if (!(x >= min)) throw BadValue("value '" + val + "' must be >= " + std::to_string(min));
}

void requirePositive(real_t x, const std::string& val) {
    if (!(x > 0) || !std::isfinite(x)) throw BadValue("value '" + val + "' must be finite and positive");
}

using Setter = std::function<void(Config&, const std::string&)>;

// Every key (including its aliases) maps to the setter of its Config field.
const std::unordered_map<std::string, Setter>& setters() {
    static const std::unordered_map<std::string, Setter> table = [] {
        std::unordered_map<std::string, Setter> t;
        auto add = [&t](std::initializer_list<const char*> keys, Setter s) {
            for (const char* k : keys) t.emplace(k, s);
        };
        auto real = [](real_t Config::* f) { return [f](Config& c, const std::string& v) { c.*f = parseReal(v); }; };
        auto realMin = [](real_t Config::* f, real_t min) {
            return [f, min](Config& c, const std::string& v) {
                const real_t x = parseReal(v);
                requireAtLeast(x, min, v);
                c.*f = x;
            };
        };
        auto realPos = [](real_t Config::* f) {
            return [f](Config& c, const std::string& v) {
                const real_t x = parseReal(v);
                requirePositive(x, v);
                c.*f = x;
            };
        };
        auto intMin = [](int Config::* f, int min) {
            return [f, min](Config& c, const std::string& v) {
                const int x = parseInt(v);
                requireAtLeast(x, min, v);
                c.*f = x;
            };
        };
        auto integer = [](int Config::* f) { return [f](Config& c, const std::string& v) { c.*f = parseInt(v); }; };
        auto boolean = [](bool Config::* f) { return [f](Config& c, const std::string& v) { c.*f = parseBool(v); }; };
        auto text = [](std::string Config::* f) { return [f](Config& c, const std::string& v) { c.*f = v; }; };
        auto choice = [](std::string Config::* f, std::initializer_list<const char*> allowed) {
            const std::vector<const char*> list(allowed);
            return [f, list](Config& c, const std::string& v) {
                const std::string l = lower(v);
                std::string all;
                for (const char* a : list) {
                    if (l == a) {
                        c.*f = l;
                        return;
                    }
                    all += (all.empty() ? "" : ", ") + std::string(a);
                }
                throw BadValue("expected one of [" + all + "], got '" + v + "'");
            };
        };

        // Projected solvers (PJ, PGS, condensed); `solver.` is accepted as prefix as well.
        add({"pj.max_iterations", "solver.max_iterations"}, intMin(&Config::pj_max_iterations, 1));
        add({"pj.tol_abs", "solver.tol_abs"}, realPos(&Config::pj_tol_abs));
        add({"pj.tol_rel", "solver.tol_rel"}, realMin(&Config::pj_tol_rel, 0));
        add({"pj.criterion", "solver.criterion"}, [](Config& c, const std::string& v) {
            c.pj_criterion = parseChoice(v, {"velocity_rms", "natural_map"}) == "natural_map" ? PjTerminationCriterion::NaturalMap : PjTerminationCriterion::VelocityRms;
        });
        add({"pj.final_check", "solver.final_check"}, [](Config& c, const std::string& v) {
            c.pj_final_check = parseChoice(v, {"none", "natural_map"}) == "natural_map" ? PjFinalCheck::NaturalMap : PjFinalCheck::None;
        });
        add({"pj.relaxation", "solver.relaxation"}, realPos(&Config::pj_relaxation));
        add({"pj.alpha", "solver.alpha"}, realPos(&Config::pj_alpha));
        add({"pj.nesterov", "solver.nesterov"}, boolean(&Config::pj_nesterov));
        add({"pj.nesterov_beta_threshold", "solver.nesterov_beta_threshold"}, real(&Config::pj_nesterov_beta_threshold));
        add({"pj.nesterov_restart_limit", "solver.nesterov_restart_limit"}, intMin(&Config::pj_nesterov_restart_limit, 0));
        add({"pj.chebyshev", "solver.chebyshev"}, boolean(&Config::pj_chebyshev));
        add({"pj.anderson", "solver.anderson"}, boolean(&Config::pj_anderson));
        add({"pj.anderson_m", "solver.anderson_m"}, intMin(&Config::pj_anderson_m, 1));
        add({"pj.warmstart", "solver.warmstart"}, boolean(&Config::pj_warmstart));
        add({"pj.rdiag_true_delassus", "solver.rdiag_true_delassus"}, boolean(&Config::pj_rdiag_true_delassus));
        add({"debug.rb"}, boolean(&Config::debug_rb));
        add({"debug.pj"}, boolean(&Config::debug_pj));

        add({"solver.name", "integrator"}, [](Config& c, const std::string& v) {
            parseChoice(v, {"moreau", "first_order", "theta"});
            c.integrator = IntegratorType::Moreau;
        });
        add({"solver.type"}, [](Config& c, const std::string& v) {
            const std::string s = parseChoice(v, {"projected_jacobi", "pj", "condensed", "projected_gauss_seidel", "pgs", "naive_schur", "qoco", "clarabel", "conicxx", "scs", "mosek"});
            if (s == "projected_jacobi" || s == "pj")
                c.solver = SolverType::ProjectedJacobi;
            else if (s == "condensed")
                c.solver = SolverType::Condensed;
            else if (s == "projected_gauss_seidel" || s == "pgs")
                c.solver = SolverType::ProjectedGaussSeidel;
            else if (s == "naive_schur")
                c.solver = SolverType::NaiveSchur;
            else if (s == "qoco")
                c.solver = SolverType::Qoco;
            else if (s == "clarabel")
                c.solver = SolverType::Clarabel;
            else if (s == "conicxx")
                c.solver = SolverType::Conicxx;
            else if (s == "scs")
                c.solver = SolverType::Scs;
            else
                c.solver = SolverType::Mosek;
        });

        // Condensed solver
        add({"condensed.sweep_mode"}, choice(&Config::condensed_sweep_mode, {"gauss_seidel", "colored"}));
        add({"condensed.local_solve"}, choice(&Config::condensed_local_solve, {"projection", "newton"}));
        add({"condensed.num_threads"}, intMin(&Config::condensed_num_threads, 0));
        add({"condensed.newton_max_iters"}, intMin(&Config::condensed_newton_max_iters, 1));
        add({"condensed.newton_tol"}, realPos(&Config::condensed_newton_tol));
        add({"condensed.newton_rho_strategy"}, choice(&Config::condensed_newton_rho_strategy, {"split", "full"}));
        add({"condensed.projection_rho_strategy"}, choice(&Config::condensed_projection_rho_strategy, {"diagonal", "split", "full"}));
        add({"condensed.true_schur"}, boolean(&Config::condensed_true_schur));

        // Interior-point backends (QOCO, Clarabel, ConicXX)
        add({"solver.kkt_static_reg", "ip.kkt_static_reg", "qoco.kkt_static_reg", "clarabel.static_regularization_constant"}, realMin(&Config::ip_kkt_static_reg, 0));
        add({"solver.kkt_dynamic_reg", "ip.kkt_dynamic_reg", "qoco.kkt_dynamic_reg", "clarabel.dynamic_regularization_delta"}, realMin(&Config::ip_kkt_dynamic_reg, 0));
        add({"ip.tol_abs"}, realPos(&Config::ip_tol_abs));
        add({"ip.tol_rel"}, realMin(&Config::ip_tol_rel, 0));
        add({"solver.iter_ref_iters", "ip.iter_ref_iters", "qoco.iter_ref_iters", "clarabel.iterative_refinement_max_iter"}, intMin(&Config::ip_iter_ref_iters, 0));
        add({"conicxx.warm_start"}, boolean(&Config::conicxx_warm_start));
        add({"conicxx.warm_lambda"}, real(&Config::conicxx_warm_lambda));
        add({"conicxx.warm_start_after_setup"}, boolean(&Config::conicxx_warm_start_after_setup));
        add({"conicxx.compliant_rows"}, boolean(&Config::conicxx_compliant_rows));
        add({"conicxx.verbose"}, intMin(&Config::conicxx_verbose, 0));
        add({"conicxx.dynamic_reg_eps"}, real(&Config::conicxx_dynamic_reg_eps));
        add({"conicxx.refine_tol"}, real(&Config::conicxx_refine_tol));
        add({"conicxx.max_step_fraction"}, real(&Config::conicxx_max_step_fraction));
        add({"conicxx.equilibrate"}, boolean(&Config::conicxx_equilibrate));
        add({"conicxx.equilibrate_max_iter"}, intMin(&Config::conicxx_equilibrate_max_iter, 0));
        add({"conicxx.equilibrate_min_scale"}, real(&Config::conicxx_equilibrate_min_scale));
        add({"conicxx.equilibrate_max_scale"}, real(&Config::conicxx_equilibrate_max_scale));
        add({"conicxx.validate_inputs"}, boolean(&Config::conicxx_validate_inputs));

        // MOSEK backend (negative values and "default" keep MOSEK's own defaults)
        add({"mosek.num_threads"}, intMin(&Config::mosek_num_threads, 0));
        add({"mosek.tol_rel_gap"}, real(&Config::mosek_tol_rel_gap));
        add({"mosek.tol_pfeas"}, real(&Config::mosek_tol_pfeas));
        add({"mosek.tol_dfeas"}, real(&Config::mosek_tol_dfeas));
        add({"mosek.max_iterations"}, integer(&Config::mosek_max_iterations));
        add({"mosek.max_time"}, real(&Config::mosek_max_time));
        add({"mosek.presolve"}, choice(&Config::mosek_presolve, {"default", "on", "off", "free"}));
        add({"mosek.scaling"}, boolean(&Config::mosek_scaling));

        // SCS backend (negative values keep SCS's own defaults)
        add({"scs.eps_abs"}, real(&Config::scs_eps_abs));
        add({"scs.eps_rel"}, real(&Config::scs_eps_rel));
        add({"scs.eps_infeas"}, real(&Config::scs_eps_infeas));
        add({"scs.time_limit_secs"}, real(&Config::scs_time_limit_secs));
        add({"scs.alpha"}, real(&Config::scs_alpha));
        add({"scs.scale"}, real(&Config::scs_scale));
        add({"scs.rho_x"}, real(&Config::scs_rho_x));
        add({"scs.max_iters"}, integer(&Config::scs_max_iters));
        add({"scs.acceleration_lookback"}, integer(&Config::scs_acceleration_lookback));
        add({"scs.acceleration_interval"}, integer(&Config::scs_acceleration_interval));
        add({"scs.adaptive_scale"}, boolean(&Config::scs_adaptive_scale));
        add({"scs.normalize"}, boolean(&Config::scs_normalize));
        add({"scs.warm_start"}, boolean(&Config::scs_warm_start));
        add({"scs.carry_scale"}, boolean(&Config::scs_carry_scale));
        add({"scs.stats_csv"}, text(&Config::scs_stats_csv));
        add({"scs.dump_dir"}, text(&Config::scs_dump_dir));
        add({"scs.dump_every"}, intMin(&Config::scs_dump_every, 1));

        // Time integration and simulation
        add({"moreau.theta", "solver.theta"}, [](Config& c, const std::string& v) {
            const real_t x = parseReal(v);
            if (!(x > 0 && x <= 1)) throw BadValue("moreau.theta must lie in (0, 1], got '" + v + "'");
            c.moreau_theta = x;
        });
        add({"moreau.implicit_gyroscopy", "solver.implicit_gyroscopy"}, boolean(&Config::moreau_implicit_gyroscopy));
        add({"sim.T"}, realPos(&Config::sim_T));
        add({"sim.dt"}, realPos(&Config::sim_dt));
        add({"sim.gravity"}, [](Config& c, const std::string& v) {
            std::string s = v;
            std::replace(s.begin(), s.end(), ',', ' ');  // commas or spaces as separators
            std::istringstream ss(s);
            double g[3];
            std::string rest;
            if (!(ss >> g[0] >> g[1] >> g[2]) || (ss >> rest)) throw BadValue("expected three numbers, got '" + v + "'");
            c.sim_gravity = Vector3r((real_t)g[0], (real_t)g[1], (real_t)g[2]);
        });

        // Output
        add({"output.write_contact_manifolds"}, boolean(&Config::output_write_contact_manifolds));
        add({"output.interval_steps"}, intMin(&Config::output_interval_steps, 0));
        add({"output.folder"}, text(&Config::output_folder));
        add({"output.filename_prefix"}, text(&Config::output_filename_prefix));

        // Collision
        add({"collision.broadphase"}, choice(&Config::collision_broadphase, {"dynamic_aabb", "daabb", "default", "dynamic_aabb_array", "daabb_array", "naive", "bruteforce", "sap", "ssap"}));
        add({"collision.disable_all"}, boolean(&Config::collision_disable_all));
        add({"collision.max_raw_contacts"}, intMin(&Config::collision_max_raw_contacts, 1));
        add({"collision.max_patches"}, [](Config& c, const std::string& v) {
            const int x = parseInt(v);
            requireAtLeast(x, 0, v);
            c.collision_max_patches = static_cast<std::size_t>(x);
        });
        add({"collision.patch_num_samples"}, [](Config& c, const std::string& v) {
            const int x = parseInt(v);
            requireAtLeast(x, 4, v);
            c.collision_patch_num_samples = static_cast<std::size_t>(x);
        });
        add({"collision.patch_tolerance"}, realPos(&Config::collision_patch_tolerance));
        add({"collision.use_patch_vertices"}, boolean(&Config::collision_use_patch_vertices));
        add({"collision.match_max_dist"}, realMin(&Config::collision_match_max_dist, 0));
        add({"collision.min_pair_contact_distance"}, realMin(&Config::collision_min_pair_contact_distance, 0));
        add({"collision.security_margin"}, realMin(&Config::collision_security_margin, 0));

        // Friction and restitution
        add({"friction.enable"}, boolean(&Config::friction_enable));
        add({"friction.combine"}, choice(&Config::friction_combine, {"min", "arith", "geom"}));
        add({"friction.default_mu"}, realMin(&Config::friction_default_mu, 0));
        add({"restitution.default_normal"}, realMin(&Config::restitution_default_normal, 0));
        add({"restitution.default_tangential"}, real(&Config::restitution_default_tangential));

        add({"scene.name"}, text(&Config::scene_name));
        return t;
    }();
    return table;
}

}  // namespace

real_t Config::sceneParam(const std::string& name, real_t fallback) const {
    const auto it = scene_params.find(name);
    if (it == scene_params.end()) return fallback;
    try {
        return parseReal(it->second);
    } catch (const BadValue& e) {
        throw std::runtime_error("scene." + name + ": " + e.what());
    }
}

Config ConfigReader::fromFile(const std::string& path) {
    Config cfg;  // start with defaults from header
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Cannot open config file '" + path + "'");
    const auto& table = setters();
    std::string line;
    int lineNo = 0;
    while (std::getline(f, line)) {
        ++lineNo;
        // Strip inline comments starting with '#'
        std::string s = line.substr(0, line.find('#'));
        s = trim(s);
        if (s.empty()) continue;
        const auto where = path + ":" + std::to_string(lineNo);
        const auto pos = s.find('=');
        if (pos == std::string::npos) throw std::runtime_error(where + ": expected 'key = value', got '" + s + "'");
        const std::string key = trim(s.substr(0, pos));
        const std::string val = trim(s.substr(pos + 1));

        if (const auto it = table.find(key); it != table.end()) {
            try {
                it->second(cfg, val);
            } catch (const BadValue& e) {
                throw std::runtime_error(where + ": " + key + ": " + e.what());
            }
        } else if (key.rfind("scene.", 0) == 0) {
            cfg.scene_params[key.substr(6)] = val;
        } else {
            std::cerr << "Warning: " << where << ": unrecognized config key '" << key << "' is ignored." << std::endl;
        }
    }
    return cfg;
}

}  // namespace cardillo::config
