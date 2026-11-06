// Config parser: values and aliases are read, invalid values and missing files are errors.
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include <cardillo/config/config.hpp>

using namespace cardillo;

namespace {

int failures = 0;

void check(bool ok, const std::string& message) {
    if (!ok) {
        ++failures;
        std::cerr << "FAILED: " << message << '\n';
    }
}

const std::string kPath = "test_config.tmp.config";

config::Config parse(const std::string& text) {
    std::ofstream(kPath) << text;
    return config::ConfigReader::fromFile(kPath);
}

bool throws(const std::string& text) {
    try {
        parse(text);
    } catch (const std::runtime_error&) {
        return true;
    }
    return false;
}

}  // namespace

int main() {
    const auto cfg = parse(
        "# comment\n"
        "solver.type = condensed   # inline comment\n"
        "solver.max_iterations = 42\n"
        "moreau.theta=0.75\n"
        "pj.nesterov = off\n"
        "sim.gravity = 0, 0, -1.5\n"
        "condensed.sweep_mode = Colored\n"
        "output.interval_steps = 0\n"
        "scene.height = 0.25\n");
    check(cfg.solver == config::SolverType::Condensed, "solver.type");
    check(cfg.pj_max_iterations == 42, "solver. alias of pj.max_iterations");
    check(cfg.moreau_theta == (real_t)0.75, "moreau.theta");
    check(!cfg.pj_nesterov, "boolean off");
    check(cfg.sim_gravity.isApprox(Vector3r(0, 0, -1.5)), "sim.gravity");
    check(cfg.condensed_sweep_mode == "colored", "choices are case-insensitive");
    check(cfg.output_interval_steps == 0, "output.interval_steps = 0 disables VTK output");
    check(cfg.sceneParam("height", 1) == (real_t)0.25 && cfg.sceneParam("missing", 1) == 1, "scene parameters");

    check(throws("pj.max_iterations = 10x\n"), "trailing garbage in a number");
    check(throws("pj.nesterov = maybe\n"), "invalid boolean");
    check(throws("solver.type = cg\n"), "unknown solver type");
    check(throws("moreau.theta = 1.5\n"), "theta outside (0, 1]");
    check(throws("friction.combine = max\n"), "invalid choice");
    check(throws("just text\n"), "line without '='");
    check(!throws("some.unknown_key = 1\n"), "unknown keys only warn");
    bool missing = false;
    try {
        config::ConfigReader::fromFile("does/not/exist.config");
    } catch (const std::runtime_error&) {
        missing = true;
    }
    check(missing, "missing config file");
    bool badScene = false;
    try {
        parse("scene.height = high\n").sceneParam("height", 0);
    } catch (const std::runtime_error&) {
        badScene = true;
    }
    check(badScene, "non-numeric scene parameter");

    std::remove(kPath.c_str());
    if (failures == 0) std::cout << "test_config: all checks passed" << std::endl;
    return failures == 0 ? 0 : 1;
}
