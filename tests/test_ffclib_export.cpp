// ffclib export (io/ffclib/): one-step problems written as ffclib files (fclib's global problem with the
// compliance C in /fclib_global/C; reading fclib's files and fclib reading ours is tested in ffclib itself).
//  1. Round trip: ffclib::read() returns the captured problem, the solver's solution and its warm start bit for
//     bit.
//  2. Errors: the velocity of the stored percussions is the stored velocity (balance), the bilateral rows from
//     the file reproduce the engine's final check, and for the Coulomb solvers ffclib's Coulomb error agrees
//     with the engine's natural map (same law; ffclib uses the Euclidean norm, the engine max(|x_N|, ||x_T||),
//     so the two differ by at most a factor sqrt(2)); this catches sign, scaling and frame errors.
//  3. ffclib's cone program of the file equals the engine's ConeProgramAssembler, except for the frozen
//     De Saxce shift on the normal rows that the engine adds and ffclib's convexified law leaves out.
//  4. Determinism, and the exporter's step selection, file names, titles and manifest.
// Scenes: a sliding ball, a box stack with restitution (with and without friction), and two spheres on the
// ground tied by a perfect distance constraint and a compliant spring-damper.

#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include <ffclib/error.hpp>
#include <ffclib/io.hpp>
#include <ffclib/solver/conic/cone_program.hpp>

#include <cardillo/config/config.hpp>
#include <cardillo/io/ffclib/ffclib_capture.hpp>
#include <cardillo/io/ffclib/ffclib_exporter.hpp>
#include <cardillo/physics/api/physics_engine.hpp>
#include <cardillo/physics/assembly/dynamics_assembler.hpp>
#include <cardillo/physics/ecs_types.hpp>
#include <cardillo/physics/integration/integration_base.hpp>
#include <cardillo/physics/solver/natural_map.hpp>
#include <cardillo/physics/solver/solver_base.hpp>

#if defined(CARDILLO_HAVE_CLARABEL) || defined(CARDILLO_HAVE_SCS)
#define FFCLIB_TEST_ENGINE_CONE_PROGRAM 1
#include <cardillo/physics/solver/conic/cone_program_assembler.hpp>
#endif

using namespace cardillo;
namespace fx = cardillo::io::ffclib_export;

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAILED: " << message << '\n';
    }
}

bool sameSparse(const CscMatrix& a, const CscMatrix& b) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) return false;
    const CscMatrix d = a - b;
    for (Eigen::Index j = 0; j < d.outerSize(); ++j)
        for (CscMatrix::InnerIterator it(d, j); it; ++it)
            if (it.value() != 0) return false;
    return true;
}

bool sameVec(const VectorXr& a, const VectorXr& b) {
    return a.size() == b.size() && (a.array() == b.array()).all();
}

std::string sci(real_t x) {
    char buf[32];
    std::snprintf(buf, sizeof buf, "%.3e", x);
    return buf;
}

real_t relDiff(real_t a, real_t b) {
    return std::abs(a - b) / std::max<real_t>(1e-300, std::max(std::abs(a), std::abs(b)));
}

std::filesystem::path tempDir(const std::string& name) {
    const auto dir = std::filesystem::temp_directory_path() / ("cardillo_test_ffclib_" + name);
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);
    return dir;
}

// ---------- scenes ----------

using Scene = std::function<void(physics::PhysicsEngine&)>;

void ballScene(physics::PhysicsEngine& e) {
    e.addStaticBody(physics::PlaneShape(Vector3r(0, 0, 1), Vector3r(0, 1, 0), 5, 5), physics::RigidState{});
    e.addRigidBody(physics::SphereShape(0.5), physics::RigidState(Vector3r(0, 0, 0.5 - 1e-6), Vector3r(1, 0, 0)), physics::RigidProps(1));
}

// Box stack: patch contacts (several points per face, redundant), the top box with restitution, captured
// at its impact in step 0. Also run with friction disabled (frictionless contacts: one row each).
void boxScene(physics::PhysicsEngine& e) {
    e.addStaticBody(physics::PlaneShape(Vector3r(0, 0, 1), Vector3r(0, 1, 0), 5, 5), physics::RigidState{});
    const real_t a = 0.1;
    e.addRigidBody(physics::CubeShape(Vector3r(a, a, a)), physics::RigidState(Vector3r(0, 0, a - 1e-5), Vector3r(0.2, 0, 0)), physics::RigidProps(2));
    const auto top =
        e.addRigidBody(physics::CubeShape(Vector3r(0.6 * a, 0.6 * a, 0.6 * a)), physics::RigidState(Vector3r(0.01, 0.02, 2 * a + 0.6 * a - 2e-5), Vector3r(0, 0, -0.05)), physics::RigidProps(0.5));
    e.ecs().emplace_or_replace<C_Restitution>(top, C_Restitution{0.3, 0});
}

// Two spheres on the ground: a perfect distance constraint to a static anchor and a compliant
// spring-damper between them (spring and damper rows).
void jointScene(physics::PhysicsEngine& e) {
    e.addStaticBody(physics::PlaneShape(Vector3r(0, 0, 1), Vector3r(0, 1, 0), 5, 5), physics::RigidState{});
    const auto anchor = e.addStaticBody(physics::SphereShape(0.01), physics::RigidState(Vector3r(-0.5, 0, 0.5)));
    const auto s1 = e.addRigidBody(physics::SphereShape(0.1), physics::RigidState(Vector3r(0, 0, 0.1 - 1e-6), Vector3r(0.3, 0.1, 0)), physics::RigidProps(1));
    const auto s2 = e.addRigidBody(physics::SphereShape(0.1), physics::RigidState(Vector3r(0.3, 0, 0.1 - 1e-6), Vector3r(-0.2, 0, 0)), physics::RigidProps(1.5));
    e.addLinearDistanceConstraint(anchor, s1);                                                  // perfect
    e.addLinearDistanceConstraint(s1, s2, Vector3r::Zero(), Vector3r::Zero(), 2e3, 5.0, 0.28);  // compliant spring and damper
}

config::Config baseConfig(config::SolverType solver) {
    config::Config cfg;
    cfg.solver = solver;
    cfg.sim_dt = 1e-3;
    cfg.sim_T = 0;  // stepped by the test (no progress bar)
    cfg.sim_gravity = Vector3r(0, 0, -9.81);
    cfg.output_interval_steps = 0;
    cfg.friction_default_mu = 0.4;
    cfg.moreau_theta = 0.6;
    cfg.pj_max_iterations = 200000;
    cfg.pj_tol_abs = 1e-12;
    cfg.pj_tol_rel = 1e-12;
    cfg.solver_criterion = config::TerminationCriterion::NaturalMap;
    cfg.condensed_num_threads = 1;
    cfg.ip_tol_abs = 1e-10;
    cfg.ip_tol_rel = 1e-10;
    cfg.scene_name = "test";
    return cfg;
}

// Captures the problem of one step like the exporter, writes it, and keeps the engine's view of the step.
struct Capture : integration::StepObserver {
    int target_step{0};
    std::string path;
    int step{-1};
    fx::WarmStart warm;
    ffclib::Record rec;
    solver::StepResidual engine;
    bool coulomb{true};  // the solver's law
    bool captured{false};
    VectorXr Cdiag_scaled;  // compliances of the engine's bilateral rows (zero: perfect constraint)
#ifdef FFCLIB_TEST_ENGINE_CONE_PROGRAM
    CscMatrix P, A;
    VectorXr q, b;
#endif

    void beforeSolve(const physics::DynamicsAssembler& dyn) override {
        if (++step == target_step) warm = fx::warmStart(dyn);
    }
    void afterSolve(physics::DynamicsAssembler& dyn, const VectorXr& u, real_t dt, real_t theta, const solver::SolverBase& s, const std::optional<solver::StepResidual>&, real_t) override {
        if (step != target_step) return;
        engine = solver::stepResidual(dyn, u, s.contactLaw(), dt, theta);
        coulomb = s.contactLaw() == solver::ContactLaw::Coulomb;
        rec = fx::capture(dyn, u, warm, dt, theta);
        rec.info.title = "test step " + std::to_string(step);
        ffclib::write(path, rec);
#ifdef FFCLIB_TEST_ENGINE_CONE_PROGRAM
        physics::assembly::ConeProgramAssembler cpa(dyn);
        P = cpa.P(dt, theta);
        q = cpa.q(dt, theta);
        A = cpa.A(dt, theta);
        b = cpa.b(dt, theta);
#endif
        Cdiag_scaled = VectorXr(dyn.numSprings() + dyn.numDampers());
        if (dyn.numSprings() > 0) Cdiag_scaled.head(dyn.numSprings()) = physics::DynamicsAssembler::springComplianceScale(dt, theta) * dyn.Cdiag();
        if (dyn.numDampers() > 0) Cdiag_scaled.tail(dyn.numDampers()) = physics::DynamicsAssembler::damperComplianceScale(dt, theta) * dyn.Adiag();
        captured = true;
    }
};

Capture runAndCapture(const Scene& scene, const config::Config& cfg, int targetStep, const std::string& path) {
    Capture cap;
    cap.target_step = targetStep;
    cap.path = path;
    physics::PhysicsEngine engine(cfg);
    scene(engine);
    engine.setStepObserver(&cap);
    for (int k = 0; k <= targetStep; ++k) engine.step();
    engine.setStepObserver(nullptr);
    return cap;
}

bool sameSolution(const ffclib::Solution& a, const ffclib::Solution& b) {
    return sameVec(a.v, b.v) && sameVec(a.r, b.r) && sameVec(a.lambda, b.lambda);
}

// The largest contact velocity of the motion without forces, max_c ||H_c^T M^-1 f + w_c||: the scale of the
// velocity errors.
real_t velocityScale(const ffclib::Problem& pb) {
    const VectorXr u = pb.H.transpose() * ffclib::MassInverse(pb.M)(pb.f) + pb.w;
    real_t scale = 0;
    for (int c = 0; c < pb.nc(); ++c) scale = std::max(scale, u.segment<3>(3 * c).norm());
    return scale > 0 ? scale : 1;
}

// ---------- checks ----------

void checkRoundTrip(const std::string& label, const Capture& cap) {
    const ffclib::Record back = ffclib::read(cap.path);
    const ffclib::Problem& a = cap.rec.problem;
    const ffclib::Problem& b = back.problem;
    check(sameSparse(a.M, b.M) && sameSparse(a.H, b.H) && sameSparse(a.G, b.G) && sameSparse(a.C, b.C), label + ": M, H, G, C round trip");
    check(sameVec(a.f, b.f) && sameVec(a.w, b.w) && sameVec(a.b, b.b) && sameVec(a.mu, b.mu), label + ": f, w, b, mu round trip");
    check(back.solution && sameSolution(*back.solution, *cap.rec.solution), label + ": solution round trip");
    check(back.guesses.size() == 1 && sameSolution(back.guesses[0], cap.rec.guesses[0]), label + ": warm start round trip");
    check(back.info.title == cap.rec.info.title && back.info.math_info == cap.rec.info.math_info, label + ": info round trip");
    // Only compliant rows have entries in C, on its diagonal.
    int compliant = 0;
    for (int k = 0; k < (int)cap.Cdiag_scaled.size(); ++k) compliant += cap.Cdiag_scaled[k] > 0;
    check(b.C.nonZeros() == compliant && sameVec(VectorXr(b.C.diagonal()), cap.Cdiag_scaled), label + ": C = diag(C_g) with the compliant rows only");
}

void checkErrors(const std::string& label, const Capture& cap, real_t balanceTol) {
    const ffclib::Record back = ffclib::read(cap.path);
    const ffclib::Problem& pb = back.problem;
    const ffclib::Solution& sol = *back.solution;
    const real_t scale = velocityScale(pb);
    const VectorXr v = ffclib::velocity(pb, sol.r, sol.lambda);
    const real_t balance = (v - sol.v).lpNorm<Eigen::Infinity>();
    check(balance <= balanceTol * scale, label + ": momentum balance " + sci(balance));

    // ffclib's errors at the stored solution: the bilateral rows are the engine's final check; the contact error
    // is the engine's natural map up to the norm (Euclidean per contact in ffclib).
    const ffclib::Errors e = ffclib::errors(pb, sol);
    check(relDiff(e.bilateral, cap.engine.bilateral) < 1e-9 || std::abs(e.bilateral - cap.engine.bilateral) < 1e-13,
          label + ": bilateral rows from the file " + sci(e.bilateral) + " vs engine " + sci(cap.engine.bilateral));

    if (!cap.coulomb) return;  // the engine's conic backends solve the convexified law with the frozen shift
    const real_t coulomb = e.coulomb, roundOff = 1e-14 * scale;
    check(coulomb >= cap.engine.coulomb - roundOff && coulomb <= std::sqrt(2.0) * cap.engine.coulomb + roundOff,
          label + ": ffclib's Coulomb error " + sci(coulomb) + " vs engine " + sci(cap.engine.coulomb));
    check(coulomb < 1e-5 * scale, label + ": the solution solves Coulomb's law (" + sci(coulomb) + ")");
}

#ifdef FFCLIB_TEST_ENGINE_CONE_PROGRAM
// ffclib's cone program of the file vs ConeProgramAssembler, whose perfect constraints have zero columns and
// whose normal rows carry the frozen De Saxce shift (>= 0, zero without friction).
void checkConeProgram(const std::string& label, const Capture& cap) {
    const ffclib::Record back = ffclib::read(cap.path);
    const ffclib::ConeProgram cp = ffclib::coneProgram(back.problem);
    const int n = back.problem.nv();
    std::vector<int> keep;
    for (int j = 0; j < n; ++j) keep.push_back(j);
    for (int k = 0; k < (int)cap.Cdiag_scaled.size(); ++k)
        if (cap.Cdiag_scaled[k] > 0) keep.push_back(n + k);
    auto selectCols = [&](const CscMatrix& M, bool rowsToo) {
        std::vector<Eigen::Triplet<real_t>> t;
        std::vector<int> map(M.cols(), -1);
        for (std::size_t i = 0; i < keep.size(); ++i) map[keep[i]] = (int)i;
        for (Eigen::Index j = 0; j < M.outerSize(); ++j) {
            if (map[j] < 0) continue;
            for (CscMatrix::InnerIterator it(M, j); it; ++it) {
                if (it.value() == 0) continue;
                const int row = rowsToo ? map[it.row()] : (int)it.row();
                if (row >= 0) t.emplace_back(row, map[j], it.value());
            }
        }
        CscMatrix out(rowsToo ? (Eigen::Index)keep.size() : M.rows(), (Eigen::Index)keep.size());
        out.setFromTriplets(t.begin(), t.end());
        return out;
    };
    const CscMatrix Pe = selectCols(cap.P, true);
    const CscMatrix Ae = selectCols(cap.A, false);
    VectorXr qe(keep.size());
    for (std::size_t i = 0; i < keep.size(); ++i) qe[(Eigen::Index)i] = cap.q[keep[i]];
    auto relErr = [](const CscMatrix& x, const CscMatrix& y) {
        const CscMatrix d = x - y;
        real_t m = 0, s = 0;
        for (Eigen::Index j = 0; j < d.outerSize(); ++j)
            for (CscMatrix::InnerIterator it(d, j); it; ++it) m = std::max(m, std::abs(it.value()));
        for (Eigen::Index j = 0; j < y.outerSize(); ++j)
            for (CscMatrix::InnerIterator it(y, j); it; ++it) s = std::max(s, std::abs(it.value()));
        return m / std::max<real_t>(s, 1e-300);
    };
    check(Pe.rows() == cp.P.rows() && relErr(cp.P, Pe) <= 1e-15, label + ": P equals ConeProgramAssembler");
    check(Ae.rows() == cp.A.rows() && Ae.cols() == cp.A.cols() && relErr(cp.A, Ae) <= 1e-15, label + ": A equals ConeProgramAssembler");
    check(qe.size() == cp.q.size() && (qe - cp.q).lpNorm<Eigen::Infinity>() <= 1e-15 * std::max<real_t>(1, qe.lpNorm<Eigen::Infinity>()), label + ": q equals ConeProgramAssembler");
    // b: equal except the shift on the normal rows (the first row of each cone; the nonnegative rows of the
    // frictionless contacts have none).
    bool bOk = cap.b.size() == cp.b.size();
    const real_t tol = 1e-14 * std::max<real_t>(1, cap.b.lpNorm<Eigen::Infinity>());
    int row = 0;
    for (; bOk && row < cp.dims.zero + cp.dims.nonneg; ++row) bOk = std::abs(cap.b[row] - cp.b[row]) <= tol;
    for (std::size_t k = 0; bOk && k < cp.dims.soc.size(); ++k) {
        bOk = cap.b[row] - cp.b[row] >= -tol;  // the shift on the normal row
        for (int j = 1; bOk && j < cp.dims.soc[k]; ++j) bOk = std::abs(cap.b[row + j] - cp.b[row + j]) <= tol;
        row += cp.dims.soc[k];
    }
    check(bOk, label + ": b equals ConeProgramAssembler up to the frozen shift");
}
#endif

void testScene(const std::string& name, const Scene& scene, int step, bool friction = true) {
    struct Variant {
        config::SolverType solver;
        const char* label;
        real_t balanceTol;
    };
    std::vector<Variant> variants = {{config::SolverType::ProjectedJacobi, "pj", 1e-12}, {config::SolverType::Condensed, "condensed", 1e-10}};
#ifdef CARDILLO_HAVE_CLARABEL
    variants.push_back({config::SolverType::Clarabel, "clarabel", 1e-7});
#endif
    const auto dir = tempDir(name);
    for (const auto& var : variants) {
        const std::string label = name + "/" + var.label;
        config::Config cfg = baseConfig(var.solver);
        if (var.solver == config::SolverType::Clarabel) cfg.solver_criterion = config::TerminationCriterion::VelocityRms;
        if (var.solver == config::SolverType::Condensed) cfg.condensed_true_schur = true;
        cfg.friction_enable = friction;
        const std::string path = (dir / (std::string(var.label) + ".hdf5")).string();
        Capture cap = runAndCapture(scene, cfg, step, path);
        check(cap.captured && cap.rec.problem.nc() > 0, label + ": a step with contacts was captured");
        if (!cap.captured || cap.rec.problem.nc() == 0) continue;
        checkRoundTrip(label, cap);
        checkErrors(label, cap, var.balanceTol);
#ifdef FFCLIB_TEST_ENGINE_CONE_PROGRAM
        checkConeProgram(label, cap);
#endif
        std::cout << label << ": nv " << cap.rec.problem.nv() << ", contacts " << cap.rec.problem.nc() << " (frictionless " << (cap.rec.problem.mu.array() == 0).count() << "), bilateral "
                  << cap.rec.problem.nb() << " (compliant " << cap.rec.problem.C.nonZeros() << "), engine natural map " << cap.engine.contact << '\n';
    }
}

void testDeterminism() {
    const auto dir = tempDir("determinism");
    const config::Config cfg = baseConfig(config::SolverType::ProjectedJacobi);
    const Capture a = runAndCapture(jointScene, cfg, 30, (dir / "a.hdf5").string());
    const Capture b = runAndCapture(jointScene, cfg, 30, (dir / "b.hdf5").string());
    const ffclib::Record ra = ffclib::read(a.path), rb = ffclib::read(b.path);
    check(sameSparse(ra.problem.H, rb.problem.H) && sameVec(ra.problem.f, rb.problem.f) && sameVec(ra.problem.w, rb.problem.w) && sameVec(ra.solution->r, rb.solution->r),
          "determinism: two runs export identical arrays");
}

void testExporter() {
    const auto dir = tempDir("exporter");
    config::Config cfg = baseConfig(config::SolverType::ProjectedJacobi);
    cfg.ffclib_folder = dir.string();
    cfg.ffclib_interval_steps = 3;
    cfg.ffclib_first_step = 2;
    cfg.ffclib_last_step = 8;
    cfg.ffclib_steps = {3, 11};
    cfg.output_filename_prefix = "ball";
    physics::PhysicsEngine engine(cfg);
    ballScene(engine);
    for (int k = 0; k < 13; ++k) engine.step();
    const std::vector<int> expected = {2, 3, 5, 8, 11};
    for (int s : expected) {
        char name[64];
        std::snprintf(name, sizeof name, "ball_%06d.hdf5", s);
        check(std::filesystem::exists(dir / name), std::string("exporter wrote ") + name);
    }
    int files = 0;
    for (const auto& entry : std::filesystem::directory_iterator(dir)) files += entry.path().extension() == ".hdf5";
    check(files == (int)expected.size(), "exporter wrote only the selected steps");
    std::ifstream manifest(dir / "manifest.csv");
    int lines = 0;
    for (std::string line; std::getline(manifest, line);) ++lines;
    check(lines == 1 + (int)expected.size(), "manifest has a header and one line per file");
    const ffclib::Record r3 = ffclib::read((dir / "ball_000003.hdf5").string());
    check(r3.info.title == "test step 3" && r3.info.description.find("h = 0.001 s") != std::string::npos && r3.solution && r3.guesses.size() == 1,
          "title, description, solution and guess of an exported step");
}

}  // namespace

int main() {
    testScene("ball", ballScene, 5);
    testScene("boxes", boxScene, 0);
    testScene("boxes_frictionless", boxScene, 0, false);
    testScene("joints", jointScene, 25);
    testDeterminism();
    testExporter();
    if (failures > 0) {
        std::cerr << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "test_ffclib_export: all checks passed\n";
    return 0;
}
