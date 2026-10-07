#include "mosek_solver.hpp"

#include <iostream>
#include <stdexcept>

namespace cardillo::solver {

namespace {

mosek::Options makeOptions(const config::Config& cfg) {
    mosek::Options o;
    o.num_threads = cfg.mosek_num_threads;
    o.tol_rel_gap = cfg.mosek_tol_rel_gap;
    o.tol_pfeas = cfg.mosek_tol_pfeas;
    o.tol_dfeas = cfg.mosek_tol_dfeas;
    o.max_iterations = cfg.mosek_max_iterations;
    o.max_time = cfg.mosek_max_time;
    if (cfg.mosek_presolve == "off")
        o.presolve = mosek::Options::Presolve::Off;
    else if (cfg.mosek_presolve == "on")
        o.presolve = mosek::Options::Presolve::On;
    else if (cfg.mosek_presolve == "free")
        o.presolve = mosek::Options::Presolve::Free;
    o.scaling = cfg.mosek_scaling;
    o.verbose = cfg.debug_pj;
    o.equality_duals = false;  // Cardillo reads only x and the contact (cone) duals
    return o;
}

}  // namespace

MosekSolver::MosekSolver(physics::DynamicsAssembler& dyn, const config::Config& cfg) : m_dyn(dyn), m_assembler(m_dyn), m_cfg(cfg), m_solver(makeOptions(cfg)) {}

MosekSolver::~MosekSolver() = default;

VectorXr MosekSolver::solve(real_t dt, real_t theta) {
    const SparseMatrix<Eigen::ColMajor>* P = nullptr;
    VectorXr* q = nullptr;
    const SparseMatrix<Eigen::ColMajor>* A = nullptr;
    VectorXr* b = nullptr;
    {
        auto sc = m_dyn.timings()->scope(misc::TimingManager::TimerId::MosekAssembly);
        P = &m_assembler.P(dt, theta);
        q = &m_assembler.q(dt, theta);
        A = &m_assembler.A(dt, theta);
        b = &m_assembler.b(dt, theta);
        m_dims.zero = static_cast<int>(m_dyn.numSprings() + m_dyn.numDampers());
        m_dims.nonneg = static_cast<int>(m_dyn.numFrictionlessContacts());
        m_dims.soc.assign(static_cast<std::size_t>(m_dyn.numFrictionalContacts()), 3);
    }

    if (m_solver.structureMatches(*P, *A, m_dims)) {
        auto sc = m_dyn.timings()->scope(misc::TimingManager::TimerId::MosekUpdate);
        m_solver.update(*P, *q, *A, *b);
    } else {
        auto sc = m_dyn.timings()->scope(misc::TimingManager::TimerId::MosekSetup);
        m_solver.setup(*P, *q, *A, *b, m_dims);
        if (m_first_setup) {
            const auto& o = m_solver.options();
            std::cout << "Initializing MOSEK solver with settings:\n";
            std::cout << "  Threads: " << (o.num_threads > 0 ? std::to_string(o.num_threads) : std::string("auto")) << "\n";
            std::cout << "  Tolerances (rel gap / pfeas / dfeas): " << (o.tol_rel_gap > 0 ? std::to_string(o.tol_rel_gap) : std::string("default")) << " / "
                      << (o.tol_pfeas > 0 ? std::to_string(o.tol_pfeas) : std::string("default")) << " / " << (o.tol_dfeas > 0 ? std::to_string(o.tol_dfeas) : std::string("default"))
                      << "\n";
            std::cout << "  Cones: zero=" << m_dims.zero << " nonneg=" << m_dims.nonneg << " soc=" << m_dims.soc.size() << "\n";
            m_first_setup = false;
        }
    }

    mosek::Status status;
    {
        auto sc = m_dyn.timings()->scope(misc::TimingManager::TimerId::MosekSolve);
        status = m_solver.optimize();
    }
    m_last_iters = m_solver.iterations();

    if (status != mosek::Status::Optimal) {
        std::cerr << "\nError solving MOSEK problem: " << mosek::toString(status) << "\n";
        std::cerr << "  Iterations: " << m_solver.iterations() << "\n";
        std::cerr << "  Termination code: " << static_cast<int>(m_solver.terminationCode()) << "\n";
        if (!mosek::isAcceptable(status)) throw std::runtime_error(std::string("MOSEK solver failed to solve the problem: ") + mosek::toString(status));
    }

    auto sc = m_dyn.timings()->scope(misc::TimingManager::TimerId::MosekExtract);
    m_solver.extract();

    const auto z = m_solver.coneDuals();
    const auto x = m_solver.x();
    VectorXr Smu = m_assembler.computeSmu();
    VectorXr impulse = z.cwiseProduct(Smu);
    solver::WarmstartProvider::storeImpulse(impulse, m_dyn);

    const int nV = static_cast<int>(m_dyn.numV());
    const int nSprings = static_cast<int>(m_dyn.numSprings());
    const int nDampers = static_cast<int>(m_dyn.numDampers());
    if (nSprings > 0) m_dyn.setLambda_g(x.segment(nV, nSprings));
    if (nDampers > 0) m_dyn.setLambda_gamma(x.segment(nV + nSprings, nDampers));

    return x.head(nV);
}

}  // namespace cardillo::solver
