#include "moreau.hpp"
#include "../solver/projected_jacobi.hpp"
#include <vector>
#include <algorithm>
#include <limits>
#include <cstdlib>
#include <iostream>

namespace cardillo::integration {

void MoreauSolver::step(real_t dt)
{
    // 1) Get current state vectors
    m_dyn.refreshState();
    // const auto& qn = m_dyn.qVec();
    const auto& vn = m_dyn.vVec();
    const auto& fn_ext = m_dyn.fVecExternal();    // gravity + applied external forces
    const auto& fn_gyro = m_dyn.fVecGyroscopic(); // gyroscopic forces from current state

    const bool implicitGyro = m_config.moreau_implicit_gyroscopy;
    const bool lambdaTheta = m_config.moreau_lambda_theta;

    // 2) Inplace midpoint position update
    explicitPositionUpdate(m_world, (1.0 - m_theta) * dt);

    // 3) Evaluate contacts at midpoint positions
    m_dyn.refreshCollisionsAndSprings(dt, m_theta, implicitGyro, lambdaTheta);

    // 4) Build the full extended RHS and solve the extended system S * x = b
    const auto& Wg = m_dyn.Wg().asSparse();
    const auto& Wgamma = m_dyn.Wgamma().asSparse();
    const auto& M_diag = m_dyn.MDiag();
    const int totalV = (m_dyn.bodyVelOffsets().empty() ? 0 : m_dyn.bodyVelOffsets().back());
    const int nSprings = (int)m_dyn.Cdiag().size();
    const int nDampers = (int)m_dyn.Adiag().size();
    const int extV = totalV + nSprings + nDampers;
    const auto& Cdiag = m_dyn.Cdiag();

    // Lambda vectors may be uninitialized on first step; copy locally and ensure correct sizes
    VectorXr Lambda_g = m_dyn.Lambda_g();
    if ((int)Lambda_g.size() != nSprings) Lambda_g = VectorXr::Zero(nSprings);
    VectorXr Lambda_gamma = m_dyn.Lambda_gamma();
    if ((int)Lambda_gamma.size() != nDampers) Lambda_gamma = VectorXr::Zero(nDampers);

    // RHS: M*vn + dt*f_ext (+ dt*f_gyro if treated explicitly)
    VectorXr rhs = VectorXr::Zero((index_t)extV);
    rhs.segment(0, totalV) = M_diag.cwiseProduct(vn) + dt * fn_ext;
    if (!implicitGyro) rhs.segment(0, totalV) += dt * fn_gyro;
    if (lambdaTheta && (nSprings > 0 || nDampers > 0)) {
        VectorXr corr = VectorXr::Zero(totalV);
        if (nSprings > 0) corr.noalias() += Wg.transpose() * Lambda_g;
        if (nDampers > 0) corr.noalias() += Wgamma.transpose() * Lambda_gamma;
        rhs.segment(0, totalV).noalias() -= (1.0 - m_theta) * corr;
    }
    if (nSprings > 0) rhs.segment(totalV, nSprings) = -(1.0 / (m_theta * dt * dt)) * Cdiag.cwiseProduct(Lambda_g) - (1.0 - m_theta) / m_theta * Wg * vn;
    if (nDampers > 0) rhs.segment(totalV + nSprings, nDampers) = -(1.0 - m_theta) / m_theta * Wgamma * vn;

    // 5) Solve extended system
    auto xnp1 = m_solver.solve(rhs, m_config.pj_tol_abs);
    // auto xnp1 = m_dyn.solveS(rhs);                             // No contacts
    VectorXr vnp1 = xnp1.segment(0, totalV);
    if (nSprings > 0) m_dyn.setLambda_g(xnp1.segment(totalV, nSprings)); else m_dyn.setLambda_g(VectorXr(0));
    if (nDampers > 0) m_dyn.setLambda_gamma(xnp1.segment(totalV + nSprings, nDampers)); else m_dyn.setLambda_gamma(VectorXr(0));

    // 6) Write final velocity to ECS
    m_dyn.writeVelocityToSystem(vnp1);

    // 7) Inplace final position update
    explicitPositionUpdate(m_world, m_theta * dt);
}

}