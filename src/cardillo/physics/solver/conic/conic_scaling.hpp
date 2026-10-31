#pragma once

#include <cmath>

#include "../../../misc/types.hpp"
#include "../../assembly/dynamics_assembler.hpp"

namespace cardillo::physics::assembly {

// Row scaling S of the contact rows for the conic (interior-point) backends. Each frictional contact
// is scaled by S = diag(1, mu, mu) (Acary et al.), which maps the Coulomb cone
// {||lambda_T|| <= mu lambda_N} onto the unit second-order cone of the scaled variable z with
// lambda = S z. Frictionless rows keep the scale 1.
inline VectorXr frictionConeScaling(const DynamicsAssembler& dyn) {
    VectorXr S = VectorXr::Ones(dyn.numContactRows());
    const VectorXr& mu = dyn.muVec();
    for (int i = dyn.numFrictionlessContacts(); i < dyn.numContactRows(); i += 3) {
        S[i + 1] = mu[i];
        S[i + 2] = mu[i];
    }
    return S;
}

// Contact part of the conic offset: S * biasImpulse plus the De Saxce shift on the normal rows. The
// slack is then s = (xi_N + mu ||xi_T(v^n)||, mu xi_T), and s in the unit second-order cone is the
// dual Coulomb cone {mu ||xi_T|| <= xi_N} shifted by mu ||xi_T(v^n)||.
inline VectorXr scaledContactOffset(const DynamicsAssembler& dyn, const VectorXr& S, const VectorXr& biasImpulse, const VectorXr& gamma_old) {
    VectorXr h = S.cwiseProduct(biasImpulse);
    const VectorXr& mu = dyn.muVec();
    for (int i = dyn.numFrictionlessContacts(); i < dyn.numContactRows(); i += 3) {
        const real_t y_1 = gamma_old[i + 1];
        const real_t y_2 = gamma_old[i + 2];
        h[i] += mu[i] * std::sqrt(y_1 * y_1 + y_2 * y_2);
    }
    return h;
}

}  // namespace cardillo::physics::assembly
