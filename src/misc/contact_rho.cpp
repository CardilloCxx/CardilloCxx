#include "contact_rho.hpp"

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>

namespace cardillo::misc {

Vector3r computeContactRho(const Matrix33r& G, ContactRhoStrategy strategy, real_t eps) {
    if (strategy == ContactRhoStrategy::FullSpectral) {
        // see Acary2018 equation (110): https://doi.org/10.1007/978-3-319-75972-2_10
        // G need not be exactly symmetric (an implicit-gyroscopic body makes Minv, and hence Gii,
        // generally non-symmetric) -- symmetrizing first, rather than using a general
        // (possibly-complex-eigenvalue) solver, mirrors what Siconos's own closed-form routine
        // assumes and keeps this a real scalar step size in every case.
        const Matrix33r Gsym = (real_t)0.5 * (G + G.transpose());
        Eigen::SelfAdjointEigenSolver<Matrix33r> es(Gsym, Eigen::EigenvaluesOnly);
        const real_t lambdaMax = es.eigenvalues()(2);  // ascending order -- largest is last
        const real_t rhoAll = (lambdaMax > eps) ? (real_t)1 / lambdaMax : (real_t)-1;
        return Vector3r(rhoAll, rhoAll, rhoAll);
    }

    // see Acary2018 equation (111): https://doi.org/10.1007/978-3-319-75972-2_10
    if (G(0, 0) <= eps) return Vector3r((real_t)-1, (real_t)-1, (real_t)-1);
    const real_t rhoN = (real_t)1 / G(0, 0);

    // Closed-form (quadratic-formula) largest eigenvalue of the symmetrized 2x2 tangential block
    // -- already exact for a symmetric input, so routing this through a generic
    // SelfAdjointEigenSolver would add setup overhead for no accuracy/robustness gain. `b` is
    // averaged from both off-diagonal entries (not just G(1,2)) for the same non-symmetry reason
    // as above -- using only G(1,2) would silently ignore G(2,1).
    const real_t a = G(1, 1), d = G(2, 2), b = (real_t)0.5 * (G(1, 2) + G(2, 1));
    const real_t disc = std::max((real_t)0, (a - d) * (a - d) + (real_t)4 * b * b);
    const real_t lambdaMaxTT = (real_t)0.5 * ((a + d) + std::sqrt(disc));
    const real_t rhoT = (lambdaMaxTT > eps) ? (real_t)1 / lambdaMaxTT : (real_t)-1;

    return Vector3r(rhoN, rhoT, rhoT);
}

}  // namespace cardillo::misc
