#pragma once

#include "../../misc/types.hpp"

namespace cardillo::physics {
class DynamicsAssembler;
}

namespace cardillo::solver {

// Contact law of a solver's one-step problem. The conic backends solve the convexified law, in
// which the De Saxce shift mu ||xi_T|| is frozen at the pre-step slip mu ||gamma_T^-||
// (conic/conic_scaling.hpp); all other solvers solve Coulomb's law.
enum class ContactLaw { Coulomb, Convexified };

// Natural map of the contact law, the accuracy measure of the final check after every step, the
// same for every solver (solver.final_check = natural_map). For contact b with percussion lambda_b and
// velocity xi_b (rows of DynamicsAssembler::W(), restitution included):
//
//   r_b = (lambda_b - proj_K(lambda_b - rho_b (xi_b + s_b e_N))) / rho_b,
//
// with K = { ||lambda_T|| <= mu lambda_N } (frictionless contacts: lambda_N >= 0), the De Saxce
// shift s_b = mu ||xi_T|| or a given (frozen) shift, and one step rho_b = 1 / lambda_max(G_b) per
// contact from its Delassus block G_b = W_b M^-1 W_b^T with the diagonal mass. r_b = 0 iff the
// contact satisfies its law; for Coulomb's law this zero set is the one of the Alart-Curnier
// projection that the projected solvers iterate. r_b is a velocity, for a sticking contact the
// velocity error itself, and does not depend on the masses. Contacts with G_b = 0 (no dynamic body)
// are skipped. The projection solvers stop on the natural map with their own step instead
// (misc::naturalMapError2), which they get with every sweep.
class ContactNaturalMap {
   public:
    // W: contact rows (frictionless rows first, then (n, t1, t2) per frictional contact).
    ContactNaturalMap(const CsrMatrix& W, const VectorXr& MinvDiag, const VectorXr& mu, int numFrictionlessRows);

    // max_b max(|r_N|, ||r_T||), NaN if the input is not finite.
    // Coulomb's law: De Saxce shift mu ||xi_T||.
    real_t residual(const VectorXr& lambda, const VectorXr& xi) const;
    // Convexified law: the given (frozen) shift per row, read on the normal rows of frictional contacts.
    real_t residual(const VectorXr& lambda, const VectorXr& xi, const VectorXr& frozenShift) const;

   private:
    // Both laws; frozenShift == nullptr: Coulomb's law.
    real_t maxResidual(const VectorXr& lambda, const VectorXr& xi, const VectorXr* frozenShift) const;

    int m_numFrictionlessRows{0};
    VectorXr m_mu;      // per row
    VectorXr m_rho;     // per row, the same on the rows of a contact; 0: contact skipped
    VectorXr m_invRho;  // 1 / m_rho (lambda_max(G_b)), 0 for skipped contacts
};

// Final check of a step, the same for every solver (solver.final_check = natural_map): the natural map
// of the stored contact percussions (DynamicsAssembler::contactPercussions()) at the velocities u
// returned by the solver, and the residual W_i^T u + b_i + c_i Lambda_i of the spring and damper
// rows (DynamicsAssembler::springBias()). All values are raw residuals in velocity units (m/s or
// rad/s), independent of the solver's tolerances. Momentum balance between u and the percussions
// is not checked: all solvers but the conic backends compute u from the percussions. Must be called
// after the solve and before u is written to the system.
struct StepResidual {
    real_t contact{0};    // natural map of the solver's contact law
    real_t coulomb{0};    // natural map of Coulomb's law; differs from `contact` for the convexified law only
    real_t bilateral{0};  // largest spring/damper row residual
};
StepResidual stepResidual(physics::DynamicsAssembler& dyn, const VectorXr& u, ContactLaw law, real_t dt, real_t theta);

}  // namespace cardillo::solver
