#pragma once

#include <optional>
#include <vector>
#include "../../collision/types.hpp"
#include "../assembly/dynamics_assembler.hpp"
#include "misc/types.hpp"

namespace cardillo::solver {

// WarmstartProvider: abstract interface for warmstart strategies.
class WarmstartProvider {
   public:
    WarmstartProvider() = default;
    ~WarmstartProvider() = default;

    // Apply a warmstart hint using per-contact stored impulses in `contactsAll`.
    //
    // Note: the normal-impulse component is seeded as-is, without clamping to a sign here.
    // PGS uses the convention normal impulse <= 0, while Projected Jacobi uses >= 0; each
    // solver's own projection step re-enforces its convention on the very first iteration, so
    // clamping here would only ever discard a correctly-signed warmstart for one of the two
    // solvers instead of improving convergence.
    static void applyWarmstart(cardillo::VectorXr& p, cardillo::physics::DynamicsAssembler& dyn) {
        const auto& contacts = dyn.contacts();
        for (const auto& c : contacts) {
            if (c.impulse_base_index < 0) continue;
            const int base = c.impulse_base_index;
            if (base >= p.size()) continue;

            p[base] = c.last_impulse(0);
            if (c.impulse_size > 1 && base + 1 < p.size()) p[base + 1] = c.last_impulse(1);
            if (c.impulse_size > 2 && base + 2 < p.size()) p[base + 2] = c.last_impulse(2);
        }
    }

    // Extract per-contact impulses and store them into contactsAll (by reference to underlying
    // contacts)
    static void storeImpulse(const cardillo::VectorXr& p, cardillo::physics::DynamicsAssembler& dyn) {
        // Iterate contacts and write last_impulse directly into the assembler's contact storage
        const auto& contacts = dyn.contacts();
        for (const auto& c : contacts) {
            if (c.impulse_base_index < 0) continue;
            const int base = c.impulse_base_index;
            if (base >= p.size()) continue;
            cardillo::Vector3r imp = cardillo::Vector3r::Zero();
            imp(0) = p[base];
            if (c.impulse_size > 1 && base + 1 < p.size()) imp(1) = p[base + 1];
            if (c.impulse_size > 2 && base + 2 < p.size()) imp(2) = p[base + 2];

            // store into assembler via API
            dyn.setContactLastImpulse(c.global_out_index, imp);
        }
    }
};

}  // namespace cardillo::solver
