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

    // Canonical sign convention for `Contact::last_impulse`'s normal (index 0) component:
    // non-negative, i.e. its magnitude is the compressive impulse. This matches Projected Jacobi,
    // QOCO and Clarabel (whose cone-program duals are non-negative by construction), and is what
    // VTK contact-pressure/force-vector output already assumes when consuming `last_impulse`
    // (see io/vtk_writer.cpp). Projected Gauss-Seidel is the sole exception: its internal
    // `project()` clamps the normal impulse <= 0. `invertNormalSign` lets PGS convert at this
    // boundary without touching its own internal algebra or the canonical storage convention.
    //
    // Apply a warmstart hint using per-contact stored impulses in `contactsAll`.
    static void applyWarmstart(VectorXr& p, physics::DynamicsAssembler& dyn, bool invertNormalSign = false) {
        const auto& contacts = dyn.contacts();
        for (const auto& c : contacts) {
            if (c.impulse_base_index < 0) continue;
            const int base = c.impulse_base_index;
            if (base >= p.size()) continue;

            p[base] = invertNormalSign ? -c.last_impulse(0) : c.last_impulse(0);
            if (c.impulse_size > 1 && base + 1 < p.size()) p[base + 1] = c.last_impulse(1);
            if (c.impulse_size > 2 && base + 2 < p.size()) p[base + 2] = c.last_impulse(2);
        }
    }

    // Extract per-contact impulses and store them into contactsAll (by reference to underlying
    // contacts)
    static void storeImpulse(const VectorXr& p, physics::DynamicsAssembler& dyn, bool invertNormalSign = false) {
        // Iterate contacts and write last_impulse directly into the assembler's contact storage
        const auto& contacts = dyn.contacts();
        for (const auto& c : contacts) {
            if (c.impulse_base_index < 0) continue;
            const int base = c.impulse_base_index;
            if (base >= p.size()) continue;
            Vector3r imp = Vector3r::Zero();
            imp(0) = invertNormalSign ? -p[base] : p[base];
            if (c.impulse_size > 1 && base + 1 < p.size()) imp(1) = p[base + 1];
            if (c.impulse_size > 2 && base + 2 < p.size()) imp(2) = p[base + 2];

            // store into assembler via API
            dyn.setContactLastImpulse(c.global_out_index, imp);
        }
    }
};

}  // namespace cardillo::solver
