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

    // Canonical convention for `Contact::last_impulse` = (pn, pt1, pt2): the percussion acting on
    // body B of the contact is pn * normal + pt1 * tangent1 + pt2 * tangent2 (normal from A to B),
    // so pn >= 0 is the compressive impulse. This matches Projected Jacobi, QOCO, Clarabel and
    // ConicXX, and is what the VTK contact output assumes (see io/vtk_writer.cpp). Projected
    // Gauss-Seidel and the condensed solver work internally with the negated percussion
    // Lambda = -pi in all three components (normal impulse <= 0); `invertSign` converts the whole
    // vector at this boundary.
    //
    // Apply a warmstart hint using per-contact stored impulses in `contactsAll`.
    static void applyWarmstart(VectorXr& p, physics::DynamicsAssembler& dyn, bool invertSign = false) {
        const auto& contacts = dyn.contacts();
        for (const auto& c : contacts) {
            if (c.impulse_base_index < 0) continue;
            const int base = c.impulse_base_index;
            if (base >= p.size()) continue;

            const real_t sign = invertSign ? (real_t)-1 : (real_t)1;
            p[base] = sign * c.last_impulse(0);
            if (c.impulse_size > 1 && base + 1 < p.size()) p[base + 1] = sign * c.last_impulse(1);
            if (c.impulse_size > 2 && base + 2 < p.size()) p[base + 2] = sign * c.last_impulse(2);
        }
    }

    // Extract per-contact impulses and store them into contactsAll (by reference to underlying
    // contacts)
    static void storeImpulse(const VectorXr& p, physics::DynamicsAssembler& dyn, bool invertSign = false) {
        // Iterate contacts and write last_impulse directly into the assembler's contact storage
        const auto& contacts = dyn.contacts();
        for (const auto& c : contacts) {
            if (c.impulse_base_index < 0) continue;
            const int base = c.impulse_base_index;
            if (base >= p.size()) continue;
            const real_t sign = invertSign ? (real_t)-1 : (real_t)1;
            Vector3r imp = Vector3r::Zero();
            imp(0) = sign * p[base];
            if (c.impulse_size > 1 && base + 1 < p.size()) imp(1) = sign * p[base + 1];
            if (c.impulse_size > 2 && base + 2 < p.size()) imp(2) = sign * p[base + 2];

            // store into assembler via API
            dyn.setContactLastImpulse(c.global_out_index, imp);
        }
    }
};

}  // namespace cardillo::solver
