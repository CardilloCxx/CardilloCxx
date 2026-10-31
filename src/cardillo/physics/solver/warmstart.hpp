#pragma once

#include <optional>
#include <vector>
#include "../../collision/types.hpp"
#include "../assembly/dynamics_assembler.hpp"
#include "cardillo/misc/types.hpp"

namespace cardillo::solver {

// WarmstartProvider: abstract interface for warmstart strategies.
class WarmstartProvider {
   public:
    WarmstartProvider() = default;
    ~WarmstartProvider() = default;

    // Canonical convention for `Contact::last_impulse` = (pn, pt1, pt2): the percussion acting on
    // body B of the contact is pn * normal + pt1 * tangent1 + pt2 * tangent2 (normal from A to B),
    // so pn >= 0 is the compressive impulse. Every contact solver (Projected Jacobi, Projected
    // Gauss-Seidel, Condensed, QOCO, Clarabel, ConicXX) works with this convention, and it is what
    // the VTK contact output assumes (see io/vtk_writer.cpp).
    //
    // Apply a warmstart hint using per-contact stored impulses in `contactsAll`.
    static void applyWarmstart(VectorXr& p, physics::DynamicsAssembler& dyn) {
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
    static void storeImpulse(const VectorXr& p, physics::DynamicsAssembler& dyn) {
        // Iterate contacts and write last_impulse directly into the assembler's contact storage
        const auto& contacts = dyn.contacts();
        for (const auto& c : contacts) {
            if (c.impulse_base_index < 0) continue;
            const int base = c.impulse_base_index;
            if (base >= p.size()) continue;
            Vector3r imp = Vector3r::Zero();
            imp(0) = p[base];
            if (c.impulse_size > 1 && base + 1 < p.size()) imp(1) = p[base + 1];
            if (c.impulse_size > 2 && base + 2 < p.size()) imp(2) = p[base + 2];

            // store into assembler via API
            dyn.setContactLastImpulse(c.global_out_index, imp);
        }
    }
};

}  // namespace cardillo::solver
