#pragma once

#include <vector>
#include <optional>
#include "misc/types.hpp"
#include "../../collision/types.hpp"
#include "../assembly/dynamics_assembler.hpp"

namespace cardillo::solver {

// WarmstartProvider: abstract interface for warmstart strategies.
class WarmstartProvider {
public:
    WarmstartProvider() = default;
    virtual ~WarmstartProvider() = default;

    // Apply a warmstart hint using per-contact stored impulses in `contactsAll`.
    virtual void applyWarmstart(cardillo::VectorXr& p,
                                       cardillo::physics::DynamicsAssembler& dyn) const
    {   
        const auto& contacts  = dyn.contacts();
        for (const auto& c : contacts) {
            if (c.impulse_base_index < 0) continue;
            const int base = c.impulse_base_index;
            if (base >= p.size()) continue;

            // normal impulse (clamped >= 0)
            p[base] = std::max<real_t>(c.last_impulse(0), (real_t)0);
            if (c.impulse_size > 1 && base + 1 < p.size()) p[base + 1] = c.last_impulse(1);
            if (c.impulse_size > 2 && base + 2 < p.size()) p[base + 2] = c.last_impulse(2);
        }
    }

    // Extract per-contact impulses and store them into contactsAll (by reference to underlying contacts)
    virtual void storeImpulse(const cardillo::VectorXr& p,
                                    cardillo::physics::DynamicsAssembler& dyn)
    {
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

    virtual void clear() {}
};

} // namespace cardillo::solver
