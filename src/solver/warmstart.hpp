#pragma once

#include <vector>
#include <optional>
#include "misc/types.hpp"
#include "../collision/types.hpp"

namespace cardillo::solver {

// Simple per-contact impulse container (normal + two tangents)
struct ContactImpulse {
    real_t pn{0};  // normal impulse
    real_t pt1{0}; // tangent-1 impulse
    real_t pt2{0}; // tangent-2 impulse
};

// WarmstartProvider: abstract interface for warmstart strategies.
class WarmstartProvider {
public:
    WarmstartProvider() = default;
    virtual ~WarmstartProvider() = default;

    // Store last-frame impulses aligned with last-frame contacts (flattened order).
    virtual void store(const std::vector<collision::Contact>& lastContacts,
                       const std::vector<ContactImpulse>& lastImpulses) = 0;

    // Build a per-contact warmstart hint for current contacts using prev_global_out_index mapping.
    virtual std::vector<ContactImpulse> makeHint(const std::vector<collision::Contact>& currContacts) const = 0;

    // Direct query
    virtual std::optional<ContactImpulse> get(int prev_global_index) const = 0;

    virtual void clear() = 0;
};

// A tiny cache to support warmstarting across time steps and solvers.
// Stores last frame's per-contact impulses in an array indexed by the
// previous frame contacts' global_out_index. To generate a hint for the
// current frame, use each contact.prev_global_out_index to fetch.
class WarmstartCache : public WarmstartProvider {
public:
    WarmstartCache() = default;

    // Store last-frame impulses aligned with last-frame contacts (flattened order).
    // Each contact must have its global_out_index set (as produced by the collision pipeline).
    void store(const std::vector<collision::Contact>& lastContacts,
               const std::vector<ContactImpulse>& lastImpulses) override {
        if (lastContacts.size() != lastImpulses.size()) return;
        // Determine required size
        std::size_t maxIdx = 0;
        for (const auto& c : lastContacts) {
            if (c.global_out_index >= 0) maxIdx = std::max<std::size_t>(maxIdx, (std::size_t)c.global_out_index);
        }
        m_last.resize(maxIdx + 1, ContactImpulse{});
        // Write by global_out_index
        for (std::size_t i = 0; i < lastContacts.size(); ++i) {
            int idx = lastContacts[i].global_out_index;
            if (idx >= 0) m_last[(std::size_t)idx] = lastImpulses[i];
        }
    }

    // Build a per-contact warmstart hint for current contacts using prev_global_out_index mapping.
    // If no previous mapping exists, the hint entry is zero-initialized.
    std::vector<ContactImpulse> makeHint(const std::vector<collision::Contact>& currContacts) const override {
        std::vector<ContactImpulse> out(currContacts.size());
        for (std::size_t i = 0; i < currContacts.size(); ++i) {
            int prevIdx = currContacts[i].prev_global_out_index;
            if (prevIdx >= 0 && (std::size_t)prevIdx < m_last.size()) {
                out[i] = m_last[(std::size_t)prevIdx];
            } // else leaves zeroes
        }
        return out;
    }

    // Direct query
    std::optional<ContactImpulse> get(int prev_global_index) const override {
        if (prev_global_index >= 0 && (std::size_t)prev_global_index < m_last.size())
            return m_last[(std::size_t)prev_global_index];
        return std::nullopt;
    }

    void clear() override { m_last.clear(); }

private:
    std::vector<ContactImpulse> m_last;
};

} // namespace cardillo::solver
