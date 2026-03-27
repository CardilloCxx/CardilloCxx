#pragma once

#include <vector>
#include <optional>
#include "misc/types.hpp"
#include "../../collision/types.hpp"

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
    // Default implementation writes into an internal cache indexed by
    // contact.global_out_index. Providers can override for custom storage.
    virtual void store(const std::vector<collision::Contact>& lastContacts,
                       const std::vector<ContactImpulse>& lastImpulses)
    {
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

    // Apply a warmstart hint directly into a solver row-vector `p` using
    // contact row grouping and dynamic->original mapping supplied by the
    // solver. Default implementation delegates to `makeHint()` and
    // expands per-contact impulses into `p`. Concrete providers may
    // override this if they have a different storage layout.
    virtual void applyHintToRowVector(cardillo::VectorXr& p,
                                      const std::vector<std::vector<int>>& contactRowGroups,
                                      const std::vector<int>& dynamicContactToOriginalAll,
                                      const std::vector<collision::Contact>& contactsAll) const
    {
        const std::vector<ContactImpulse> hints = makeHint(contactsAll);
        for (const auto& rows : contactRowGroups) {
            if (rows.empty()) continue;
            const int rn = rows[0];
            const int origIdx = (rn >= 0 && rn < (int)dynamicContactToOriginalAll.size()) ? dynamicContactToOriginalAll[(size_t)rn] : -1;
            if (origIdx < 0 || origIdx >= (int)hints.size()) continue;
            const ContactImpulse& hi = hints[(size_t)origIdx];
            p[rn] = std::max<real_t>(hi.pn, (real_t)0);
            if (rows.size() > 1) p[rows[1]] = hi.pt1;
            if (rows.size() > 2) p[rows[2]] = hi.pt2;
        }
    }

    // Extract per-contact impulses from a solver row-vector `p` and store
    // them in the provider's cache for the next step. Default implementation
    // assembles a per-contact vector and calls `store()`; providers may
    // override to avoid intermediate allocations.
    virtual void storeFromRowVector(const cardillo::VectorXr& p,
                                    const std::vector<std::vector<int>>& contactRowGroups,
                                    const std::vector<int>& dynamicContactToOriginalAll,
                                    const std::vector<collision::Contact>& contactsAll)
    {
        std::vector<ContactImpulse> lastImp(contactsAll.size());
        for (const auto& rows : contactRowGroups) {
            if (rows.empty()) continue;
            const int rn = rows[0];
            const int origIdx = (rn >= 0 && rn < (int)dynamicContactToOriginalAll.size()) ? dynamicContactToOriginalAll[(size_t)rn] : -1;
            if (origIdx < 0 || origIdx >= (int)lastImp.size()) continue;
            ContactImpulse ci{};
            ci.pn = p[rn];
            if (rows.size() > 1) ci.pt1 = p[rows[1]]; else ci.pt1 = (real_t)0;
            if (rows.size() > 2) ci.pt2 = p[rows[2]]; else ci.pt2 = (real_t)0;
            lastImp[(size_t)origIdx] = ci;
        }
        store(contactsAll, lastImp);
    }

    // Build a per-contact warmstart hint for current contacts using prev_global_out_index mapping.
    virtual std::vector<ContactImpulse> makeHint(const std::vector<collision::Contact>& currContacts) const
    {
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
    virtual std::optional<ContactImpulse> get(int prev_global_index) const
    {
        if (prev_global_index >= 0 && (std::size_t)prev_global_index < m_last.size())
            return m_last[(std::size_t)prev_global_index];
        return std::nullopt;
    }

    virtual void clear() { m_last.clear(); }

protected:
    // internal cache for default provider implementation
    std::vector<ContactImpulse> m_last;
};

} // namespace cardillo::solver
