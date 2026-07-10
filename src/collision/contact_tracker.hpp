#ifndef CONTACTTRACKER_H
#define CONTACTTRACKER_H

#pragma once

#include "../config/config.hpp"
#include "../misc/timings/TimingManager.hpp"
#include "types.hpp"

namespace cardillo {
namespace physics {
class DynamicsAssembler;
}
}  // namespace cardillo

namespace cardillo::collision {
class ContactTracker {
   public:
    ContactTracker(config::Config& cfg, cardillo::misc::TimingManager* timings) : m_cfg(cfg), m_timings(timings) {}
    ~ContactTracker();

    void registerNextContacts(ContactMap& curr);

    /**
     * @brief Propagates warm-start impulses from the previous step's flattened contact buffer to
     * this step's flattened buffer.
     * @param currFlattened This step's flattened contacts (Contact::global_out_index indexes into
     * this vector); receives the propagated impulses.
     * @param prevFlattened The previous step's flattened, post-solve contacts
     * (Contact::last_impulse is authoritative there).
     * @param currMatched This step's ContactMap, after registerNextContacts() has run on it, so
     * each entry's `prev_global_out_index` is populated. Passed explicitly (rather than read from
     * internal state) so this function's correctness does not depend on being called immediately
     * after registerNextContacts() in the same detectAll() pass.
     */
    void applyPrevImpulses(std::vector<Contact>& currFlattened, const std::vector<Contact>& prevFlattened, const ContactMap& currMatched);

    /**
     * @brief Moves @p curr into the previous-contacts cache used for next step's matching.
     * @warning Must be called only after applyPrevImpulses() has finished reading @p curr (which
     * registerNextContacts() populated with matched prev_global_out_index fields) -- once moved
     * from, @p curr is left in an unspecified (likely empty) state. Kept as a separate call rather
     * than folded into registerNextContacts() so that step, unlike this one, never mutates its
     * argument's contents, and so the (relatively large, ~350B/contact) buffer is moved rather than
     * deep-copied once nothing else needs the old `curr`.
     */
    void commitContacts(ContactMap&& curr) { m_prevContacts = std::move(curr); }

   private:
    config::Config& m_cfg;
    cardillo::misc::TimingManager* m_timings;

    mutable ContactMap m_prevContacts;
    static void matchContactsForPair(const ContactList& prev, ContactList& curr, real_t maxDist);
};
}  // namespace cardillo::collision

#endif