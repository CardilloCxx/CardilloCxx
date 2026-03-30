#ifndef CONTACTTRACKER_H
#define CONTACTTRACKER_H

#pragma once

#include "types.hpp"
#include "../config/config.hpp"
#include "../misc/timings/TimingManager.hpp"

namespace cardillo { namespace physics { class DynamicsAssembler; } }

namespace cardillo::collision {
    class ContactTracker {
    public:
        ContactTracker(config::Config& cfg, cardillo::misc::TimingManager* timings) : m_cfg(cfg), m_timings(timings) {}
        ~ContactTracker();

        void registerNextContacts(ContactMap& curr);
        void applyPrevImpulses(std::vector<Contact>& currFlattened,  const std::vector<Contact>& prevFlattened);

    private:
        config::Config& m_cfg;
        cardillo::misc::TimingManager* m_timings;

        mutable ContactMap m_prevContacts;
        static void matchContactsForPair(const ContactList& prev, ContactList& curr, real_t maxDist);

    };
}

#endif