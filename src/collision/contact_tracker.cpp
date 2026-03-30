#include "contact_tracker.hpp"

namespace cardillo::collision {

ContactTracker::~ContactTracker()
{

}

void ContactTracker::registerNextContacts(ContactMap& curr) const
{
    // Timings
    auto sc = m_timings->scope(cardillo::misc::TimingManager::TimerId::CollisionMatching);

    // For each pair in the current map, find the corresponding pair in the previous map and match contacts
    for (auto& [pairKey, currList] : curr) {
        auto it = m_prevContacts.find(pairKey);
        if (it != m_prevContacts.end()) {
            const ContactList& prevList = it->second;
            matchContactsForPair(prevList, currList, m_cfg.collision_match_max_dist);
        } else {
            // No previous contacts for this pair; mark all as unmatched
            for (Contact& c : currList) {
                c.prev_global_out_index = -1;
            }
        }
    }
    m_prevContacts = std::move(curr);
}


void ContactTracker::matchContactsForPair(const ContactList& prev, ContactList& curr, real_t maxDist) {
    if (curr.empty()) return;
    const real_t max2 = maxDist * maxDist;
    for (std::size_t i = 0; i < curr.size(); ++i) {
        int best = -1;
        real_t best2 = std::numeric_limits<real_t>::infinity();
        const Vector3r& p = curr[i].point;
        for (std::size_t j = 0; j < prev.size(); ++j) {
            const Vector3r d = p - prev[j].point;
            const real_t d2 = d.squaredNorm();
            if (d2 < best2) { best2 = d2; best = (int)j; }
        }
        if (best >= 0 && best2 <= max2) {
            curr[i].prev_global_out_index = prev[best].global_out_index;

            // Copy last impulse for warmstarting if matched
            curr[i].last_impulse = prev[best].last_impulse;
        } else {
            curr[i].prev_global_out_index = -1;
        }
    }
}

} // namespace cardillo::collision


