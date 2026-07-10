#include "contact_tracker.hpp"

#include <algorithm>
#include <utility>
#include <vector>

namespace cardillo::collision {

ContactTracker::~ContactTracker() {}

void ContactTracker::applyPrevImpulses(std::vector<Contact>& currFlattened, const std::vector<Contact>& prevFlattened, const ContactMap& currMatched) {
    for (const auto& [pairKey, list] : currMatched) {
        for (const Contact& c : list) {
            if (c.global_out_index < 0 || (std::size_t)c.global_out_index >= currFlattened.size()) continue;
            if (c.prev_global_out_index >= 0 && (std::size_t)c.prev_global_out_index < prevFlattened.size()) {
                currFlattened[(std::size_t)c.global_out_index].last_impulse = prevFlattened[(std::size_t)c.prev_global_out_index].last_impulse;
            }
        }
    }
}

void ContactTracker::registerNextContacts(ContactMap& curr) {
    // Timings
    auto sc = m_timings->scope(cardillo::misc::TimingManager::TimerId::CollisionMatching);

    // For each pair in the current map, find the corresponding pair in the previous map and match
    // contacts
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
}

void ContactTracker::matchContactsForPair(const ContactList& prev, ContactList& curr, real_t maxDist) {
    if (curr.empty()) return;
    const real_t max2 = maxDist * maxDist;

    // Greedy nearest-neighbor matching, closest pairs first, each previous contact consumed at
    // most once. Without the exclusion, two current contacts could both match the same previous
    // contact (duplicating its warm-start impulse) while a third current contact that should have
    // matched it instead gets no warm start at all.
    std::vector<bool> prevUsed(prev.size(), false);
    std::vector<std::pair<real_t, std::pair<std::size_t, std::size_t>>> candidates;
    candidates.reserve(curr.size() * prev.size());
    for (std::size_t i = 0; i < curr.size(); ++i) {
        curr[i].prev_global_out_index = -1;
        for (std::size_t j = 0; j < prev.size(); ++j) {
            const real_t d2 = (curr[i].point - prev[j].point).squaredNorm();
            if (d2 <= max2) candidates.emplace_back(d2, std::make_pair(i, j));
        }
    }
    std::sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

    std::vector<bool> currMatched(curr.size(), false);
    for (const auto& [d2, ij] : candidates) {
        const auto [i, j] = ij;
        if (currMatched[i] || prevUsed[j]) continue;
        currMatched[i] = true;
        prevUsed[j] = true;
        curr[i].prev_global_out_index = prev[j].global_out_index;
        curr[i].last_impulse = prev[j].last_impulse;
    }
}

}  // namespace cardillo::collision
