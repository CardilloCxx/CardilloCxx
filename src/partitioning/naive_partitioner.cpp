#include "naive_partitioner.hpp"
#include <algorithm>
#include <cmath>

namespace cardillo::partitioning {

PartitionerResult NaivePartitioner::build(int Nb, const std::vector<ContactEdge>& contacts, bool buildContactAdjacency) const {
    PartitionerResult res;

    // MPI world info
    MPI_Comm_rank(MPI_COMM_WORLD, &res.worldRank);
    MPI_Comm_size(MPI_COMM_WORLD, &res.worldSize);

    // Even split bodies per rank
    const int bodiesPerRank = (Nb + res.worldSize - 1) / res.worldSize; // ceil
    res.bodyStart = res.worldRank * bodiesPerRank;
    res.bodyEnd   = std::min(Nb, res.bodyStart + bodiesPerRank);

    const int localCount = std::max(0, res.bodyEnd - res.bodyStart);
    res.bodyToContacts.assign((size_t)localCount, {});

    // Classify contacts
    res.allContacts.clear();
    res.boundaryContacts.clear();
    res.halo.clear();
    res.allContacts.reserve(contacts.size());
    res.boundaryContacts.reserve(contacts.size()/2);
    res.halo.reserve(contacts.size()/2);

    for (int cid = 0; cid < (int)contacts.size(); ++cid) {
        const auto& c = contacts[cid];
        const int a = c.bodyA;
        const int b = c.bodyB;

        // Owners
        const bool aValid = (a >= 0);
        const bool bValid = (b >= 0);
        const int ownerA = aValid ? ownerOf(a, res.worldSize, bodiesPerRank) : -1;
        const int ownerB = bValid ? ownerOf(b, res.worldSize, bodiesPerRank) : -1;

        const bool aLocal = aValid && (a >= res.bodyStart && a < res.bodyEnd);
        const bool bLocal = bValid && (b >= res.bodyStart && b < res.bodyEnd);

        if (aLocal) {
            res.allContacts.push_back(cid);
            res.bodyToContacts[(size_t)(a - res.bodyStart)].push_back(cid);
        }
        if (bLocal && (!aLocal || a != b)) {
            res.allContacts.push_back(cid);
            res.bodyToContacts[(size_t)(b - res.bodyStart)].push_back(cid);
        }

        // boundary if exactly one endpoint is local and the other endpoint is a valid remote body
        if ((aLocal ^ bLocal) && ((aLocal && bValid) || (bLocal && aValid))) {
            res.boundaryContacts.push_back(cid);
            int localBody = aLocal ? a : b;
            int remoteBody = aLocal ? b : a;
            int remoteRank = aLocal ? ownerB : ownerA;
            if (remoteRank >= 0) res.halo.push_back(HaloEntry{localBody, remoteBody, remoteRank, cid});
        }
    }

    // Deduplicate allContacts to ensure each contact is processed once per iteration
    std::sort(res.allContacts.begin(), res.allContacts.end());
    res.allContacts.erase(std::unique(res.allContacts.begin(), res.allContacts.end()), res.allContacts.end());

    // Build neighbor plans to avoid per-iteration hashing in communication
    {
        // Discover neighbors and collect raw lists
        std::unordered_map<int, std::vector<int>> sendBodiesMap;
        std::unordered_map<int, std::vector<int>> recvBodiesMap;
        std::unordered_map<int, std::vector<int>> pSendCidsMap;
        std::unordered_map<int, std::vector<int>> pRecvCidsMap;

        // Compose per-halo info
        for (const auto& h : res.halo) {
            // Body velocities: owner of localBody sends to neighbor; expect to receive remoteBody from neighbor
            if (h.localBody >= res.bodyStart && h.localBody < res.bodyEnd)
                sendBodiesMap[h.remoteRank].push_back(h.localBody);
            recvBodiesMap[h.remoteRank].push_back(h.remoteBody);
        }

        // For p ownership: owner is bodyA's owner. We need contacts list per neighbor for send/recv.
        // Build a quick lookup for owner rank of a body
        const int bodiesPerRank = (Nb + res.worldSize - 1) / res.worldSize;
        auto ownerOfBody = [&](int b) { return (b >= 0) ? ownerOf(b, res.worldSize, bodiesPerRank) : -1; };
        for (int cid : res.boundaryContacts) {
            // Identify the remote rank involved in this contact
            // Pick either endpoint that is remote relative to our local ownership
            int a = contacts[(size_t)cid].bodyA;
            int b = contacts[(size_t)cid].bodyB;
            int ownerA = ownerOfBody(a);
            int ownerB = ownerOfBody(b);
            // If neither endpoint valid continue
            if (ownerA == -1 && ownerB == -1) continue;
            // Skip if both bodies local (shouldn't be boundary)
            if ((a >= res.bodyStart && a < res.bodyEnd) && (b >= res.bodyStart && b < res.bodyEnd)) continue;
            // Determine neighbor rank for this boundary contact
            int nbr = -1;
            if (a >= res.bodyStart && a < res.bodyEnd) nbr = ownerB;
            else if (b >= res.bodyStart && b < res.bodyEnd) nbr = ownerA;
            if (nbr < 0) continue;

            // Who owns p? Owner of bodyA
            bool weOwnP = (ownerA == res.worldRank);
            if (weOwnP) pSendCidsMap[nbr].push_back(cid);
            else pRecvCidsMap[nbr].push_back(cid);
        }

        // Move to vectors with sorted unique entries
        res.neighbors.clear();
        res.bodiesSendPerNeighbor.clear();
        res.bodiesRecvPerNeighbor.clear();
        res.pSendCidsPerNeighbor.clear();
        res.pRecvCidsPerNeighbor.clear();

        // Collect union of neighbor keys
        std::unordered_map<int, bool> seen;
        for (auto& kv : sendBodiesMap) seen[kv.first] = true;
        for (auto& kv : recvBodiesMap) seen[kv.first] = true;
        for (auto& kv : pSendCidsMap) seen[kv.first] = true;
        for (auto& kv : pRecvCidsMap) seen[kv.first] = true;
        for (auto& kv : seen) res.neighbors.push_back(kv.first);
        std::sort(res.neighbors.begin(), res.neighbors.end());

        auto uniq_sort = [](std::vector<int>& v){ std::sort(v.begin(), v.end()); v.erase(std::unique(v.begin(), v.end()), v.end()); };

        for (int nbr : res.neighbors) {
            auto sb = sendBodiesMap[nbr]; uniq_sort(sb);
            auto rb = recvBodiesMap[nbr]; uniq_sort(rb);
            auto ps = pSendCidsMap[nbr]; uniq_sort(ps);
            auto pr = pRecvCidsMap[nbr]; uniq_sort(pr);
            res.bodiesSendPerNeighbor.push_back(std::move(sb));
            res.bodiesRecvPerNeighbor.push_back(std::move(rb));
            res.pSendCidsPerNeighbor.push_back(std::move(ps));
            res.pRecvCidsPerNeighbor.push_back(std::move(pr));
        }
    }

    if (buildContactAdjacency) {
        // Build set of relevant contacts (local + boundary)
    std::vector<int> relevant = res.allContacts;
        const int R = (int)relevant.size();
        // Map contactId -> 0..R-1 compact index
        std::unordered_map<int,int> cmap; cmap.reserve(R*2);
        for (int i = 0; i < R; ++i) cmap[relevant[i]] = i;

        res.contactNeighbors.assign((size_t)R, {});

        // For each local body, all incident contacts are neighbors of each other
        for (int lb = 0; lb < localCount; ++lb) {
            const auto& lst = res.bodyToContacts[(size_t)lb];
            // Add undirected edges among all contacts in lst if both are in cmap
            for (size_t i = 0; i < lst.size(); ++i) {
                auto it_i = cmap.find(lst[i]); if (it_i == cmap.end()) continue;
                for (size_t j = i + 1; j < lst.size(); ++j) {
                    auto it_j = cmap.find(lst[j]); if (it_j == cmap.end()) continue;
                    int ci = it_i->second, cj = it_j->second;
                    res.contactNeighbors[(size_t)ci].push_back(cj);
                    res.contactNeighbors[(size_t)cj].push_back(ci);
                }
            }
        }
        // Optionally: deduplicate neighbors
        for (auto& nbrs : res.contactNeighbors) {
            std::sort(nbrs.begin(), nbrs.end());
            nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
        }
    }

    return res;
}

}
