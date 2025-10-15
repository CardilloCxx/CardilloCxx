#pragma once

#include <vector>
#include <unordered_map>
#include <utility>
#include <mpi.h>
#include "../misc/types.hpp"

namespace cardillo::partitioning {

struct ContactEdge {
    int bodyA;
    int bodyB;
    MatrixXXr WblockA;    // dofPxdofV block for bodyA
    MatrixXXr WblockB;    // dofPxdofV block for bodyB
    MatrixXXr WblockMinv; // dofVxdofV block for Minv 
};

struct HaloEntry {
    int localBody;
    int remoteBody;
    int remoteRank;
    int contactId; // index into contacts vector provided to build()
};

struct PartitionerResult {
    int worldRank = 0;
    int worldSize = 1;
    int bodyStart = 0; // inclusive
    int bodyEnd = 0;   // exclusive

    // Contacts classified for this rank
    // allContacts: contacts that touch at least one local body (includes boundary)
    // boundaryContacts: subset of allContacts where the other endpoint is a valid remote body (not -1)
    std::vector<int> allContacts;      // indices into input contacts (local + boundary)
    std::vector<int> boundaryContacts; // indices into input contacts touching a remote rank
    std::vector<HaloEntry> halo;       // halo communication descriptors

    // Body -> contacts adjacency for local bodies [bodyStart, bodyEnd)
    std::vector<std::vector<int>> bodyToContacts; // local body index -> list of contact indices

    // Optional: contact adjacency (each contact neighbors others sharing a body)
    std::vector<std::vector<int>> contactNeighbors; // only for local+boundary contacts

    // Precomputed neighbor plans (to avoid per-iteration hashing in communication):
    // Distinct neighbor ranks we communicate with (sorted ascending)
    std::vector<int> neighbors;
    // For each neighbors[i], list of local body indices (global ids) this rank sends to neighbor
    std::vector<std::vector<int>> bodiesSendPerNeighbor;
    // For each neighbors[i], list of remote body indices (global ids) this rank expects to receive from neighbor
    std::vector<std::vector<int>> bodiesRecvPerNeighbor;
    // For each neighbors[i], list of contact ids for which this rank owns p (owner = bodyA owner) and must send to neighbor
    std::vector<std::vector<int>> pSendCidsPerNeighbor;
    // For each neighbors[i], list of contact ids this rank expects to receive p values for (owned by neighbor)
    std::vector<std::vector<int>> pRecvCidsPerNeighbor;
};

class NaivePartitioner {
public:
    // Build a naive per-timestep partitioning for bodies 0..Nb-1 and contacts array.
    // Does not store any global W; consumes per-contact blocks.
    PartitionerResult build(int Nb, const std::vector<ContactEdge>& contacts, bool buildContactAdjacency = false) const;

    // Public helper: compute owner rank for a body given worldSize and bodies-per-rank (ceil split)
    static inline int ownerOf(int body, int worldSize, int bodiesPerRank) {
        int r = body / bodiesPerRank;
        if (r >= worldSize) r = worldSize - 1;
        if (r < 0) r = 0;
        return r;
    }
};

}
