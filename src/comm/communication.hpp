#pragma once

#include <vector>
#include <unordered_map>
#include <mpi.h>
#include "../misc/types.hpp"
#include "../partitioning/naive_partitioner.hpp"

namespace cardillo::comm {

class Communication {
public:
    // Owner-push body velocities u to neighbors that reference them in halo
    static void exchangeBodyVelocitiesOwnerPush(
        std::vector<VectorXr>& u,
        const cardillo::partitioning::PartitionerResult& res);

    // Owner-push per-contact percussions p: ownership defined by bodyA owner (first endpoint)
    static void exchangePercussionsOwnerPush(
        VectorXr& p,
        const cardillo::partitioning::PartitionerResult& res);

    // Replicate all body velocities to every rank via Allreduce (sum of owner contributions)
    static void replicateAllBodyVelocities(
        std::vector<VectorXr>& u,
        const cardillo::partitioning::PartitionerResult& res);

    // Concatenated variants using body offsets
    static void exchangeBodyVelocitiesOwnerPushConcat(
        VectorXr& u_concat,
        const cardillo::partitioning::PartitionerResult& res,
        const std::vector<int>& bodyOffsets);

    static void replicateAllBodyVelocitiesConcat(
        VectorXr& u_concat,
        const cardillo::partitioning::PartitionerResult& res,
        const std::vector<int>& bodyOffsets);
};

}
