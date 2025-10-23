#include "communication.hpp"
#include <algorithm>

namespace cardillo::comm {

void Communication::exchangeBodyVelocitiesOwnerPush(
    std::vector<VectorXr>& u,
    const cardillo::partitioning::PartitionerResult& res)
{
    // Use precomputed neighbor plans
    const auto& neighbors = res.neighbors;
    const auto& sendLists = res.bodiesSendPerNeighbor;
    const auto& recvLists = res.bodiesRecvPerNeighbor;

    // Exchange counts and payloads per neighbor
    for (size_t i = 0; i < neighbors.size(); ++i) {
        int nbr = neighbors[i];
        const auto& bodies = sendLists[i];
        int sendCount = (int)bodies.size();
        const auto& recvBodies = recvLists[i];
        int recvCount = (int)recvBodies.size();

        // Pack send buffer: concat vectors in order
        int dof = (sendCount > 0) ? (int)u[(size_t)bodies[0]].size() : 0;
        static thread_local std::vector<real_t> sendBuf;
        static thread_local std::vector<real_t> recvBuf;
        sendBuf.resize((size_t)dof * (size_t)sendCount);
        for (int j = 0; j < sendCount; ++j) {
            const auto& vec = u[(size_t)bodies[(size_t)j]];
            std::copy(vec.data(), vec.data() + vec.size(), sendBuf.data() + j * dof);
        }
        // Probe recv dof if any
        int recvDof = dof;
        // Exchange dof
        MPI_Sendrecv(&dof, 1, MPI_INT, nbr, 102,
                     &recvDof, 1, MPI_INT, nbr, 102,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        recvBuf.resize((size_t)recvDof * (size_t)recvCount);
        MPI_Sendrecv(sendBuf.data(), (int)sendBuf.size(), MPI_DOUBLE, nbr, 103,
                     recvBuf.data(), (int)recvBuf.size(), MPI_DOUBLE, nbr, 103,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Unpack to u for received bodies
        for (int j = 0; j < recvCount; ++j) {
            int b = recvBodies[(size_t)j];
            if (b >= 0 && b < (int)u.size()) {
                if ((int)u[(size_t)b].size() != recvDof) u[(size_t)b].resize(recvDof);
                std::copy(recvBuf.data() + j * recvDof, recvBuf.data() + (j+1) * recvDof, u[(size_t)b].data());
            }
        }
    }
}

void Communication::exchangePercussionsOwnerPush(
    VectorXr& p,
    const cardillo::partitioning::PartitionerResult& res)
{
    // Use precomputed neighbor plans for percussions
    const auto& neighbors = res.neighbors;
    const auto& pSendLists = res.pSendCidsPerNeighbor;
    const auto& pRecvLists = res.pRecvCidsPerNeighbor;

    for (size_t i = 0; i < neighbors.size(); ++i) {
        int nbr = neighbors[i];
        const auto& sendIds = pSendLists[i];
        const auto& recvIds = pRecvLists[i];

        int sendCount = (int)sendIds.size();
        int recvCount = (int)recvIds.size();

        // Build and exchange values
        static thread_local std::vector<real_t> sendVals;
        static thread_local std::vector<real_t> recvVals;
        sendVals.resize((size_t)sendCount);
        for (int j = 0; j < sendCount; ++j) sendVals[(size_t)j] = p[sendIds[(size_t)j]];
        recvVals.resize((size_t)recvCount);
        MPI_Sendrecv(sendVals.data(), sendCount, MPI_DOUBLE, nbr, 202,
                     recvVals.data(), recvCount, MPI_DOUBLE, nbr, 202,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int j = 0; j < recvCount; ++j) {
            int cid = recvIds[(size_t)j];
            if (cid >= 0 && cid < p.size()) p[(size_t)cid] = recvVals[(size_t)j];
        }
    }
}

void Communication::replicateAllBodyVelocities(
    std::vector<VectorXr>& u,
    const cardillo::partitioning::PartitionerResult& res)
{
    const int Nb = (int)u.size();
    // Build contiguous layout
    std::vector<int> sizes(Nb, 0), offsets(Nb, 0);
    int total = 0;
    for (int b = 0; b < Nb; ++b) { sizes[b] = (int)u[(size_t)b].size(); offsets[b] = total; total += sizes[b]; }
    std::vector<real_t> sendbuf((size_t)total, (real_t)0), recvbuf((size_t)total, (real_t)0);
    // Pack only local bodies; others remain zero so Allreduce sum yields owners' values
    for (int b = res.bodyStart; b < res.bodyEnd; ++b) {
        int n = sizes[(size_t)b]; if (n <= 0) continue;
        auto* dst = sendbuf.data() + offsets[(size_t)b];
        const auto& ub = u[(size_t)b];
        std::copy(ub.data(), ub.data() + n, dst);
    }
    // Allreduce (sum) to combine owner contributions into a full vector on every rank
    MPI_Allreduce(sendbuf.data(), recvbuf.data(), total, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // Unpack into u
    for (int b = 0; b < Nb; ++b) {
        int n = sizes[(size_t)b]; if (n <= 0) continue;
        u[(size_t)b] = Eigen::Map<const VectorXr>(recvbuf.data() + offsets[(size_t)b], n);
    }
}

void Communication::exchangeBodyVelocitiesOwnerPushConcat(
    VectorXr& u_concat,
    const cardillo::partitioning::PartitionerResult& res,
    const std::vector<int>& bodyOffsets)
{
    const auto& neighbors = res.neighbors;
    const auto& sendLists = res.bodiesSendPerNeighbor;
    const auto& recvLists = res.bodiesRecvPerNeighbor;

    for (size_t i = 0; i < neighbors.size(); ++i) {
        int nbr = neighbors[i];
        const auto& bodies = sendLists[i];
        int sendCount = (int)bodies.size();
        const auto& recvBodies = recvLists[i];
        int recvCount = (int)recvBodies.size();

        // Determine dof for this neighbor exchange (assumed uniform as before)
        int dof = 0;
        if (sendCount > 0) {
            int b0 = bodies[0];
            dof = bodyOffsets[(size_t)b0 + 1] - bodyOffsets[(size_t)b0];
        } else if (recvCount > 0) {
            int b0 = recvBodies[0];
            dof = bodyOffsets[(size_t)b0 + 1] - bodyOffsets[(size_t)b0];
        }
        static thread_local std::vector<real_t> sendBuf;
        static thread_local std::vector<real_t> recvBuf;
        sendBuf.resize((size_t)dof * (size_t)sendCount);
        for (int j = 0; j < sendCount; ++j) {
            int b = bodies[(size_t)j];
            int off = bodyOffsets[(size_t)b];
            std::copy(u_concat.data() + off, u_concat.data() + off + dof, sendBuf.data() + j * dof);
        }
        int recvDof = dof;
        MPI_Sendrecv(&dof, 1, MPI_INT, nbr, 112,
                     &recvDof, 1, MPI_INT, nbr, 112,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        recvBuf.resize((size_t)recvDof * (size_t)recvCount);
        MPI_Sendrecv(sendBuf.data(), (int)sendBuf.size(), MPI_DOUBLE, nbr, 113,
                     recvBuf.data(), (int)recvBuf.size(), MPI_DOUBLE, nbr, 113,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int j = 0; j < recvCount; ++j) {
            int b = recvBodies[(size_t)j];
            int off = bodyOffsets[(size_t)b];
            std::copy(recvBuf.data() + j * recvDof, recvBuf.data() + (j+1) * recvDof, u_concat.data() + off);
        }
    }
}

void Communication::replicateAllBodyVelocitiesConcat(
    VectorXr& u_concat,
    const cardillo::partitioning::PartitionerResult& res,
    const std::vector<int>& bodyOffsets)
{
    const int Nb = (int)bodyOffsets.size() - 1;
    // total size is last offset
    int total = (Nb >= 0) ? bodyOffsets[(size_t)Nb] : 0;
    std::vector<real_t> sendbuf((size_t)total, (real_t)0), recvbuf((size_t)total, (real_t)0);
    // Pack local bodies only
    for (int b = res.bodyStart; b < res.bodyEnd; ++b) {
        int off = bodyOffsets[(size_t)b];
        int n = bodyOffsets[(size_t)b+1] - off;
        if (n > 0) std::copy(u_concat.data() + off, u_concat.data() + off + n, sendbuf.data() + off);
    }
    MPI_Allreduce(sendbuf.data(), recvbuf.data(), total, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // Unpack
    if ((int)u_concat.size() != total) u_concat.resize(total);
    std::copy(recvbuf.data(), recvbuf.data() + total, u_concat.data());
}

}
