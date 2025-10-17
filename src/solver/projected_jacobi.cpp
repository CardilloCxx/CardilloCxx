#include "projected_jacobi.hpp"
#include <iostream>
#include <mpi.h>
#include <unordered_map>
#include <algorithm>
#include <limits>
#include <cmath>
#include <cstring>
#include "../partitioning/naive_partitioner.hpp"
#include "../comm/communication.hpp"

namespace cardillo::solver {

// Small helpers
static inline std::vector<real_t> precompute_Rdiag(const std::vector<partitioning::ContactEdge>& edges,
												   const std::vector<int>& relevantContacts,
												   const std::vector<MatrixXXr>& MinvBlocks,
												   real_t alpha,
												   real_t compliance) {
	const int C = (int)edges.size();
	std::vector<real_t> Rdiag_local(C, 0);
	for (int cid : relevantContacts) {
		const auto& e = edges[cid];
		real_t Dii = 0;
		if (e.bodyA >= 0 && e.WblockA.rows() > 0 && e.WblockA.cols() > 0) {
			const auto& MinvA = MinvBlocks[(size_t)e.bodyA];
			if (MinvA.rows() == e.WblockA.cols()) {
				VectorXr wA = e.WblockA.row(0).transpose();
				Dii += wA.transpose() * (MinvA * wA);
			}
		}
		if (e.bodyB >= 0 && e.WblockB.rows() > 0 && e.WblockB.cols() > 0) {
			const auto& MinvB = MinvBlocks[(size_t)e.bodyB];
			if (MinvB.rows() == e.WblockB.cols()) {
				VectorXr wB = e.WblockB.row(0).transpose();
				Dii += wB.transpose() * (MinvB * wB);
			}
		}
		// Add compliance to stabilize near-singular Dii
		Dii += compliance;
		Rdiag_local[cid] = (Dii > (real_t)0) ? (alpha / Dii) : (real_t)0;
	}
	return Rdiag_local;
}

std::vector<VectorXr> ProjectedJacobiSolver::iterateWithPreliminaryVelocity(const std::vector<VectorXr>& v_pre_blocks, real_t tol) {
	// Use block data from the assembler (SoA)
	const auto& Wblocks = m_dyn.WBlocks();           // size may be <= 2*C if static endpoints are omitted
	const auto& MinvBlocksVec = m_dyn.MinvBlocks();
	const auto& bodyToBlocks = m_dyn.WBlocksFromBodyAll();       // size Nb, lists of block indices
	const auto& blockToBody = m_dyn.WBlockToBodyAll();           // size 2*C, W-block index -> body
	const auto& blockToContact = m_dyn.WBlockToContactAll();     // size 2*C, W-block index -> contact id
	const auto& contactToBlocks = m_dyn.WBlocksFromContactAll(); // size C, pairs of block indices
	
	const int Nb = (int)v_pre_blocks.size();
	const int C = (int)contactToBlocks.size();
	if (C == 0) return v_pre_blocks; // no contacts, return preliminary velocity as-is

	int worldSize = 1, worldRank = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

	// Convert to partitioner edges per contact using mappings
	std::vector<partitioning::ContactEdge> partEdges(C);
	for (int i = 0; i < C; ++i) {
		auto [idxA, idxB] = contactToBlocks[(size_t)i];
		partEdges[i].bodyA = (idxA >= 0 ? blockToBody[(size_t)idxA] : -1);
		partEdges[i].bodyB = (idxB >= 0 ? blockToBody[(size_t)idxB] : -1);
		partEdges[i].WblockA = (idxA >= 0 ? Wblocks[(size_t)idxA] : MatrixXXr());
		partEdges[i].WblockB = (idxB >= 0 ? Wblocks[(size_t)idxB] : MatrixXXr());
	}

	partitioning::NaivePartitioner part;
	auto res = part.build(Nb, partEdges, false);

	// Precompute R = alpha / D per relevant contact (allContacts includes boundary)
	std::vector<real_t> Rdiag_local = precompute_Rdiag(partEdges, res.allContacts, MinvBlocksVec, alpha(), m_compliance);

	// Initialize state
	std::vector<VectorXr> u = v_pre_blocks;           // current body velocities
	std::vector<VectorXr> u_prev = u;                 // previous iteration
	const std::vector<VectorXr> u_free = v_pre_blocks;// preliminary velocities
	VectorXr p = VectorXr::Zero(C);
	if (m_warmStart && m_last_p.has_value() && m_last_p->size() == C) p = m_last_p.value();

	// Hoisted temporaries: per-body accumulator for W^T * p
	std::vector<VectorXr> sumW_per_body((size_t)Nb);
	for (int b = 0; b < Nb; ++b) sumW_per_body[(size_t)b].resize(u[(size_t)b].size());

	// Main iteration loop
	real_t err_global = std::numeric_limits<real_t>::max(); index_t iteration = 0;
	int stagnant = 0; // backtracking patience counter
	while (iteration < m_maxIterations) {
		// Save previous u for error on local bodies
		for (int b = res.bodyStart; b < res.bodyEnd; ++b) u_prev[(size_t)b] = u[(size_t)b];

		// p update: p_np1 = -Prox(R * (W * u) - p)
		for (int cid : res.allContacts) {
			auto [idxA, idxB] = contactToBlocks[(size_t)cid];
			real_t y = 0;
			if (idxA >= 0) y += Wblocks[(size_t)idxA].row(0) * u[blockToBody[(size_t)idxA]];
			if (idxB >= 0) y += Wblocks[(size_t)idxB].row(0) * u[blockToBody[(size_t)idxB]];
			real_t p_new = -std::min<real_t>(0, Rdiag_local[cid] * y - p[cid]);
			// Relaxation blending p
			p[cid] = (real_t)1.0 * ((real_t)1.0 - m_relax) * p[cid] + m_relax * p_new;
		}

		// Exchange boundary p values: owner (bodyA owner) pushes p to neighbors
		cardillo::comm::Communication::exchangePercussionsOwnerPush(p, res, contactToBlocks, blockToBody);

		// u update: u_np1 = u_free + Minv * W^T * p (uses updated p)
		for (int b = res.bodyStart; b < res.bodyEnd; ++b) {
			auto& sumW = sumW_per_body[(size_t)b];
			sumW.setZero();
			for (int idx : bodyToBlocks[(size_t)b]) sumW.noalias() += Wblocks[(size_t)idx].transpose() * p[blockToContact[(size_t)idx]];
			u[(size_t)b] = u_free[(size_t)b] + MinvBlocksVec[(size_t)b] * sumW;
		}

		// Exchange boundary u values: owner (body owner) pushes u to neighbors
		cardillo::comm::Communication::exchangeBodyVelocitiesOwnerPush(u, res);

		// Convergence: global norm across local bodies
		real_t loc_sum = 0; for (int b = res.bodyStart; b < res.bodyEnd; ++b) loc_sum += (u[(size_t)b] - u_prev[(size_t)b]).squaredNorm();
		MPI_Allreduce(&loc_sum, &err_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); err_global = std::sqrt((double)err_global);

		++iteration;
		if (err_global <= tol) break;

		if (iteration % 1000 == 0 && worldRank == 0) {
			std::cout << "[ProjectedJacobi] Iteration " << iteration << ", error = " << err_global << std::endl;
		}
	}

	// Store p for warm start next call
	m_last_p = p;

	// Final synchronization: replicate all body velocities to every rank
	cardillo::comm::Communication::replicateAllBodyVelocities(u, res);

	return u;
}

} // namespace cardillo::solver