#include "projected_jacobi.hpp"
#include <iostream>
#include <mpi.h>
#include <algorithm>
#include <limits>
#include <cmath>
#include <unordered_map>
#include <vector>
#include "../partitioning/naive_partitioner.hpp"
#include "../comm/communication.hpp"

namespace cardillo::solver {

// Precompute R = alpha / diag(G + compliance I) for relevant contacts
// G = W * Minv * W^T, with Minv diagonal (given as MinvDiag)
static inline VectorXr precompute_Rdiag_sparse(
	const Eigen::SparseMatrix<real_t, Eigen::RowMajor>& W,
	const VectorXr& MinvDiag,
	const std::vector<int>& relevantContacts,
	real_t alpha,
	real_t compliance)
{
	const int C = (int)W.rows();
	VectorXr R = VectorXr::Zero(C);
	for (int cid : relevantContacts) {
		if (cid < 0 || cid >= C) continue;
		real_t Dii = 0;
		for (Eigen::SparseMatrix<real_t, Eigen::RowMajor>::InnerIterator it(W, cid); it; ++it) {
			const int col = it.col();
			const real_t w = it.value();
			Dii += w * w * MinvDiag[col];
		}
		Dii += compliance;
		R[cid] = (Dii > (real_t)0) ? (alpha / Dii) : (real_t)0;
	}
	return R;
}

// Small context to avoid re-passing many references for each sweep
namespace {
struct PJIterContext {
	const Eigen::SparseMatrix<real_t, Eigen::RowMajor>& W;
	const VectorXr& MinvDiag;
	std::vector<int> bodyOffsets;
	partitioning::PartitionerResult res;
	VectorXr Rdiag;
	// Grouping of contact rows by original contact index (for friction cone projection)
	std::vector<std::vector<int>> contactRowGroups; // each group has 1 (normal only) or 3 rows (n,t1,t2) sorted by row id
	std::vector<real_t> groupMu; // same size as contactRowGroups: per-contact friction mu
	// Iteration runtime state and buffers
	VectorXr y_contact;   // contact-space buffer (size = #contacts)
	VectorXr tmp_concat;  // velocity-space buffer (size = concat dofs)
	// Controls and diagnostics for the loop
	int worldRank{0};
	int maxIterations{1000000};
	real_t tol{(real_t)1e-5};
	real_t relax{(real_t)1};
	bool debug{false};
	// Outputs
	real_t err_global{std::numeric_limits<real_t>::max()};
	index_t iteration{0};
};
}

// Compact holder for grouped contact rows and friction coefficients
struct ContactGroups {
	std::vector<std::vector<int>> rows; // 1 or 3 rows per contact
	std::vector<real_t> mu;             // same length as rows
};

static inline ContactGroups build_contact_groups(const cardillo::physics::DynamicsAssembler& dyn)
{
	ContactGroups cg;
	const auto& W = dyn.W();
	const int C = (int)W.rows();
	const auto& dynToOrig = dyn.dynamicContactToOriginalAll();
	const auto& contactsAll = dyn.contacts();
	std::unordered_map<int, std::vector<int>> groups;
	groups.reserve((size_t)C);
	for (int r = 0; r < C; ++r) {
		int origIdx = (r >= 0 && r < (int)dynToOrig.size()) ? dynToOrig[(size_t)r] : -1;
		if (origIdx < 0 || origIdx >= (int)contactsAll.size()) continue;
		groups[origIdx].push_back(r);
	}
	cg.rows.reserve(groups.size());
	cg.mu.reserve(groups.size());
	for (auto& kv : groups) {
		auto& rows = kv.second;
		std::sort(rows.begin(), rows.end());
		cg.rows.push_back(rows);
		cg.mu.push_back(contactsAll[(size_t)kv.first].friction_mu);
	}
	return cg;
}

static inline PJIterContext build_context(
	cardillo::physics::DynamicsAssembler& dyn,
	real_t alpha,
	real_t compliance)
{
	const auto& W = dyn.W();
	const auto& MinvDiagVec = dyn.MinvDiag();
	const auto& bodyOffsets = dyn.bodyVelOffsets();
	const int C = (int)W.rows();

	std::vector<partitioning::ContactEdge> partEdges(C);
	const auto& dynToOrig = dyn.dynamicContactToOriginalAll();
	const auto& reg = dyn.system().ecs();
	const auto& contactsAll = dyn.contacts();
	for (int i = 0; i < C; ++i) {
		int origIdx = (i >= 0 && i < (int)dynToOrig.size()) ? dynToOrig[(size_t)i] : -1;
		int ba = -1, bb = -1;
		if (origIdx >= 0 && origIdx < (int)contactsAll.size()) {
			const auto& c = contactsAll[(size_t)origIdx];
			if (reg.any_of<cardillo::PhysicsSystem::C_BodyIndex>(c.a)) ba = reg.get<cardillo::PhysicsSystem::C_BodyIndex>(c.a).b;
			if (reg.any_of<cardillo::PhysicsSystem::C_BodyIndex>(c.b)) bb = reg.get<cardillo::PhysicsSystem::C_BodyIndex>(c.b).b;
		}
		partEdges[(size_t)i].bodyA = ba;
		partEdges[(size_t)i].bodyB = bb;
	}

	partitioning::NaivePartitioner part;
	int Nb = (int)bodyOffsets.size() - 1;
	auto res = part.build(Nb, partEdges, false);

	PJIterContext ctx{dyn.W(), dyn.MinvDiag(), std::vector<int>{}, std::move(res), VectorXr()};
	ctx.bodyOffsets = dyn.bodyVelOffsets();
	ctx.Rdiag = precompute_Rdiag_sparse(W, MinvDiagVec, ctx.res.allContacts, alpha, compliance);
	// Build contact row groups and per-contact mu
	const ContactGroups cg = build_contact_groups(dyn);
	ctx.contactRowGroups = cg.rows;
	ctx.groupMu = cg.mu;
	// allocate buffers
	const int dofConcat = (ctx.bodyOffsets.empty() ? 0 : ctx.bodyOffsets.back());
	ctx.y_contact = VectorXr::Zero(C);
	ctx.tmp_concat = VectorXr::Zero(dofConcat);
	return ctx;
}

// Compute provisional gradient step z = p_old - R * (W * u_in)
static inline VectorXr provisional_step(PJIterContext& ctx, const VectorXr& u_in, const VectorXr& p_old)
{
	const auto& W = ctx.W; const auto& Rdiag = ctx.Rdiag; const auto& res = ctx.res;
	ctx.y_contact.noalias() = W * u_in;
	VectorXr z = p_old; // start from old p
	for (int cid : res.allContacts) {
		z[cid] = p_old[cid] - Rdiag[cid] * ctx.y_contact[cid];
	}
	return z;
}

// Project per-contact groups: 1-row -> unilateral; 3-rows -> Coulomb cone (disk in tangential plane)
static inline void project_groups(const PJIterContext& ctx, const VectorXr& z, const VectorXr& p_old, VectorXr& p_proj)
{
	p_proj = p_old; // copy, then overwrite affected entries
	const auto& groups = ctx.contactRowGroups; const auto& mus = ctx.groupMu;
	for (size_t gi = 0; gi < groups.size(); ++gi) {
		const auto& rows = groups[gi];
		if (rows.empty()) continue;
		if (rows.size() == 1) {
			const int rn = rows[0];
			p_proj[rn] = std::max<real_t>(z[rn], (real_t)0);
		} else {
			const int rn = rows[0];
			const int rt1 = rows.size() > 1 ? rows[1] : rows[0];
			const int rt2 = rows.size() > 2 ? rows[2] : rows[0];
			const real_t mu = mus[gi];
			const real_t pn = std::max<real_t>(z[rn], (real_t)0);
			const real_t t1 = z[rt1], t2 = z[rt2];
			const real_t tnorm = std::sqrt((double)(t1*t1 + t2*t2));
			const real_t s = (tnorm > (real_t)0) ? std::min<real_t>((real_t)1, (mu * pn) / (real_t)tnorm) : (real_t)1;
			p_proj[rn] = pn;
			if (rows.size() > 1) p_proj[rt1] = s * t1;
			if (rows.size() > 2) p_proj[rt2] = s * t2;
		}
	}
}

// Apply relaxation and exchange p; then update u and exchange u
static inline void relax_exchange_update(PJIterContext& ctx, const VectorXr& u_free,
										 const VectorXr& p_old, const VectorXr& p_proj,
										 VectorXr& p_inout, VectorXr& u_out)
{
	const auto& W = ctx.W; const auto& Minv = ctx.MinvDiag; const auto& bodyOffsets = ctx.bodyOffsets; const auto& res = ctx.res;
	for (int cid : res.allContacts) {
		const real_t p_new = p_proj[cid];
		p_inout[cid] = ((real_t)1 - ctx.relax) * p_old[cid] + ctx.relax * p_new;
	}
	cardillo::comm::Communication::exchangePercussionsOwnerPush(p_inout, res);
	ctx.tmp_concat.noalias() = W.transpose() * p_inout;
	u_out = u_free + Minv.cwiseProduct(ctx.tmp_concat);
	cardillo::comm::Communication::exchangeBodyVelocitiesOwnerPushConcat(u_out, res, bodyOffsets);
}

// One Projected-Jacobi sweep: update p with current u, then compute new u
static inline void pj_sweep(PJIterContext& ctx, const VectorXr& u_in, const VectorXr& u_free,
							VectorXr& p_inout, VectorXr& u_out)
{
	const VectorXr z = provisional_step(ctx, u_in, p_inout);
	VectorXr p_proj;
	project_groups(ctx, z, p_inout, p_proj);
	relax_exchange_update(ctx, u_free, p_inout, p_proj, p_inout, u_out);
}

// Compute global L2 norm of difference across local body segments
static inline real_t global_segment_norm(const PJIterContext& ctx, const VectorXr& a, const VectorXr& b) {
	const auto& bo = ctx.bodyOffsets; const auto& res = ctx.res;
	double loc_sum = 0.0;
	for (int br = res.bodyStart; br < res.bodyEnd; ++br) {
		int off = bo[(size_t)br]; int n = bo[(size_t)br+1] - off;
		if (n > 0) {
			Eigen::Map<const VectorXr> va(a.data()+off, n), vb(b.data()+off, n);
			loc_sum += (va - vb).squaredNorm();
		}
	}
	double gsum = 0.0; MPI_Allreduce(&loc_sum, &gsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return (real_t)std::sqrt(gsum);
}

static inline double global_segment_dot(const PJIterContext& ctx, const VectorXr& a, const VectorXr& b) {
	const auto& bo = ctx.bodyOffsets; const auto& res = ctx.res;
	double loc = 0.0;
	for (int br = res.bodyStart; br < res.bodyEnd; ++br) {
		int off = bo[(size_t)br]; int n = bo[(size_t)br+1] - off;
		if (n > 0) {
			Eigen::Map<const VectorXr> va(a.data()+off, n), vb(b.data()+off, n);
			loc += va.dot(vb);
		}
	}
	double g = 0.0; MPI_Allreduce(&loc, &g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); return g;
}

namespace {

// Standard fixed-point iteration loop
static inline void standard_loop(PJIterContext& ctx,
								VectorXr& u, const VectorXr& u_free, VectorXr& p) {
	ctx.err_global = std::numeric_limits<real_t>::max();
	ctx.iteration = 0;
	VectorXr u_prev = u;
	while (ctx.iteration < ctx.maxIterations) {
		u_prev = u;
		pj_sweep(ctx, u_prev, u_free, p, u);
		ctx.err_global = global_segment_norm(ctx, u, u_prev);
		++ctx.iteration;
		if (ctx.err_global <= ctx.tol) break;
		if (ctx.debug && ctx.worldRank == 0 && (ctx.iteration % 1000 == 0)) {
			std::cout << "[ProjectedJacobi] Iteration " << ctx.iteration << ", error = " << ctx.err_global << std::endl;
		}
	}
}

// Nesterov-accelerated loop
static inline void nesterov_loop(PJIterContext& ctx,
								VectorXr& u, const VectorXr& u_free, VectorXr& p,
								double beta_threshold, int restart_limit) {
	ctx.err_global = std::numeric_limits<real_t>::max();
	ctx.iteration = 0;
	VectorXr xk = u;
	VectorXr yk = xk;
	VectorXr xk1 = xk;
	double thk = 1.0;
	real_t err_prev = std::numeric_limits<real_t>::infinity();
	int restarts = 0;
	bool momentum_disabled = false;
	while (ctx.iteration < ctx.maxIterations) {
		pj_sweep(ctx, yk, u_free, p, xk1);
		ctx.err_global = global_segment_norm(ctx, xk1, yk);
		++ctx.iteration;
		if (ctx.err_global <= ctx.tol) { xk = xk1; break; }

		bool restart = false;
		if (!std::isfinite((double)ctx.err_global)) restart = true;
		else if (err_prev < std::numeric_limits<real_t>::infinity() && ctx.err_global > (real_t)1.05 * err_prev) restart = true;

		double thk1 = 0.5 * (1.0 + std::sqrt(4.0 * thk * thk + 1.0));
		double betak1 = (thk - 1.0) / thk1;
		if (betak1 > beta_threshold) restart = true;
		{
			VectorXr a = yk - xk1; VectorXr b = xk1 - xk;
			double d = global_segment_dot(ctx, a, b);
			if (d > 0.0) restart = true;
		}
		if (!std::isfinite(betak1) || betak1 < 0.0 || betak1 > 1.0) restart = true;

		if (restart) {
			yk = xk1;
			thk = 1.0;
			xk = xk1;
			if (++restarts >= restart_limit) momentum_disabled = true;
		} else {
			if (momentum_disabled) {
				yk = xk1;
				thk = 1.0;
			} else {
				yk = xk1 + (real_t)betak1 * (xk1 - xk);
				thk = thk1;
			}
			xk = xk1;
		}
		err_prev = ctx.err_global;

		if (ctx.debug && ctx.worldRank == 0 && (ctx.iteration % 1000 == 0)) {
			std::cout << "[ProjectedJacobi+Nesterov] Iter " << ctx.iteration << ", error = " << ctx.err_global
					  << ", beta = " << betak1 << std::endl;
		}
	}
	u = xk;
}

}

VectorXr ProjectedJacobiSolver::iterateWithPreliminaryVelocity(const VectorXr& v_pre, real_t tol) {
	// Use sparse W and block-diagonal Minv
	const auto& Wref = m_dyn.W();
	const int C = (int)Wref.rows();
	if (C == 0 || Wref.nonZeros() == 0) return v_pre; // no contacts, return preliminary velocity as-is

	int worldRank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

	// Build iteration context and precomputations
	PJIterContext ctx = build_context(m_dyn, alpha(), m_compliance);
	const auto& bodyOffsets = ctx.bodyOffsets;

	// Initialize state
	VectorXr u = v_pre;                  // current stacked velocities
	const VectorXr u_free = v_pre;       // preliminary velocities
	VectorXr p = VectorXr::Zero(C);
	if (m_warmStart && m_wsProvider != nullptr) {
		const auto& contactsAll = m_dyn.contacts();
		std::vector<ContactImpulse> hints = m_wsProvider->makeHint(contactsAll);
		const auto& dynToOrig = m_dyn.dynamicContactToOriginalAll();
		// Expand per-contact hints into row vector according to grouped rows
		for (const auto& rows : ctx.contactRowGroups) {
			if (rows.empty()) continue;
			const int rn = rows[0];
			const int origIdx = (rn >= 0 && rn < (int)dynToOrig.size()) ? dynToOrig[(size_t)rn] : -1;
			if (origIdx < 0 || origIdx >= (int)hints.size()) continue;
			const ContactImpulse& hi = hints[(size_t)origIdx];
			p[rn] = std::max<real_t>(hi.pn, (real_t)0);
			if (rows.size() > 1) p[rows[1]] = hi.pt1;
			if (rows.size() > 2) p[rows[2]] = hi.pt2;
		}
	}

	// Configure iteration context
	ctx.worldRank = worldRank;
	ctx.maxIterations = m_maxIterations;
	ctx.tol = tol;
	ctx.relax = m_relax;
	ctx.debug = m_debug;

	// Initial ghost exchanges to ensure consistent data across ranks before first sweep
	cardillo::comm::Communication::exchangeBodyVelocitiesOwnerPushConcat(u, ctx.res, ctx.bodyOffsets);
	if (m_warmStart) {
		cardillo::comm::Communication::exchangePercussionsOwnerPush(p, ctx.res);
	}

	if (!m_useNesterov) {
		standard_loop(ctx, u, u_free, p);
		if (m_debug && worldRank == 0) {
			std::cout << "[ProjectedJacobi] Converged in " << ctx.iteration << " iterations, final error = " << ctx.err_global << std::endl;
		}
	} else {
		nesterov_loop(ctx, u, u_free, p, m_nest_beta_threshold, m_nest_restart_limit);
		if (m_debug && worldRank == 0) {
			std::cout << "[ProjectedJacobi+Nesterov] Converged in " << ctx.iteration << " iterations, final error = " << ctx.err_global << std::endl;
		}
	}

	// Store last impulses into cache for next step warmstarting
	if (m_wsProvider != nullptr) {
		const auto& contactsAll = m_dyn.contacts();
		const auto& dynToOrig = m_dyn.dynamicContactToOriginalAll();
		std::vector<ContactImpulse> lastImp(contactsAll.size());
		for (const auto& rows : ctx.contactRowGroups) {
			if (rows.empty()) continue;
			const int rn = rows[0];
			const int origIdx = (rn >= 0 && rn < (int)dynToOrig.size()) ? dynToOrig[(size_t)rn] : -1;
			if (origIdx < 0 || origIdx >= (int)lastImp.size()) continue;
			ContactImpulse ci{};
			ci.pn = p[rn];
			if (rows.size() > 1) ci.pt1 = p[rows[1]]; else ci.pt1 = (real_t)0;
			if (rows.size() > 2) ci.pt2 = p[rows[2]]; else ci.pt2 = (real_t)0;
			lastImp[(size_t)origIdx] = ci;
		}
		m_wsProvider->store(contactsAll, lastImp);
	}

	// Final synchronization: replicate all body velocities to every rank
	cardillo::comm::Communication::replicateAllBodyVelocitiesConcat(u, ctx.res, ctx.bodyOffsets);

	return u;
}

} // namespace cardillo::solver