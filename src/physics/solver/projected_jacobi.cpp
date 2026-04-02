#include "projected_jacobi.hpp"
#include <iostream>
#include <algorithm>
#include <limits>
#include <cmath>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <filesystem>

namespace cardillo::solver {

// Precompute R = alpha / diag(G + compliance I) for relevant contacts
// G = W * Minv * W^T, with Minv diagonal (given as MinvDiag)
static inline VectorXr precompute_Rdiag_sparse(
	const Eigen::SparseMatrix<real_t, Eigen::RowMajor>& W,
	const VectorXr& MinvDiag,
	real_t alpha)
{
	const int C = (int)W.rows();
	VectorXr R = VectorXr::Zero(C);

	for (int cid = 0; cid < C; ++cid) {
		if (cid < 0 || cid >= C) continue;
		real_t Dii = 0;
		for (Eigen::SparseMatrix<real_t, Eigen::RowMajor>::InnerIterator it(W, cid); it; ++it) {
			const int col = it.col();
			const real_t w = it.value();
			Dii += w * w * MinvDiag[col];
		}
		R[cid] = (Dii > (real_t)0) ? (alpha / Dii) : (real_t)0;
	}
	return R;
}

// Small context to avoid re-passing many references for each sweep
struct PJIterContext {
	const Eigen::SparseMatrix<real_t, Eigen::RowMajor>& W;
	cardillo::physics::DynamicsAssembler* dyn; // pointer to assembler (for S-solve etc.)
	cardillo::physics::assembly::PjAssembler* assembler; // pointer to PJAssembler (for S-solve etc.)
	std::vector<int> bodyOffsets;
	VectorXr Rdiag;
	// Grouping of contact rows by original contact index (for friction cone projection)
	std::vector<std::vector<int>> contactRowGroups; // each group has 1 (normal only) or 3 rows (n,t1,t2) sorted by row id
	std::vector<real_t> groupMu; // same size as contactRowGroups: per-contact friction mu
	// Iteration runtime state and buffers
	VectorXr y_contact;   // contact-space buffer (size = #contacts)
	VectorXr tmp_concat;  // velocity-space buffer (size = concat dofs)

	// Dimensions
	int nV{0};  // velocity dofs
	int nP{0};  // percussion dofs
	int nS{0};  // spring dofs (lagrange multipliers)
	int nD{0};  // damper dofs (gamma rows)

	// Controls and diagnostics for the loop
	int maxIterations{1000000};
	real_t tol{(real_t)1e-5};
	real_t relax{(real_t)1};
	bool debug{false};
	real_t eps_rel{(real_t)0};
	// Outputs
	real_t err_global{std::numeric_limits<real_t>::max()};
	index_t iteration{0};
};

// Compact holder for grouped contact rows and friction coefficients
struct ContactGroups {
	std::vector<std::vector<int>> rows; // 1 or 3 rows per contact
	std::vector<real_t> mu;             // same length as rows
};

static inline ContactGroups build_contact_groups(const cardillo::physics::DynamicsAssembler& dyn)
{
	ContactGroups cg;
	const auto& contactsAll = dyn.contacts();
	cg.rows.reserve(contactsAll.size());
	cg.mu.reserve(contactsAll.size());
	for (const auto& c : contactsAll) {
		if (c.impulse_base_index < 0) continue;
		std::vector<int> rows;
		rows.reserve((size_t)c.impulse_size);
		rows.push_back(c.impulse_base_index);
		if (c.impulse_size > 1) rows.push_back(c.impulse_base_index + 1);
		if (c.impulse_size > 2) rows.push_back(c.impulse_base_index + 2);
		std::sort(rows.begin(), rows.end());
		cg.rows.push_back(std::move(rows));
		cg.mu.push_back(c.friction_mu);
	}
	return cg;
}

static inline PJIterContext build_context(
	cardillo::physics::DynamicsAssembler& dyn,
	cardillo::physics::assembly::PjAssembler& m_assembler,
	real_t alpha)
{
	const auto& W = dyn.W().asSparseRowMajor();
	const int C = (int)W.rows();

	PJIterContext ctx{W, &dyn, &m_assembler, std::vector<int>{}, VectorXr()};
	ctx.bodyOffsets = dyn.bodyVelOffsets();
	// Use Minv diagonal from DynamicsAssembler when computing diagonal of G = W * Minv * W^T
	ctx.Rdiag = precompute_Rdiag_sparse(W, dyn.MinvDiag(), alpha);
	// Populate dimension helpers: velocity dofs, percussion (contact) rows, and spring dofs
	const int dofConcat = (ctx.bodyOffsets.empty() ? 0 : ctx.bodyOffsets.back());
	ctx.nV = dofConcat;
	ctx.nP = C;
	ctx.nS = (int)dyn.Cdiag().size();
	ctx.nD = (int)dyn.Adiag().size();
	// Build contact row groups and per-contact mu
	const ContactGroups cg = build_contact_groups(dyn);
	ctx.contactRowGroups = cg.rows;
	ctx.groupMu = cg.mu;
	// allocate buffers
	ctx.y_contact = VectorXr::Zero(C);
	ctx.tmp_concat = VectorXr::Zero(dofConcat);
	return ctx;
}

// Compute provisional gradient step z = p_in - R * (W * u_in)
static inline VectorXr provisional_step(PJIterContext& ctx, const VectorXr& u_in, const VectorXr& p_in)
{
	const auto& W = ctx.W; const auto& Rdiag = ctx.Rdiag;
	ctx.y_contact.noalias() = W * u_in;
	VectorXr z = VectorXr::Zero(p_in.size());
	for (int cid = 0; cid < p_in.size(); ++cid) {
		z[cid] = p_in[cid] - Rdiag[cid] * ctx.y_contact[cid];
	}
	return z;
}

// Project per-contact groups: 1-row -> unilateral; 3-rows -> Coulomb cone (disk in tangential plane)
static inline void project_groups(const PJIterContext& ctx, const VectorXr& z, const VectorXr& p_old, VectorXr& p_proj)
{
	// TODO: is this valid?
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
static inline void relax_exchange_update(PJIterContext& ctx, const VectorXr& x_old, 
										 const VectorXr& p_old, VectorXr& x_new, VectorXr& p_new)
{
	const auto& W = ctx.W; cardillo::physics::DynamicsAssembler* dyn = ctx.dyn;

	// Use effective-mass solve (S^{-1}) to update u
	// Build extended RHS [tmp_concat; 0; 0] and solve full S * delta_x = rhs_ext, then apply only velocity part
	VectorXr rhs_ext = VectorXr::Zero(ctx.nV + ctx.nS + ctx.nD);
	rhs_ext.segment(0,  ctx.nV) = W.transpose() * (p_new - p_old);
	x_new = x_old + ctx.assembler->solveS(rhs_ext);
}

// One Projected-Jacobi sweep: update p with current u, then compute new u
static inline void pj_sweep(PJIterContext& ctx, const VectorXr& x_in, 
						    const VectorXr& p_in, VectorXr& x_out,
							VectorXr& p_out)
{
	const VectorXr z = provisional_step(ctx, x_in.segment(0, ctx.nV), p_in);
	project_groups(ctx, z, p_in, p_out);
	relax_exchange_update(ctx, x_in, p_in, x_out, p_out);
}

// Compute global L2 norm of difference across local body segments
static inline real_t global_segment_norm(const PJIterContext& ctx, const VectorXr& a, const VectorXr& b, real_t eps_abs) {
	const auto& bo = ctx.bodyOffsets;
	double loc_sum = 0.0;
	for (int br = 0; br + 1 < (int)bo.size(); ++br) {
		int off = bo[(size_t)br]; int n = bo[(size_t)br+1] - off;
		if (n > 0) {
			Eigen::Map<const VectorXr> va(a.data()+off, n), vb(b.data()+off, n);
			if (ctx.eps_rel > (real_t)0) {
				const VectorXr scale = VectorXr::Constant(n, eps_abs)
					+ ctx.eps_rel * va.cwiseAbs().cwiseMax(vb.cwiseAbs());
				const VectorXr diff = (va - vb).cwiseQuotient(scale);
				loc_sum += diff.squaredNorm();
			} else {
				const VectorXr diff = (va - vb) / eps_abs;
				loc_sum += diff.squaredNorm();
			}
		}
	}

	if (std::isnan(loc_sum)) {
		std::cerr << "[ProjectedJacobi] Warning: NaN detected in global_segment_norm" << std::endl;
	}

	const real_t denom = (ctx.nV > 0) ? (real_t)std::sqrt((double)ctx.nV) : (real_t)1;
	return (real_t)std::sqrt(loc_sum) / denom;
}

static inline real_t global_contact_norm(const PJIterContext& ctx, const VectorXr& a, const VectorXr& b) {
	(void)ctx;
	double loc_sum = 0.0;
	for (int cid = 0; cid < a.size(); ++cid) {
		real_t diff = a[cid] - b[cid];
		loc_sum += static_cast<double>(diff) * static_cast<double>(diff);
	}
	return (real_t)std::sqrt(loc_sum);
}

static inline real_t global_velocity_l2_norm(const PJIterContext& ctx, const VectorXr& a, const VectorXr& b) {
	const auto& bo = ctx.bodyOffsets;
	double loc_sum = 0.0;
	for (int br = 0; br + 1 < (int)bo.size(); ++br) {
		int off = bo[(size_t)br]; int n = bo[(size_t)br+1] - off;
		if (n > 0) {
			Eigen::Map<const VectorXr> va(a.data()+off, n), vb(b.data()+off, n);
			const VectorXr diff = (va - vb);
			loc_sum += diff.squaredNorm();
		}
	}
	return (real_t)std::sqrt(loc_sum);
}

static inline double global_segment_dot(const PJIterContext& ctx, const VectorXr& a, const VectorXr& b) {
	const auto& bo = ctx.bodyOffsets;
	double loc = 0.0;
	for (int br = 0; br + 1 < (int)bo.size(); ++br) {
		int off = bo[(size_t)br]; int n = bo[(size_t)br+1] - off;
		if (n > 0) {
			Eigen::Map<const VectorXr> va(a.data()+off, n), vb(b.data()+off, n);
			loc += va.dot(vb);
		}
	}
	return loc;
}

// Standard fixed-point iteration loop
static inline void standard_loop(PJIterContext& ctx,
				 VectorXr& u, VectorXr& p) {
	ctx.err_global = std::numeric_limits<real_t>::max();
	ctx.iteration = 0;
	VectorXr u_prev = u;
	VectorXr p_prev = p;
	while (ctx.iteration < ctx.maxIterations) {
		pj_sweep(ctx, u_prev, p_prev, u, p);
		ctx.err_global = global_segment_norm(ctx, u, u_prev, ctx.tol);
		++ctx.iteration;
		if (ctx.err_global <= (real_t)1) break;
		if (ctx.debug && (ctx.iteration % 1000 == 0)) {
			std::cout << "[ProjectedJacobi] Iteration " << ctx.iteration << ", error = " << ctx.err_global << std::endl;
		}
		u_prev = u;
		p_prev = p;
	}
}

// Nesterov-accelerated loop
static inline void nesterov_loop(PJIterContext& ctx,
				VectorXr& u, VectorXr& p,
				double beta_threshold, int restart_limit) {
	ctx.err_global = std::numeric_limits<real_t>::max();
	ctx.iteration = 0;
	VectorXr xuk = u;
	VectorXr xpk = p;
	VectorXr xuk1 = u;
	VectorXr xpk1 = p;
	VectorXr yuk = u;
	VectorXr ypk = p;
	double thk = 1.0;
	real_t err_prev = std::numeric_limits<real_t>::infinity();
	int restarts = 0;
	bool momentum_disabled = false;
	
	while (ctx.iteration < ctx.maxIterations) {
		pj_sweep(ctx, yuk, ypk, xuk1, xpk1);

		ctx.err_global = global_segment_norm(ctx, xuk1, yuk, ctx.tol);
		++ctx.iteration;
		double thk1 = 0.5 * (1.0 + std::sqrt(4.0 * thk * thk + 1.0));
		double betak1 = (thk - 1.0) / thk1;

		bool restart = false;
		std::string resetReason;
		if (!std::isfinite((double)ctx.err_global)) { restart = true; resetReason = "nan_error"; }
		else if (err_prev < std::numeric_limits<real_t>::infinity() && ctx.err_global > (real_t)1.05 * err_prev) { restart = true; resetReason = "error_increase"; }
		if (!restart && betak1 > beta_threshold) { restart = true; resetReason = "beta_threshold"; }
		if (!restart) {
			VectorXr a = yuk - xuk1; VectorXr b = xuk1 - xuk;
			double d = global_segment_dot(ctx, a, b);
			if (d > 0.0) { restart = true; resetReason = "direction_conflict"; }
		}
		if (!restart && (!std::isfinite(betak1) || betak1 < 0.0 || betak1 > 1.0)) { restart = true; resetReason = "beta_invalid"; }
		if (ctx.err_global <= (real_t)1) break;

		if (restart) {
			yuk = xuk1;
			ypk = xpk1;
			xuk = xuk1;
			xpk = xpk1;
			thk = 1.0;
			if (++restarts >= restart_limit) momentum_disabled = true;
		} else {
			if (momentum_disabled) {
				yuk = xuk1;
				ypk = xpk1;
				thk = 1.0;
			} else {
				yuk = xuk1 + (real_t)betak1 * (xuk1 - xuk);
				ypk = xpk1 + (real_t)betak1 * (xpk1 - xpk);
				thk = thk1;
			}
			xuk = xuk1;
			xpk = xpk1;
		}
		err_prev = ctx.err_global;

		if (ctx.debug && (ctx.iteration % 1000 == 0)) {
			std::cout << "[ProjectedJacobi+Nesterov] Iter " << ctx.iteration << ", error = " << ctx.err_global
					  << ", beta = " << betak1 << std::endl;
		}
	}

	u = xuk1;
	p = xpk1;
}

VectorXr ProjectedJacobiSolver::solve(real_t dt, real_t theta)  {
	m_dyn.updateStateDependentTerms();
	m_assembler.buildAndFactorS(dt, theta);

	auto sc_solve = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::ProjectedJacobi);

	// Use sparse W and the effective-mass S (assembled & factorized in DynamicsAssembler)
	const auto& Wref = m_dyn.W().asSparseRowMajor();
	const int C = (int)Wref.rows();
	const int Nv = (Wref.cols() > 0) ? Wref.cols() : 0;

	// Build iteration context and precomputations
	PJIterContext ctx = build_context(m_dyn, m_assembler, alpha());

	// Initialize impulses
	VectorXr p = VectorXr::Zero(C);
	if (m_warmStart) {
		cardillo::solver::WarmstartProvider::applyWarmstart(p, m_dyn);
	}

	// Initialize state vector x consisting of velocities and lagrange multipliers (springs)
	VectorXr rhs = m_assembler.rhs(dt, theta);
	rhs.segment(0, Nv).noalias() += Wref.transpose() * p;

	VectorXr x = m_assembler.solveS(rhs);
	
	if (C == 0 || Wref.nonZeros() == 0) {
		m_last_iters= 0;
		m_lastError = (real_t)0;
		
		// Track spring and damper forces, return velocity.
		const int nSprings = (int)m_dyn.Cdiag().size();
		const int nDampers = (int)m_dyn.Adiag().size();
		if (nSprings > 0) m_dyn.setLambda_g(x.segment(Nv, nSprings)); 
		if (nDampers > 0) m_dyn.setLambda_gamma(x.segment(Nv + nSprings, nDampers));
		return x.segment(0, Nv);
	}

	// Configure iteration context
	ctx.maxIterations = m_maxIterations;
	ctx.tol = m_cfg.pj_tol_abs;
	ctx.relax = m_relax;
	ctx.debug = m_debug;
	ctx.eps_rel = m_epsRel;

	if (!m_useNesterov) {
		standard_loop(ctx, x, p);
		if (m_debug) {
			std::cout << "[ProjectedJacobi] Converged in " << ctx.iteration << " iterations, final error = " << ctx.err_global << std::endl;
		}
	} else {
		nesterov_loop(ctx, x, p, m_nest_beta_threshold, m_nest_restart_limit);
		if (m_debug) {
			std::cout << "[ProjectedJacobi+Nesterov] Converged in " << ctx.iteration << " iterations, final error = " << ctx.err_global << std::endl;
		}
	}

	m_last_iters = static_cast<int>(ctx.iteration);
	m_lastError = ctx.err_global;

	// Store last impulses into cache for next step
	cardillo::solver::WarmstartProvider::storeImpulse(p, m_dyn);

	// Track spring and damper forces, return velocity.
	const int nSprings = (int)m_dyn.Cdiag().size();
	const int nDampers = (int)m_dyn.Adiag().size();
	if (nSprings > 0) m_dyn.setLambda_g(x.segment(Nv, nSprings)); 
    if (nDampers > 0) m_dyn.setLambda_gamma(x.segment(Nv + nSprings, nDampers));
	return x.segment(0, Nv);
}

} // namespace cardillo::solver