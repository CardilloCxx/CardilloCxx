#pragma once

#include <cstdint>
#include <unordered_map>
#include <vector>
#include <entt/entt.hpp>
#include "../../config/config.hpp"
#include "../../misc/block_sparse_ldlt.hpp"
#include "../../misc/types.hpp"
#include "dynamics_assembler.hpp"

namespace cardillo::physics::assembly {

// One row-block of the condensed system: a spring, a damper, a frictionless contact (dim==1), or
// a frictional contact (dim==3, normal+2 tangential). Everything here is a small dense
// per-body/per-block quantity -- no global sparse matrix is ever built or referenced.
struct RowBlock {
    enum class Kind : uint8_t { Spring, Damper, ContactFrictionless, ContactFrictional };

    Kind kind{Kind::Spring};
    int dim{0};     // rows in this block: 1..6 (Spring/Damper), 1 (frictionless), 3 (frictional)
    int offset{0};  // row offset in the packed [springs | dampers | frictionless | frictional] layout

    int bodyIndexA{-1}, bodyIndexB{-1};  // -1 if static/absent
    int aOff{0}, aDof{0};
    int bOff{0}, bDof{0};

    MatrixXXr Ja;  // dim x aDof (empty if aDof==0)
    MatrixXXr Jb;  // dim x bDof (empty if bDof==0)

    // dim x dim mass-only Delassus block: Ja*MinvA*Ja^T + Jb*MinvB*Jb^T. MinvA/MinvB are normally
    // diag(MinvDiag) -- except for a body with an active implicit-gyroscopic override
    // (CondensedTopology::gyroMinvBlocks), where it's that body's (generally non-symmetric) full
    // block instead, making Gii itself generally non-symmetric too.
    MatrixXXr Gii;
    VectorXr complianceDiag;   // size dim; 0 for contacts, Crow/(theta*dt^2) or Arow/(theta*dt) for springs/dampers
    MatrixXXr GiiInv;          // invertSmallSpd(Gii + diag(complianceDiag))

    real_t mu{0};          // ContactFrictional only
    int contactIndex{-1};  // index into DynamicsAssembler::contacts(); -1 for Spring/Damper
};

struct CondensedTopology {
    std::vector<RowBlock> blocks;
    std::vector<std::vector<int>> blocksOfBody;  // per body: incident block indices (as A or B)
    int numLambda{0};
    int springRows{0}, damperRows{0}, frictionlessRows{0}, frictionalRows{0};
    // Blocks are packed [springs | dampers | frictionless contacts | frictional contacts], so the
    // bilateral (spring+damper) blocks are exactly blocks[0:numBilateralBlocks) and the contact
    // blocks are blocks[numBilateralBlocks:] -- no separate index list needed. Distinct from
    // springRows/damperRows, which count *rows* (a block can have dim>1), not *blocks*.
    int numBilateralBlocks{0};

    // Per-body override of the (block-diagonal) inverse mass, keyed by body index -- populated only
    // for rigid bodies with an active implicit-gyroscopic correction (moreau.implicit_gyroscopy=true
    // AND the body actually has rotational dof), EMPTY otherwise (the overwhelming common case).
    // Every place that would otherwise read `MinvDiag.segment(off,dof)` as a plain diagonal must
    // check here first: a body present here has a genuinely non-diagonal, generally non-symmetric
    // effective inverse mass (translational 3x3 unchanged/diagonal, rotational 3x3 =
    // invertSmallGeneral(I - dt*Grot), see CondensedAssembler::buildTopology()'s implementation and
    // PjAssembler::buildAndFactorS() for the reference formula this mirrors). Absent from this map
    // == "use MinvDiag as before", so an empty map is exactly behaviorally equivalent to this
    // feature not existing.
    std::unordered_map<int, Matrixr<6, 6>> gyroMinvBlocks;
};

// Builds and evaluates the condensed (block-sparse, matrix-free) system in exactly the same
// mathematical convention as PgsAssembler (see docs/chapters/solvers/projected_gauss_seidel.rst):
// same rhs()/compliance formulas, same sign convention (contact normal <= 0 internally). Every
// per-block dense piece is built directly from DynamicsAssembler::constraintResults()/contacts(),
// never through Eigen::SparseMatrix.
class CondensedAssembler {
   public:
    CondensedAssembler(cardillo::physics::DynamicsAssembler& dyn, const cardillo::config::Config& cfg) : m_dyn(&dyn), m_cfg(cfg) {}

    // Builds the RowBlock list + body incidence from the current constraintResults()/contacts(),
    // and (see CondensedTopology::gyroMinvBlocks) each implicit-gyroscopic body's effective inverse
    // mass block. Independent of theta; `dt` is only needed for the gyroscopic correction (M_eff =
    // M - dt*Grot depends on it) -- geometry/mass otherwise wouldn't need it at all. Call once per
    // solve() (topology can change every step as contacts appear/disappear).
    CondensedTopology buildTopology(real_t dt) const;

    // Fills complianceDiag/GiiInv for every block (dt/theta-dependent). Mutates topo in place.
    void updateCompliance(CondensedTopology& topo, real_t dt, real_t theta) const;

    // Identical formula to PgsAssembler::ufree(), generalized for implicit-gyroscopic bodies via
    // topo.gyroMinvBlocks (see CondensedTopology's doc comment).
    VectorXr ufree(real_t dt, real_t theta, const CondensedTopology& topo) const;

    // Per-block target vector, packed in topo's row layout. A de-sparsified port of
    // PgsAssembler::rhs() -- same formulas/signs, transcribed line by line, no
    // Eigen::SparseMatrix/TripletMatrix::asSparse() calls anywhere in this class.
    VectorXr rhs(const CondensedTopology& topo, real_t dt, real_t theta, const VectorXr& u_free) const;

    // Builds and factors the block-sparse LDLT of Sbb = the bilateral-only (springs+dampers)
    // sub-block of S = W*Minv*W^T + diag(C), restricted to bilateral RowBlocks
    // (blocks[0:numBilateralBlocks)) and their shared-body coupling. Used by CondensedSolver's
    // condensed.true_schur=true path to eliminate the compliant chain exactly every outer
    // iteration instead of relying on the Gauss-Seidel/Jacobi/colored sweep to slowly diffuse
    // through it -- see docs/chapters/solvers/condensed.rst. Call after updateCompliance() (needs
    // each bilateral block's Gii+complianceDiag). Natural order (block index order) is used as the
    // elimination order: zero fill-in for every current example scene's compliant network except
    // `hangbridge`'s branching tripod/deck topology, where correctness still holds, just with
    // some fill-in. Builds with BlockSparseLDLT's symmetric=false (block-LU) mode -- roughly double
    // the Schur-update cost -- only if topo.gyroMinvBlocks makes some bilateral-graph coupling
    // genuinely non-symmetric; otherwise (every scene without an active implicit-gyroscopic body in
    // its compliant chain) uses the default symmetric=true mode, unchanged in cost from before this
    // capability existed.
    cardillo::misc::BlockSparseLDLT buildBilateralFactorization(const CondensedTopology& topo) const;

   private:
    cardillo::physics::DynamicsAssembler* m_dyn;
    cardillo::config::Config m_cfg;

    // Cache of the last elimination order computed by buildBilateralFactorization(), keyed on the
    // bilateral graph's *structure* (dims + edge node pairs -- not the numeric Gii/complianceDiag
    // values, which change every call regardless). The bilateral (spring+damper) topology is
    // static for a scene's lifetime in every current example (constraints aren't created/destroyed
    // at runtime), so this cache hits every call after the first and turns the O(n^2) minimum-
    // degree symbolic pass into a one-time cost instead of a per-timestep one. `mutable`: this is
    // memoization behind an otherwise-const interface, not an observable behavior change -- see
    // buildBilateralFactorization()'s doc comment for the correctness argument (a cache miss always
    // falls back to recomputing from scratch, so a stale/wrong cache can only cost performance,
    // never correctness).
    mutable std::vector<int> m_cachedBilateralDims;
    mutable std::vector<std::array<int, 2>> m_cachedBilateralEdgeNodes;
    mutable std::vector<int> m_cachedBilateralOrder;
    mutable bool m_hasCachedBilateralOrder{false};
};

}  // namespace cardillo::physics::assembly
