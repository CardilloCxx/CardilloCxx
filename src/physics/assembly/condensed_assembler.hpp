#pragma once

#include <cstdint>
#include <vector>
#include <entt/entt.hpp>
#include "../../config/config.hpp"
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

    MatrixXXr Gii;             // dim x dim mass-only Delassus block: Ja*diag(MinvA)*Ja^T + Jb*diag(MinvB)*Jb^T
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
};

// Builds and evaluates the condensed (block-sparse, matrix-free) system in exactly the same
// mathematical convention as PgsAssembler (see docs/chapters/solvers/projected_gauss_seidel.rst):
// same rhs()/compliance formulas, same sign convention (contact normal <= 0 internally). Every
// per-block dense piece is built directly from DynamicsAssembler::constraintResults()/contacts(),
// never through Eigen::SparseMatrix.
class CondensedAssembler {
   public:
    CondensedAssembler(cardillo::physics::DynamicsAssembler& dyn, const cardillo::config::Config& cfg) : m_dyn(&dyn), m_cfg(cfg) {}

    // Builds the RowBlock list + body incidence from the current constraintResults()/contacts().
    // Pure geometry+mass, independent of dt/theta; call once per solve() (topology can change
    // every step as contacts appear/disappear).
    CondensedTopology buildTopology() const;

    // Fills complianceDiag/GiiInv for every block (dt/theta-dependent). Mutates topo in place.
    void updateCompliance(CondensedTopology& topo, real_t dt, real_t theta) const;

    // Identical formula to PgsAssembler::ufree().
    VectorXr ufree(real_t dt, real_t theta) const;

    // Per-block target vector, packed in topo's row layout. A de-sparsified port of
    // PgsAssembler::rhs() -- same formulas/signs, transcribed line by line, no
    // Eigen::SparseMatrix/TripletMatrix::asSparse() calls anywhere in this class.
    VectorXr rhs(const CondensedTopology& topo, real_t dt, real_t theta, const VectorXr& u_free) const;

   private:
    cardillo::physics::DynamicsAssembler* m_dyn;
    cardillo::config::Config m_cfg;
};

}  // namespace cardillo::physics::assembly
