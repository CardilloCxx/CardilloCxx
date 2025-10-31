// Implementation of generalized spring constraint deformation and Jacobians
#include "force_interaction.hpp"
#include <limits>
#include <cmath>

namespace cardillo {
namespace physics {

// Manager methods
size_t ForceInteractionManager::addConstraint(const SpringConstraint& c)
{
    m_constraints.push_back(c);
    return m_constraints.size() - 1;
}

void ForceInteractionManager::removeConstraint(size_t idx)
{
    if (idx < m_constraints.size()) {
        m_constraints.erase(m_constraints.begin() + static_cast<std::ptrdiff_t>(idx));
    }
}

// Per-constraint computations
void SpringConstraint::recomputeDeformations()
{
    // Registry must be provided at construction (Option B)
    entt::registry* reg = registry;
    // Default world values (static world)
    Vector3r posA = Vector3r::Zero();
    Vector3r posB = Vector3r::Zero();
    Quaternion4r qa = Quaternion4r::Identity();
    Quaternion4r qb = Quaternion4r::Identity();
    Vector3r vA = Vector3r::Zero();
    Vector3r vB = Vector3r::Zero();
    Vector3r wA = Vector3r::Zero();
    Vector3r wB = Vector3r::Zero();

    if (reg && bodyA != entt::null) {
        if (reg->all_of<cardillo::PhysicsSystem::C_Position3>(bodyA)) posA = reg->get<cardillo::PhysicsSystem::C_Position3>(bodyA).value;
        if (reg->all_of<cardillo::PhysicsSystem::C_Orientation>(bodyA)) qa = reg->get<cardillo::PhysicsSystem::C_Orientation>(bodyA).q;
        if (reg->all_of<cardillo::PhysicsSystem::C_LinearVelocity3>(bodyA)) vA = reg->get<cardillo::PhysicsSystem::C_LinearVelocity3>(bodyA).value;
        if (reg->all_of<cardillo::PhysicsSystem::C_AngularVelocity3>(bodyA)) wA = reg->get<cardillo::PhysicsSystem::C_AngularVelocity3>(bodyA).value;
    }
    if (reg && bodyB != entt::null) {
        if (reg->all_of<cardillo::PhysicsSystem::C_Position3>(bodyB)) posB = reg->get<cardillo::PhysicsSystem::C_Position3>(bodyB).value;
        if (reg->all_of<cardillo::PhysicsSystem::C_Orientation>(bodyB)) qb = reg->get<cardillo::PhysicsSystem::C_Orientation>(bodyB).q;
        if (reg->all_of<cardillo::PhysicsSystem::C_LinearVelocity3>(bodyB)) vB = reg->get<cardillo::PhysicsSystem::C_LinearVelocity3>(bodyB).value;
        if (reg->all_of<cardillo::PhysicsSystem::C_AngularVelocity3>(bodyB)) wB = reg->get<cardillo::PhysicsSystem::C_AngularVelocity3>(bodyB).value;
    }

    // World rotations
    RA = qa.toRotationMatrix();
    RB = qb.toRotationMatrix();

    // World attachment points
    xA = posA + RA * rA_local;
    xB = posB + RB * rB_local;

    // Velocities at attachments (world-frame)
    cardillo::Vector3r wA_world = RA * wA;
    cardillo::Vector3r wB_world = RB * wB;
    cardillo::Vector3r vA_attach = vA + wA_world.cross(RA * rA_local);
    cardillo::Vector3r vB_attach = vB + wB_world.cross(RB * rB_local);


    // Transform local (body-A) K/D into world frame so assemblers can
    // read K/D directly. The 6x6 transform is block-diagonal with RA
    // for both translational and rotational parts.
    Matrix6r T = Matrix6r::Zero();
    T.block<3,3>(0,0) = RA;
    T.block<3,3>(3,3) = RA;
    K = T * K_local * T.transpose();
    D = T * D_local * T.transpose();
}

// Constructor: capture registry pointer and initialize rest pose from current ECS state
SpringConstraint::SpringConstraint(entt::registry& reg,
                                   entt::entity a,
                                   entt::entity b,
                                   const Vector3r& rA_loc,
                                   const Vector3r& rB_loc)
    : bodyA(a), bodyB(b), rA_local(rA_loc), rB_local(rB_loc), registry(&reg)
{
    // Default world values
    Vector3r posA = Vector3r::Zero();
    Vector3r posB = Vector3r::Zero();
    Quaternion4r qa = Quaternion4r::Identity();
    Quaternion4r qb = Quaternion4r::Identity();

    if (bodyA != entt::null) {
        if (reg.all_of<cardillo::PhysicsSystem::C_Position3>(bodyA)) posA = reg.get<cardillo::PhysicsSystem::C_Position3>(bodyA).value;
        if (reg.all_of<cardillo::PhysicsSystem::C_Orientation>(bodyA)) qa = reg.get<cardillo::PhysicsSystem::C_Orientation>(bodyA).q;
    }
    if (bodyB != entt::null) {
        if (reg.all_of<cardillo::PhysicsSystem::C_Position3>(bodyB)) posB = reg.get<cardillo::PhysicsSystem::C_Position3>(bodyB).value;
        if (reg.all_of<cardillo::PhysicsSystem::C_Orientation>(bodyB)) qb = reg.get<cardillo::PhysicsSystem::C_Orientation>(bodyB).q;
    }

    // Initialize cached transforms/attachments
    RA = qa.toRotationMatrix();
    RB = qb.toRotationMatrix();
    xA = posA + RA * rA_local;
    xB = posB + RB * rB_local;

    // Rest pose measured from A->B in world frame
    rest_translation = xB - xA;
    cardillo::Quaternion4r qrel = qa.conjugate() * qb;
    qrel.normalize();
    cardillo::AngleAxis3r aa(qrel);
    rest_rotation = RA * (aa.axis() * aa.angle()); // <-- rotate axis from A-local to world
}

void SpringConstraint::computeWBlocks(Matrix6r & WgA_out, Matrix6r & WgB_out,
                                      Matrix6r & WgammaA_out, Matrix6r & WgammaB_out) const
{
    // Zero outputs first
    WgA_out.setZero();
    WgB_out.setZero();
    WgammaA_out.setZero();
    WgammaB_out.setZero();

    // Precompute world lever arms
    cardillo::Vector3r rA_w = RA * rA_local;
    cardillo::Vector3r rB_w = RB * rB_local;
    cardillo::Matrix33r Sra = cardillo::skew_from_vector(rA_w);
    cardillo::Matrix33r Srb = cardillo::skew_from_vector(rB_w);

    // WgA: columns [vA (3) | wA (3)]
    WgA_out.block<3,3>(0,0) = -cardillo::Matrix33r::Identity(); // d g_pos / d vA
    WgA_out.block<3,3>(0,3) =  Sra;                             // d g_pos / d wA
    WgA_out.block<3,3>(3,3) = -cardillo::Matrix33r::Identity(); // d g_rot / d wA

    // WgB: columns [vB (3) | wB (3)]
    WgB_out.block<3,3>(0,0) =  cardillo::Matrix33r::Identity(); // d g_pos / d vB
    WgB_out.block<3,3>(0,3) = -Srb;                              // d g_pos / d wB
    WgB_out.block<3,3>(3,3) =  cardillo::Matrix33r::Identity(); // d g_rot / d wB

    // Wgamma: for viscous damping that acts on gdot, the Jacobian is identical to Wg
    WgammaA_out = WgA_out;
    WgammaB_out = WgB_out;
}

void SpringConstraint::addTranslationalSpring(const cardillo::Vector3r& axis_world, real_t k, real_t d)
{
    // Store axis in body-A local coordinates so it rotates with body A.
    cardillo::Vector3r n_world = axis_world.normalized();
    cardillo::Vector3r n_local = RA.transpose() * n_world;
    cardillo::Matrix33r P_local = n_local * n_local.transpose();
    K_local.block<3,3>(0,0) += k * P_local;
    D_local.block<3,3>(0,0) += d * P_local;
}

void SpringConstraint::addRotationalSpring(const cardillo::Vector3r& axis_world, real_t k, real_t d)
{
    // Rotational spring axis expressed in body-A local frame
    cardillo::Vector3r n_world = axis_world.normalized();
    cardillo::Vector3r n_local = RA.transpose() * n_world;
    cardillo::Matrix33r P_local = n_local * n_local.transpose();
    K_local.block<3,3>(3,3) += k * P_local;
    D_local.block<3,3>(3,3) += d * P_local;
}

void SpringConstraint::fixTranslation()
{
    real_t inf_like = (real_t)1e12; // large finite stiffness to approximate infinity
    K_local.block<3,3>(0,0).setConstant(inf_like);
    D_local.block<3,3>(0,0).setZero();
}

void SpringConstraint::fixRotation()
{
    real_t inf_like = (real_t)1e12;
    K_local.block<3,3>(3,3).setConstant(inf_like);
    D_local.block<3,3>(3,3).setZero();
}

} // namespace physics
} // namespace cardillo
