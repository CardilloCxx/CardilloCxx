#include "physics_system.hpp"

namespace cardillo {

// Construction and configuration
PhysicsSystem::PhysicsSystem() {
    PetscPrintf(PETSC_COMM_WORLD, "Hello from PETSc + Eigen!\n");
    m_gravity = Vector3r(0, 0, -9.81);
}

void PhysicsSystem::setGravity(const Vector3r& g) { m_gravity = g; m_forces_dirty = true; }

// Entity creation
index_t PhysicsSystem::addPointMass(real_t mass, const Vector3r& x0, const Vector3r& v0, real_t radius) {
    auto e = m_reg.create();
    m_reg.emplace<C_PhysicsObject>(e);
    m_reg.emplace<C_PointMassTag>(e);
    m_reg.emplace<C_Collidable>(e);
    m_reg.emplace<C_VisualObject>(e);
    m_reg.emplace<C_PointVisualTag>(e);
    m_reg.emplace<C_Mass>(e, C_Mass{mass});
    m_reg.emplace<C_Position3>(e, C_Position3{x0});
    m_reg.emplace<C_LinearVelocity3>(e, C_LinearVelocity3{v0});
    m_reg.emplace<C_Radius>(e, C_Radius{radius});
    m_structure_dirty = true;
    return static_cast<index_t>(entt::to_integral(e));
}

entt::entity PhysicsSystem::createRigidVisualEntity_(const Vector3r& center) {
    auto e = m_reg.create();
    m_reg.emplace<C_RigidBodyTag>(e);
    m_reg.emplace<C_Collidable>(e);
    m_reg.emplace<C_VisualObject>(e);
    m_reg.emplace<C_Position3>(e, C_Position3{center});
    return e;
}

index_t PhysicsSystem::addRigidBody(const Plane& p) {
    auto e = createRigidVisualEntity_(p.center);
    m_reg.emplace<C_PlaneVisualTag>(e);
    m_reg.emplace<C_Plane>(e, C_Plane{p.normal, p.up, p.sizeX, p.sizeY});
    return static_cast<index_t>(entt::to_integral(e));
}

index_t PhysicsSystem::addRigidBody(const Cube& c) {
    auto e = createRigidVisualEntity_(c.center);
    m_reg.emplace<C_CubeVisualTag>(e);
    m_reg.emplace<C_Cube>(e, C_Cube{c.halfExtents, c.R});
    return static_cast<index_t>(entt::to_integral(e));
}

} // namespace cardillo

// ---------- Dynamics getters implementations ----------
namespace cardillo {

MatrixXXr PhysicsSystem::getMass(entt::entity e) const {
    // For point masses: 3x3 diagonal with mass on the diagonal
    if (m_reg.any_of<C_PointMassTag, C_Mass>(e)) {
        const auto& m = m_reg.get<C_Mass>(e).m;
        MatrixXXr M = MatrixXXr::Zero(3,3);
        M.diagonal().setConstant(m);
        return M;
    }
    // Default: zero-sized (no DOFs)
    return MatrixXXr(0,0);
}

MatrixXXr PhysicsSystem::getMassInverse(entt::entity e) const {
    if (m_reg.any_of<C_PointMassTag, C_Mass>(e)) {
        const auto& m = m_reg.get<C_Mass>(e).m;
        MatrixXXr Minv = MatrixXXr::Zero(3,3);
        const real_t invm = (m != (real_t)0) ? (real_t)1 / m : (real_t)0;
        Minv.diagonal().setConstant(invm);
        return Minv;
    }
    return MatrixXXr(0,0);
}

VectorXr PhysicsSystem::getPosition(entt::entity e) const {
    if (m_reg.any_of<C_PointMassTag, C_Position3>(e)) {
        const auto& x = m_reg.get<C_Position3>(e).value;
        VectorXr q(3); q << x[0], x[1], x[2];
        return q;
    }
    return VectorXr(0);
}

VectorXr PhysicsSystem::getVelocity(entt::entity e) const {
    if (m_reg.any_of<C_PointMassTag, C_LinearVelocity3>(e)) {
        const auto& v = m_reg.get<C_LinearVelocity3>(e).value;
        VectorXr out(3); out << v[0], v[1], v[2];
        return out;
    }
    return VectorXr(0);
}

VectorXr PhysicsSystem::getForce(entt::entity e) const {
    if (m_reg.any_of<C_PointMassTag, C_Mass>(e)) {
        const auto& m = m_reg.get<C_Mass>(e).m;
        Vector3r fg = m * m_gravity;
        VectorXr out(3); out << fg[0], fg[1], fg[2];
        return out;
    }
    return VectorXr(0);
}

int PhysicsSystem::numBodies() const {
    int count = 0;
    auto view = m_reg.view<C_BodyIndex, C_PhysicsObject>();
    for (auto e : view) (void)e, ++count;
    return count;
}

} // namespace cardillo