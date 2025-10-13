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
