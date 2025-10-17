#include "physics_system.hpp"

namespace cardillo {

namespace {
// Compute diagonal inertia for a box with half-extents (body frame), mass m
inline Vector3r boxInertiaDiag(real_t m, const Vector3r& he) {
    // For a box with full extents 2*he, Ixx = 1/12 m ((2hy)^2 + (2hz)^2) = 1/3 m (hy^2 + hz^2)
    Vector3r I;
    I.x() = (real_t)1.0/3.0 * m * (he.y()*he.y() + he.z()*he.z());
    I.y() = (real_t)1.0/3.0 * m * (he.x()*he.x() + he.z()*he.z());
    I.z() = (real_t)1.0/3.0 * m * (he.x()*he.x() + he.y()*he.y());
    return I;
}
}

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

// Helper: create a purely visual/collidable rigid-like entity (no physics)
entt::entity PhysicsSystem::createRigidVisualEntity_(const Vector3r& center) {
    auto e = m_reg.create();
    m_reg.emplace<C_Collidable>(e);
    m_reg.emplace<C_VisualObject>(e);
    m_reg.emplace<C_Position3>(e, C_Position3{center});
    return e;
}

// Obstacles (static visuals)
index_t PhysicsSystem::addObstacleBody(const Plane& p) {
    auto e = createRigidVisualEntity_(p.center);
    m_reg.emplace<C_PlaneVisualTag>(e);
    m_reg.emplace<C_Plane>(e, C_Plane{p.normal, p.up, p.sizeX, p.sizeY});
    return static_cast<index_t>(entt::to_integral(e));
}

index_t PhysicsSystem::addObstacleBody(const Cube& c) {
    auto e = createRigidVisualEntity_(c.center);
    m_reg.emplace<C_CubeVisualTag>(e);
    m_reg.emplace<C_Cube>(e, C_Cube{c.halfExtents});
    m_reg.emplace<C_Orientation>(e, C_Orientation{c.q});
    return static_cast<index_t>(entt::to_integral(e));
}

// Dynamic rigid body
index_t PhysicsSystem::addRigidBody(real_t mass,
                                    const Vector3r& position,
                                    const Quaternion4r& orientation,
                                    const Vector3r& linearVelocity,
                                    const Vector3r& angularVelocity,
                                    const Cube& shape) {
    auto e = m_reg.create();
    m_reg.emplace<C_PhysicsObject>(e);
    m_reg.emplace<C_RigidBodyTag>(e);
    m_reg.emplace<C_Collidable>(e);
    m_reg.emplace<C_VisualObject>(e);
    m_reg.emplace<C_Position3>(e, C_Position3{position});
    m_reg.emplace<C_Orientation>(e, C_Orientation{orientation});
    m_reg.emplace<C_LinearVelocity3>(e, C_LinearVelocity3{linearVelocity});
    m_reg.emplace<C_AngularVelocity3>(e, C_AngularVelocity3{angularVelocity});
    m_reg.emplace<C_Mass>(e, C_Mass{mass});
    m_reg.emplace<C_CubeVisualTag>(e);
    m_reg.emplace<C_Cube>(e, C_Cube{shape.halfExtents});
    m_structure_dirty = true;
    return static_cast<index_t>(entt::to_integral(e));
}

// ---------- Dynamics getters ----------

MatrixXXr PhysicsSystem::getMass(entt::entity e) const {
    // Rigid body first
    if (m_reg.all_of<C_RigidBodyTag, C_PhysicsObject, C_Mass, C_Cube>(e)) {
        const real_t m = m_reg.get<C_Mass>(e).m;
        const Vector3r he = m_reg.get<C_Cube>(e).halfExtents;
        const Vector3r Idiag = boxInertiaDiag(m, he);
        MatrixXXr M = MatrixXXr::Zero(6,6);
        M(0,0) = m; M(1,1) = m; M(2,2) = m;
        M(3,3) = Idiag.x();
        M(4,4) = Idiag.y();
        M(5,5) = Idiag.z();
        return M;
    }
    // Point mass
    if (m_reg.all_of<C_PointMassTag, C_PhysicsObject, C_Mass>(e)) {
        const real_t m = m_reg.get<C_Mass>(e).m;
        MatrixXXr M = MatrixXXr::Zero(3,3);
        M.diagonal().setConstant(m);
        return M;
    }
    return MatrixXXr(0,0);
}

MatrixXXr PhysicsSystem::getMassInverse(entt::entity e) const {
    // Rigid body first
    if (m_reg.all_of<C_RigidBodyTag, C_PhysicsObject, C_Mass, C_Cube>(e)) {
        const real_t m = m_reg.get<C_Mass>(e).m;
        const Vector3r he = m_reg.get<C_Cube>(e).halfExtents;
        const Vector3r Idiag = boxInertiaDiag(m, he);
        MatrixXXr Minv = MatrixXXr::Zero(6,6);
        const real_t invm = (m != (real_t)0) ? (real_t)1 / m : (real_t)0;
        Minv(0,0) = invm; Minv(1,1) = invm; Minv(2,2) = invm;
        Minv(3,3) = (Idiag.x() != (real_t)0) ? (real_t)1 / Idiag.x() : (real_t)0;
        Minv(4,4) = (Idiag.y() != (real_t)0) ? (real_t)1 / Idiag.y() : (real_t)0;
        Minv(5,5) = (Idiag.z() != (real_t)0) ? (real_t)1 / Idiag.z() : (real_t)0;
        return Minv;
    }
    // Point mass
    if (m_reg.all_of<C_PointMassTag, C_PhysicsObject, C_Mass>(e)) {
        const real_t m = m_reg.get<C_Mass>(e).m;
        MatrixXXr Minv = MatrixXXr::Zero(3,3);
        const real_t invm = (m != (real_t)0) ? (real_t)1 / m : (real_t)0;
        Minv.diagonal().setConstant(invm);
        return Minv;
    }
    return MatrixXXr(0,0);
}

VectorXr PhysicsSystem::getPosition(entt::entity e) const {
    // Prefer rigid body
    if (m_reg.all_of<C_RigidBodyTag, C_PhysicsObject, C_Position3, C_Orientation>(e)) {
        const auto& x = m_reg.get<C_Position3>(e).value;
        Quaternion4r qn = m_reg.get<C_Orientation>(e).q;
        qn.normalize();
        VectorXr q(7);
        q[0] = x[0]; q[1] = x[1]; q[2] = x[2];
        q[3] = qn.w(); q[4] = qn.x(); q[5] = qn.y(); q[6] = qn.z();
        return q;
    }
    // Point mass
    if (m_reg.all_of<C_PointMassTag, C_PhysicsObject, C_Position3>(e)) {
        const auto& x = m_reg.get<C_Position3>(e).value;
        VectorXr q(3); q << x[0], x[1], x[2];
        return q;
    }
    return VectorXr(0);
}

VectorXr PhysicsSystem::getVelocity(entt::entity e) const {
    // Prefer rigid body
    if (m_reg.all_of<C_RigidBodyTag, C_PhysicsObject, C_LinearVelocity3, C_AngularVelocity3>(e)) {
        const auto& v = m_reg.get<C_LinearVelocity3>(e).value;
        const auto& w = m_reg.get<C_AngularVelocity3>(e).value; // body-frame angular velocity
        VectorXr out(6);
        out << v[0], v[1], v[2], w[0], w[1], w[2];
        return out;
    }
    // Point mass
    if (m_reg.all_of<C_PointMassTag, C_PhysicsObject, C_LinearVelocity3>(e)) {
        const auto& v = m_reg.get<C_LinearVelocity3>(e).value;
        VectorXr out(3); out << v[0], v[1], v[2];
        return out;
    }
    return VectorXr(0);
}

VectorXr PhysicsSystem::getForce(entt::entity e) const {
    // Prefer rigid body
    if (m_reg.all_of<C_RigidBodyTag, C_PhysicsObject, C_Mass>(e)) {
        const real_t m = m_reg.get<C_Mass>(e).m;
        Vector3r fg = m * m_gravity;
        Vector3r tau = Vector3r::Zero();
        if (m_reg.all_of<C_Cube, C_AngularVelocity3>(e)) {
            const Vector3r he = m_reg.get<C_Cube>(e).halfExtents;
            const Vector3r Idiag = boxInertiaDiag(m, he);
            const Vector3r w = m_reg.get<C_AngularVelocity3>(e).value; // body-frame
            const Vector3r Iw(Idiag.x()*w.x(), Idiag.y()*w.y(), Idiag.z()*w.z());
            tau = - w.cross(Iw);
        }
        VectorXr out(6); out.setZero();
        out[0] = fg[0]; out[1] = fg[1]; out[2] = fg[2];
        out[3] = tau.x(); out[4] = tau.y(); out[5] = tau.z();
        return out;
    }
    // Point mass
    if (m_reg.all_of<C_PointMassTag, C_PhysicsObject, C_Mass>(e)) {
        const real_t m = m_reg.get<C_Mass>(e).m;
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