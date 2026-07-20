#include "world.hpp"
#include <fstream>
#include <iomanip>
#include <optional>
#include <sstream>
#include "../collision/collision_coal.hpp"
#include "assets/assets.hpp"
#include "constraints/constraint_factory.hpp"
#include "constraints/constraints.hpp"
#include "solver/warmstart.hpp"
#include "synchronization/derived_entity_sync.hpp"

namespace cardillo {

World::~World() = default;

World::World(const config::Config& cfg) {
    m_cfg = cfg;
    setGravity(m_cfg.sim_gravity);
}

void World::setGravity(const Vector3r& g) {
    m_gravity = g;
    m_forces_dirty = true;
}

// Assets access ----------------------------------

PhysicsAssets& World::assets() {
    if (!m_assets) m_assets = std::make_shared<PhysicsAssets>();
    return *m_assets;
}
const PhysicsAssets& World::assets() const {
    if (!m_assets) const_cast<World*>(this)->m_assets = std::make_shared<PhysicsAssets>();
    return *m_assets;
}

void World::track(entt::entity e, const std::string& name) {
    if (!m_reg.valid(e)) return;
    if (!m_reg.any_of<C_TrackTag>(e)) {
        m_reg.emplace<C_TrackTag>(e, C_TrackTag{name});
    } else {
        m_reg.get<C_TrackTag>(e).name = name;
    }
}

// ---------- Dynamics getters ----------

VectorXr World::getMassDiag(entt::entity e) const {
    // Prefer rigid body
    if (m_reg.all_of<C_RigidBodyTag, C_PhysicsObject, C_Mass, C_InertiaDiag>(e)) {
        VectorXr diag(6);
        diag.head<3>() = Vector3r::Constant(m_reg.get<C_Mass>(e).m);
        diag.tail<3>() = m_reg.get<C_InertiaDiag>(e).I;
        return diag;
    }

    // Point mass
    if (m_reg.all_of<C_PointMassTag, C_PhysicsObject, C_Mass>(e)) {
        return VectorXr::Constant(3, m_reg.get<C_Mass>(e).m);
    }

    throw std::runtime_error("World::getMassDiag() called on entity without known mass.");
}

MatrixXXr World::getMass(entt::entity e) const {
    return getMassDiag(e).asDiagonal();
}

VectorXr World::getMassInverseDiag(entt::entity e) const {
    return getMassDiag(e).cwiseInverse();
}

VectorXr World::getPosition(entt::entity e) const {
    // Prefer rigid body
    if (m_reg.all_of<C_RigidBodyTag, C_PhysicsObject, C_Position3, C_Orientation>(e)) {
        VectorXr q(7);
        q.head<3>() = m_reg.get<C_Position3>(e).value;
        // store normalized quaternion coefficients
        q.tail<4>() = m_reg.get<C_Orientation>(e).value.coeffs().normalized();
        return q;
    }
    // Point mass
    if (m_reg.all_of<C_PointMassTag, C_PhysicsObject, C_Position3>(e)) {
        return m_reg.get<C_Position3>(e).value;
    }
    throw std::runtime_error("World::getPosition() called on entity without position.");
}

VectorXr World::getVelocity(entt::entity e) const {
    // Prefer rigid body
    if (m_reg.all_of<C_RigidBodyTag, C_PhysicsObject, C_LinearVelocity3, C_AngularVelocity3>(e)) {
        VectorXr v(6);
        v.head<3>() = m_reg.get<C_LinearVelocity3>(e).value; // linear velocity in inertial basis
        v.tail<3>() = m_reg.get<C_AngularVelocity3>(e).value; // angular velocity in body-fixed basis
        return v;
    }
    // Point mass
    if (m_reg.all_of<C_PointMassTag, C_PhysicsObject, C_LinearVelocity3>(e)) {
        return m_reg.get<C_LinearVelocity3>(e).value; // linear velocity in inertial basis
    }

    // throw std::runtime_error("World::getVelocity() called on entity without velocity.");

    // TODO: Where is this required?
    auto state = RigidBody::getState(m_reg, e);
    VectorXr v(6);
    v.head<3>() = state.linearVelocity;
    v.tail<3>() = state.angularVelocity;
    return v;
}

VectorXr World::getForceExternal(entt::entity e) const {
    // Prefer rigid body
    if (m_reg.all_of<C_RigidBodyTag, C_PhysicsObject, C_Mass>(e)) {
        // gravity
        const real_t m = m_reg.get<C_Mass>(e).m;
        Vector3r fg = m * m_gravity;

        // external forces
        if (m_reg.any_of<C_ExternalForce>(e)) {
            fg += m_reg.get<C_ExternalForce>(e).f;
        }

        // external torques
        Vector3r tau = Vector3r::Zero();
        if (m_reg.any_of<C_ExternalTorque>(e)) {
            tau += m_reg.get<C_ExternalTorque>(e).tau;
        }

        VectorXr out(6);
        out.head<3>() = fg;
        out.tail<3>() = tau;
        return out;
    }
    // Point mass
    if (m_reg.all_of<C_PointMassTag, C_PhysicsObject, C_Mass>(e)) {
        // gravity
        const real_t m = m_reg.get<C_Mass>(e).m;
        Vector3r fg = m * m_gravity;

        // external forces
        if (m_reg.any_of<C_ExternalForce>(e)) {
            fg += m_reg.get<C_ExternalForce>(e).f;
        }

        return fg;
    }
    throw std::runtime_error("World::getForceExternal() called on entity without external force.");
}

VectorXr World::getForceGyroscopic(entt::entity e) const {
    // Gyroscopic forces only apply to rigid bodies
    if (m_reg.all_of<C_RigidBodyTag, C_PhysicsObject, C_Mass>(e)) {
        VectorXr out(6);
        out.setZero();
        if (m_reg.all_of<C_InertiaDiag>(e)) {
            const Vector3r& w = m_reg.get<C_AngularVelocity3>(e).value;
            const Vector3r& Idiag = m_reg.get<C_InertiaDiag>(e).I;
            out.tail<3>() = -w.cross(Idiag.cwiseProduct(w));
        }
        return out;
    }

    // Point masses have no gyroscopic contribution; returning a zero torque 
    // is still required for the current interface.
    if (m_reg.all_of<C_PointMassTag, C_PhysicsObject, C_Mass>(e)) {
        return VectorXr::Zero(3);
    }

    throw std::runtime_error("World::getForceGyroscopic() called on entity without gyroscopy.");
    // return VectorXr(0);
}

VectorXr World::getForce(entt::entity e) const {
    // Combined: external + gyroscopic
    return getForceExternal(e) + getForceGyroscopic(e);
}

Vector3r World::getInertiaDiag(entt::entity e) const {
    // If explicitly stored (mesh or custom bodies)
    if (m_reg.all_of<C_InertiaDiag>(e)) {
        return m_reg.get<C_InertiaDiag>(e).I;
    }

    throw std::runtime_error("World::getInertiaDiag() called on entity without inertia.");
}

real_t World::getKineticEnergy(entt::entity e) const {
    real_t KE = (real_t)0;
    
    // Rigid body
    if (m_reg.all_of<C_RigidBodyTag, C_PhysicsObject, C_Mass, C_InertiaDiag, C_LinearVelocity3, C_AngularVelocity3>(e)) {
        const real_t m = m_reg.get<C_Mass>(e).m;
        const Vector3r& Idiag = m_reg.get<C_InertiaDiag>(e).I;
        const Vector3r& v = m_reg.get<C_LinearVelocity3>(e).value;
        const Vector3r& w = m_reg.get<C_AngularVelocity3>(e).value;
        KE += (real_t)0.5 * m * v.squaredNorm(); // linear velocity
        KE += (real_t)0.5 * w.dot(Idiag.cwiseProduct(w)); // angular velocity
        return KE;
    }

    // Point mass
    if (m_reg.all_of<C_PointMassTag, C_PhysicsObject, C_Mass, C_LinearVelocity3>(e)) {
        const real_t m = m_reg.get<C_Mass>(e).m;
        const Vector3r& v = m_reg.get<C_LinearVelocity3>(e).value;
        KE += (real_t)0.5 * m * v.squaredNorm();
        return KE;
    }

    return KE;
}

int World::numBodies() const {
    if (m_num_bodies_dirty) {
        int count = 0;
        // Prefer assembler-assigned body indices when available
        auto viewIndexed = m_reg.view<C_BodyIndex, C_PhysicsObject>();
        for (auto e : viewIndexed) (void)e, ++count;
        if (count == 0) {
            // Fallback: count dynamic physics objects directly before assembler runs
            auto viewDyn = m_reg.view<C_PhysicsObject>();
            for (auto e : viewDyn) (void)e, ++count;
        }
        m_num_bodies_cached = count;
        m_num_bodies_dirty = false;
    }
    return m_num_bodies_cached;
}

void World::applyForce(entt::entity e, const Vector3r& force_world, const Vector3r& torque_body) {
    if (!m_reg.valid(e)) return;
    if (force_world.allFinite() && !force_world.isZero()) {
        if (m_reg.any_of<C_ExternalForce>(e))
            m_reg.get<C_ExternalForce>(e).f += force_world;
        else
            m_reg.emplace<C_ExternalForce>(e, C_ExternalForce{force_world});
    }
    if (torque_body.allFinite() && !torque_body.isZero()) {
        if (m_reg.any_of<C_ExternalTorque>(e))
            m_reg.get<C_ExternalTorque>(e).tau += torque_body;
        else
            m_reg.emplace<C_ExternalTorque>(e, C_ExternalTorque{torque_body});
    }
    m_forces_dirty = true;
}

// Apply a pure moment expressed in world coordinates.
void World::applyInertialTorque(entt::entity e, const Vector3r& torque_world) {
    if (!m_reg.valid(e)) return;
    if (!torque_world.allFinite() || torque_world.isZero()) return;
    if (!m_reg.any_of<C_Orientation>(e)) return;
    const Quaternion4r& q = m_reg.get<C_Orientation>(e).value;
    const Vector3r torque_body = q.conjugate() * torque_world;  // R^T * tau_world
    applyForce(e, Vector3r::Zero(), torque_body);
}

void World::makeStatic(entt::entity e) {
    if (!m_reg.valid(e)) return;
    // Remove dynamic physics components; keep visuals/collidable and shape tags
    if (m_reg.any_of<C_PhysicsObject>(e)) m_reg.remove<C_PhysicsObject>(e);
    if (m_reg.any_of<C_RigidBodyTag>(e)) m_reg.remove<C_RigidBodyTag>(e);
    if (m_reg.any_of<C_PointMassTag>(e)) m_reg.remove<C_PointMassTag>(e);
    if (m_reg.any_of<C_BodyIndex>(e)) m_reg.remove<C_BodyIndex>(e);
    if (m_reg.any_of<C_Mass>(e)) m_reg.remove<C_Mass>(e);
    if (m_reg.any_of<C_InertiaDiag>(e)) m_reg.remove<C_InertiaDiag>(e);
    if (m_reg.any_of<C_ExternalForce>(e)) m_reg.remove<C_ExternalForce>(e);
    if (m_reg.any_of<C_ExternalTorque>(e)) m_reg.remove<C_ExternalTorque>(e);
    markStructureDirty();
}

// ---------- Minimal setters ----------

void World::setPosition(entt::entity e, const Vector3r& p) {
    if (!m_reg.valid(e)) return;
    if (m_reg.any_of<C_Position3>(e))
        m_reg.get<C_Position3>(e).value = p;
    else
        m_reg.emplace<C_Position3>(e, C_Position3{p});
    markStateDirty();
}

void World::setOrientation(entt::entity e, const Quaternion4r& q_in) {
    if (!m_reg.valid(e)) return;
    if (m_reg.any_of<C_Orientation>(e)) {
        const Quaternion4r& q_ref = m_reg.get<C_Orientation>(e).value;
        m_reg.get<C_Orientation>(e).setValue(MathHelper::alignQuaternionTo(q_in, q_ref));
    } else {
        Quaternion4r q = q_in;
        q.normalize();
        m_reg.emplace<C_Orientation>(e, C_Orientation::fromQuaternion(q));
    }
    markStateDirty();
}

void World::setLinearVelocity(entt::entity e, const Vector3r& v) {
    if (!m_reg.valid(e)) return;
    if (m_reg.any_of<C_LinearVelocity3>(e))
        m_reg.get<C_LinearVelocity3>(e).value = v;
    else
        m_reg.emplace<C_LinearVelocity3>(e, C_LinearVelocity3{v});
    if (m_reg.any_of<C_LinearAcceleration3>(e))
        m_reg.get<C_LinearAcceleration3>(e).value = Vector3r::Zero();
    else
        m_reg.emplace<C_LinearAcceleration3>(e, C_LinearAcceleration3{Vector3r::Zero()});
    markStateDirty();
}

void World::setAngularVelocity(entt::entity e, const Vector3r& w) {
    if (!m_reg.valid(e)) return;
    if (m_reg.any_of<C_AngularVelocity3>(e))
        m_reg.get<C_AngularVelocity3>(e).value = w;
    else
        m_reg.emplace<C_AngularVelocity3>(e, C_AngularVelocity3{w});
    if (m_reg.any_of<C_AngularAcceleration3>(e))
        m_reg.get<C_AngularAcceleration3>(e).value = Vector3r::Zero();
    else
        m_reg.emplace<C_AngularAcceleration3>(e, C_AngularAcceleration3{Vector3r::Zero()});
    markStateDirty();
}

void World::setVelocityByForce(entt::entity e, const Vector3r& v, const Vector3r& w) {
    if (!m_reg.valid(e)) return;

    Vector3r currentV = Vector3r::Zero();
    Vector3r currentW = Vector3r::Zero();
    if (m_reg.any_of<C_LinearVelocity3>(e)) currentV = m_reg.get<C_LinearVelocity3>(e).value;
    if (m_reg.any_of<C_AngularVelocity3>(e)) currentW = m_reg.get<C_AngularVelocity3>(e).value;

    Vector3r deltaV = currentV - v;
    Vector3r deltaW = currentW - w;
    const real_t dt = m_cfg.sim_dt;

    if (!deltaV.isZero() || !deltaW.isZero()) {
        real_t mass = (m_reg.any_of<C_Mass>(e) ? m_reg.get<C_Mass>(e).m : (real_t)0);
        Vector3r inertia = (m_reg.any_of<C_InertiaDiag>(e) ? m_reg.get<C_InertiaDiag>(e).I : Vector3r::Zero());
        applyForce(e, -deltaV * mass / dt, -deltaW.cwiseProduct(inertia) / dt);
    }
    m_forces_dirty = true;
}

void World::setTrajectory(entt::entity e, std::optional<std::function<TrajectoryPose(real_t)>> positionFunc, std::optional<std::function<TrajectoryTwist(real_t)>> velocityFunc) {
    if (!m_reg.valid(e)) return;

    if (!RigidBody::isStatic(m_reg, e)) {
        std::cout << "Warning: setTrajectory() called on non-static entity; making static.\n";
        makeStatic(e);
    }

    if (!positionFunc.has_value() && !velocityFunc.has_value()) {
        removeTrajectory(e);
        return;
    }

    C_StaticTrajectory traj;
    traj.positionFunc = std::move(positionFunc);
    traj.velocityFunc = std::move(velocityFunc);
    traj.elapsed = (real_t)0;
    traj.initialized = false;
    traj.previousPosition = std::nullopt;

    m_reg.emplace_or_replace<C_StaticTrajectory>(e, std::move(traj));
    markStateDirty();
}

void World::removeTrajectory(entt::entity e) {
    if (!m_reg.valid(e)) return;
    if (m_reg.any_of<C_StaticTrajectory>(e)) {
        m_reg.remove<C_StaticTrajectory>(e);
        markStateDirty();
    }
}

// ---------- Mesh asset access ----------

const ::cardillo::MeshAsset& World::getMeshAsset(entt::entity e) const {
    if (!m_reg.any_of<C_Mesh>(e)) throw std::runtime_error("getMeshAsset(entity): entity has no C_Mesh");
    const auto& cm = m_reg.get<C_Mesh>(e);
    const bool isDynamic = m_reg.any_of<C_PhysicsObject>(e);
    return assets().getMesh(cm.path, cm.scale, /*normalized*/ isDynamic);
}
}  // namespace cardillo