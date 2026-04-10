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

// namespace {

// Vector3r worldPointFromLocal(const entt::registry& reg, entt::entity e, const Vector3r& r_local)
// {
//     Vector3r p = r_local;
//     if (reg.valid(e) && reg.all_of<C_Position3>(e)) {
//         p = reg.get<C_Position3>(e).value;
//         if (reg.all_of<C_Orientation>(e)) {
//             p += reg.get<C_Orientation>(e).value.toRotationMatrix() * r_local;
//         } else {
//             p += r_local;
//         }
//     }
//     return p;
// }

// } // namespace

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

// const config::Config& World::config() const { return m_cfg; }
// config::Config& World::config() { return m_cfg; }

void World::track(entt::entity e, const std::string& name) {
    if (!m_reg.valid(e)) return;
    if (!m_reg.any_of<C_TrackTag>(e)) {
        m_reg.emplace<C_TrackTag>(e, C_TrackTag{name});
    } else {
        m_reg.get<C_TrackTag>(e).name = name;
    }
}

// ---------- Dynamics getters ----------

MatrixXXr World::getMass(entt::entity e) const {
    // Rigid body
    if (m_reg.all_of<C_RigidBodyTag, C_PhysicsObject, C_Mass, C_InertiaDiag>(e)) {
        const real_t m = m_reg.get<C_Mass>(e).m;
        const Vector3r Idiag = m_reg.get<C_InertiaDiag>(e).I;
        MatrixXXr M = MatrixXXr::Zero(6, 6);

        M(0, 0) = m;
        M(1, 1) = m;
        M(2, 2) = m;
        M(3, 3) = Idiag.x();
        M(4, 4) = Idiag.y();
        M(5, 5) = Idiag.z();
        return M;
    }

    // Point mass
    if (m_reg.all_of<C_PointMassTag, C_PhysicsObject, C_Mass>(e)) {
        const real_t m = m_reg.get<C_Mass>(e).m;
        MatrixXXr M = MatrixXXr::Zero(3, 3);
        M.diagonal().setConstant(m);
        return M;
    }
    std::cout << "Warning: getMass() called on entity without known mass.\n";
    return MatrixXXr(0, 0);
}

VectorXr World::getMassInverseDiag(entt::entity e) const {
    // Compute from getMass() to keep a single source of truth for layout
    const MatrixXXr M = getMass(e);
    const index_t n = (index_t)M.rows();
    if (n == 0) return VectorXr(0);
    VectorXr d = M.diagonal();
    VectorXr inv = VectorXr::Zero(n);
    for (index_t i = 0; i < n; ++i) {
        const real_t v = d[i];
        if (v > (real_t)0) inv[i] = (real_t)1 / v;
    }
    return inv;
}

VectorXr World::getPosition(entt::entity e) const {
    // Prefer rigid body
    if (m_reg.all_of<C_RigidBodyTag, C_PhysicsObject, C_Position3, C_Orientation>(e)) {
        const auto& x = m_reg.get<C_Position3>(e).value;
        Quaternion4r qn = m_reg.get<C_Orientation>(e).value;
        qn.normalize();
        VectorXr q(7);
        q.head<3>() = x;
        q.tail<4>() = qn.coeffs();
        return q;
    }
    // Point mass
    if (m_reg.all_of<C_PointMassTag, C_PhysicsObject, C_Position3>(e)) {
        const auto& x = m_reg.get<C_Position3>(e).value;
        VectorXr q(3);
        q << x[0], x[1], x[2];
        return q;
    }
    return VectorXr(0);
}

VectorXr World::getVelocity(entt::entity e) const {
    // Prefer rigid body
    if (m_reg.all_of<C_RigidBodyTag, C_PhysicsObject, C_LinearVelocity3, C_AngularVelocity3>(e)) {
        const auto& v = m_reg.get<C_LinearVelocity3>(e).value;
        const auto& w = m_reg.get<C_AngularVelocity3>(e).value;  // body-frame angular velocity
        VectorXr out(6);
        out << v[0], v[1], v[2], w[0], w[1], w[2];
        return out;
    }
    // Point mass
    if (m_reg.all_of<C_PointMassTag, C_PhysicsObject, C_LinearVelocity3>(e)) {
        const auto& v = m_reg.get<C_LinearVelocity3>(e).value;
        VectorXr out(3);
        out << v[0], v[1], v[2];
        return out;
    }
    return VectorXr(0);
}

VectorXr World::getForceExternal(entt::entity e) const {
    // Prefer rigid body
    if (m_reg.all_of<C_RigidBodyTag, C_PhysicsObject, C_Mass>(e)) {
        const real_t m = m_reg.get<C_Mass>(e).m;
        Vector3r fg = m * m_gravity;
        Vector3r tau = Vector3r::Zero();
        // Add any external one-shot forces/torques
        if (m_reg.any_of<C_ExternalForce>(e)) {
            fg += m_reg.get<C_ExternalForce>(e).f;
        }
        if (m_reg.any_of<C_ExternalTorque>(e)) {
            tau += m_reg.get<C_ExternalTorque>(e).tau;
        }
        VectorXr out(6);
        out.setZero();
        out[0] = fg[0];
        out[1] = fg[1];
        out[2] = fg[2];
        out[3] = tau[0];
        out[4] = tau[1];
        out[5] = tau[2];
        return out;
    }
    // Point mass
    if (m_reg.all_of<C_PointMassTag, C_PhysicsObject, C_Mass>(e)) {
        const real_t m = m_reg.get<C_Mass>(e).m;
        Vector3r fg = m * m_gravity;
        if (m_reg.any_of<C_ExternalForce>(e)) {
            fg += m_reg.get<C_ExternalForce>(e).f;
        }
        VectorXr out(3);
        out << fg[0], fg[1], fg[2];
        return out;
    }
    return VectorXr(0);
}

VectorXr World::getForceGyroscopic(entt::entity e) const {
    // Gyroscopic forces only apply to rigid bodies
    if (m_reg.all_of<C_RigidBodyTag, C_PhysicsObject, C_Mass>(e)) {
        Vector3r tau = Vector3r::Zero();
        if (m_reg.all_of<C_InertiaDiag>(e)) {
            const Vector3r w = m_reg.get<C_AngularVelocity3>(e).value;  // body-frame
            Vector3r Idiag = m_reg.get<C_InertiaDiag>(e).I;
            const Vector3r Iw = Idiag.cwiseProduct(w);
            tau = -w.cross(Iw);
        }
        VectorXr out(6);
        out.setZero();
        out.tail<3>() = tau;
        return out;
    }
    // Point masses have no gyroscopic contribution
    if (m_reg.all_of<C_PointMassTag, C_PhysicsObject, C_Mass>(e)) {
        return VectorXr::Zero(3);
    }
    return VectorXr(0);
}

VectorXr World::getForce(entt::entity e) const {
    // Combined: external + gyroscopic
    VectorXr f_ext = getForceExternal(e);
    VectorXr f_gyro = getForceGyroscopic(e);
    return f_ext + f_gyro;
}

Vector3r World::getInertiaDiag(entt::entity e) const {
    // If explicitly stored (mesh or custom bodies)
    if (m_reg.all_of<C_InertiaDiag>(e)) {
        return m_reg.get<C_InertiaDiag>(e).I;
    }

    // Default: zero
    return Vector3r::Zero();
}

real_t World::getKineticEnergy(entt::entity e) const {
    real_t KE = (real_t)0;
    // Rigid body
    if (m_reg.all_of<C_RigidBodyTag, C_PhysicsObject, C_Mass, C_InertiaDiag, C_LinearVelocity3, C_AngularVelocity3>(e)) {
        const real_t m = m_reg.get<C_Mass>(e).m;
        const Vector3r Idiag = m_reg.get<C_InertiaDiag>(e).I;
        const Vector3r v = m_reg.get<C_LinearVelocity3>(e).value;
        const Vector3r w = m_reg.get<C_AngularVelocity3>(e).value;  // body-frame
        KE += (real_t)0.5 * m * v.squaredNorm();
        Vector3r Iw = Idiag.cwiseProduct(w);
        KE += (real_t)0.5 * w.dot(Iw);
    }
    // Point mass
    else if (m_reg.all_of<C_PointMassTag, C_PhysicsObject, C_Mass, C_LinearVelocity3>(e)) {
        const real_t m = m_reg.get<C_Mass>(e).m;
        const Vector3r v = m_reg.get<C_LinearVelocity3>(e).value;
        KE += (real_t)0.5 * m * v.squaredNorm();
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
    const Quaternion4r q = m_reg.get<C_Orientation>(e).value;
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
    if (m_reg.any_of<C_LinearVelocity3>(e)) m_reg.remove<C_LinearVelocity3>(e);
    if (m_reg.any_of<C_LinearAcceleration3>(e)) m_reg.remove<C_LinearAcceleration3>(e);
    if (m_reg.any_of<C_AngularVelocity3>(e)) m_reg.remove<C_AngularVelocity3>(e);
    if (m_reg.any_of<C_AngularAcceleration3>(e)) m_reg.remove<C_AngularAcceleration3>(e);
    if (m_reg.any_of<C_ExternalForce>(e)) m_reg.remove<C_ExternalForce>(e);
    if (m_reg.any_of<C_ExternalTorque>(e)) m_reg.remove<C_ExternalTorque>(e);
    // Ensure pose remains; visuals/collidable components are left intact
    if (!m_reg.any_of<C_Position3>(e)) m_reg.emplace<C_Position3>(e, C_Position3{Vector3r::Zero()});
    if (!m_reg.any_of<C_Orientation>(e)) m_reg.emplace<C_Orientation>(e, C_Orientation{Quaternion4r::Identity()});
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
        const Quaternion4r q_ref = m_reg.get<C_Orientation>(e).value;
        m_reg.get<C_Orientation>(e).value = MathHelper::alignQuaternionTo(q_in, q_ref);
    } else {
        Quaternion4r q = q_in;
        q.normalize();
        m_reg.emplace<C_Orientation>(e, C_Orientation{q});
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

// ---------- Mesh / HeightField asset access ----------

const ::cardillo::MeshAsset& World::getMeshAsset(entt::entity e) const {
    if (!m_reg.any_of<C_Mesh>(e)) throw std::runtime_error("getMeshAsset(entity): entity has no C_Mesh");
    const auto& cm = m_reg.get<C_Mesh>(e);
    const bool isDynamic = m_reg.any_of<C_PhysicsObject>(e);
    return assets().getMesh(cm.path, cm.scale, /*normalized*/ isDynamic);
}
const ::cardillo::HeightFieldAsset& World::getHeightFieldAsset(entt::entity e) const {
    if (!m_reg.any_of<C_HeightField>(e)) throw std::runtime_error("getHeightFieldAsset(entity): entity has no C_HeightField");
    const auto& ch = m_reg.get<C_HeightField>(e);
    return assets().getHeightField(ch.path, ch.x_dim, ch.y_dim, ch.z_scale, ch.min_height);
}
}  // namespace cardillo