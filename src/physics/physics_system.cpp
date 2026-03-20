#include <optional>
#include "physics_system.hpp"
#include "assets.hpp"
#include "constraints.hpp"
#include "../collision/collision_coal.hpp"
#include "../solver/warmstart.hpp"
#include "../io/csv_writer.hpp"
#include <sstream>
#include <iomanip>
#include <fstream>

namespace cardillo {

// Construction and configuration
PhysicsSystem::PhysicsSystem() 
{
    m_gravity = Vector3r(0, 0, -9.81);
    // Default warmstart provider: simple global-index cache
    m_warmstart_provider = std::make_unique<cardillo::solver::WarmstartCache>();
}

PhysicsSystem::~PhysicsSystem() = default;

PhysicsSystem::PhysicsSystem(const config::Config& cfg) : PhysicsSystem() {
    setConfig(cfg);
}

void PhysicsSystem::setGravity(const Vector3r& g) { m_gravity = g; m_forces_dirty = true; }

// Assets access ----------------------------------
PhysicsAssets& PhysicsSystem::assets() {
    if (!m_assets) m_assets = std::make_shared<PhysicsAssets>();
    return *m_assets;
}
const PhysicsAssets& PhysicsSystem::assets() const {
    if (!m_assets) const_cast<PhysicsSystem*>(this)->m_assets = std::make_shared<PhysicsAssets>();
    return *m_assets;
}

void PhysicsSystem::track(entt::entity e, const std::string& name)
{
    if (!m_reg.valid(e)) return;
    if (!m_reg.any_of<C_TrackTag>(e)) {
        m_reg.emplace<C_TrackTag>(e, C_TrackTag{name});
    } else {
        m_reg.get<C_TrackTag>(e).name = name;
    }
}


void PhysicsSystem::writeTrackedStateToCsv(real_t t)
{
    static cardillo::io::CsvWriter writer;
    static bool initialized = false;
    if (!initialized) {
        std::string path = m_cfg.output_folder + "/" + m_cfg.output_filename_prefix + "_tracked.csv";
        std::vector<std::string> header = {
            "t",
            "name",
            "px","py","pz",
            "vx","vy","vz",
            "wx","wy","wz",
            "euler_x","euler_y","euler_z"
        };
        writer.open(path, header);
        initialized = true;
    }
    if (!writer.isOpen()) return;

    auto view = m_reg.view<C_TrackTag, C_Position3, C_LinearVelocity3, C_AngularVelocity3, C_Orientation>();
    for (auto e : view) {
        const auto& tag   = view.get<C_TrackTag>(e);
        const auto& pos   = view.get<C_Position3>(e).value;
        const auto& v     = view.get<C_LinearVelocity3>(e).value;
        const auto& w     = view.get<C_AngularVelocity3>(e).value;
        const auto& euler = view.get<C_Orientation>(e).value.toRotationMatrix().eulerAngles(0, 1, 2);
        writer.writeRow(
            t,
            tag.name,
            pos.x(), pos.y(), pos.z(),
            v.x(), v.y(), v.z(),
            w.x(), w.y(), w.z(),
            euler.x(), euler.y(), euler.z()
        );
    }
}


// Add a general constraint pattern (ownership transferred)
size_t PhysicsSystem::addConstraint(std::unique_ptr<cardillo::physics::ConstraintPattern> pattern)
{
    if (pattern) {
    m_constraints_new.emplace_back(std::move(pattern));
    markStructureDirty();
        return m_constraints_new.size() - 1;
    }
    return static_cast<size_t>(-1);
}

// ---------- Dynamics getters ----------

MatrixXXr PhysicsSystem::getMass(entt::entity e) const {
    // Rigid body
    if (m_reg.all_of<C_RigidBodyTag, C_PhysicsObject, C_Mass, C_InertiaDiag>(e)) {
        const real_t m = m_reg.get<C_Mass>(e).m;
        const Vector3r Idiag = m_reg.get<C_InertiaDiag>(e).I;
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
    std::cout << "Warning: getMass() called on entity without known mass.\n";
    return MatrixXXr(0,0);
}

VectorXr PhysicsSystem::getMassInverseDiag(entt::entity e) const {
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

VectorXr PhysicsSystem::getPosition(entt::entity e) const {
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

VectorXr PhysicsSystem::getForceExternal(entt::entity e) const {
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
        VectorXr out(6); out.setZero();
        out[0] = fg[0]; out[1] = fg[1]; out[2] = fg[2];
        out[3] = tau[0]; out[4] = tau[1]; out[5] = tau[2];
        return out;
    }
    // Point mass
    if (m_reg.all_of<C_PointMassTag, C_PhysicsObject, C_Mass>(e)) {
        const real_t m = m_reg.get<C_Mass>(e).m;
        Vector3r fg = m * m_gravity;
        if (m_reg.any_of<C_ExternalForce>(e)) {
            fg += m_reg.get<C_ExternalForce>(e).f;
        }
        VectorXr out(3); out << fg[0], fg[1], fg[2];
        return out;
    }
    return VectorXr(0);
}

VectorXr PhysicsSystem::getForceGyroscopic(entt::entity e) const {
    // Gyroscopic forces only apply to rigid bodies
    if (m_reg.all_of<C_RigidBodyTag, C_PhysicsObject, C_Mass>(e)) {
        Vector3r tau = Vector3r::Zero();
        if (m_reg.all_of<C_InertiaDiag>(e)) {
            const Vector3r w = m_reg.get<C_AngularVelocity3>(e).value; // body-frame
            Vector3r Idiag = m_reg.get<C_InertiaDiag>(e).I;
            const Vector3r Iw = Idiag.cwiseProduct(w);
            tau = -w.cross(Iw);
        }
        VectorXr out(6); out.setZero();
        out.tail<3>() = tau;
        return out;
    }
    // Point masses have no gyroscopic contribution
    if (m_reg.all_of<C_PointMassTag, C_PhysicsObject, C_Mass>(e)) {
        return VectorXr::Zero(3);
    }
    return VectorXr(0);
}

VectorXr PhysicsSystem::getForce(entt::entity e) const {
    // Combined: external + gyroscopic
    VectorXr f_ext = getForceExternal(e);
    VectorXr f_gyro = getForceGyroscopic(e);
    return f_ext + f_gyro;
}

Vector3r PhysicsSystem::getInertiaDiag(entt::entity e) const {
    // If explicitly stored (mesh or custom bodies)
    if (m_reg.all_of<C_InertiaDiag>(e)) {
        return m_reg.get<C_InertiaDiag>(e).I;
    }
    // Cube fallback
    if (m_reg.all_of<C_Cube, C_Mass>(e)) {
        const real_t m = m_reg.get<C_Mass>(e).m;
        const Vector3r he = m_reg.get<C_Cube>(e).halfExtents;
        Vector3r I;
        I.x() = (real_t)1.0 / 3.0 * m * (he.y() * he.y() + he.z() * he.z());
        I.y() = (real_t)1.0 / 3.0 * m * (he.x() * he.x() + he.z() * he.z());
        I.z() = (real_t)1.0 / 3.0 * m * (he.x() * he.x() + he.y() * he.y());
        return I;
    }
    // Default: zero
    return Vector3r::Zero();
}

real_t PhysicsSystem::getKineticEnergy(entt::entity e) const {
    real_t KE = (real_t)0;
    // Rigid body
    if (m_reg.all_of<C_RigidBodyTag, C_PhysicsObject, C_Mass, C_InertiaDiag, C_LinearVelocity3, C_AngularVelocity3>(e)) {
        const real_t m = m_reg.get<C_Mass>(e).m;
        const Vector3r Idiag = m_reg.get<C_InertiaDiag>(e).I;
        const Vector3r v = m_reg.get<C_LinearVelocity3>(e).value;
        const Vector3r w = m_reg.get<C_AngularVelocity3>(e).value; // body-frame
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

int PhysicsSystem::numBodies() const {
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

void PhysicsSystem::applyForce(entt::entity e, const Vector3r& force_world, const Vector3r& torque_body) {
    if (!m_reg.valid(e)) return;
    if (force_world.allFinite() && !force_world.isZero()) {
        if (m_reg.any_of<C_ExternalForce>(e)) m_reg.get<C_ExternalForce>(e).f += force_world;
        else m_reg.emplace<C_ExternalForce>(e, C_ExternalForce{force_world});
    }
    if (torque_body.allFinite() && !torque_body.isZero()) {
        if (m_reg.any_of<C_ExternalTorque>(e)) m_reg.get<C_ExternalTorque>(e).tau += torque_body;
        else m_reg.emplace<C_ExternalTorque>(e, C_ExternalTorque{torque_body});
    }
    m_forces_dirty = true;
}

// Apply a pure moment expressed in world coordinates.
void PhysicsSystem::applyInertialTorque(entt::entity e, const Vector3r& torque_world) {
    if (!m_reg.valid(e)) return;
    if (!torque_world.allFinite() || torque_world.isZero()) return;
    if (!m_reg.any_of<C_Orientation>(e)) return;
    const Quaternion4r q = m_reg.get<C_Orientation>(e).value;
    const Vector3r torque_body = q.conjugate() * torque_world; // R^T * tau_world
    applyForce(e, Vector3r::Zero(), torque_body);
}

void PhysicsSystem::makeStatic(entt::entity e) {
    if (!m_reg.valid(e)) return;
    // Remove dynamic physics components; keep visuals/collidable and shape tags
    if (m_reg.any_of<C_PhysicsObject>(e)) m_reg.remove<C_PhysicsObject>(e);
    if (m_reg.any_of<C_RigidBodyTag>(e)) m_reg.remove<C_RigidBodyTag>(e);
    if (m_reg.any_of<C_PointMassTag>(e)) m_reg.remove<C_PointMassTag>(e);
    if (m_reg.any_of<C_BodyIndex>(e)) m_reg.remove<C_BodyIndex>(e);
    if (m_reg.any_of<C_Mass>(e)) m_reg.remove<C_Mass>(e);
    if (m_reg.any_of<C_InertiaDiag>(e)) m_reg.remove<C_InertiaDiag>(e);
    if (m_reg.any_of<C_LinearVelocity3>(e)) m_reg.remove<C_LinearVelocity3>(e);
    if (m_reg.any_of<C_AngularVelocity3>(e)) m_reg.remove<C_AngularVelocity3>(e);
    if (m_reg.any_of<C_ExternalForce>(e)) m_reg.remove<C_ExternalForce>(e);
    if (m_reg.any_of<C_ExternalTorque>(e)) m_reg.remove<C_ExternalTorque>(e);
    // Ensure pose remains; visuals/collidable components are left intact
    if (!m_reg.any_of<C_Position3>(e)) m_reg.emplace<C_Position3>(e, C_Position3{Vector3r::Zero()});
    if (!m_reg.any_of<C_Orientation>(e)) m_reg.emplace<C_Orientation>(e, C_Orientation{Quaternion4r::Identity()});
    markStructureDirty();
}

// ---------- Minimal setters ----------

void PhysicsSystem::setPosition(entt::entity e, const Vector3r& p) {
    if (!m_reg.valid(e)) return;
    if (m_reg.any_of<C_Position3>(e)) m_reg.get<C_Position3>(e).value = p; else m_reg.emplace<C_Position3>(e, C_Position3{p});
    markStateDirty();
}

void PhysicsSystem::setOrientation(entt::entity e, const Quaternion4r& q_in) {
    if (!m_reg.valid(e)) return;
    if (m_reg.any_of<C_Orientation>(e)) {
        const Quaternion4r q_ref = m_reg.get<C_Orientation>(e).value;
        m_reg.get<C_Orientation>(e).value = PhysicsSystem::alignQuaternionTo(q_in, q_ref);
    } else {
        Quaternion4r q = q_in; q.normalize();
        m_reg.emplace<C_Orientation>(e, C_Orientation{q});
    }
    markStateDirty();
}

void PhysicsSystem::setLinearVelocity(entt::entity e, const Vector3r& v) {
    if (!m_reg.valid(e)) return;
    if (m_reg.any_of<C_LinearVelocity3>(e)) m_reg.get<C_LinearVelocity3>(e).value = v; else m_reg.emplace<C_LinearVelocity3>(e, C_LinearVelocity3{v});
    markStateDirty();
}

void PhysicsSystem::setAngularVelocity(entt::entity e, const Vector3r& w) {
    if (!m_reg.valid(e)) return;
    if (m_reg.any_of<C_AngularVelocity3>(e)) m_reg.get<C_AngularVelocity3>(e).value = w; else m_reg.emplace<C_AngularVelocity3>(e, C_AngularVelocity3{w});
    markStateDirty();
}

void PhysicsSystem::setVelocityByForce(entt::entity e, const Vector3r& v, const Vector3r& w)
{
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

// ---------- Subsystems ----------

cardillo::collision::CollisionCoal& PhysicsSystem::collisionManager() {
    if (!m_collision_mgr) {
        m_collision_mgr = std::make_unique<cardillo::collision::CollisionCoal>();
        m_collision_mgr->registerSystem(this);
    }
    return *m_collision_mgr;
}

const cardillo::collision::CollisionCoal& PhysicsSystem::collisionManager() const {
    // const-correct lazy init: cast away const for creation then return const ref
    if (!m_collision_mgr) {
        auto* self = const_cast<PhysicsSystem*>(this);
        self->m_collision_mgr = std::make_unique<cardillo::collision::CollisionCoal>();
        self->m_collision_mgr->registerSystem(self);
    }
    return *m_collision_mgr;
}

cardillo::misc::TimingManager& PhysicsSystem::timings() {
    if (!m_timings) m_timings = std::make_unique<cardillo::misc::TimingManager>();
    return *m_timings;
}
const cardillo::misc::TimingManager& PhysicsSystem::timings() const {
    auto* self = const_cast<PhysicsSystem*>(this);
    if (!self->m_timings) self->m_timings = std::make_unique<cardillo::misc::TimingManager>();
    return *self->m_timings;
}

// ---------- Collision pair control ----------

void PhysicsSystem::disableCollisionBetween(entt::entity a, entt::entity b) {
    collisionManager().disablePair(a, b);
}

// ---------- Mesh / HeightField asset access ----------

const ::cardillo::MeshAsset& PhysicsSystem::getMeshAsset(entt::entity e) const {
    if (!m_reg.any_of<C_Mesh>(e)) throw std::runtime_error("getMeshAsset(entity): entity has no C_Mesh");
    const auto& cm = m_reg.get<C_Mesh>(e);
    const bool isDynamic = m_reg.any_of<C_PhysicsObject>(e);
    return assets().getMesh(cm.path, cm.scale, /*normalized*/ isDynamic);
}
const ::cardillo::HeightFieldAsset& PhysicsSystem::getHeightFieldAsset(entt::entity e) const {
    if (!m_reg.any_of<C_HeightField>(e)) throw std::runtime_error("getHeightFieldAsset(entity): entity has no C_HeightField");
    const auto& ch = m_reg.get<C_HeightField>(e);
    return assets().getHeightField(ch.path, ch.x_dim, ch.y_dim, ch.z_scale, ch.min_height);
}

void PhysicsSystem::explicitPositionUpdate(real_t h) {
    auto _sc = timings().scope(cardillo::misc::TimingManager::TimerId::Integration);

    // position update
    auto position_view = m_reg.view<C_Position3, const C_LinearVelocity3>();
    for (auto [e, pos, vel] : position_view.each()) {
        pos.value += h * vel.value;
    }

    // auto director_orientation_view = m_reg.view<C_DirectorTriad, const C_DirectorTriadVelocity>();
    // for (auto [e, pos, vel] : position_view.each()) {
    //     pos.value += h * vel.value;
    // }

    // Matrix33r A_IK = VectorXr::Zero(9).reshaped<3,3>();
    // Vector3r d1 = A_IK.col<0>();
    // Vector3r d2 = A_IK.col<1>();
    // Vector3r d3 = A_IK.col<2>();
    // VectorXr q = VectorXr::Zero(9);
    // Vector3r d1 = q.head<3>();
    // Vector3r d2 = q.segment<3>(3);
    // Vector3r d3 = q.tail<3>();

    // orientations
    auto orientation_view = m_reg.view<C_Orientation, const C_AngularVelocity3>();
    for (auto [e, orientation, angularVel] : orientation_view.each()) {
        const Vector3r& omega = angularVel.value;
        const Quaternion4r q_prev = orientation.value;

        Vector4r& P = orientation.value.coeffs();
        real_t w = P(3);
        P(3) -= h * 0.5 * P.head<3>().dot(omega);
        P.head<3>() += h * 0.5 * (w * omega + P.head<3>().cross(omega));

        // re-normalize to avoid quaternion drift and keep hemisphere consistent with previous
        orientation.value = PhysicsSystem::alignQuaternionTo(orientation.value, q_prev);
    }
    updateEntities();
}

void PhysicsSystem::linearImplicitPositionUpdate(real_t h) {
    auto _sc = timings().scope(cardillo::misc::TimingManager::TimerId::Integration);

    // position update
    auto position_view = m_reg.view<C_Position3, const C_LinearVelocity3>();
    for (auto [e, pos, vel] : position_view.each()) {
        pos.value += h * vel.value;
    }

    // orientations
    auto orientation_view = m_reg.view<C_Orientation, const C_AngularVelocity3>();
    for (auto [e, orientation, angularVel] : orientation_view.each()) {
        const Vector3r& omega = angularVel.value;
        const Quaternion4r q_prev = orientation.value;

        // skew matrix q_dot = 0.5 * D * P
        Matrix44r D = Matrix44r::Zero();
        D.block<3, 3>(0, 0) = skew_from_vector(-omega);
        D.block<3, 1>(0, 3) = omega;
        D.block<1, 3>(3, 0) = -omega.transpose();

        // iteration matrix (I - h / 2 * D) P_{n+1} = P_{n+1/2}
        const Matrix44r A = Matrix44r::Identity() - 0.5 * h * D;

        // quaternion coefficients
        Vector4r& P = orientation.value.coeffs();

        // linear implicit update
        P = A.inverse() * P;

        // re-normalize to avoid quaternion drift and keep hemisphere consistent with previous
        orientation.value = PhysicsSystem::alignQuaternionTo(orientation.value, q_prev);
    }
    // Update beam element lengths after full pose (position + orientation) update
    updateEntities();
}

void PhysicsSystem::updateBeamElementEntity(entt::entity e) {
    if (!m_reg.valid(e) || !m_reg.any_of<C_BeamElement, C_Position3>(e)) return;
    auto& be = m_reg.get<C_BeamElement>(e);
    const auto& pos = m_reg.get<C_Position3>(e).value;
    real_t newLen = be.l0;

    auto getDesiredLengthBetween = [&](entt::entity a, entt::entity b)->real_t {
        if (!(m_reg.valid(a) && m_reg.valid(b))) {
            std::cout << "Beam length: invalid entities a=" << (int)a << " b=" << (int)b << "\eaen";
        }
        if (!m_reg.any_of<C_Position3>(a) || !m_reg.any_of<C_Position3>(b)) {
            std::cout << "Beam length: missing positions a=" << m_reg.any_of<C_Position3>(a) << " b=" << m_reg.any_of<C_Position3>(b) << " for entities a=" << (int)a << " b=" << (int)b << "\n";
        }
        if (!m_reg.any_of<C_Orientation>(a) || !m_reg.any_of<C_Orientation>(b)) {
            std::cout << "Beam length: missing orientations a=" << m_reg.any_of<C_Orientation>(a) << " b=" << m_reg.any_of<C_Orientation>(b) << " for entities a=" << (int)a << " b=" << (int)b << "\n";
        }
        if (m_reg.valid(a) && m_reg.valid(b) &&
            m_reg.any_of<C_Position3>(a) && m_reg.any_of<C_Position3>(b) &&
            m_reg.any_of<C_Orientation>(a) && m_reg.any_of<C_Orientation>(b)) {
            const auto& pa = m_reg.get<C_Position3>(a).value;
            const auto& pb = m_reg.get<C_Position3>(b).value;
            auto r_AB = pb - pa;
            auto R_A = m_reg.get<C_Orientation>(a).value.toRotationMatrix();
            auto R_B = m_reg.get<C_Orientation>(b).value.toRotationMatrix();
            
            index_t x_col_A = (m_reg.any_of<C_Capsule>(a) || m_reg.any_of<C_Cylinder>(a)) ? 2 : 0; 
            index_t x_col_B = (m_reg.any_of<C_Capsule>(b) || m_reg.any_of<C_Cylinder>(b)) ? 2 : 0;

            auto e_Ax = R_A.col(x_col_A);
            auto e_Bx = R_B.col(x_col_B);
            real_t base = r_AB.norm();
            if (base > (real_t)0) {
                auto r_hat = r_AB.normalized();
                real_t d1 = std::clamp(e_Ax.dot(r_hat), (real_t)-1, (real_t)1);
                real_t d2 = std::clamp(e_Bx.dot(r_hat), (real_t)-1, (real_t)1);
                // Circular mean of angles: atan2(sum sin, sum cos)
                // Since alpha in [0,pi], sin is non-negative; compute from cos using sin=√(1-cos^2)
                real_t s1 = std::sqrt(std::max((real_t)0, (real_t)1 - d1*d1));
                real_t s2 = std::sqrt(std::max((real_t)0, (real_t)1 - d2*d2));
                real_t base_angle = std::atan2(s1 + s2, d1 + d2);
                real_t c = std::cos(base_angle);
                if (std::abs(c) > (real_t)1e-12) return base / c;
            }
        }
        return (real_t) 0;
    };

    index_t contributions = 0;
    real_t totalLen = (real_t)0;

    if (be.prev.has_value()) {
        totalLen += getDesiredLengthBetween(be.prev.value(), e);
        contributions++;
    }

    if (be.next.has_value()) {
        totalLen += getDesiredLengthBetween(e, be.next.value());
        contributions++;
    }

    if (contributions > 0) {
        real_t avgLen = totalLen / (real_t)contributions;
        newLen = avgLen;
    }

    // Write back current length
    const real_t prevLen = be.l;
    be.l = newLen;

    // If length changed, update collider/visual geometry to reflect new x-extent
    const real_t eps = (real_t)1e-8;
    if (std::abs(be.l - prevLen) > eps) {
        bool shapeChanged = false;
        // Visual cube/capsule/cylinder
        if (m_reg.any_of<C_Cube>(e)) {
            auto& cb = m_reg.get<C_Cube>(e);
            const real_t newHalfX = be.l * (real_t)0.5;
            if (std::abs(cb.halfExtents.x() - newHalfX) > eps) { cb.halfExtents.x() = newHalfX; shapeChanged = true; }
        }
        if (m_reg.any_of<C_Capsule>(e)) {
            auto& cap = m_reg.get<C_Capsule>(e);
            const real_t newHalf = be.l * (real_t)0.5;
            if (std::abs(cap.halfLength - newHalf) > eps) { cap.halfLength = newHalf; shapeChanged = true; }
        }
        if (m_reg.any_of<C_Cylinder>(e)) {
            auto& cyl = m_reg.get<C_Cylinder>(e);
            const real_t newHalf = be.l * (real_t)0.5;
            if (std::abs(cyl.halfLength - newHalf) > eps) { cyl.halfLength = newHalf; shapeChanged = true; }
        }
        // Collider cube/capsule
        if (m_reg.any_of<C_RB_Cube>(e)) {
            auto& cb = m_reg.get<C_RB_Cube>(e);
            const real_t newHalfX = be.l * (real_t)0.5;
            if (std::abs(cb.halfExtents.x() - newHalfX) > eps) { cb.halfExtents.x() = newHalfX; shapeChanged = true; }
        }
        if (m_reg.any_of<C_RB_Capsule>(e)) {
            auto& cap = m_reg.get<C_RB_Capsule>(e);
            const real_t newHalf = be.l * (real_t)0.5;
            if (std::abs(cap.halfLength - newHalf) > eps) { cap.halfLength = newHalf; shapeChanged = true; }
        }
        if (m_reg.any_of<C_RB_Cylinder>(e)) {
            auto& cyl = m_reg.get<C_RB_Cylinder>(e);
            const real_t newHalf = be.l * (real_t)0.5;
            if (std::abs(cyl.halfLength - newHalf) > eps) { cyl.halfLength = newHalf; shapeChanged = true; }
        }
        if (shapeChanged) {
            markStructureDirty();
        }
    }
}

void PhysicsSystem::updateEntities() {
    auto view = m_reg.view<C_BeamElement, const C_Position3>();
    for (auto [e, be, pos] : view.each()) {
        (void)be; (void)pos;
        updateBeamElementEntity(e);
    }
}

} // namespace cardillo::