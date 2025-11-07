#include "physics_system.hpp"
#include "assets.hpp"
#include "constraints.hpp"
#include "../misc/spline.hpp"
#include "../collision/collision_coal.hpp"
#include "../solver/warmstart.hpp"
#include "../io/softbody_loader.hpp"
#include <coal/mesh_loader/loader.h>
#include <coal/BVH/BVH_model.h>
#include <sstream>
#include <iomanip>
#include <fstream>
#include "../io/heightmap_loader.hpp"
#include <coal/hfield.h>
#include <coal/shape/geometric_shapes.h>

namespace cardillo {

// Compute diagonal inertia for a box with half-extents (body frame), mass m
inline Vector3r boxInertiaDiag(real_t m, const Vector3r& he) {
    // For a box with full extents 2*he, Ixx = 1/12 m ((2hy)^2 + (2hz)^2) = 1/3 m (hy^2 + hz^2)
    Vector3r I;
    I.x() = (real_t)1.0/3.0 * m * (he.y()*he.y() + he.z()*he.z());
    I.y() = (real_t)1.0/3.0 * m * (he.x()*he.x() + he.z()*he.z());
    I.z() = (real_t)1.0/3.0 * m * (he.x()*he.x() + he.y()*he.y());
    return I;
}
// Compute diagonal inertia for a solid sphere of radius r and mass m (body frame)
inline Vector3r sphereInertiaDiag(real_t m, real_t r) {
    const real_t c = (real_t)2.0/5.0 * m * r * r;
    return Vector3r(c, c, c);
}

inline Vector3r capsuleInertiaDiag(real_t m, real_t radius, real_t halfLength) {
    coal::Capsule capsule((coal::CoalScalar)radius, (coal::CoalScalar)(halfLength * 2));
    const real_t volume = static_cast<real_t>(capsule.computeVolume());
    if (volume <= (real_t)0) return Vector3r::Zero();
    const auto Iunit = capsule.computeMomentofInertia();
    const real_t scale = m / volume;
    return Vector3r(static_cast<real_t>(Iunit(0,0)),
                    static_cast<real_t>(Iunit(1,1)),
                    static_cast<real_t>(Iunit(2,2))) * scale;
}

// Helper to populate common rigid body components
inline void emplaceRigidBodyCommon(entt::registry& reg,
                                   entt::entity e,
                                   real_t mass,
                                   const Vector3r& position,
                                   const Quaternion4r& orientation,
                                   const Vector3r& linearVelocity,
                                   const Vector3r& angularVelocity,
                                   real_t friction_default_mu)
{
    reg.emplace<PhysicsSystem::C_PhysicsObject>(e);
    reg.emplace<PhysicsSystem::C_RigidBodyTag>(e);
    reg.emplace<PhysicsSystem::C_Collidable>(e);
    reg.emplace<PhysicsSystem::C_VisualObject>(e);
    reg.emplace<PhysicsSystem::C_Position3>(e, PhysicsSystem::C_Position3{position});
    reg.emplace<PhysicsSystem::C_Orientation>(e, PhysicsSystem::C_Orientation{orientation});
    reg.emplace<PhysicsSystem::C_LinearVelocity3>(e, PhysicsSystem::C_LinearVelocity3{linearVelocity});
    reg.emplace<PhysicsSystem::C_AngularVelocity3>(e, PhysicsSystem::C_AngularVelocity3{angularVelocity});
    reg.emplace<PhysicsSystem::C_Mass>(e, PhysicsSystem::C_Mass{mass});
    if (!reg.any_of<PhysicsSystem::C_Friction>(e)) reg.emplace<PhysicsSystem::C_Friction>(e, PhysicsSystem::C_Friction{friction_default_mu});
}

// Construction and configuration
PhysicsSystem::PhysicsSystem() {
    PetscPrintf(PETSC_COMM_WORLD, "Hello from PETSc + Eigen!\n");
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

// Unified rigid body creation ------------------------------------------------
entt::entity PhysicsSystem::addRigidBody(const RigidShape& shape,
                                         const RigidState& state,
                                         const RigidProps& props) {
    auto e = m_reg.create();

    // Always pose components
    m_reg.emplace<C_Position3>(e, C_Position3{state.position});
    m_reg.emplace<C_Orientation>(e, C_Orientation{Quaternion4r(state.orientation).normalized()});
    m_reg.emplace<C_LinearVelocity3>(e, C_LinearVelocity3{state.linearVelocity});
    m_reg.emplace<C_AngularVelocity3>(e, C_AngularVelocity3{state.angularVelocity});

    if (props.visual)    m_reg.emplace<C_VisualObject>(e);
    if (props.collidable)m_reg.emplace<C_Collidable>(e);

    // Determine mass from props.mass or density + shape volume
    std::optional<real_t> massOpt = props.mass;
    real_t computedMass = (real_t)0;
    auto computeVolumeCube = [](const Vector3r& he){ return (real_t)8 * he.x()*he.y()*he.z(); };
    auto computeVolumeCapsule = [](real_t r, real_t h){ return (real_t)M_PI * r*r * (2*h + (real_t)4.0/3.0 * r); };
    auto computeVolumeSphere = [](real_t r){ return (real_t)4.0/3.0 * (real_t)M_PI * r*r*r; };
    real_t densityUsed = props.density.value_or((real_t)0);

    // Pre-extract shape volume if density is provided and mass absent
    if (!massOpt.has_value() && props.density.has_value()) {
        std::visit([&](auto&& s){ 
            using T = std::decay_t<decltype(s)>; 
            if constexpr (std::is_same_v<T, CubeShape>) computedMass = densityUsed * computeVolumeCube(s.halfExtents); 
            else if constexpr (std::is_same_v<T, CapsuleShape>) computedMass = densityUsed * computeVolumeCapsule(s.radius, s.halfLength); 
            else if constexpr (std::is_same_v<T, SphereShape>) computedMass = densityUsed * computeVolumeSphere(s.radius); 
            else if constexpr (std::is_same_v<T, PlaneShape>) computedMass = 0; 
          else if constexpr (std::is_same_v<T, MeshShape>) { const ::cardillo::MeshAsset& ma = assets().getMesh(s.path, s.scale, true); 
                 if (ma.volume > (real_t)0) computedMass = densityUsed * ma.volume; } }, shape);
        if (computedMass > (real_t)0) massOpt = computedMass;
    }

    const real_t mass = std::max((real_t)0, massOpt.value_or((real_t)0));
    

    // Friction
    real_t mu = (props.friction >= (real_t)0) ? props.friction : m_cfg.friction_default_mu;
    m_reg.emplace<C_Friction>(e, C_Friction{mu});

    // Populate according to shape variant
    std::visit([&](auto&& s){
        using T = std::decay_t<decltype(s)>;
        if constexpr (std::is_same_v<T, CubeShape>) {
            if (props.visual)    m_reg.emplace<C_CubeVisualTag>(e);
            if (props.collidable)m_reg.emplace<C_RB_Cube>(e, C_RB_Cube{s.halfExtents});
            m_reg.emplace<C_Cube>(e, C_Cube{s.halfExtents});
            if (mass > 0) {
                m_reg.emplace<C_PhysicsObject>(e); m_reg.emplace<C_RigidBodyTag>(e); m_reg.emplace<C_Mass>(e, C_Mass{mass});
                m_reg.emplace<C_InertiaDiag>(e, C_InertiaDiag{boxInertiaDiag(mass, s.halfExtents)});
            }
        } else if constexpr (std::is_same_v<T, SphereShape>) {
            if (props.visual)    m_reg.emplace<C_PointVisualTag>(e);
            if (props.collidable)m_reg.emplace<C_RB_Sphere>(e);
            m_reg.emplace<C_Radius>(e, C_Radius{s.radius});
            if (mass > 0) {
                m_reg.emplace<C_PhysicsObject>(e); m_reg.emplace<C_RigidBodyTag>(e); m_reg.emplace<C_Mass>(e, C_Mass{mass});
                m_reg.emplace<C_InertiaDiag>(e, C_InertiaDiag{sphereInertiaDiag(mass, s.radius)});
            }
        } else if constexpr (std::is_same_v<T, CapsuleShape>) {
            if (props.visual)    m_reg.emplace<C_CapsuleVisualTag>(e);
            if (props.collidable)m_reg.emplace<C_RB_Capsule>(e, C_RB_Capsule{s.radius, s.halfLength});
            m_reg.emplace<C_Capsule>(e, C_Capsule{s.radius, s.halfLength});
            if (mass > 0) {
                m_reg.emplace<C_PhysicsObject>(e); m_reg.emplace<C_RigidBodyTag>(e); m_reg.emplace<C_Mass>(e, C_Mass{mass});
                m_reg.emplace<C_InertiaDiag>(e, C_InertiaDiag{capsuleInertiaDiag(mass, s.radius, s.halfLength)});
            }
        } else if constexpr (std::is_same_v<T, PlaneShape>) {
            if (props.visual)    m_reg.emplace<C_PlaneVisualTag>(e);
            if (props.collidable)m_reg.emplace<C_RB_Plane>(e, C_RB_Plane{s.normal, s.up, s.sizeX, s.sizeY});
            m_reg.emplace<C_Plane>(e, C_Plane{s.normal, s.up, s.sizeX, s.sizeY});

        } else if constexpr (std::is_same_v<T, MeshShape>) {
            if (props.visual)    m_reg.emplace<C_MeshVisualTag>(e);
            if (props.collidable)m_reg.emplace<C_RB_Mesh>(e);
            m_reg.emplace<C_Mesh>(e, C_Mesh{s.path, s.scale});
            const bool dynamic = mass > 0;
            const ::cardillo::MeshAsset& asset = assets().getMesh(s.path, s.scale, dynamic);
            if (dynamic) {
                // Adjust pose by principal axes & COM
                Quaternion4r q_rpa(asset.Rpa);
                Quaternion4r q_new = state.orientation * q_rpa;
                Vector3r pos_new = state.position + (state.orientation * asset.com);
                m_reg.get<C_Position3>(e).value = pos_new;
                m_reg.get<C_Orientation>(e).q = q_new;
                m_reg.emplace<C_PhysicsObject>(e); m_reg.emplace<C_RigidBodyTag>(e); m_reg.emplace<C_Mass>(e, C_Mass{mass});
                if (asset.volume > (real_t)0) {
                    const real_t rho = mass / asset.volume;
                    Vector3r Idiag = rho * asset.inertia_diag_unit.cwiseMax(Vector3r::Zero());
                    m_reg.emplace<C_InertiaDiag>(e, C_InertiaDiag{Idiag});
                }
            }
        }
    }, shape);

    markStructureDirty();
    return e;
}

entt::entity PhysicsSystem::addStaticBody(const RigidShape& shape,
                                          const RigidState& state) {
    RigidProps staticProps; // neither mass nor density set => static
    return addRigidBody(shape, state, staticProps);
}

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
    // Default friction coefficient from config unless already set
    if (!m_reg.any_of<C_Friction>(e)) m_reg.emplace<C_Friction>(e, C_Friction{m_cfg.friction_default_mu});
    markStructureDirty();
    return static_cast<index_t>(entt::to_integral(e));
}
// Helper: create a purely visual/collidable rigid-like entity (no physics)
entt::entity PhysicsSystem::createRigidVisualEntity_(const Vector3r& center) {
    auto e = m_reg.create();
    m_reg.emplace<C_Collidable>(e);
    m_reg.emplace<C_VisualObject>(e);
    m_reg.emplace<C_Position3>(e, C_Position3{center});
    markStructureDirty();
    return e;
}


// Static HeightField obstacle (terrain) from EXR
index_t PhysicsSystem::addObstacleHeightField(const Vector3r& position,
                                              const Quaternion4r& orientation,
                                              const std::string& exrPath,
                                              real_t x_dim,
                                              real_t y_dim,
                                              real_t z_scale,
                                              real_t min_height) {
    auto e = m_reg.create();
    m_reg.emplace<C_Collidable>(e);
    m_reg.emplace<C_VisualObject>(e);
    m_reg.emplace<C_Position3>(e, C_Position3{position});
    m_reg.emplace<C_Orientation>(e, C_Orientation{orientation});
    // Visual + collider markers
    m_reg.emplace<C_HeightFieldVisualTag>(e);
    m_reg.emplace<C_HeightField>(e, C_HeightField{exrPath, x_dim, y_dim, z_scale, min_height});
    m_reg.emplace<C_RB_HeightField>(e);
    if (!m_reg.any_of<C_Friction>(e)) m_reg.emplace<C_Friction>(e, C_Friction{m_cfg.friction_default_mu});
    markStructureDirty();
    // Touch the asset so failures surface early
    (void)getHeightFieldAsset(e);
    return static_cast<index_t>(entt::to_integral(e));
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

// Create a mass-spring soft body from an OBJ: one point-mass per vertex, springs along triangle edges
// Overload with pose + velocities + total mass (uniform per-vertex mass)
std::vector<entt::entity> PhysicsSystem::addSoftBody(const std::string& objPath,
                                                     real_t stiffness,
                                                     real_t damping,
                                                     const Vector3r& position,
                                                     const Quaternion4r& orientation,
                                                     const Vector3r& linearVelocity,
                                                     const Vector3r& angularVelocity,
                                                     real_t totalMass)
{
    std::vector<entt::entity> nodes;
    cardillo::io::SoftBodyMesh sb;
    if (!cardillo::io::load_obj_softbody(objPath, sb)) {
        std::printf("[SoftBody] Failed to load OBJ: %s\n", objPath.c_str());
        return nodes;
    }

    const size_t N = sb.positions.size();
    if (N == 0) return nodes;

    // Uniform mass per vertex (fallback to small default if non-positive totalMass)
    real_t nodeMass = (totalMass > (real_t)0) ? (totalMass / (real_t)N) : (real_t)0.02;
    const real_t nodeRadius = (real_t)0.02;

    Matrix33r R = orientation.toRotationMatrix();
    nodes.reserve(N);
    for (const auto& p0 : sb.positions) {
        Vector3r pw = position + R * p0;
        Vector3r vw = linearVelocity + angularVelocity.cross(pw - position);
        index_t id = addPointMass(nodeMass, pw, vw, nodeRadius);
        nodes.push_back(entt::entity(static_cast<uint32_t>(id)));
    }

    // Edge springs (create distance constraints between edge vertices)
    for (const auto& e : sb.edges) {
        int i = e.first;
        int j = e.second;
        if (i >= 0 && j >= 0 && (size_t)i < nodes.size() && (size_t)j < nodes.size()) {
            entt::entity A = nodes[(size_t)i];
            entt::entity B = nodes[(size_t)j];
            addConstraint<physics::LinearDistanceConstraint>(m_reg, A, B, Vector3r::Zero(), Vector3r::Zero(), stiffness, damping);
        }
    }

    // Create a visual surface entity to represent the softbody as a mesh
    if (!sb.triangles.empty()) {
        entt::entity surf = m_reg.create();
        m_reg.emplace<C_VisualObject>(surf);
        m_reg.emplace<C_SoftBodyVisualTag>(surf);
        C_SoftBodySurface surfComp;
        surfComp.triangles = sb.triangles;
        surfComp.nodes = nodes;
        m_reg.emplace<C_SoftBodySurface>(surf, std::move(surfComp));
    }

    markStructureDirty();
    return nodes;
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

VectorXr PhysicsSystem::getForce(entt::entity e) const {
    // Prefer rigid body
    if (m_reg.all_of<C_RigidBodyTag, C_PhysicsObject, C_Mass>(e)) {
        const real_t m = m_reg.get<C_Mass>(e).m;
        Vector3r fg = m * m_gravity;
        Vector3r tau = Vector3r::Zero();
        if (m_reg.all_of<C_InertiaDiag>(e)) {
            const Vector3r w = m_reg.get<C_AngularVelocity3>(e).value; // body-frame
            Vector3r Idiag = m_reg.get<C_InertiaDiag>(e).I;
            const Vector3r Iw = Idiag.cwiseProduct(w);
            tau = -w.cross(Iw);
        }
        // Add any external one-shot forces/torques (world-frame)
        if (m_reg.any_of<C_ExternalForce>(e)) {
            fg += m_reg.get<C_ExternalForce>(e).f;
        }
        if (m_reg.any_of<C_ExternalTorque>(e)) {
            tau += m_reg.get<C_ExternalTorque>(e).tau;
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
        if (m_reg.any_of<C_ExternalForce>(e)) {
            fg += m_reg.get<C_ExternalForce>(e).f;
        }
        VectorXr out(3); out << fg[0], fg[1], fg[2];
        return out;
    }
    return VectorXr(0);
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
        return boxInertiaDiag(m, he);
    }
    // Default: zero
    return Vector3r::Zero();
}

int PhysicsSystem::numBodies() const {
    if (m_num_bodies_dirty) {
        int count = 0;
        auto view = m_reg.view<C_BodyIndex, C_PhysicsObject>();
        for (auto e : view) (void)e, ++count;
        m_num_bodies_cached = count;
        m_num_bodies_dirty = false;
    }
    return m_num_bodies_cached;
}

void PhysicsSystem::applyForce(entt::entity e, const Vector3r& force_world, const Vector3r& torque_world) {
    if (!m_reg.valid(e)) return;
    if (force_world.allFinite() && !force_world.isZero()) {
        if (m_reg.any_of<C_ExternalForce>(e)) m_reg.get<C_ExternalForce>(e).f += force_world;
        else m_reg.emplace<C_ExternalForce>(e, C_ExternalForce{force_world});
    }
    if (torque_world.allFinite() && !torque_world.isZero()) {
        if (m_reg.any_of<C_ExternalTorque>(e)) m_reg.get<C_ExternalTorque>(e).tau += torque_world;
        else m_reg.emplace<C_ExternalTorque>(e, C_ExternalTorque{torque_world});
    }
    m_forces_dirty = true;
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
    Quaternion4r q = q_in; q.normalize();
    if (m_reg.any_of<C_Orientation>(e)) m_reg.get<C_Orientation>(e).value = q; else m_reg.emplace<C_Orientation>(e, C_Orientation{q});
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

std::pair<entt::entity, entt::entity> PhysicsSystem::createBeam(const misc::SplinePattern& spline,
                                                                size_t segments,
                                                                real_t width,
                                                                real_t height,
                                                                real_t density,
                                                                real_t E,
                                                                real_t nu) {
    return createBeam(spline, segments, width, height, density, E, nu, 1,1,1,1,1, Vector3r::Zero(), Vector3r::Zero());
}

// Extended overload with stiffness scaling and damping
std::pair<entt::entity, entt::entity> PhysicsSystem::createBeam(const misc::SplinePattern& spline,
                                                                size_t segments,
                                                                real_t width,
                                                                real_t height,
                                                                real_t density,
                                                                real_t E,
                                                                real_t nu,
                                                                real_t axialScale,
                                                                real_t shearScale,
                                                                real_t torsionScale,
                                                                real_t bendYScale,
                                                                real_t bendZScale,
                                                                const Vector3r& Dg,
                                                                const Vector3r& Df) {
    using namespace physics;
    const real_t totalLen = spline.totalLength();
    const real_t segLen = totalLen / (real_t)segments;

    // Material and section properties
    const real_t G = E / (2 * (1 + nu));
    const real_t A = width * height;
    const real_t Iy = width * std::pow(height, 3) / 12;
    const real_t Iz = std::pow(width, 3) * height / 12;
    const real_t Ip = Iy + Iz;

    Vector3r Ke(E * A / segLen , G * A/ segLen, G * A / segLen);
    Vector3r Kf(G * Ip / segLen, E * Iy / segLen, E * Iz / segLen);

    // Apply user scaling to decouple rope-like properties
    Ke.x() *= axialScale;        // stretch
    Ke.y() *= shearScale;        // shear Y
    Ke.z() *= shearScale;        // shear Z
    Kf.x() *= torsionScale;      // torsion (GJ)
    Kf.y() *= bendYScale;        // bending about local Y (EI_y)
    Kf.z() *= bendZScale;        // bending about local Z (EI_z)

    // Rigid body shape
    PhysicsSystem::CubeShape shape(Vector3r(segLen / 2, width / 2, height / 2));

    // Mass per segment from density and volume
    const real_t mass = density * (A * segLen);

    Vector3r vlin = Vector3r::Zero();
    Vector3r vang = Vector3r::Zero();

    // Root
    misc::SplineSample s0 = spline.sample((real_t)0);
    Matrix33r R0; R0.col(0)=s0.tangent; R0.col(1)=s0.normal; R0.col(2)=s0.binormal;
    Quaternion4r q0(R0); q0.normalize();
    entt::entity a = addRigidBody(shape, RigidState(s0.position, vlin, q0, vang), RigidProps(mass));

    entt::entity root = a;
    entt::entity end = a;

    index_t segments_to_place = spline.isLoop() ? segments - 1 : segments;
    for (size_t i = 1; i <= segments_to_place; ++i) {
        real_t alpha = (real_t)i / (real_t)segments;
        misc::SplineSample si = spline.sample(alpha);
        Matrix33r Ri; Ri.col(0)=si.tangent; Ri.col(1)=si.normal; Ri.col(2)=si.binormal;
        Quaternion4r qi(Ri); qi.normalize();
        entt::entity b = addRigidBody(shape, RigidState(si.position, vlin, qi, vang), RigidProps(mass));
        addConstraint<BeamConstraint>(ecs(), a, b, Ke, Kf, Dg, Df);
        disableCollisionBetween(a, b);
        a = b;
        end = b;
    }

    if (spline.isLoop()) {
        addConstraint<BeamConstraint>(ecs(), a, root, Ke, Kf, Dg, Df);
        disableCollisionBetween(a, root);
    }

    return {root, end};
}

void PhysicsSystem::explicitPositionUpdate(real_t h) {
    // position update
    auto position_view = m_reg.view<C_Position3, C_LinearVelocity3>();
    for (auto [e, pos, vel] : position_view.each()) {
        pos.value += h * vel.value;
    }

    // orientations
    auto orientation_view = m_reg.view<C_Orientation, C_AngularVelocity3>();
    for (auto [e, orientation, angularVel] : orientation_view.each()) {
        const Vector3r& omega = angularVel.value;
        Vector4r& P = orientation.value.coeffs();
        real_t w = P(3);
        P(3) -= h * 0.5 * P.head<3>().dot(omega);
        P.head<3>() += h * 0.5 * (w * omega + P.head<3>().cross(omega));
    }
}

void PhysicsSystem::linearImplicitPositionUpdate(real_t h) {
    // TODO: This has to be implemented (do explicit update instead)
    explicitPositionUpdate(h);
}

} // namespace cardillo::