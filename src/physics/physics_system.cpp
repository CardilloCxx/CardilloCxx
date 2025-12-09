#include "physics_system.hpp"
#include "assets.hpp"
#include "constraints.hpp"
#include "../misc/spline.hpp"
#include "../collision/collision_coal.hpp"
#include "../solver/warmstart.hpp"
#include "../io/softbody_loader.hpp"
#include "../io/csv_writer.hpp"
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
    if (props.collidable && !m_cfg.collision_disable_all) m_reg.emplace<C_Collidable>(e);

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
            if (props.collidable && !m_cfg.collision_disable_all) m_reg.emplace<C_RB_Cube>(e, C_RB_Cube{s.halfExtents});
            m_reg.emplace<C_Cube>(e, C_Cube{s.halfExtents});
            if (mass > 0) {
                m_reg.emplace<C_PhysicsObject>(e); m_reg.emplace<C_RigidBodyTag>(e); m_reg.emplace<C_Mass>(e, C_Mass{mass});
                m_reg.emplace<C_InertiaDiag>(e, C_InertiaDiag{boxInertiaDiag(mass, s.halfExtents)});
            }
        } else if constexpr (std::is_same_v<T, SphereShape>) {
            if (props.visual)    m_reg.emplace<C_PointVisualTag>(e);
            if (props.collidable && !m_cfg.collision_disable_all) m_reg.emplace<C_RB_Sphere>(e);
            m_reg.emplace<C_Radius>(e, C_Radius{s.radius});
            if (mass > 0) {
                m_reg.emplace<C_PhysicsObject>(e); m_reg.emplace<C_RigidBodyTag>(e); m_reg.emplace<C_Mass>(e, C_Mass{mass});
                m_reg.emplace<C_InertiaDiag>(e, C_InertiaDiag{sphereInertiaDiag(mass, s.radius)});
            }
        } else if constexpr (std::is_same_v<T, CapsuleShape>) {
            if (props.visual)    m_reg.emplace<C_CapsuleVisualTag>(e);
            if (props.collidable && !m_cfg.collision_disable_all) m_reg.emplace<C_RB_Capsule>(e, C_RB_Capsule{s.radius, s.halfLength});
            m_reg.emplace<C_Capsule>(e, C_Capsule{s.radius, s.halfLength});
            if (mass > 0) {
                m_reg.emplace<C_PhysicsObject>(e); m_reg.emplace<C_RigidBodyTag>(e); m_reg.emplace<C_Mass>(e, C_Mass{mass});
                m_reg.emplace<C_InertiaDiag>(e, C_InertiaDiag{capsuleInertiaDiag(mass, s.radius, s.halfLength)});
            }
        } else if constexpr (std::is_same_v<T, PlaneShape>) {
            if (props.visual)    m_reg.emplace<C_PlaneVisualTag>(e);
            if (props.collidable && !m_cfg.collision_disable_all) m_reg.emplace<C_RB_Plane>(e, C_RB_Plane{s.normal, s.up, s.sizeX, s.sizeY});
            m_reg.emplace<C_Plane>(e, C_Plane{s.normal, s.up, s.sizeX, s.sizeY});

        } else if constexpr (std::is_same_v<T, MeshShape>) {
            if (props.visual)    m_reg.emplace<C_MeshVisualTag>(e);
            if (props.collidable && !m_cfg.collision_disable_all) m_reg.emplace<C_RB_Mesh>(e);
            m_reg.emplace<C_Mesh>(e, C_Mesh{s.path, s.scale});
            const bool dynamic = mass > 0;
            const ::cardillo::MeshAsset& asset = assets().getMesh(s.path, s.scale, dynamic);
            if (dynamic) {
                // Adjust pose by principal axes & COM
                Quaternion4r q_rpa(asset.Rpa);
                Quaternion4r q_new = state.orientation * q_rpa;
                Vector3r pos_new = state.position + (state.orientation * asset.com);
                m_reg.get<C_Position3>(e).value = pos_new;
                m_reg.get<C_Orientation>(e).value = q_new;
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

void PhysicsSystem::track(entt::entity e, const std::string& name)
{
    if (!m_reg.valid(e)) return;
    if (!m_reg.any_of<C_TrackTag>(e)) {
        m_reg.emplace<C_TrackTag>(e, C_TrackTag{name});
    } else {
        m_reg.get<C_TrackTag>(e).name = name;
    }
}

// Entity creation
index_t PhysicsSystem::addPointMass(real_t mass, const Vector3r& x0, const Vector3r& v0, real_t radius) {
    auto e = m_reg.create();
    m_reg.emplace<C_PhysicsObject>(e);
    m_reg.emplace<C_PointMassTag>(e);
    if (!m_cfg.collision_disable_all) m_reg.emplace<C_Collidable>(e);
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
    if (!m_cfg.collision_disable_all) m_reg.emplace<C_Collidable>(e);
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
    if (!m_cfg.collision_disable_all) m_reg.emplace<C_Collidable>(e);
    m_reg.emplace<C_VisualObject>(e);
    m_reg.emplace<C_Position3>(e, C_Position3{position});
    m_reg.emplace<C_Orientation>(e, C_Orientation{orientation});
    // Visual + collider markers
    m_reg.emplace<C_HeightFieldVisualTag>(e);
    m_reg.emplace<C_HeightField>(e, C_HeightField{exrPath, x_dim, y_dim, z_scale, min_height});
    if (!m_cfg.collision_disable_all) m_reg.emplace<C_RB_HeightField>(e);
    if (!m_reg.any_of<C_Friction>(e)) m_reg.emplace<C_Friction>(e, C_Friction{m_cfg.friction_default_mu});
    markStructureDirty();
    // Touch the asset so failures surface early
    (void)getHeightFieldAsset(e);
    return static_cast<index_t>(entt::to_integral(e));
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
            "wx","wy","wz"
        };
        writer.open(path, header);
        initialized = true;
    }
    if (!writer.isOpen()) return;

    auto view = m_reg.view<C_TrackTag, C_Position3, C_LinearVelocity3, C_AngularVelocity3>();
    for (auto e : view) {
        const auto& tag   = view.get<C_TrackTag>(e);
        const auto& pos   = view.get<C_Position3>(e).value;
        const auto& v     = view.get<C_LinearVelocity3>(e).value;
        const auto& w     = view.get<C_AngularVelocity3>(e).value;
        writer.writeRow(
            t,
            tag.name,
            pos.x(), pos.y(), pos.z(),
            v.x(), v.y(), v.z(),
            w.x(), w.y(), w.z()
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
        // Add any external one-shot forces/torques
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

std::pair<entt::entity, entt::entity> PhysicsSystem::createBeam(const misc::SplinePattern& spline,
                                                               const BeamCrossSection& section,
                                                               const BeamSpringParams& springs,
                                                               const RigidState& stateDefaults,
                                                               const RigidProps& propsDefaults,
                                                               size_t segments) {
    using namespace physics;
    const real_t totalLen = spline.totalLength();
    const real_t segLen   = totalLen / (real_t)segments;

    // Mass per segment
    real_t massPerSegment = 0;
    if (propsDefaults.mass.has_value()) massPerSegment = *propsDefaults.mass / (real_t)segments;
    else if (propsDefaults.density.has_value()) massPerSegment = *propsDefaults.density * (section.area() * segLen);

    // Shape
    RigidShape shape;
    Matrix33r Rshape = Matrix33r::Identity();
    if (section.type == BeamBodyType::Cube) {
        shape = CubeShape(Vector3r(segLen * (real_t)0.5, section.width*(real_t)0.5, section.height*(real_t)0.5));
    } else {
        real_t r = std::min(section.width, section.height) * (real_t)0.5;
        shape = CapsuleShape(r, segLen * (real_t)0.5);
        Rshape = Quaternion4r::FromTwoVectors(Vector3r::UnitZ(), Vector3r::UnitX()).toRotationMatrix();
    }

    RigidProps segProps = propsDefaults;
    segProps.mass = (massPerSegment > 0 ? std::optional<real_t>(massPerSegment) : std::nullopt);

    // Per-segment stiffness (material or direct overrides)
    const Vector3r Ke = springs.Ke(segLen, section);
    const Vector3r Kf = springs.Kf(segLen, section);
    const Vector3r De = springs.De;
    const Vector3r Df = springs.Df;

    entt::entity root = entt::null;
    entt::entity prev = entt::null;
    entt::entity end  = entt::null;

    const bool loop = spline.isLoop();
    // Compute spline COM in world coordinates and compose with state
    Vector3r splineCOMWorld = spline.centerOfMass();
    Matrix33r Rbody = stateDefaults.orientation.toRotationMatrix();
    // Rotating about COM keeps COM itself unchanged; state translation shifts COM
    Vector3r worldCOM = splineCOMWorld + stateDefaults.position;
    // Precompute state (body-frame) velocities expressed in world frame
    Vector3r v_body_world = Rbody * stateDefaults.linearVelocity;
    Vector3r w_body_world = Rbody * stateDefaults.angularVelocity; // body-frame to world

    Quaternion4r q_prev = Quaternion4r::Identity();

    for (size_t i = 0; i < segments; ++i) {
        real_t alpha = (real_t)i / (real_t)segments; // alpha in [0, 1)
        misc::SplineSample si = spline.sample(alpha);
        Matrix33r Rlocal; Rlocal.col(0)=si.tangent; Rlocal.col(1)=si.normal; Rlocal.col(2)=si.binormal;
        Matrix33r Rworld =  Rbody * Rlocal * Rshape;
        Quaternion4r qworld(Rworld); qworld.normalize();
        
        // Ensure quaternion continuity
        if (q_prev.dot(qworld) < (real_t)0) qworld.coeffs() = -qworld.coeffs();
        q_prev = qworld;

        // Rotate about spline COM, then add COM shift + state translation
        Vector3r worldPos = splineCOMWorld + stateDefaults.position + Rbody * (si.position - splineCOMWorld);

        // v = v_body_world + w_body_world x (worldPos - worldCOM)
        Vector3r v_world = v_body_world + w_body_world.cross(worldPos - worldCOM);
        RigidState segState;
        segState.position = worldPos;
        segState.orientation = qworld;
        segState.linearVelocity = v_world;
        // Convert provided body-frame angular velocity (state frame) to final body frame (state*spline)
        segState.angularVelocity = Rlocal.transpose() * stateDefaults.angularVelocity;
        entt::entity cur = addRigidBody(shape, segState, segProps);

        // Beam element component and neighbor setup (data lives on entity)
        C_BeamElement be_cur;
        be_cur.l0 = segLen;
        be_cur.l  = segLen;
        if (prev != entt::null) {
            be_cur.prev = prev;
            // Ensure prev has a component and set its next
            if (!m_reg.any_of<C_BeamElement>(prev)) {
                C_BeamElement be_prev;
                be_prev.l0 = segLen;
                be_prev.l  = segLen;
                be_prev.next = cur;
                m_reg.emplace<C_BeamElement>(prev, be_prev);
            } else {
                m_reg.get<C_BeamElement>(prev).next = cur;
            }
        }
        m_reg.emplace<C_BeamElement>(cur, be_cur);

        if (prev != entt::null) {
            addConstraint<BeamConstraint>(ecs(), prev, cur, springs, section);
            disableCollisionBetween(prev, cur);
        }
        if (root == entt::null) root = cur;
        prev = cur;
        end = cur;
    }
    if (loop && root != entt::null && end != entt::null && end != root) {
        addConstraint<BeamConstraint>(ecs(), end, root, springs, section);
        disableCollisionBetween(end, root);
        // Close neighbor links for looped beam
        if (m_reg.any_of<C_BeamElement>(end)) m_reg.get<C_BeamElement>(end).next = root;
        if (m_reg.any_of<C_BeamElement>(root)) m_reg.get<C_BeamElement>(root).prev = end;
    }
    return {root, end};
}
std::pair<entt::entity, entt::entity> PhysicsSystem::createBeams(const std::vector<const misc::SplinePattern*>& splines,
                                                                  const BeamCrossSection& section,
                                                                  const BeamSpringParams& springs,
                                                                  const RigidState& stateDefaults,
                                                                  const RigidProps& propsDefaults,
                                                                  size_t segments) {
    real_t totalLen = 0;
    for (const auto* sp : splines) {
        if (sp) totalLen += sp->totalLength();
    }

    entt::entity first = entt::null;;
    entt::entity second = entt::null;
    entt::entity prevEnd = entt::null;
    for (size_t i = 0; i < splines.size(); ++i) {

        auto pair = createBeam(*splines[i], section, springs, stateDefaults, propsDefaults, (size_t) (segments * (splines[i]->totalLength() / totalLen)));
        if (first == entt::null) first = pair.first;
        if (prevEnd != entt::null && pair.first != entt::null) {
            addConstraint<physics::RigidConstraint>(ecs(), prevEnd, pair.first);
            disableCollisionBetween(prevEnd, pair.first);
        }
        prevEnd = pair.second;   
    }
    second = prevEnd;
    return {first, second};
}
void PhysicsSystem::explicitPositionUpdate(real_t h) {
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

        Vector4r& P = orientation.value.coeffs();
        real_t w = P(3);
        P(3) -= h * 0.5 * P.head<3>().dot(omega);
        P.head<3>() += h * 0.5 * (w * omega + P.head<3>().cross(omega));

        // re-normalize to avoid quaternion drift
        orientation.value.normalize();
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

        // re-normalize to avoid quaternion drift
        orientation.value.normalize();
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
            
            index_t x_col_A = m_reg.any_of<C_Capsule>(a) ? 2 : 0; 
            index_t x_col_B = m_reg.any_of<C_Capsule>(b) ? 2 : 0;

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
    else std::cout << "Warning: beam element has no neighbors to compute length from.\n";

    // Write back current length
    const real_t prevLen = be.l;
    be.l = newLen;

    // If length changed, update collider/visual geometry to reflect new x-extent
    const real_t eps = (real_t)1e-8;
    if (std::abs(be.l - prevLen) > eps) {
        bool shapeChanged = false;
        // Visual cube/capsule
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