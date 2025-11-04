#include "physics_system.hpp"
#include "constraints.hpp"
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

// Obstacles (static visuals)
index_t PhysicsSystem::addObstacleBody(const Plane& p) {
    auto e = createRigidVisualEntity_(p.center);
    m_reg.emplace<C_PlaneVisualTag>(e);
    m_reg.emplace<C_Plane>(e, C_Plane{p.normal, p.up, p.sizeX, p.sizeY});
    // Rigid-body default kinematics
    m_reg.emplace<C_Orientation>(e, C_Orientation{Quaternion4r::Identity()});
    m_reg.emplace<C_LinearVelocity3>(e, C_LinearVelocity3{Vector3r::Zero()});
    m_reg.emplace<C_AngularVelocity3>(e, C_AngularVelocity3{Vector3r::Zero()});
    // Attach unified rigid-body plane tag
    m_reg.emplace<C_RB_Plane>(e, C_RB_Plane{p.normal, p.up, p.sizeX, p.sizeY});
    if (!m_reg.any_of<C_Friction>(e)) m_reg.emplace<C_Friction>(e, C_Friction{m_cfg.friction_default_mu});
    return static_cast<index_t>(entt::to_integral(e));
}

index_t PhysicsSystem::addObstacleBody(const Cube& c) {
    auto e = createRigidVisualEntity_(c.center);
    m_reg.emplace<C_CubeVisualTag>(e);
    m_reg.emplace<C_Cube>(e, C_Cube{c.halfExtents});
    m_reg.emplace<C_Orientation>(e, C_Orientation{c.q});
    // Rigid-body default vels
    m_reg.emplace<C_LinearVelocity3>(e, C_LinearVelocity3{Vector3r::Zero()});
    m_reg.emplace<C_AngularVelocity3>(e, C_AngularVelocity3{Vector3r::Zero()});
    // Attach unified rigid-body cube tag
    m_reg.emplace<C_RB_Cube>(e, C_RB_Cube{c.halfExtents});
    // Default friction for obstacles as well
    if (!m_reg.any_of<C_Friction>(e)) m_reg.emplace<C_Friction>(e, C_Friction{m_cfg.friction_default_mu});
    return static_cast<index_t>(entt::to_integral(e));
}

// Static mesh obstacle: visual + collider only (no physics object)
index_t PhysicsSystem::addObstacleMesh(const Vector3r& position,
                                       const Quaternion4r& orientation,
                                       const std::string& meshPath,
                                       const Vector3r& scale) {
    auto e = m_reg.create();
    m_reg.emplace<C_Collidable>(e);
    m_reg.emplace<C_VisualObject>(e);
    m_reg.emplace<C_Position3>(e, C_Position3{position});
    m_reg.emplace<C_Orientation>(e, C_Orientation{orientation});
    // Visual + collider markers
    m_reg.emplace<C_MeshVisualTag>(e);
    m_reg.emplace<C_Mesh>(e, C_Mesh{meshPath, scale});
    m_reg.emplace<C_RB_Mesh>(e);
    if (!m_reg.any_of<C_Friction>(e)) m_reg.emplace<C_Friction>(e, C_Friction{m_cfg.friction_default_mu});
    markStructureDirty();
    return static_cast<index_t>(entt::to_integral(e));
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

// Dynamic rigid body
entt::entity PhysicsSystem::addRigidBody(real_t mass,
                                    const Vector3r& position,
                                    const Quaternion4r& orientation,
                                    const Vector3r& linearVelocity,
                                    const Vector3r& angularVelocity,
                                    const Cube& shape) {
    auto e = m_reg.create();
    emplaceRigidBodyCommon(m_reg, e, mass, position, orientation, linearVelocity, angularVelocity, m_cfg.friction_default_mu);
    m_reg.emplace<C_CubeVisualTag>(e);
    m_reg.emplace<C_Cube>(e, C_Cube{shape.halfExtents});
    m_reg.emplace<C_RB_Cube>(e, C_RB_Cube{shape.halfExtents});
    m_reg.emplace<C_InertiaDiag>(e, C_InertiaDiag{boxInertiaDiag(mass, shape.halfExtents)});
    markStructureDirty();
    return e;
}


// Mesh-based rigid body
entt::entity PhysicsSystem::addRigidBodyMesh(real_t mass,
                                        const Vector3r& position,
                                        const Quaternion4r& orientation,
                                        const Vector3r& linearVelocity,
                                        const Vector3r& angularVelocity,
                                        const std::string& meshPath,
                                        const Vector3r& scale) {
    auto e = m_reg.create();
    // Visual + collider markers
    m_reg.emplace<C_MeshVisualTag>(e);
    m_reg.emplace<C_Mesh>(e, C_Mesh{meshPath, scale});
    m_reg.emplace<C_RB_Mesh>(e);
    // Compute and store diagonal inertia using normalized asset metadata if available
    
    emplaceRigidBodyCommon(m_reg, e, mass, position, orientation, linearVelocity, angularVelocity, m_cfg.friction_default_mu);

    const MeshAsset& asset = getMeshAsset(e);

    // Revert to original orientation/position 
    const Matrix33r Rpa = asset.Rpa;
    const Vector3r com = asset.com;
    Quaternion4r q_rpa(Rpa);
    Quaternion4r q_new = orientation * q_rpa;
    Vector3r pos_new = position + (orientation * com);

    // replace position and orientation with adjusted versions
    m_reg.get<C_Position3>(e).value = pos_new;
    m_reg.get<C_Orientation>(e).q = q_new;

    if (asset.bvh && asset.volume != (real_t)0.0) {
        const real_t rho = mass / asset.volume; // density implied by mass
        Vector3r Idiag = rho * asset.inertia_diag_unit.cwiseMax(Vector3r::Zero());
        m_reg.emplace<C_InertiaDiag>(e, C_InertiaDiag{Idiag});
    }
    markStructureDirty();
    return e;
}

// Sphere rigid body
entt::entity PhysicsSystem::addRigidBodySphere(real_t mass,
                                         const Vector3r& position,
                                         const Quaternion4r& orientation,
                                         const Vector3r& linearVelocity,
                                         const Vector3r& angularVelocity,
                                         real_t radius) {
    auto e = m_reg.create();
    emplaceRigidBodyCommon(m_reg, e, mass, position, orientation, linearVelocity, angularVelocity, m_cfg.friction_default_mu);
    // Visualize spheres using point visual with radius; collider via sphere tag + radius
    m_reg.emplace<C_PointVisualTag>(e);
    m_reg.emplace<C_Radius>(e, C_Radius{radius});
    m_reg.emplace<C_RB_Sphere>(e);
    m_reg.emplace<C_InertiaDiag>(e, C_InertiaDiag{sphereInertiaDiag(mass, radius)});
    markStructureDirty();
    return e;
}

entt::entity PhysicsSystem::addRigidBodyCapsule(real_t mass,
                                           const Vector3r& position,
                                           const Quaternion4r& orientation,
                                           const Vector3r& linearVelocity,
                                           const Vector3r& angularVelocity,
                                           const Capsule& shape) {
    auto e = m_reg.create();
    emplaceRigidBodyCommon(m_reg, e, mass, position, orientation, linearVelocity, angularVelocity, m_cfg.friction_default_mu);
    m_reg.emplace<C_CapsuleVisualTag>(e);
    m_reg.emplace<C_Capsule>(e, C_Capsule{shape.radius, shape.halfLength});
    m_reg.emplace<C_RB_Capsule>(e, C_RB_Capsule{shape.radius, shape.halfLength});
    m_reg.emplace<C_InertiaDiag>(e, C_InertiaDiag{capsuleInertiaDiag(mass, shape.radius, shape.halfLength)});
    m_structure_dirty = true;
    return e;
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
        Quaternion4r qn = m_reg.get<C_Orientation>(e).q;
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
    if (m_reg.any_of<C_Orientation>(e)) m_reg.get<C_Orientation>(e).q = q; else m_reg.emplace<C_Orientation>(e, C_Orientation{q});
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

// ---------- Mesh asset cache ----------

namespace {
struct MeshNormalizationResult {
    coal::BVHModelPtr_t bvh;          // normalized BVH (scaled, COM-centered, principal-axes aligned)
    cardillo::Vector3r inertia_diag_unit{cardillo::Vector3r::Zero()};
    real_t volume{(real_t)0.0};
    cardillo::Matrix33r Rpa{cardillo::Matrix33r::Identity()};
    coal::Vec3s com{0,0,0};
};

// Build a normalized BVH and compute unit-density inertia diag + volume
MeshNormalizationResult normalizeMesh(const std::string& meshPath, const cardillo::Vector3r& scale) {
    MeshNormalizationResult out;
    coal::MeshLoader loader(coal::BV_AABB);
    coal::BVHModelPtr_t bvh = loader.load(meshPath);

    if (!bvh) return out;
    // 1) Scale vertices and collect triangles
    std::vector<coal::Vec3s> Vscaled;
    std::vector<coal::Triangle> F;
    if (bvh->vertices) Vscaled = *bvh->vertices;
    if (bvh->tri_indices) F = *bvh->tri_indices;
    for (auto& v : Vscaled) {
        v[0] *= (coal::CoalScalar)scale.x();
        v[1] *= (coal::CoalScalar)scale.y();
        v[2] *= (coal::CoalScalar)scale.z();
    }
    // 2) Compute COM and volume on scaled geometry
    coal::Vec3s com(0,0,0);
    if (!Vscaled.empty()) {
        auto temp0 = std::make_shared<coal::BVHModel<coal::AABB>>();
        temp0->beginModel();
        if (!F.empty()) temp0->addSubModel(Vscaled, F); else temp0->addSubModel(Vscaled);
        temp0->endModel();
        com = temp0->computeCOM();
    out.volume = (real_t)temp0->computeVolume();
    }
    out.com = com;
    // 3) Center vertices at COM
    std::vector<coal::Vec3s> Vcentered = Vscaled;
    for (auto& v : Vcentered) { v[0] -= com[0]; v[1] -= com[1]; v[2] -= com[2]; }
    // 4) Compute inertia about COM on centered geometry
    cardillo::Matrix33r Icom = cardillo::Matrix33r::Zero();
    if (!Vcentered.empty()) {
        auto temp1 = std::make_shared<coal::BVHModel<coal::AABB>>();
        temp1->beginModel();
        if (!F.empty()) temp1->addSubModel(Vcentered, F); else temp1->addSubModel(Vcentered);
        temp1->endModel();
    const auto Icoal = temp1->computeMomentofInertia();
    Icom << (real_t)Icoal(0,0), (real_t)Icoal(0,1), (real_t)Icoal(0,2),
         (real_t)Icoal(1,0), (real_t)Icoal(1,1), (real_t)Icoal(1,2),
         (real_t)Icoal(2,0), (real_t)Icoal(2,1), (real_t)Icoal(2,2);
    }
    // 5) Principal axes from Icom
    cardillo::Matrix33r Rpa = cardillo::Matrix33r::Identity();
    cardillo::Vector3r Idiag = cardillo::Vector3r::Zero();
    if (Icom.allFinite()) {
        Eigen::SelfAdjointEigenSolver<cardillo::Matrix33r> es(Icom);
        if (es.info() == Eigen::Success) {
            Rpa = es.eigenvectors();
            if (Rpa.determinant() < 0) Rpa.col(0) = -Rpa.col(0);
            Idiag = es.eigenvalues();
        }
    }
    out.Rpa = Rpa;
    out.inertia_diag_unit = Idiag;
    // 6) Rotate centered vertices to principal axes and build normalized BVH
    std::vector<coal::Vec3s> Vnorm; Vnorm.reserve(Vcentered.size());
    for (const auto& v : Vcentered) {
    cardillo::Vector3r vc((real_t)v[0], (real_t)v[1], (real_t)v[2]);
        cardillo::Vector3r vp = Rpa.transpose() * vc;
        Vnorm.emplace_back((coal::CoalScalar)vp.x(), (coal::CoalScalar)vp.y(), (coal::CoalScalar)vp.z());
    }
    auto temp = std::make_shared<coal::BVHModel<coal::AABB>>();
    temp->beginModel();
    if (!F.empty()) temp->addSubModel(Vnorm, F); else temp->addSubModel(Vnorm);
    temp->endModel(); temp->computeLocalAABB();
    out.bvh = temp;
    return out;
}

// Build a BVH scaled by 'scale' without COM-centering or principal-axes normalization
coal::BVHModelPtr_t buildScaledBVH(const std::string& meshPath, const cardillo::Vector3r& scale) {
    coal::MeshLoader loader(coal::BV_AABB);
    coal::BVHModelPtr_t bvh = loader.load(meshPath);
    if (!bvh) return nullptr;
    std::vector<coal::Vec3s> V;
    std::vector<coal::Triangle> F;
    if (bvh->vertices) V = *bvh->vertices;
    if (bvh->tri_indices) F = *bvh->tri_indices;

    // Compute Center point
    real_t min_x = std::numeric_limits<real_t>::max();;
    real_t min_y = std::numeric_limits<real_t>::max();;
    real_t min_z = std::numeric_limits<real_t>::max();;
    real_t max_x = std::numeric_limits<real_t>::lowest();;
    real_t max_y = std::numeric_limits<real_t>::lowest();;
    real_t max_z = std::numeric_limits<real_t>::lowest();;

    for (const auto& v : V) {
        if ((real_t)v[0] < min_x) min_x = (real_t)v[0];
        if ((real_t)v[1] < min_y) min_y = (real_t)v[1];
        if ((real_t)v[2] < min_z) min_z = (real_t)v[2];
        if ((real_t)v[0] > max_x) max_x = (real_t)v[0];
        if ((real_t)v[1] > max_y) max_y = (real_t)v[1];
        if ((real_t)v[2] > max_z) max_z = (real_t)v[2];
    }
    real_t cx = (min_x + max_x) / 2.0;
    real_t cy = (min_y + max_y) / 2.0;
    real_t cz = (min_z + max_z) / 2.0;

    // Scale about COM: p' = com + S*(p - com)
    std::vector<coal::Vec3s> Vscaled; Vscaled.reserve(V.size());
    for (const auto& v : V) {
        coal::CoalScalar dx = (coal::CoalScalar)scale.x() * (v[0] - cx);
        coal::CoalScalar dy = (coal::CoalScalar)scale.y() * (v[1] - cy);
        coal::CoalScalar dz = (coal::CoalScalar)scale.z() * (v[2] - cz);
        Vscaled.emplace_back(cx + dx, cy + dy, cz + dz);
    }
    auto temp = std::make_shared<coal::BVHModel<coal::AABB>>();
    temp->beginModel();
    if (!F.empty()) temp->addSubModel(Vscaled, F); else temp->addSubModel(Vscaled);
    temp->endModel(); temp->computeLocalAABB();
    return temp;
}

// Parse simple OBJ UVs once; returns (uvs, hasUV). Maps the first vt per vertex index.
std::pair<std::vector<Eigen::Vector2f>, bool> parseOBJUVs(const std::string& meshPath, std::size_t vertexCount) {
    std::vector<Eigen::Vector2f> uvs(vertexCount, Eigen::Vector2f(0.f,0.f));
    bool hasUV = false;
    std::ifstream in(meshPath);
    if (!in) return {uvs, hasUV};
    std::vector<Eigen::Vector2f> vt; vt.reserve(1024);
    std::vector<int> vtIndexPerV; vtIndexPerV.assign((int)vertexCount, -1);
    std::string line;
    while (std::getline(in, line)) {
        if (line.size() >= 2 && line[0] == 'v' && line[1] == 't') {
            std::istringstream ss(line.substr(2)); float u=0.f, v=0.f; ss >> u >> v; vt.emplace_back(u,v);
        } else if (!line.empty() && line[0] == 'f') {
            std::istringstream ss(line.substr(1)); std::string tok;
            while (ss >> tok) {
                int vi = -1, vti = -1; size_t p1 = tok.find('/');
                if (p1 == std::string::npos) { vi = std::stoi(tok); }
                else {
                    vi = std::stoi(tok.substr(0, p1));
                    size_t p2 = tok.find('/', p1+1);
                    std::string vts = (p2 == std::string::npos) ? tok.substr(p1+1) : tok.substr(p1+1, p2-(p1+1));
                    if (!vts.empty()) vti = std::stoi(vts);
                }
                if (vi <= 0) continue; int idx0 = vi - 1;
                if (vti > 0 && idx0 >= 0 && idx0 < (int)vtIndexPerV.size()) {
                    if (vtIndexPerV[idx0] == -1) vtIndexPerV[idx0] = vti - 1;
                }
            }
        }
    }
    for (size_t i = 0; i < vtIndexPerV.size(); ++i) {
        int uvi = vtIndexPerV[i];
        if (uvi >= 0 && uvi < (int)vt.size()) { uvs[(int)i] = vt[(size_t)uvi]; hasUV = true; }
    }
    return {uvs, hasUV};
}
} // anonymous namespace

static inline std::string makeMeshKey_(const std::string& path, const Vector3r& s, bool normalized) {
    std::ostringstream ss; ss.setf(std::ios::fixed); ss<<std::setprecision(6);
    ss << path << "|" << (double)s.x() << "," << (double)s.y() << "," << (double)s.z() << "|norm=" << (normalized ? 1 : 0);
    return ss.str();
}

const PhysicsSystem::MeshAsset& PhysicsSystem::getMeshAsset(entt::entity e) const {
    // Extract mesh path/scale from entity
    if (!m_reg.any_of<C_Mesh>(e)) throw std::runtime_error("getMeshAsset(entity): entity has no C_Mesh");
    const auto& cm = m_reg.get<C_Mesh>(e);
    const bool isDynamic = m_reg.any_of<C_PhysicsObject>(e);
    const std::string key = makeMeshKey_(cm.path, cm.scale, /*normalized*/ isDynamic);
    auto it = m_meshCache.find(key);
    if (it != m_meshCache.end()) return it->second;
    // Build new asset via helpers
    MeshAsset asset;
    try {
        if (!isDynamic) {
            // Static mesh: just build scaled BVH
            asset.bvh = buildScaledBVH(cm.path, cm.scale);
            asset.inertia_diag_unit = Vector3r::Zero();
            asset.volume = (real_t)0.0;
        }else {
            // Dynamic mesh: normalize (COM-center + principal axes)
            auto nm = normalizeMesh(cm.path, cm.scale);
            asset.bvh = nm.bvh;
            asset.inertia_diag_unit = nm.inertia_diag_unit;
            asset.volume = nm.volume;
            asset.Rpa = nm.Rpa;
            asset.com = Vector3r((real_t)nm.com[0], (real_t)nm.com[1], (real_t)nm.com[2]);
            asset.normalized = true;
        }

        // Optional debug printout
        if (m_cfg.debug_mesh) {
            std::printf("[Mesh] path=%s dyn=%d volume=%.6g com=(%.6g,%.6g,%.6g) Iunit=(%.6g,%.6g,%.6g)\n",
                        cm.path.c_str(), (int)isDynamic, (double)asset.volume,
                        (double)asset.com.x(), (double)asset.com.y(), (double)asset.com.z(),
                        (double)asset.inertia_diag_unit.x(), (double)asset.inertia_diag_unit.y(), (double)asset.inertia_diag_unit.z());
            if (asset.normalized) {
                std::printf("       Rpa=(%.6g,%.6g,%.6g; %.6g,%.6g,%.6g; %.6g,%.6g,%.6g)\n",
                            (double)asset.Rpa(0,0), (double)asset.Rpa(0,1), (double)asset.Rpa(0,2),
                            (double)asset.Rpa(1,0), (double)asset.Rpa(1,1), (double)asset.Rpa(1,2),
                            (double)asset.Rpa(2,0), (double)asset.Rpa(2,1), (double)asset.Rpa(2,2));
            }
            if (asset.bvh) {
                coal::Vec3s vcom = asset.bvh->computeCOM();
                real_t vol_check = (real_t)asset.bvh->computeVolume();
                auto Icheck_coal = asset.bvh->computeMomentofInertia();
                Matrix33r Icheck;
                Icheck << (real_t)Icheck_coal(0,0), (real_t)Icheck_coal(0,1), (real_t)Icheck_coal(0,2),
                           (real_t)Icheck_coal(1,0), (real_t)Icheck_coal(1,1), (real_t)Icheck_coal(1,2),
                           (real_t)Icheck_coal(2,0), (real_t)Icheck_coal(2,1), (real_t)Icheck_coal(2,2);
                std::printf("       [Check] volume=%.6g com=(%.6g,%.6g,%.6g) I=(%.6g,%.6g,%.6g; %.6g,%.6g,%.6g; %.6g,%.6g,%.6g)\n",
                            (double)vol_check, (double)vcom[0], (double)vcom[1], (double)vcom[2],
                            (double)Icheck(0,0), (double)Icheck(0,1), (double)Icheck(0,2),
                            (double)Icheck(1,0), (double)Icheck(1,1), (double)Icheck(1,2),
                            (double)Icheck(2,0), (double)Icheck(2,1), (double)Icheck(2,2));
            }
        }

        // Parse OBJ UVs once
        if (cm.path.size() >= 4 && (cm.path.rfind(".obj") == cm.path.size() - 4 || cm.path.rfind(".OBJ") == cm.path.size() - 4)) {
            std::size_t vcount = 0;
            if (asset.bvh && asset.bvh->vertices) vcount = asset.bvh->vertices->size();
            auto uvres = parseOBJUVs(cm.path, vcount);
            asset.uvs = std::move(uvres.first);
            asset.hasUV = uvres.second;
        }
    } catch (...) {
        std::printf("Warning: Failed to load mesh asset from path: %s\n", cm.path.c_str());
    }
    auto [it2, ok] = m_meshCache.emplace(key, std::move(asset));
    (void)ok;
    return it2->second;
}

static inline std::string makeHFKey_(const std::string& path, real_t xdim, real_t ydim, real_t zscale, real_t minh) {
    std::ostringstream ss; ss.setf(std::ios::fixed); ss<<std::setprecision(6);
    ss << path << "|xd=" << (double)xdim << ",yd=" << (double)ydim << ",zs=" << (double)zscale << ",min=" << (double)minh;
    return ss.str();
}

const PhysicsSystem::HeightFieldAsset& PhysicsSystem::getHeightFieldAsset(entt::entity e) const {
    if (!m_reg.any_of<C_HeightField>(e)) throw std::runtime_error("getHeightFieldAsset(entity): entity has no C_HeightField");
    const auto& ch = m_reg.get<C_HeightField>(e);
    const std::string key = makeHFKey_(ch.path, ch.x_dim, ch.y_dim, ch.z_scale, ch.min_height);
    auto it = m_hfCache.find(key);
    if (it != m_hfCache.end()) return it->second;
    // Build new asset
    HeightFieldAsset asset; asset.path = ch.path; asset.x_dim = ch.x_dim; asset.y_dim = ch.y_dim; asset.z_scale = ch.z_scale; asset.min_height = ch.min_height;
    try {
        Eigen::Matrix<coal::CoalScalar, Eigen::Dynamic, Eigen::Dynamic> H;
        int rows = 0, cols = 0;
        if (cardillo::io::load_exr_heightmap(ch.path, H, rows, cols)) {
            if (ch.z_scale != (real_t)1.0) {
                H *= (coal::CoalScalar)ch.z_scale;
            }
        } else {
            // Fallback: build a flat 2x2 patch at zero height (will be clamped to min_height)
            rows = cols = 2;
            H = Eigen::Matrix<coal::CoalScalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(rows, cols);
            std::printf("[HeightField] Using flat fallback heightfield (2x2) for path: %s\n", ch.path.c_str());
        }
        // Ensure at least 2x2
        if (rows >= 2 && cols >= 2) {
            auto hf = std::make_shared<coal::HeightField<coal::AABB>>((coal::CoalScalar)ch.x_dim, (coal::CoalScalar)ch.y_dim, H, (coal::CoalScalar)ch.min_height);
            hf->computeLocalAABB();
            asset.hf = hf;
            asset.rows = rows; asset.cols = cols;
        }
    } catch (...) {
        std::printf("Warning: Exception while building heightfield from: %s\n", ch.path.c_str());
    }
    auto [it2, ok] = m_hfCache.emplace(key, std::move(asset)); (void)ok;
    return it2->second;
}

} // namespace cardillo::