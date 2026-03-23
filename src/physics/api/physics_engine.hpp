#pragma once

#include <type_traits>
#include <utility>
#include <string>
#include <vector>
#include <limits>

#include "physics_types.hpp"

namespace cardillo { class World; }
namespace cardillo { namespace misc { struct SplinePattern; } }
namespace cardillo { namespace physics { struct JointFrame; } }

namespace cardillo {
namespace physics {

class PhysicsEngine {
public:
    explicit PhysicsEngine(cardillo::World& world) : m_world(world) {}

    // Public API
    cardillo::World& world() { return m_world; }
    const cardillo::World& world() const { return m_world; }

    cardillo::World& system() { return m_world; }
    const cardillo::World& system() const { return m_world; }

    entt::entity addRigidBody(const RigidShape& shape,
                              const RigidState& state,
                              const RigidProps& props);

    entt::entity addStaticBody(const RigidShape& shape,
                               const RigidState& state);

    index_t addPointMass(real_t mass,
                         const Vector3r& x0,
                         const Vector3r& v0,
                         real_t radius = (real_t)0.05);

    index_t addObstacleHeightField(const Vector3r& position,
                                   const Quaternion4r& orientation,
                                   const std::string& exrPath,
                                   real_t x_dim,
                                   real_t y_dim,
                                   real_t z_scale = (real_t)1.0,
                                   real_t min_height = (real_t)0.0);

    std::vector<entt::entity> addSoftBody(const std::string& objPath,
                                          real_t stiffness,
                                          real_t damping,
                                          const Vector3r& position = Vector3r::Zero(),
                                          const Quaternion4r& orientation = Quaternion4r::Identity(),
                                          const Vector3r& linearVelocity = Vector3r::Zero(),
                                          const Vector3r& angularVelocity = Vector3r::Zero(),
                                          real_t totalMass = (real_t)0.0);

    std::pair<entt::entity, entt::entity> createBeam(const cardillo::misc::SplinePattern& spline,
                                                      const BeamCrossSection& section,
                                                      const BeamSpringParams& springs,
                                                      const RigidState& stateDefaults,
                                                      const RigidProps& propsDefaults,
                                                      size_t segments);

    std::pair<entt::entity, entt::entity> createBeams(const std::vector<const cardillo::misc::SplinePattern*>& splines,
                                                       const BeamCrossSection& section,
                                                       const BeamSpringParams& springs,
                                                       const RigidState& stateDefaults,
                                                       const RigidProps& propsDefaults,
                                                       size_t segmentsPerSpline);

    size_t addLinearDistanceConstraint(entt::entity a,
                                       entt::entity b,
                                       const Vector3r& rA_local = Vector3r::Zero(),
                                       const Vector3r& rB_local = Vector3r::Zero(),
                                       real_t stiffness = std::numeric_limits<real_t>::infinity(),
                                       real_t damping = (real_t)0);

    size_t addRigidConstraint(entt::entity a, entt::entity b);

    size_t addTranslationRotationConstraint(entt::entity a,
                                            entt::entity b,
                                            const JointFrame& frame,
                                            const Vector3r& K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity()),
                                            const Vector3r& D_trans = Vector3r::Zero(),
                                            const Vector3r& K_rot = Vector3r::Zero(),
                                            const Vector3r& D_rot = Vector3r::Zero());

    size_t addTranslationalConstraint(entt::entity a,
                                      entt::entity b,
                                      const JointFrame& frame,
                                      const Vector3r& K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity()),
                                      const Vector3r& D_trans = Vector3r::Zero());

    size_t addRotationConstraint(entt::entity a,
                                 entt::entity b,
                                 const JointFrame& frame,
                                 const Vector3r& K_rot = Vector3r::Constant(std::numeric_limits<real_t>::infinity()),
                                 const Vector3r& D_rot = Vector3r::Zero());

    size_t addHingeConstraint(entt::entity a,
                              entt::entity b,
                              const JointFrame& frame,
                              real_t K_axis = (real_t)0,
                              real_t D_axis = (real_t)0,
                              const Vector3r& K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity()),
                              const Vector3r& D_trans = Vector3r::Zero());

    size_t addBeamConstraint(entt::entity a,
                             entt::entity b,
                             const BeamSpringParams& springs,
                             const BeamCrossSection& section);

    // Compatibility forwarders (public API)
    entt::registry& ecs();
    const entt::registry& ecs() const;

    void disableCollisionBetween(entt::entity a, entt::entity b);
    void makeStatic(entt::entity e);

    MatrixXXr getMass(entt::entity e) const;
    Vector3r getInertiaDiag(entt::entity e) const;
    VectorXr getPosition(entt::entity e) const;
    real_t getKineticEnergy(entt::entity e) const;

    void applyForce(entt::entity e, const Vector3r& f, const Vector3r& tau);
    void applyInertialTorque(entt::entity e, const Vector3r& tau);

    void setPosition(entt::entity e, const Vector3r& p);
    void setOrientation(entt::entity e, const Quaternion4r& q);
    void setLinearVelocity(entt::entity e, const Vector3r& v);
    void setAngularVelocity(entt::entity e, const Vector3r& w);
    void setVelocityByForce(entt::entity e, const Vector3r& v, const Vector3r& w);

    void setGravity(const Vector3r& g);
    const Vector3r& gravity() const;

    void markStructureDirty();
    void track(entt::entity e, const std::string& name);
    void writeTrackedStateToCsv(real_t t);

private:
    cardillo::World& m_world;
};

} // namespace physics
} // namespace cardillo
