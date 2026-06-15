#pragma once

#include <limits>
#include <memory>
#include <cmath>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "../../config/config.hpp"
#include "../assets/body_factory.hpp"
#include "../constraints/constraint_factory.hpp"
#include "../world.hpp"
#include "physics_types.hpp"

// Forward declarations for subsystems that will be owned by PhysicsEngine
namespace cardillo {
namespace collision {
class CollisionCoal;
}
}  // namespace cardillo
namespace cardillo {
namespace misc {
class TimingManager;
}
}  // namespace cardillo
namespace cardillo {
namespace physics {
namespace pipeline {
class PhysicsPipeline;
}
}  // namespace physics
}  // namespace cardillo

namespace cardillo {
namespace physics {

/// High-level facade for building scenes and advancing the simulation.
class PhysicsEngine {
   public:
    /// Construct an empty engine. Call initFromConfig before stepping.
    PhysicsEngine();
    /// Destroy owned subsystems.
    ~PhysicsEngine();
    /// Construct the engine and initialize it from a simulation config.
    explicit PhysicsEngine(const config::Config& cfg);
    /// Rebuild the owned world, collision manager, timings, and pipeline from a config.
    void initFromConfig(const config::Config& cfg);

    /// Create a dynamic rigid body and return its entity.
    inline entt::entity addRigidBody(const RigidShape& shape, const RigidState& state, const RigidProps& props) { return BodyFactory::addRigidBody(*m_world, shape, state, props); }

    /// Create a static rigid body and return its entity.
    inline entt::entity addStaticBody(const RigidShape& shape, const RigidState& state) { return BodyFactory::addStaticBody(*m_world, shape, state); }

    /// Create a point mass with an optional visual radius.
    inline entt::entity addPointMass(real_t mass, const Vector3r& x0, const Vector3r& v0, real_t radius = (real_t)0.05) { return BodyFactory::addPointMass(*m_world, mass, x0, v0, radius); }

    /// Spawn a static height-field obstacle from an EXR file.
    inline entt::entity addObstacleHeightField(const Vector3r& position, const Quaternion4r& orientation, const std::string& exrPath, real_t x_dim, real_t y_dim, real_t z_scale = (real_t)1.0,
                                          real_t min_height = (real_t)0.0) {
        return BodyFactory::addObstacleHeightField(*m_world, position, orientation, exrPath, x_dim, y_dim, z_scale, min_height);
    }

    /// Create a soft body from an OBJ mesh.
    inline std::vector<entt::entity> addSoftBody(const std::string& objPath, real_t stiffness, real_t damping, const Vector3r& position = Vector3r::Zero(),
                                                 const Quaternion4r& orientation = Quaternion4r::Identity(), const Vector3r& linearVelocity = Vector3r::Zero(),
                                                 const Vector3r& angularVelocity = Vector3r::Zero(), real_t totalMass = (real_t)0.0, real_t nodeRadius = (real_t)0.02) {
        return BodyFactory::addSoftBody(*m_world, objPath, stiffness, damping, position, orientation, linearVelocity, angularVelocity, totalMass, nodeRadius, m_collision_mgr.get());
    }

    /// Sample a spline into a beam chain and return the root and tip entities.
    inline std::pair<entt::entity, entt::entity> createBeam(const cardillo::misc::SplinePattern& spline, const BeamCrossSection& section, const BeamSpringParams& springs,
                                                            const RigidState& stateDefaults, const RigidProps& propsDefaults, size_t segments) {
        return BodyFactory::createBeam(*m_world, spline, section, springs, stateDefaults, propsDefaults, segments, m_collision_mgr.get());
    }

    /// Create several beams from a list of spline patterns.
    inline std::pair<entt::entity, entt::entity> createBeams(const std::vector<const cardillo::misc::SplinePattern*>& splines, const BeamCrossSection& section, const BeamSpringParams& springs,
                                                             const RigidState& stateDefaults, const RigidProps& propsDefaults, size_t segments) {
        return BodyFactory::createBeams(*m_world, splines, section, springs, stateDefaults, propsDefaults, segments, m_collision_mgr.get());
    }

    /// Advance the simulation by the engine's configured time step.
    void step();

    /// Add a linear distance constraint between two entities.
    inline size_t addLinearDistanceConstraint(entt::entity a, entt::entity b, const Vector3r& rA_local = Vector3r::Zero(), const Vector3r& rB_local = Vector3r::Zero(),
                                              real_t stiffness = std::numeric_limits<real_t>::infinity(), real_t damping = (real_t)0) {
        return ConstraintFactory::addLinearDistanceConstraint(*m_world, a, b, rA_local, rB_local, stiffness, damping);
    }

    /// Constrain two entities to the same rigid frame.
    inline size_t addRigidConstraint(entt::entity a, entt::entity b = entt::null) { return ConstraintFactory::addRigidConstraint(*m_world, a, b); }

    /// Add a six-DOF translation/rotation joint with configurable stiffness and damping.
    inline size_t addTranslationRotationConstraint(entt::entity a, entt::entity b, const JointFrame& frame, const Vector3r& K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity()),
                                                   const Vector3r& D_trans = Vector3r::Zero(), const Vector3r& K_rot = Vector3r::Zero(), const Vector3r& D_rot = Vector3r::Zero()) {
        return ConstraintFactory::addTranslationRotationConstraint(*m_world, a, b, frame, K_trans, D_trans, K_rot, D_rot);
    }

    /// Constrain only translation while leaving rotation free.
    inline size_t addTranslationalConstraint(entt::entity a, entt::entity b, const JointFrame& frame, const Vector3r& K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity()),
                                             const Vector3r& D_trans = Vector3r::Zero()) {
        return ConstraintFactory::addTranslationalConstraint(*m_world, a, b, frame, K_trans, D_trans);
    }

    /// Constrain only rotation while leaving translation free.
    inline size_t addRotationConstraint(entt::entity a, entt::entity b, const JointFrame& frame, const Vector3r& K_rot = Vector3r::Constant(std::numeric_limits<real_t>::infinity()),
                                        const Vector3r& D_rot = Vector3r::Zero()) {
        return ConstraintFactory::addRotationConstraint(*m_world, a, b, frame, K_rot, D_rot);
    }

    /// Add a hinge joint with a single rotational axis and optional translational stiffness.
    inline size_t addHingeConstraint(entt::entity a, entt::entity b, const JointFrame& frame, real_t K_axis = (real_t)0, real_t D_axis = (real_t)0,
                                     const Vector3r& K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity()), const Vector3r& D_trans = Vector3r::Zero()) {
        return ConstraintFactory::addHingeConstraint(*m_world, a, b, frame, K_axis, D_axis, K_trans, D_trans);
    }

    /// Add a beam constraint between two beam segment entities.
    inline size_t addBeamConstraint(entt::entity a, entt::entity b, const BeamSpringParams& springs, const BeamCrossSection& section) {
        return ConstraintFactory::addBeamConstraint(*m_world, a, b, springs, section);
    }

    /// Set the scalar target velocity for a constraint pattern.
    void setConstraintScalarVelocity(size_t constraintIndex, real_t v);
    /// Set the translational target velocity for a constraint pattern.
    void setConstraintLinearVelocity(size_t constraintIndex, const Vector3r& v);
    /// Set the angular target velocity for a constraint pattern.
    void setConstraintAngularVelocity(size_t constraintIndex, const Vector3r& w);

    /// Expose the underlying ECS registry for queries and custom systems.
    entt::registry& ecs() { return m_world->ecs(); }
    /// Expose the underlying ECS registry for queries and custom systems.
    const entt::registry& ecs() const { return m_world->ecs(); }

    /// Disable collision checks for a specific pair.
    void disableCollisionBetween(entt::entity a, entt::entity b);
    /// Convert a dynamic entity into a static one.
    void makeStatic(entt::entity e) { m_world->makeStatic(e); }

    /// Query the mass matrix for an entity.
    MatrixXXr getMass(entt::entity e) const { return m_world->getMass(e); }
    /// Query the body-frame inertia diagonal for an entity.
    Vector3r getInertiaDiag(entt::entity e) const { return m_world->getInertiaDiag(e); }
    /// Query the generalized position vector for an entity.
    VectorXr getPosition(entt::entity e) const { return m_world->getPosition(e); }
    /// Query the kinetic energy for an entity.
    real_t getKineticEnergy(entt::entity e) const { return m_world->getKineticEnergy(e); }

    /// Apply a world-space force and torque to an entity.
    void applyForce(entt::entity e, const Vector3r& f, const Vector3r& tau) { m_world->applyForce(e, f, tau); }
    /// Apply a pure world-space torque to an entity.
    void applyInertialTorque(entt::entity e, const Vector3r& tau) { m_world->applyInertialTorque(e, tau); }

    /// Set the position of an entity.
    void setPosition(entt::entity e, const Vector3r& p) { m_world->setPosition(e, p); }
    /// Set the orientation of an entity.
    void setOrientation(entt::entity e, const Quaternion4r& q) { m_world->setOrientation(e, q); }
    /// Set the linear velocity of an entity.
    void setLinearVelocity(entt::entity e, const Vector3r& v) { m_world->setLinearVelocity(e, v); }
    /// Set the angular velocity of an entity.
    void setAngularVelocity(entt::entity e, const Vector3r& w) { m_world->setAngularVelocity(e, w); }
    /// Set velocity directly by force-like state input.
    void setVelocityByForce(entt::entity e, const Vector3r& v, const Vector3r& w) { m_world->setVelocityByForce(e, v, w); }
    /// Attach a time-dependent trajectory callback to an entity.
    void addTrajectory(entt::entity e, std::optional<std::function<TrajectoryPose(real_t)>> positionFunc, std::optional<std::function<TrajectoryTwist(real_t)>> velocityFunc) {
        m_world->setTrajectory(e, std::move(positionFunc), std::move(velocityFunc));
    }
    /// Attach a periodic spline trajectory to an entity.
    template <class TSpline, std::enable_if_t<std::is_base_of_v<cardillo::misc::SplinePattern, TSpline>, int> = 0>
    void addTrajectory(entt::entity e, const TSpline& spline, real_t period) {
        auto splineCopy = std::make_shared<TSpline>(spline);
        const real_t safePeriod = std::max(period, (real_t)1e-8);
        m_world->setTrajectory(
            e,
            [splineCopy, safePeriod](real_t t) -> TrajectoryPose {
                const real_t phase = std::fmod(std::max(t, (real_t)0), safePeriod) / safePeriod;
                const auto sample = splineCopy->sample(phase);
                return {sample.position, Quaternion4r::Identity()};
            },
            std::nullopt);
    }
            /// Remove an entity trajectory.
    void removeTrajectory(entt::entity e) { m_world->removeTrajectory(e); }

            /// Set the world gravity vector.
    void setGravity(const Vector3r& g) { m_world->setGravity(g); }
            /// Get the world gravity vector.
    const Vector3r& gravity() const { return m_world->gravity(); }

            /// Mark the world structure dirty so the next step rebuilds derived state.
    void markStructureDirty() { m_world->markStructureDirty(); }
            /// Attach a tracking label to an entity for output.
    void track(entt::entity e, const std::string& name) { m_world->track(e, name); }

            /// Advance the simulation pipeline by dt seconds.
    void step(real_t dt);

            /// Access the collision manager owned by the engine.
    collision::CollisionCoal& collisionManager();
            /// Access the timing manager owned by the engine.
    cardillo::misc::TimingManager& timings();
            /// Query whether the pipeline reached the configured simulation end time.
    bool isFinished() const;

   private:
    std::unique_ptr<cardillo::World> m_world;
    // Owned subsystems
    config::Config m_cfg{};                                               // global config
    std::unique_ptr<cardillo::collision::CollisionCoal> m_collision_mgr;  // engine-owned
    std::unique_ptr<cardillo::misc::TimingManager> m_timings;             // engine-owned
    // Pipeline (orchestrates collision, assembly, solver, integrator)
    std::unique_ptr<cardillo::physics::pipeline::PhysicsPipeline> m_pipeline;
};

}  // namespace physics
}  // namespace cardillo
