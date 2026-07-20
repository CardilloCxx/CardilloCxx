#pragma once

#include <memory>

#include <entt/entt.hpp>

#include "../world.hpp"
#include "constraints.hpp"

namespace cardillo {
namespace physics {

class ConstraintFactory {
   public:
    static size_t addLinearDistanceConstraint(World& world, entt::entity a, entt::entity b, const Vector3r& rA_local = Vector3r::Zero(), const Vector3r& rB_local = Vector3r::Zero(),
                                              real_t stiffness = std::numeric_limits<real_t>::infinity(), real_t damping = (real_t)0, real_t targetDistance = (real_t)-1.0) {
        return insertPattern(world, std::unique_ptr<ConstraintPattern>(new LinearDistanceConstraint(world.ecs(), a, b, rA_local, rB_local, stiffness, damping, targetDistance)));
    }

    static size_t addRigidConstraint(World& world, entt::entity a, entt::entity b) {
        const auto aState = RigidBody::getState(world.ecs(), a);
        const JointFrame frame(aState.position, Matrix33r::Identity());
        const Vector3r K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity());
        const Vector3r K_rot = Vector3r::Constant(std::numeric_limits<real_t>::infinity());
        return insertPattern(world, std::unique_ptr<ConstraintPattern>(new TranslationRotationConstraint(world.ecs(), b, a, frame, K_trans, Vector3r::Zero(), K_rot, Vector3r::Zero())));
    }

    static size_t addTranslationRotationConstraint(World& world, entt::entity a, entt::entity b, const JointFrame& frame,
                                                   const Vector3r& K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity()), const Vector3r& D_trans = Vector3r::Zero(),
                                                   const Vector3r& K_rot = Vector3r::Zero(), const Vector3r& D_rot = Vector3r::Zero()) {
        return insertPattern(world, std::unique_ptr<ConstraintPattern>(new TranslationRotationConstraint(world.ecs(), a, b, frame, K_trans, D_trans, K_rot, D_rot)));
    }

    static size_t addTranslationalConstraint(World& world, entt::entity a, entt::entity b, const JointFrame& frame,
                                             const Vector3r& K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity()), const Vector3r& D_trans = Vector3r::Zero()) {
        const Vector3r K_rot = Vector3r::Zero();
        const Vector3r D_rot = Vector3r::Zero();
        return insertPattern(world, std::unique_ptr<ConstraintPattern>(new TranslationRotationConstraint(world.ecs(), a, b, frame, K_trans, D_trans, K_rot, D_rot)));
    }

    static size_t addRotationConstraint(World& world, entt::entity a, entt::entity b, const JointFrame& frame, const Vector3r& K_rot = Vector3r::Constant(std::numeric_limits<real_t>::infinity()),
                                        const Vector3r& D_rot = Vector3r::Zero()) {
        const Vector3r K_trans = Vector3r::Zero();
        const Vector3r D_trans = Vector3r::Zero();
        return insertPattern(world, std::unique_ptr<ConstraintPattern>(new TranslationRotationConstraint(world.ecs(), a, b, frame, K_trans, D_trans, K_rot, D_rot)));
    }

    static size_t addHingeConstraint(World& world, entt::entity a, entt::entity b, const JointFrame& frame, real_t K_axis = (real_t)0, real_t D_axis = (real_t)0,
                                     const Vector3r& K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity()), const Vector3r& D_trans = Vector3r::Zero()) {
        const Vector3r K_rot(K_axis, std::numeric_limits<real_t>::infinity(), std::numeric_limits<real_t>::infinity());
        const Vector3r D_rot(D_axis, (real_t)0, (real_t)0);

        return insertPattern(world, std::unique_ptr<ConstraintPattern>(new TranslationRotationConstraint(world.ecs(), a, b, frame, K_trans, D_trans, K_rot, D_rot)));
    }

    static size_t addBeamConstraint(World& world, entt::entity a, entt::entity b, const physics::BeamSpringParams& springs, const physics::BeamCrossSection& section) {
        return insertPattern(world, std::unique_ptr<ConstraintPattern>(new BeamConstraint(world.ecs(), a, b, springs, section)));
    }

    static size_t addConstraint(World& world, std::unique_ptr<ConstraintPattern> pattern) {
        return insertPattern(world, std::move(pattern));
    }

   private:
    static size_t insertPattern(World& world, std::unique_ptr<ConstraintPattern> pattern) {
        auto& patterns = world.constraintPatterns();
        patterns.push_back(std::move(pattern));
        world.markStructureDirty();
        return patterns.size() - 1;
    }
};

}  // namespace physics
}  // namespace cardillo
