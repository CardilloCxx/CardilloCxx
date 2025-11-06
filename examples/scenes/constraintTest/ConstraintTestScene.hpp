#pragma once

#include "../SceneBase.hpp"
#include "physics/constraints.hpp"
#include <Eigen/Geometry>

using namespace cardillo;

// Simple constraint test: two cubes connected at their corners via translational constraint in A's frame.
class ConstraintTestScene : public SceneBase {
public:
    ConstraintTestScene() = default;
    ~ConstraintTestScene() = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;
        // Gravity off (also set in config; duplicate here for clarity if other configs are used)
        sys.setGravity(Vector3r::Zero());

        // Floor (make static to avoid zero-mass dynamic DOFs)
        {
            PhysicsSystem::Cube shape;
            shape.halfExtents = Vector3r(15.0, 15.0, 0.1);
            real_t mass = 0.0; // static
            Vector3r pos = Vector3r(0, 0, -0.1);
            Quaternion4r q = Quaternion4r::Identity();
            Vector3r v = Vector3r::Zero();
            Vector3r w = Vector3r::Zero();
            entt::entity floor = sys.addRigidBody(mass, pos, q, v, w, shape);
            sys.makeStatic(floor);
        }

        // Translational constraint test center
        {
            // Cube shape parameters
            PhysicsSystem::Cube shape;
            shape.halfExtents = Vector3r(0.25, 0.25, 0.25); // 0.5m cube

            real_t mass = 1.0;

            // Body A
            Vector3r posA = Vector3r(3, 0, 3);
            Quaternion4r qA = Quaternion4r::Identity();
            Vector3r vA = Vector3r(1.0, 0.0, 0.0);
            Vector3r wA = Vector3r(0.0, 5.0, 0.0); // initial spin around z
            entt::entity A = sys.addRigidBody(mass, posA, qA, vA, wA, shape);

            // Body B
            const Vector3r he = shape.halfExtents;
            Vector3r posB = Vector3r(3 + 2.0 * he.x(), 2.0 * he.y(), 3 + 2.0 * he.z()); // offset so corners meet
            Quaternion4r qB = Quaternion4r::Identity();
            Vector3r vB = Vector3r(0.0, 0.0, 0.0); // initial velocity towards A
            Vector3r wB = Vector3r::Zero();
            entt::entity B = sys.addRigidBody(mass, posB, qB, vB, wB, shape);

            sys.addConstraint<physics::TranslationalConstraint>(sys.ecs(), A, B);
        }

        // Translational constraint test
        {
            // Cube shape parameters
            PhysicsSystem::Cube shape;
            shape.halfExtents = Vector3r(0.25, 0.25, 0.25); // 0.5m cube

            real_t mass = 1.0;

            // Body A
            Vector3r posA = Vector3r(0, 0, 3);
            Quaternion4r qA = Quaternion4r::Identity();
            Vector3r vA = Vector3r::Zero();
            Vector3r wA = Vector3r(0.0, 5.0, 0.0); // initial spin around z
            entt::entity A = sys.addRigidBody(mass, posA, qA, vA, wA, shape);

            // Body B
            const Vector3r he = shape.halfExtents;
            Vector3r posB = Vector3r(2.0 * he.x(), 2.0 * he.y(), 3 + 2.0 * he.z()); // offset so corners meet
            Quaternion4r qB = Quaternion4r::Identity();
            Vector3r vB = Vector3r(0.0, 0.0, 0.0); // initial velocity towards A
            Vector3r wB = Vector3r::Zero();
            entt::entity B = sys.addRigidBody(mass, posB, qB, vB, wB, shape);

            // Attachment points:
            Vector3r rA_local = Vector3r( he.x(), he.y(), he.z());
            Vector3r rB_local = Vector3r( -he.x(), -he.y(), -he.z());
            sys.addConstraint<physics::TranslationalConstraint>(sys.ecs(), A, B, rA_local, rB_local);
        }

        // Double pendulum
        {
            PhysicsSystem::Cube shape;
            shape.halfExtents = Vector3r(0.1, 0.1, 0.1); // thin vertical rectangle
            real_t mass = 1.0;

            // Body A
            Vector3r posA = Vector3r(-3.0, 0.0, 3.0);
            Quaternion4r qA = Quaternion4r::Identity();
            Vector3r vA = Vector3r::Zero();
            Vector3r wA = Vector3r::Zero();
            entt::entity A = sys.addRigidBody(mass, posA, qA, vA, wA, shape);
            sys.makeStatic(A); // fix first body

            // Body B
            Vector3r posB = Vector3r(-3.0, 0.0, 2.0);
            Quaternion4r qB = Quaternion4r::Identity();
            Vector3r vB = Vector3r(3.0, 2.0, 1.0);
            Vector3r wB = Vector3r::Zero();
            entt::entity B = sys.addRigidBody(mass, posB, qB, vB, wB, shape);

            // Body C
            Vector3r posC = Vector3r(-3.0, 0.0, 0.5);
            Quaternion4r qC = Quaternion4r::Identity();
            Vector3r vC = Vector3r(12.0, 2.0, 1.0);
            Vector3r wC = Vector3r::Zero();
            entt::entity C = sys.addRigidBody(mass, posC, qC, vC, wC, shape);

            // Constraints
            sys.addConstraint<physics::LinearDistanceConstraint>(sys.ecs(), A, B, Vector3r(0.0, 0.0, -0.1), Vector3r(0.0, 0.0, 0.1));
            sys.addConstraint<physics::LinearDistanceConstraint>(sys.ecs(), B, C, Vector3r(0.0, 0.0, -0.1), Vector3r(0.0, 0.0, 0.1));
        }
    }

    void updateScene(cardillo::PhysicsSystem& /*sys*/, real_t /*t*/, real_t /*dt*/) override {

    }
};
