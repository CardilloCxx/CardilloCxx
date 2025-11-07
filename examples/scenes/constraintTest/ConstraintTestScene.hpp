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
        // sys.setGravity(Vector3r::Zero());

        // Floor (make static to avoid zero-mass dynamic DOFs)
        PhysicsSystem::Cube shape;
        shape.halfExtents = Vector3r(15.0, 15.0, 0.1);
        real_t mass = 0.0; // static
        Vector3r pos = Vector3r(0, 0, -0.1);
        Quaternion4r q = Quaternion4r::Identity();
        Vector3r v = Vector3r::Zero();
        Vector3r w = Vector3r::Zero();
        entt::entity floor = sys.addRigidBody(mass, pos, q, v, w, shape);
        sys.makeStatic(floor);

//         // Hinge constraint test
        {
            // Tilted obstacle cube
            PhysicsSystem::Cube obstacle_shape;
            obstacle_shape.halfExtents = Vector3r(3.0, 3.0, 3.0);
            Quaternion4r q_obstacle = Quaternion4r(Eigen::AngleAxis<real_t>(M_PI / 6.0, Vector3r::UnitY()));
            Vector3r pos_obstacle = Vector3r(0, 6, 2);
            entt::entity obstacle = sys.addRigidBody(0.0, pos_obstacle, q_obstacle, Vector3r::Zero(), Vector3r::Zero(), obstacle_shape);
            sys.makeStatic(obstacle);

            Matrix33r rot_obstacle = q_obstacle.toRotationMatrix();

            // Disk to be hinged
            PhysicsSystem::Cube shape;
            shape.halfExtents = Vector3r(1.5, 1.5, 0.05); 
            real_t mass = 10.0;
            Vector3r pos = pos_obstacle + rot_obstacle * Vector3r(0.0, 0.0, 3.5);
            Vector3r v = Vector3r::Zero();
            Vector3r w = Vector3r(0.0, 0.0, 0.0);
            entt::entity disk = sys.addRigidBody(mass, pos, q_obstacle, v, w, shape);

            Vector3r a_local = Vector3r::Zero(); // hinge axis in A's local frame
            Vector3r b_local = Vector3r(0, 0, 0); // hinge axis in B's local frame
            Vector3r axis_A_frame = Vector3r::UnitZ(); 
            real_t K_hinge = 0;
            real_t D_hinge = 0;
            Vector2r K_rotation = Vector2r::Constant(std::numeric_limits<real_t>::infinity());
            Vector2r D_rotation = Vector2r::Zero();
            Vector3r K_translation = Vector3r::Constant(std::numeric_limits<real_t>::infinity());
            Vector3r D_translation = Vector3r::Zero();

            sys.addConstraint<physics::HingeConstraint>(sys.ecs(), obstacle, disk,
                                                        a_local, b_local,
                                                        axis_A_frame,
                                                        K_hinge, D_hinge,
                                                        K_rotation, D_rotation,
                                                        K_translation, D_translation);


            // Cube over the hinge
            shape.halfExtents = Vector3r(0.1, 0.1, 0.1);
            mass = 1.0;
            Vector3r pos_cube = pos_obstacle + rot_obstacle * Vector3r(1.0, 1.0, 3.61);
            Vector3r v_cube = Vector3r(0,0,0);
            Vector3r w_cube = Vector3r::Zero();
            entt::entity cube = sys.addRigidBody(mass, pos_cube, q_obstacle, v_cube, w_cube, shape);

            // sys.addConstraint<physics::TranslationalConstraint>(sys.ecs(), disk, cube, Vector3r(0.5, 0.5, 0.0), Vector3r(0.0, 0.0, 0.0), 
            //                                                 Vector3r(1e10, 1e10, 100), Vector3r::Zero());
            sys.addConstraint<physics::RigidConstraint>(sys.ecs(), disk, cube);

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
