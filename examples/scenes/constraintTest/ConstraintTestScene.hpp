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
    PhysicsSystem::CubeShape floorShape{Vector3r(15.0, 15.0, 0.1)};
    PhysicsSystem::RigidState floorState; floorState.position = Vector3r(0,0,-0.1); floorState.orientation = Quaternion4r::Identity();
    PhysicsSystem::RigidProps floorProps; (void)sys.addRigidBody(floorShape, floorState, floorProps);

//         // Hinge constraint test
        {
            // Tilted obstacle cube
            PhysicsSystem::CubeShape obstacle_shape{Vector3r(3.0, 3.0, 3.0)};
            Quaternion4r q_obstacle = Quaternion4r(Eigen::AngleAxis<real_t>(M_PI / 6.0, Vector3r::UnitY()));
            Vector3r pos_obstacle = Vector3r(0, 6, 2);
            PhysicsSystem::RigidState st_ob; st_ob.position = pos_obstacle; st_ob.orientation = q_obstacle; PhysicsSystem::RigidProps pr_ob; entt::entity obstacle = sys.addRigidBody(obstacle_shape, st_ob, pr_ob);

            Matrix33r rot_obstacle = q_obstacle.toRotationMatrix();

            // Disk to be hinged
            PhysicsSystem::CubeShape diskShape{Vector3r(1.5, 1.5, 0.05)};
            real_t mass = 10.0;
            Vector3r pos = pos_obstacle + rot_obstacle * Vector3r(0.0, 0.0, 3.5);
            PhysicsSystem::RigidState st; st.position = pos; st.orientation = q_obstacle; PhysicsSystem::RigidProps pr; pr.mass = mass; entt::entity disk = sys.addRigidBody(diskShape, st, pr);

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
            PhysicsSystem::CubeShape cubeShape{Vector3r(0.1, 0.1, 0.1)};
            mass = 1.0;
            Vector3r pos_cube = pos_obstacle + rot_obstacle * Vector3r(1.0, 1.0, 3.75);
            PhysicsSystem::RigidState stc; stc.position = pos_cube; stc.orientation = q_obstacle; PhysicsSystem::RigidProps prc; prc.mass = mass; entt::entity cube = sys.addRigidBody(cubeShape, stc, prc);

            sys.addConstraint<physics::TranslationalConstraint>(sys.ecs(), disk, cube, Vector3r(1, 1, 0.0), Vector3r(0.0, 0.0, 0.0), 
                                                            Vector3r(10, 10, 1e20), Vector3r::Zero());
            // sys.addConstraint<physics::RigidConstraint>(sys.ecs(), disk, cube);

        }

        // Translational constraint test center
        {
            // Cube shape parameters
            PhysicsSystem::CubeShape shape{Vector3r(0.25, 0.25, 0.25)}; // 0.5m cube

            real_t mass = 1.0;

            // Body A
            Vector3r posA = Vector3r(3, 0, 3);
            Quaternion4r qA = Quaternion4r::Identity();
            Vector3r vA = Vector3r(1.0, 0.0, 0.0);
            Vector3r wA = Vector3r(0.0, 5.0, 0.0); // initial spin around z
            entt::entity A; {
                PhysicsSystem::RigidState stA; stA.position = posA; stA.orientation = qA; stA.linearVelocity = vA; stA.angularVelocity = wA; PhysicsSystem::RigidProps prA; prA.mass = mass; A = sys.addRigidBody(shape, stA, prA);
            }

            // Body B
            const Vector3r he = shape.halfExtents;
            Vector3r posB = Vector3r(3 + 2.0 * he.x(), 2.0 * he.y(), 3 + 2.0 * he.z()); // offset so corners meet
            Quaternion4r qB = Quaternion4r::Identity();
            Vector3r vB = Vector3r(0.0, 0.0, 0.0); // initial velocity towards A
            Vector3r wB = Vector3r::Zero();
            entt::entity B; {
                PhysicsSystem::RigidState stB; stB.position = posB; stB.orientation = qB; stB.linearVelocity = vB; stB.angularVelocity = wB; PhysicsSystem::RigidProps prB; prB.mass = mass; B = sys.addRigidBody(shape, stB, prB);
            }

            sys.addConstraint<physics::TranslationalConstraint>(sys.ecs(), A, B);
        }

        // Translational constraint test
        {
            // Cube shape parameters
            PhysicsSystem::CubeShape shape{Vector3r(0.25, 0.25, 0.25)}; // 0.5m cube

            real_t mass = 1.0;

            // Body A
            Vector3r posA = Vector3r(0, 0, 3);
            Quaternion4r qA = Quaternion4r::Identity();
            Vector3r vA = Vector3r::Zero();
            Vector3r wA = Vector3r(0.0, 5.0, 0.0); // initial spin around z
            entt::entity A; {
                PhysicsSystem::RigidState stA; stA.position = posA; stA.orientation = qA; stA.linearVelocity = vA; stA.angularVelocity = wA; PhysicsSystem::RigidProps prA; prA.mass = mass; A = sys.addRigidBody(shape, stA, prA);
            }

            // Body B
            const Vector3r he = shape.halfExtents;
            Vector3r posB = Vector3r(2.0 * he.x(), 2.0 * he.y(), 3 + 2.0 * he.z()); // offset so corners meet
            Quaternion4r qB = Quaternion4r::Identity();
            Vector3r vB = Vector3r(0.0, 0.0, 0.0); // initial velocity towards A
            Vector3r wB = Vector3r::Zero();
            entt::entity B; {
                PhysicsSystem::RigidState stB; stB.position = posB; stB.orientation = qB; stB.linearVelocity = vB; stB.angularVelocity = wB; PhysicsSystem::RigidProps prB; prB.mass = mass; B = sys.addRigidBody(shape, stB, prB);
            }

            // Attachment points:
            Vector3r rA_local = Vector3r( he.x(), he.y(), he.z());
            Vector3r rB_local = Vector3r( -he.x(), -he.y(), -he.z());
            sys.addConstraint<physics::TranslationalConstraint>(sys.ecs(), A, B, rA_local, rB_local);
        }

        // Double pendulum
        {
            PhysicsSystem::CubeShape shape{Vector3r(0.1, 0.1, 0.1)}; // thin vertical rectangle
            real_t mass = 1.0;

            // Body A
            Vector3r posA = Vector3r(-3.0, 0.0, 3.0);
            Quaternion4r qA = Quaternion4r::Identity();
            Vector3r vA = Vector3r::Zero();
            Vector3r wA = Vector3r::Zero();
            entt::entity A; {
                PhysicsSystem::RigidState stA; stA.position = posA; stA.orientation = qA; stA.linearVelocity = vA; stA.angularVelocity = wA; PhysicsSystem::RigidProps prA; prA.mass = mass; A = sys.addRigidBody(shape, stA, prA); sys.makeStatic(A);
            }

            // Body B
            Vector3r posB = Vector3r(-3.0, 0.0, 2.0);
            Quaternion4r qB = Quaternion4r::Identity();
            Vector3r vB = Vector3r(3.0, 2.0, 1.0);
            Vector3r wB = Vector3r::Zero();
            entt::entity B; {
                PhysicsSystem::RigidState stB; stB.position = posB; stB.orientation = qB; stB.linearVelocity = vB; stB.angularVelocity = wB; PhysicsSystem::RigidProps prB; prB.mass = mass; B = sys.addRigidBody(shape, stB, prB);
            }

            // Body C
            Vector3r posC = Vector3r(-3.0, 0.0, 0.5);
            Quaternion4r qC = Quaternion4r::Identity();
            Vector3r vC = Vector3r(12.0, 2.0, 1.0);
            Vector3r wC = Vector3r::Zero();
            entt::entity C; {
                PhysicsSystem::RigidState stC; stC.position = posC; stC.orientation = qC; stC.linearVelocity = vC; stC.angularVelocity = wC; PhysicsSystem::RigidProps prC; prC.mass = mass; C = sys.addRigidBody(shape, stC, prC);
            }

            // Constraints
            sys.addConstraint<physics::LinearDistanceConstraint>(sys.ecs(), A, B, Vector3r(0.0, 0.0, -0.1), Vector3r(0.0, 0.0, 0.1));
            sys.addConstraint<physics::LinearDistanceConstraint>(sys.ecs(), B, C, Vector3r(0.0, 0.0, -0.1), Vector3r(0.0, 0.0, 0.1));
        }
    }

    void updateScene(cardillo::PhysicsSystem& /*sys*/, real_t /*t*/, real_t /*dt*/) override {

    }
};
