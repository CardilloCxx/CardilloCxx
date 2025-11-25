#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>

using namespace cardillo;

// Parcel scene.
class DiscreteRodScene : public SceneBase {
public:
    const char* sceneName() const override { return "discrete_rod"; }
    DiscreteRodScene() = default;
    ~DiscreteRodScene() = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;
        using namespace cardillo::misc;

        // Floor (static)
        sys.addStaticBody(PhysicsSystem::CubeShape(Vector3r(15.0, 15.0, 0.1)), PhysicsSystem::RigidState(Vector3r(0,0,-0.1)));

        const size_t segments = 50;
        const real_t length = M_PI;
        const real_t width = 0.05 * length;
        const real_t height = width / 2;
        const real_t density = 600;
        real_t E = 5e7;
        real_t nu = 0.3;

        // // Circle spline with circumference=length
        // const real_t radius = length / ((real_t)2 * M_PI);
        // misc::CircleSpline spline(Vector3r(0,0,4), radius, Vector3r::UnitX(), Vector3r::UnitZ());

        // PhysicsSystem::BeamCrossSection sec_circle(height, width); // width/height naming preserved
        // auto springs_circle = PhysicsSystem::BeamSpringParams::fromMaterial(E, nu);
        // (void) sys.createBeam(spline, sec_circle, springs_circle, PhysicsSystem::RigidState{}, PhysicsSystem::RigidProps::withDensity(density), segments);

        // Linear Beam
        misc::LinearSpline line(Vector3r(3,-1,4), Vector3r(3,1,4));

        PhysicsSystem::BeamCrossSection sec_line(width, height);
        auto springs_line = PhysicsSystem::BeamSpringParams::fromMaterial(E, nu);
        auto endpoints = sys.createBeam(line, sec_line, springs_line, PhysicsSystem::RigidState{}, PhysicsSystem::RigidProps::withDensity(density), segments);
        sys.makeStatic(endpoints.first);
        m_rodEnd = endpoints.second;

        // // Helix Beam
        // misc::HelixSpline helix(Vector3r(-3,0,4), Vector3r::UnitZ(), /* radius */ 0.5, /* pitch */ 1.0, /*turns*/ 2.0);
        // PhysicsSystem::BeamCrossSection sec_helix(width, height);
        // auto springs_helix = PhysicsSystem::BeamSpringParams::fromMaterial(E * 10.0, nu);
        // auto helix_endpoints = sys.createBeam(helix, sec_helix, springs_helix, PhysicsSystem::RigidState{}, PhysicsSystem::RigidProps::withDensity(density), segments);
        // sys.makeStatic(helix_endpoints.second);

        // Beam rope between two static cubes
        auto startCube = sys.addStaticBody(PhysicsSystem::CubeShape(Vector3r(0.1,0.1,0.1)), PhysicsSystem::RigidState(Vector3r(6,-0.1,4)));
        auto endCube = sys.addRigidBody(PhysicsSystem::CubeShape(Vector3r(0.1,0.1,0.1)), PhysicsSystem::RigidState(Vector3r(6,2.1,4)), PhysicsSystem::RigidProps(1e10));
        m_endCube = endCube;

        misc::LinearSpline rope1(Vector3r(6.05, 0, 4), Vector3r(6.05, 2, 4));
        {
            PhysicsSystem::BeamCrossSection sec_rope(0.03, 0.03);
            auto springs_rope = PhysicsSystem::BeamSpringParams::fromMaterial(100, 0.3,
                                                                              1e7, /*shear*/ 1e7, /*torsion*/ 1000, /*bendY*/ 0, /*bendZ*/ 0,
                                                                              Vector3r::Constant(1), Vector3r::Constant(1));
            endpoints = sys.createBeam(rope1, sec_rope, springs_rope, PhysicsSystem::RigidState{}, PhysicsSystem::RigidProps::withDensity(density), 350);
        }
        sys.addConstraint<physics::RigidConstraint>(sys.ecs(), endpoints.first, startCube);
        sys.disableCollisionBetween(endpoints.first, startCube);
        sys.addConstraint<physics::RigidConstraint>(sys.ecs(), endCube, endpoints.second);
        sys.disableCollisionBetween(endpoints.second, endCube);

        misc::LinearSpline rope2(Vector3r(5.95, 0, 4), Vector3r(5.95, 2, 4));
        {
            PhysicsSystem::BeamCrossSection sec_rope2(0.03, 0.03);
            auto springs_rope2 = PhysicsSystem::BeamSpringParams::fromMaterial(100, 0.3,
                                                                               1e7, /*shear*/ 1e7, /*torsion*/ 1000, /*bendY*/ 0, /*bendZ*/ 0,
                                                                               Vector3r::Constant(1), Vector3r::Constant(1));
            endpoints = sys.createBeam(rope2, sec_rope2, springs_rope2, PhysicsSystem::RigidState{}, PhysicsSystem::RigidProps::withDensity(density), 350);
        }
        sys.addConstraint<physics::RigidConstraint>(sys.ecs(), endpoints.first, startCube);
        sys.disableCollisionBetween(endpoints.first, startCube);
        sys.addConstraint<physics::RigidConstraint>(sys.ecs(), endCube, endpoints.second);
        sys.disableCollisionBetween(endpoints.second, endCube);

        //   // Hinge constraint test
        // {
        //     // Tilted obstacle cube
        //     PhysicsSystem::CubeShape obstacle_shape{Vector3r(3.0, 3.0, 3.0)};
        //     Quaternion4r q_obstacle = Quaternion4r(Eigen::AngleAxis<real_t>(M_PI / 6.0, Vector3r::UnitY()));
        //     Vector3r pos_obstacle = Vector3r(0, 6, 2);
        //     PhysicsSystem::RigidState st_ob; st_ob.position = pos_obstacle; st_ob.orientation = q_obstacle; PhysicsSystem::RigidProps pr_ob; entt::entity obstacle = sys.addRigidBody(obstacle_shape, st_ob, pr_ob);

        //     Matrix33r rot_obstacle = q_obstacle.toRotationMatrix();

        //     // Disk to be hinged
        //     PhysicsSystem::CubeShape diskShape{Vector3r(1.5, 1.5, 0.05)};
        //     real_t mass = 10.0;
        //     Vector3r pos = pos_obstacle + rot_obstacle * Vector3r(0.0, 0.0, 3.5);
        //     PhysicsSystem::RigidState st; st.position = pos; st.orientation = q_obstacle; PhysicsSystem::RigidProps pr; pr.mass = mass; entt::entity disk = sys.addRigidBody(diskShape, st, pr);

        //     Vector3r a_local = Vector3r::Zero(); // hinge axis in A's local frame
        //     Vector3r b_local = Vector3r(0, 0, 0); // hinge axis in B's local frame
        //     Vector3r axis_A_frame = Vector3r::UnitZ(); 
        //     real_t K_hinge = 0;
        //     real_t D_hinge = 0;
        //     Vector2r K_rotation = Vector2r::Constant(std::numeric_limits<real_t>::infinity());
        //     Vector2r D_rotation = Vector2r::Zero();
        //     Vector3r K_translation = Vector3r::Constant(std::numeric_limits<real_t>::infinity());
        //     Vector3r D_translation = Vector3r::Zero();

        //     sys.addConstraint<physics::HingeConstraint>(sys.ecs(), obstacle, disk,
        //                                                 a_local, b_local,
        //                                                 axis_A_frame,
        //                                                 K_hinge, D_hinge,
        //                                                 K_rotation, D_rotation,
        //                                                 K_translation, D_translation);


        //     // Cube over the hinge
        //     PhysicsSystem::CubeShape cubeShape{Vector3r(0.1, 0.1, 0.1)};
        //     mass = 1.0;
        //     Vector3r pos_cube = pos_obstacle + rot_obstacle * Vector3r(1.0, 1.0, 3.75);
        //     PhysicsSystem::RigidState stc; stc.position = pos_cube; stc.orientation = q_obstacle; PhysicsSystem::RigidProps prc; prc.mass = mass; entt::entity cube = sys.addRigidBody(cubeShape, stc, prc);

        //     // sys.addConstraint<physics::TranslationalConstraint>(sys.ecs(), disk, cube, Vector3r(1, 1, 0.0), Vector3r(0.0, 0.0, 0.0), 
        //     //                                                 Vector3r(10, 10, 1e20), Vector3r::Zero());
        //     sys.addConstraint<physics::RigidConstraint>(sys.ecs(), disk, cube);

        // }

        // // Translational constraint test center
        // {
        //     // Cube shape parameters
        //     PhysicsSystem::CubeShape shape{Vector3r(0.25, 0.25, 0.25)}; // 0.5m cube

        //     real_t mass = 1.0;

        //     // Body A
        //     Vector3r posA = Vector3r(3, 0, 3);
        //     Quaternion4r qA = Quaternion4r::Identity();
        //     Vector3r vA = Vector3r(1.0, 0.0, 0.0);
        //     Vector3r wA = Vector3r(0.0, 5.0, 0.0); // initial spin around z
        //     entt::entity A; {
        //         PhysicsSystem::RigidState stA; stA.position = posA; stA.orientation = qA; stA.linearVelocity = vA; stA.angularVelocity = wA; PhysicsSystem::RigidProps prA; prA.mass = mass; A = sys.addRigidBody(shape, stA, prA);
        //     }

        //     // Body B
        //     const Vector3r he = shape.halfExtents;
        //     Vector3r posB = Vector3r(3 + 2.0 * he.x(), 2.0 * he.y(), 3 + 2.0 * he.z()); // offset so corners meet
        //     Quaternion4r qB = Quaternion4r::Identity();
        //     Vector3r vB = Vector3r(0.0, 0.0, 0.0); // initial velocity towards A
        //     Vector3r wB = Vector3r::Zero();
        //     entt::entity B; {
        //         PhysicsSystem::RigidState stB; stB.position = posB; stB.orientation = qB; stB.linearVelocity = vB; stB.angularVelocity = wB; PhysicsSystem::RigidProps prB; prB.mass = mass; B = sys.addRigidBody(shape, stB, prB);
        //     }

        //     sys.addConstraint<physics::TranslationalConstraint>(sys.ecs(), A, B);
        // }

        // // Translational constraint test
        // {
        //     // Cube shape parameters
        //     PhysicsSystem::CubeShape shape{Vector3r(0.25, 0.25, 0.25)}; // 0.5m cube

        //     real_t mass = 1.0;

        //     // Body A
        //     Vector3r posA = Vector3r(9, 0, 3);
        //     Quaternion4r qA = Quaternion4r::Identity();
        //     Vector3r vA = Vector3r::Zero();
        //     Vector3r wA = Vector3r(0.0, 5.0, 0.0); // initial spin around z
        //     entt::entity A; {
        //         PhysicsSystem::RigidState stA; stA.position = posA; stA.orientation = qA; stA.linearVelocity = vA; stA.angularVelocity = wA; PhysicsSystem::RigidProps prA; prA.mass = mass; A = sys.addRigidBody(shape, stA, prA);
        //     }

        //     // Body B
        //     const Vector3r he = shape.halfExtents;
        //     Vector3r posB = Vector3r(9 + 2.0 * he.x(), 2.0 * he.y(), 3 + 2.0 * he.z()); // offset so corners meet
        //     Quaternion4r qB = Quaternion4r::Identity();
        //     Vector3r vB = Vector3r(0.0, 0.0, 0.0); // initial velocity towards A
        //     Vector3r wB = Vector3r::Zero();
        //     entt::entity B; {
        //         PhysicsSystem::RigidState stB; stB.position = posB; stB.orientation = qB; stB.linearVelocity = vB; stB.angularVelocity = wB; PhysicsSystem::RigidProps prB; prB.mass = mass; B = sys.addRigidBody(shape, stB, prB);
        //     }

        //     // Attachment points:
        //     Vector3r rA_local = Vector3r( he.x(), he.y(), he.z());
        //     Vector3r rB_local = Vector3r( -he.x(), -he.y(), -he.z());
        //     sys.addConstraint<physics::TranslationalConstraint>(sys.ecs(), A, B, rA_local, rB_local);
        // }

        // // Double pendulum
        // {
        //     PhysicsSystem::CubeShape shape{Vector3r(0.1, 0.1, 0.1)}; // thin vertical rectangle
        //     real_t mass = 1.0;

        //     // Body A
        //     Vector3r posA = Vector3r(9 + -3.0, 0.0, 3.0);
        //     Quaternion4r qA = Quaternion4r::Identity();
        //     Vector3r vA = Vector3r::Zero();
        //     Vector3r wA = Vector3r::Zero();
        //     entt::entity A; {
        //         PhysicsSystem::RigidState stA; stA.position = posA; stA.orientation = qA; stA.linearVelocity = vA; stA.angularVelocity = wA; PhysicsSystem::RigidProps prA; prA.mass = mass; A = sys.addRigidBody(shape, stA, prA); sys.makeStatic(A);
        //     }

        //     // Body B
        //     Vector3r posB = Vector3r(9 + -3.0, 0.0, 2.0);
        //     Quaternion4r qB = Quaternion4r::Identity();
        //     Vector3r vB = Vector3r(3.0, 2.0, 1.0);
        //     Vector3r wB = Vector3r::Zero();
        //     entt::entity B; {
        //         PhysicsSystem::RigidState stB; stB.position = posB; stB.orientation = qB; stB.linearVelocity = vB; stB.angularVelocity = wB; PhysicsSystem::RigidProps prB; prB.mass = mass; B = sys.addRigidBody(shape, stB, prB);
        //     }

        //     // Body C
        //     Vector3r posC = Vector3r(9 + -3.0, 0.0, 0.5);
        //     Quaternion4r qC = Quaternion4r::Identity();
        //     Vector3r vC = Vector3r(12.0, 2.0, 1.0);
        //     Vector3r wC = Vector3r::Zero();
        //     entt::entity C; {
        //         PhysicsSystem::RigidState stC; stC.position = posC; stC.orientation = qC; stC.linearVelocity = vC; stC.angularVelocity = wC; PhysicsSystem::RigidProps prC; prC.mass = mass; C = sys.addRigidBody(shape, stC, prC);
        //     }

        //     // Constraints
        //     sys.addConstraint<physics::LinearDistanceConstraint>(sys.ecs(), A, B, Vector3r(0.0, 0.0, -0.1), Vector3r(0.0, 0.0, 0.1));
        //     sys.addConstraint<physics::LinearDistanceConstraint>(sys.ecs(), B, C, Vector3r(0.0, 0.0, -0.1), Vector3r(0.0, 0.0, 0.1));
        // }
    }

    void updateScene(cardillo::PhysicsSystem& sys, real_t t, real_t /*dt*/) override {
        // Apply a twisting moment at the rod end
        // real_t torque_magnitude = 0.05;
        real_t torque_max = 1e3;
        real_t t1 = 1.5;
        real_t torque_magnitude = torque_max * std::min((t - 0.0) / (t1 - 0.0), 1.0);
        sys.applyForce(m_rodEnd, Vector3r::Zero(), Vector3r(0, -torque_magnitude, 0));
        sys.applyForce(m_ropeEnd, Vector3r::Zero(), Vector3r(0.1, 0, 0));

        // if (t < 0.05) {
        //     sys.setLinearVelocity(m_endCube, Vector3r(0, -1, 0));
        //     sys.setAngularVelocity(m_endCube, Vector3r::Zero());
        // }else {
        //     sys.setLinearVelocity(m_endCube, Vector3r::Zero());
        //     sys.setAngularVelocity(m_endCube, Vector3r(0, 2, 0));
        // }
        sys.setLinearVelocity(m_endCube, Vector3r::Zero());
        sys.setAngularVelocity(m_endCube, Vector3r(0, 2, 0));
        sys.applyForce(m_endCube, -sys.gravity() * sys.getMass(m_endCube).diagonal(), Vector3r::Zero());
    }

    private:
        entt::entity m_rodEnd{entt::null};
        entt::entity m_ropeEnd{entt::null};
        entt::entity m_endCube{entt::null};
};
