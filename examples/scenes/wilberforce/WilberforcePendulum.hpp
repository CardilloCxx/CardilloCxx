#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <memory>
#include <cmath>

using namespace cardillo;

// Wilberforce Pendulum:
// A mass-spring system where torsional and vertical modes are tuned to exchange energy.
// We model the spring as a helical beam (capsule cross-section) anchored at the top,
// with a rigid mass at the bottom. Initial conditions excite vertical mode and small twist.
// see https://onlinelibrary.wiley.com/doi/epdf/10.1002/pamm.202100110
class WilberforcePendulumScene : public SceneBase {
public:
    const char* sceneName() const override { return "wilberforce"; }
    WilberforcePendulumScene() = default;
    ~WilberforcePendulumScene() = default;

    void populate(physics::PhysicsEngine& engine) override {
        using namespace cardillo;
        using namespace misc;

        // Spring parameters
        const real_t wireDiameter = 0.001;               // wire diameter 1 mm
        const real_t wireRadius = 0.5 * wireDiameter;    // wire radius 0.5 mm
        const real_t coilDiameter = 0.032;               // spring mean radius 32 mm
        const real_t coilRadius = 0.5 * coilDiameter;    // spring mean radius 16 mm
        // const size_t turns = 20;                         // helical turns count
        const size_t turns = 12;                          // helical turns count
        const size_t segmentsPerTurn = 30;               // recommended minimum segments per turn (24)
        const size_t segments = turns * segmentsPerTurn;
        
        // Material properties
        const real_t density = 7850.0;               // kg/m^3
        const real_t E = 206e9;                      // Young's modulus (GPa)
        const real_t nu = 0.2638;                    // Poisson ratio (G ~ 81.5 GPa)
        const real_t G = E / (2.0 * (1.0 + nu));     // Shear modulus

        engine.setGravity(Vector3r(0, 0, -9.81)); // no gravity

        // Beam cross-section (capsule) used by createBeam
        physics::BeamCrossSection sec(wireDiameter, wireDiameter, physics::BeamBodyType::Capsule);
        auto springs = physics::BeamSpringParams::fromMaterial(E, nu);
        // springs.setDampingFromFactor(0.001); // set damping factor

        const real_t pitch = wireDiameter;
        const real_t freeLength = static_cast<real_t>(turns) * pitch;
        // Build spring via CompoundSpline (single helix segment now, extensible for future pieces)
        auto helix = misc::HelixSpline(Vector3r::Zero(), -Vector3r::UnitZ(), coilRadius, pitch, static_cast<real_t>(turns), Vector3r::UnitX());

        // Build sequence of splines; create beams per spline and connect with rigid constraints.
        std::vector<const misc::SplinePattern*> parts{&helix};
        auto endpoints = engine.createBeams(parts, sec, springs, physics::RigidState{}, physics::RigidProps::withDensity(density), segments);
        m_top = endpoints.first;
        // cube_constraint = engine.addRigidConstraint(m_top);
        // TODO: I think `getPosition` should be named `getPose`
        VectorXr pose0 = engine.getPosition(m_top);
        Vector3r position0 = pose0.head<3>();
        Quaternion4r orientation0 = Quaternion4r(pose0.tail<4>().data());
        engine.makeStatic(m_top);
        engine.addTrajectory(m_top,
                            [position0, orientation0](real_t t) {
                                TrajectoryPose pose;
                                pose.first = position0 + Vector3r(0.0, 0.05 * std::sin(t), 0.0);
                                pose.second = (
                                    Quaternion4r(Eigen::AngleAxis<real_t>(0.25 * M_PI * std::sin(t), Vector3r::UnitX()))
                                    * orientation0
                                ).normalized();
                                return pose;
                            },
                            std::nullopt
                            );
                            //  std::nullopt, 
                            //  [](real_t t) {
                            //      TrajectoryTwist twist;
                            //      twist.first = {0, 0.0 * std::sin(10 * t), 0.0};
                            //      twist.second = Vector3r::Zero();
                            //      return twist;
                            //  });

        m_bottom = endpoints.second;

        real_t tunedMass = 0.5240569245269475;

        m_bob = engine.addRigidBody(physics::MeshShape(std::string(PROJECT_SOURCE_DIR) + "/res/meshes/bob3.obj"), physics::RigidState(Vector3r(0, 0, -0.025 - static_cast<real_t>(turns) * pitch - wireDiameter), Vector3r(0, 0, 0)),
                                    physics::RigidProps(tunedMass));

        // working:
        // Bob inertia diag: Ixx = 0.000105331, Iyy = 0.000124083, Izz = 0.000124083
        // Bob inertia world frame:
        //  0.000124083 -4.08802e-12  1.32047e-10
        // -4.08802e-12  0.000124083 -1.34392e-10
        //  1.32047e-10 -1.34392e-10  0.000105331
        // Bob mass: 0.524057 kg

        // Set inertia
        engine.ecs().get<C_InertiaDiag>(m_bob).I = Vector3r(0.000124083, 0.000124083, 0.000105331);

        // track bob trajectory in csv file
        engine.track(m_bob, "bob");
        // Print Inertia and Mass for verification
        auto Idiag = engine.ecs().get<C_InertiaDiag>(m_bob).I;
        std::cout << "Bob inertia diag: Ixx = " << Idiag.x() << ", Iyy = " << Idiag.y() << ", Izz = " << Idiag.z() << std::endl;

        auto RotMat = engine.ecs().get<C_Orientation>(m_bob).rotation;
        std::cout << "Bob inertia world frame:\n" << RotMat * Idiag.asDiagonal() * RotMat.transpose() << std::endl;
        std::cout << "Bob mass: " << engine.getMass(m_bob).col(0).row(0) << " kg" << std::endl;

        // Constraint the bob to only allow vertical and torsional motion
        // const real_t inf = std::numeric_limits<real_t>::infinity();
        // engine.addTranslationRotationConstraint(m_bob, m_top, physics::JointFrame{},
        //     Vector3r(inf, inf, 0), Vector3r::Zero(),
        //     Vector3r(inf, inf, 0), Vector3r::Zero());

        // Pin bottom endpoint to bob using a rigid constraint
        engine.addRigidConstraint(m_bottom, m_bob);
    }

    void updateScene(physics::PhysicsEngine& engine, real_t t, real_t /*dt*/) override 
    {
        // Pull bob downward slowly to start vertical oscillation
        const real_t vz0 = -0.5;
        const real_t t0 = 0.3;
        if (t < t0) {
            engine.setVelocityByForce(m_bob, Vector3r(0,0,vz0 * t / t0), Vector3r(0,0,0));
        }

        // engine.setConstraintLinearVelocity(cube_constraint, Vector3r(0, 0, 0.5 * std::sin(t)));
    }

private:
    entt::entity m_top{entt::null};
    entt::entity m_bottom{entt::null};
    entt::entity m_bob{entt::null};
    bool m_reset{false};
    index_t cube_constraint{-1};
};