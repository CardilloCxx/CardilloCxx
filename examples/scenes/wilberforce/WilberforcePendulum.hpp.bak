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

    void populate(cardillo::physics::PhysicsEngine& engine) override {
        auto& sys = engine.world();
        using namespace cardillo;
        using namespace cardillo::misc;

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

        sys.setGravity(Vector3r(0, 0, -9.81)); // no gravity

        // real_t custom_mass = 0.5240569245269475;
        // // real_t Ixx = 0.00013237;
        // // real_t Iyy = Ixx; 
        // // real_t Izz = 0.00016377;
        // real_t Ixx = 0.00016377;
        // real_t Iyy = 0.00013237;
        // real_t Izz = Iyy; 
        // // const real_t sizeZ = std::sqrt(6 / custom_mass * (Ixx + Iyy - Izz));
        // // const real_t sizeZ = std::sqrt(6 / custom_mass * (Iyy + Izz - Ixx));
        // const real_t sizeX = std::sqrt(6 / custom_mass * (Iyy + Izz - Ixx));
        // const real_t sizeY = std::sqrt(6 / custom_mass * (Ixx + Izz - Iyy));
        // const real_t sizeZ = std::sqrt(6 / custom_mass * (Ixx + Iyy - Izz));

        // Beam cross-section (capsule) used by createBeam
        physics::BeamCrossSection sec(wireDiameter, wireDiameter, physics::BeamBodyType::Capsule);
        auto springs = physics::BeamSpringParams::fromMaterial(E, nu);
        // springs.setDampingFromFactor(0.001); // set damping factor

        const real_t pitch = wireDiameter;
        const real_t freeLength = static_cast<real_t>(turns) * pitch;
        // Build spring via CompoundSpline (single helix segment now, extensible for future pieces)
        auto helix = misc::HelixSpline(Vector3r::Zero(), -Vector3r::UnitZ(), coilRadius, pitch, static_cast<real_t>(turns), Vector3r::UnitX());
        auto linear = misc::LinearSpline(Vector3r(coilRadius, 0, -freeLength), Vector3r(0, 0, -freeLength));
        // auto linear2 = misc::LinearSpline( Vector3r(0, 0, -freeLength), Vector3r(0, 0, -freeLength -sizeX));
        auto linear2 = misc::LinearSpline( Vector3r(0, 0, -freeLength), Vector3r(0, 0, -freeLength - 0.01));

        // Build sequence of splines; create beams per spline and connect with rigid constraints.
        // std::vector<const misc::SplinePattern*> parts{&helix, &linear, &linear2};
        // std::vector<const misc::SplinePattern*> parts{&helix, &linear};
        std::vector<const misc::SplinePattern*> parts{&helix};
        auto endpoints = engine.createBeams(parts, sec, springs, physics::RigidState{}, physics::RigidProps::withDensity(density), segments);
        m_top = endpoints.first; sys.makeStatic(m_top);
        m_bottom = endpoints.second;
    
        // // NEW: Modify mass and inertia of the bottom node
        // sys.ecs().get<World::C_Mass>(m_bottom).m *= 0.5;
        // sys.ecs().get<World::C_InertiaDiag>(m_bottom).I *= 0.5;
    
        // // NEW: Modify mass and inertia of the bottom node
        // sys.ecs().get<World::C_Mass>(m_bottom).m = custom_mass;
        
        // // Set custom diagonal inertia for anisotropic behavior
        // sys.ecs().get<World::C_InertiaDiag>(m_bottom).I = Vector3r(Ixx, Iyy, Izz);

        // sys.ecs().remove<World::C_Capsule>(m_bottom);
        // sys.ecs().remove<World::C_CapsuleVisualTag>(m_bottom);

        // sys.ecs().emplace<World::C_CubeVisualTag>(m_bottom);
        // Vector3r halfExtents(sizeX/2.0, sizeY/2.0, sizeZ/2.0);
        // Quaternion4r orientation = Quaternion4r(Eigen::AngleAxis<real_t>(M_PI / 2, Vector3r::UnitY()));
        // // sys.ecs().emplace<World::C_Cube>(m_bottom, World::C_Cube{Vector3r(0, 0, -freeLength -sizeX/2), halfExtents, orientation});
        // sys.ecs().emplace<World::C_Cube>(m_bottom, World::C_Cube{Vector3r(0, 0, 0), halfExtents, orientation});

        // // // Add cube visual tag and component with appropriate dimensions
        // // // (adjust halfExtents to match your beam's cross-section)
        // // Vector3r halfExtents(beamLength/2.0, beamWidth/2.0, beamHeight/2.0);
        // // sys.ecs().emplace<C_CubeVisualTag>(bottomBeamEntity);
        // // sys.ecs().emplace<C_Cube>(bottomBeamEntity, 
        // //     C_Cube{Vector3r::Zero(), halfExtents, Quaternion4r::Identity()});

        // // track bob trajectory in csv file
        // sys.ecs().emplace<World::C_TrackTag>(m_bottom, World::C_TrackTag{});

        // // we want K / m = lambda / Iz
        // const real_t d = wireDiameter;      // wire diameter
        // const real_t D = coilRadius * 2.0;  // mean coil diameter
        // const real_t n = static_cast<real_t>(turns); // number of active coils
        // const real_t K =      (G * std::pow(d, 4)) / (8.0 * std::pow(D, 3) * n); // helix axial stiffness
        // const real_t lambda = (G * std::pow(d, 4) )/ (32.0 * n * D);        // helix torsional stiffness
        // const real_t tunedMass = 0.469;
        // const real_t tunedMass = 0.95;
        // const real_t Iz = lambda / (K / tunedMass);
        //////////////////////////////////////////////////////////////////
        // parameter marco
        // mass = 0.5240569245269475
        // Inertia = array([[0.00013237, 0.        , 0.        ],
        //    [0.        , 0.00013237, 0.        ],
        //    [0.        , 0.        , 0.00016377]])
        //////////////////////////////////////////////////////////////////
        real_t tunedMass = 0.5240569245269475;
        const real_t Ix = 0.00013237;
        const real_t Iy = Ix;
        const real_t Iz = 0.00016377;
        // const real_t Iz = 0.0002;
        // const real_t Iz = 0.00025;
        const real_t sizeX = std::sqrt(6 / tunedMass * (Iy + Iz - Ix));
        const real_t sizeY = std::sqrt(6 / tunedMass * (Ix + Iz - Iy));
        const real_t sizeZ = std::sqrt(6 / tunedMass * (Ix + Iy - Iz));

        // Ixx = 0.000104211, Iyy = 0.000122214, Izz = 0.000122214
        // const real_t tunedSize = std::sqrt(6 * Iz / tunedMass); // Iz = 1 / 12 m (a^2 + a^2) = 1/6 m a^2 => a = sqrt(6 Iz / m)

        // std::cout << "Tuning info: K = " << K << " N/m, lambda = " << lambda << " Nm/rad, m = " << tunedMass
        //           << " kg, Iz = " << Iz << " kg m^2, cube size = " << tunedSize << " m" << "Expected frequencies: f_z = " << std::sqrt(K / tunedMass) / (2 * M_PI) << " Hz"
        //           << " f_theta = " << std::sqrt(lambda / Iz) / (2 * M_PI) << " Hz, Spring mass: " << sys.getMass(m_bottom).col(0).row(0) * segments << " kg" << std::endl;

        std::cout << "sizeX: " << sizeX << ", sizeY: " << sizeY << ", sizeZ: " << sizeZ << std::endl;

        // m_bob = engine.addRigidBody(physics::CubeShape(Vector3r(tunedSize,tunedSize,tunedSize)),
        // m_bob = engine.addRigidBody(physics::CubeShape(Vector3r(sizeX/2, sizeY/2, sizeZ/2)),
        m_bob = engine.addRigidBody(physics::MeshShape("res/meshes/bob2.obj"),
            // physics::RigidState(Vector3r(0,0,-tunedSize) + sys.getPosition(m_bottom).head<3>()),
            // physics::RigidState(Vector3r(0,0,-sizeZ) + sys.getPosition(m_bottom).head<3>()),
            physics::RigidState(Vector3r(0,0,-0.025 - static_cast<real_t>(turns) * pitch - wireDiameter), Vector3r(0, 0, 0)),
            physics::RigidProps(tunedMass));

        

        // Set inertia to match target Iz for cube
        sys.ecs().get<World::C_InertiaDiag>(m_bob).I *= 0.75;

        // track bob trajectory in csv file
        sys.track(m_bob, "bob");
        // Print Inertia for verification
        auto Idiag = sys.ecs().get<World::C_InertiaDiag>(m_bob).I;
        std::cout << "Bob inertia diag: Ixx = " << Idiag.x() << ", Iyy = " << Idiag.y() << ", Izz = " << Idiag.z() << std::endl;

        auto RotMat = sys.ecs().get<World::C_Orientation>(m_bob).value.toRotationMatrix();
        std::cout << "Bob inertia world frame:\n" << RotMat * Idiag.asDiagonal() * RotMat.transpose() << std::endl;
        std::cout << "Bob mass: " << sys.getMass(m_bob).col(0).row(0) << " kg" << std::endl;


        // Constraint the bob to only allow vertical and torsional motion
        // const real_t inf = std::numeric_limits<real_t>::infinity();
        // sys.addTranslationRotationConstraint(m_bob, m_top, physics::JointFrame{},
        //     Vector3r(inf, inf, 0), Vector3r::Zero(),
        //     Vector3r(inf, inf, 0), Vector3r::Zero());
        
 
        // Pin bottom endpoint to bob using a rigid constraint
        // sys.addRigidConstraint(m_bottom, m_bob, Vector3r::Zero(), Vector3r(0,0,tunedSize));
        sys.addRigidConstraint(m_bottom, m_bob);
        
        // const real_t vz0 = 0.0; // m/s downward
        // const real_t wz0 = 0.0;  // rad/s small spin around Z to couple torsion & vertical modes

        // sys.setLinearVelocity(m_bob, Vector3r(0,0,vz0));
        // sys.setAngularVelocity(m_bob, Vector3r(0,0,wz0));
    }

    void updateScene(cardillo::physics::PhysicsEngine& engine, real_t t, real_t /*dt*/) override 
    {
        auto& sys = engine.world();
        // Pull bob downward slowly to start vertical oscillation
        const real_t vz0 = -0.5;
        const real_t t0 = 0.3;
        if (t < t0) {
            sys.setVelocityByForce(m_bob, Vector3r(0,0,vz0 * t / t0), Vector3r(0,0,0));
        }

    }

private:
    entt::entity m_top{entt::null};
    entt::entity m_bottom{entt::null};
    entt::entity m_bob{entt::null};
    bool m_reset{false};
};