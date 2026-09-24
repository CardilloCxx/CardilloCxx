#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace cardillo;

// vibrations of a cantilever beam
class CantileverScene : public SceneBase {
public:
    const char* sceneName() const override { return "cantilever"; }
    CantileverScene() = default;
    ~CantileverScene() = default;

    void populate(physics::PhysicsEngine& engine) override {
        using namespace cardillo;
        using namespace misc;

        engine.setGravity(Vector3r(0, 0, -9.81)); // no gravity for this scene

        // Material and geometry
        const real_t L  = (real_t)0.61415; // 31 cm
        const real_t r  = (real_t)0.01;    // 10 mm radius
        const real_t d  = (real_t)2 * r;   // 20 mm diameter
        const real_t E  = (real_t)7e7;     // 7 GPa
        const real_t nu = (real_t)0.75;    // Poisson ratio
        // const real_t rho = (real_t)0.4 / (M_PI * r * r * L);  // kg/m^3 (mass 0.4 kg)
        const real_t rho = (real_t)1350.0;  // kg/m^3 (typical plastic)

        // Beam spring params from material (no extra damping by default).
        auto springs = physics::BeamSpringParams::fromMaterial(E, nu);

        // Build a straight beam along +X of length L, elevated above ground
        const Vector3r p0(0, 0, 0);
        const Vector3r p1 = p0 + Vector3r(L, 0, 0);
        misc::LinearSpline spline(p0, p1);

        // Number of beam segments
        const size_t segments = 5;

        // Default state: identity orientation; default density props
        physics::RigidState stateDefaults(Vector3r::Zero(), Vector3r::Zero(), Quaternion4r::Identity());
        physics::RigidProps props = physics::RigidProps::withDensity(rho);

        // const auto section = physics::BeamCrossSection::triangle(d,d);
        // const auto section = physics::BeamCrossSection::round(d);
        const auto section = physics::BeamCrossSection::rectangular(d * 2.0, d);
        auto beam_ends = engine.createBeam(spline, section, springs, stateDefaults, props, segments, physics::BeamColliderMode::InterSegmentHull);
        engine.makeStatic(beam_ends.first);

        const Vector3r initalEndPos = engine.getPosition(beam_ends.second).head<3>();
        const Vector4r initialEndRot = engine.getPosition(beam_ends.second).tail<4>();

        std::cout << "Initial end position: " << initalEndPos.transpose() << std::endl;

    //        void addTrajectory(entt::entity e, std::optional<std::function<TrajectoryPose(real_t)>> positionFunc, std::optional<std::function<TrajectoryTwist(real_t)>> velocityFunc) {
    //     m_world->setTrajectory(e, std::move(positionFunc), std::move(velocityFunc));
    // }

        engine.addTrajectory(beam_ends.second, [L, initialEndRot, initalEndPos](real_t t) {
            TrajectoryPose pose;
            
            const real_t duration_per_phase = 3.0; 
            const real_t current_phase = std::floor(t / duration_per_phase);
            const real_t tau = (t - current_phase * duration_per_phase) / duration_per_phase;
            const real_t smooth_cycle = 0.5 * (1.0 - std::cos(2.0 * M_PI * tau));
            const real_t translation_dist = 0.1;               // Distance to move away and back
            const real_t max_angle = M_PI / 4.0;                // 45 degrees tilt in radians

            Vector3r position_offset(0.0, 0.0, 0.0);
            Vector3r rotation_axis(0.0, 0.0, 0.0);
            real_t angle = 0.0;

            // Determine motion based on active phase
            int phase_idx = static_cast<int>(current_phase) % 6;
            switch (phase_idx) {
                case 0:
                    position_offset.x() = translation_dist * smooth_cycle;
                    break;
                case 1: 
                    position_offset.y() = translation_dist * smooth_cycle;
                    break;
                case 2: 
                    position_offset.z() = 0.5 * translation_dist * smooth_cycle;
                    break;
                case 3:
                    rotation_axis = Vector3r(0.0, 0.0, 1.0);
                    angle = max_angle * smooth_cycle;
                    break;
                case 4: 
                    rotation_axis = Vector3r(0.0, 1.0, 0.0);
                    angle = 0.5 * max_angle * smooth_cycle;
                    break;
                case 5:
                    rotation_axis = Vector3r(1.0, 0.0, 0.0);
                    angle = 0.5 * max_angle * smooth_cycle;
                    break;
                   
            }

            // Set final position
            pose.first = Vector3r(initalEndPos) + position_offset;

            // Apply rotation perturbation to initial orientation
            Quaternion4r delta_rot(Eigen::AngleAxis<real_t>(angle, rotation_axis));
            pose.second = Quaternion4r(initialEndRot.data()) * delta_rot; 

            return pose;
        }, std::nullopt);

        // engine.addTranslationalConstraint(beam_ends.second, entt::null, physics::JointFrame(beam_ends.second), );


//         // Create a second beam
//         const auto section2 = physics::BeamCrossSection::rectangular(0.75* d, 0.1 * d);
//         stateDefaults.setOrientation(Quaternion4r(Eigen::AngleAxis<real_t>(M_PI_2 / 2.0, Vector3r::UnitZ())));
// 
//         stateDefaults.angularVelocity = Vector3r(0, 0, 10.0);
//         misc::LinearSpline spline2(Vector3r(L * 0.5, -L * 0.5, 0.1), Vector3r(L * 0.5, L * 0.5, 0.1));
//         auto springs2 = physics::BeamSpringParams::fromMaterial(E * 0.4, nu);
//         beam_ends = engine.createBeam(spline2, section2, springs2, stateDefaults, props, segments * 2, physics::BeamColliderMode::InterSegmentHull);

        // spheres falling on the beam
        // engine.addRigidBody(physics::SphereShape(0.01), physics::RigidState(Vector3r(L * 0.5, 0.1, 0.06)), physics::RigidProps::withDensity(20000.0));

        // Cubes sitting on the beam lined up in axial direction
        const real_t cubeSize = 0.02;
        const size_t numCubes = 15;
        for (size_t i = 0; i < numCubes; ++i) {
            const real_t x = L * 0.15 + (real_t)i * cubeSize * 1.5;
            engine.addRigidBody(physics::CubeShape(Vector3r(cubeSize * 0.5, cubeSize * 0.5, cubeSize * 0.5)), physics::RigidState(Vector3r(x, 0, cubeSize * 0.51 + r)), physics::RigidProps::withDensity(200.0));
        }

        m_Kf = springs.Kf(L / segments, section);
        m_L = L;
        m_segments = segments;

    }

    void updateScene(physics::PhysicsEngine& engine, real_t t, real_t /*dt*/) override {
        real_t t1 = 3.0;
        if (t < t1 * 10.0) {
            // engine.applyForce(m_beamRightEnd, Vector3r(-0.5 * std::max(t1, t), -0.5 * std::max(t1, t), 1.5 * std::max(t1, t)), Vector3r::Zero());
            
            real_t time = std::min(t, t1);
            Vector3r force(time / t1, time / t1, 0.0);
            // Vector3r moment(0.0, -2 * M_PI * m_Kf(1) * m_L / m_segments / m_L * t / t1 / 2, 0.0);
            Vector3r torsion(5.0 * m_Kf(2) * m_L / m_segments / m_L * time / t1, 0.0, 0.0);
            // engine.applyForce(m_beamRightEnd1, force * 1.0, torsion * 0.0);

            // engine.applyForce(m_beamRightEnd3, force, moment);
            // engine.applyForce(m_beamRightEnd, Vector3r(0, -0.5 * std::max(t1, t), 0), Vector3r(0, -1.0 * std::max(t1, t), 0));
            // engine.applyForce(m_beamRightEnd, Vector3r::Zero(), Vector3r(-1.0 * std::max(t1, t), 0, 0));
            // engine.applyForce(m_beamRightEnd, Vector3r::Zero(), Vector3r(0, -1.0 * std::max(t1, t), 0));
            // engine.applyForce(m_beamRightEnd, Vector3r::Zero(), Vector3r(0, 0, -1.0 * std::max(t1, t)));
        }

        // if (t > 0.1) 
        // {
        //     // Make every entity static except the sphere:
        //     auto rb_view = engine.ecs().view<C_PhysicsObject>();
        //     for (auto [entity] : rb_view.each()) {
        //         if (entity == m_beamRightEnd1 || entity == m_beamRightEnd2 || entity == m_beamRightEnd3) continue;
        //         engine.makeStatic(entity);
        //     }
        // }
    }

    private:
     entt::entity m_beamRightEnd1{entt::null};
     Vector3r m_Kf{Vector3r::Zero()};
     real_t m_L{0.0};
     size_t m_segments{0};
};