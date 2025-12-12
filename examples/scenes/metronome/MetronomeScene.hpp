#pragma once

#include "../SceneBase.hpp"
#include <random>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cstdio>

// Grid of metronomes resting on a shared sheet supported by two rollers.
// Each metronome is a mesh body (optionally box-collided via config), with a lever
// attached through a hinge with rotational spring/damper. Collisions between body
// and lever are disabled. Levers receive staggered initial kicks in pseudo-random order.
class MetronomeScene : public SceneBase {
public:
    const char* sceneName() const override { return "metronome"; }
    MetronomeScene() = default;
    ~MetronomeScene() override = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;

        // Static ground
        auto floor = sys.addStaticBody(PhysicsSystem::CubeShape(Vector3r(5.0, 5.0, 0.1)), PhysicsSystem::RigidState(Vector3r(0, 0, -0.1)));

        const int nx = 1;
        const int ny = 3;
        const real_t spacing = (real_t)0.12;

        // Two supporting rollers (dynamic, aligned along +X)
        const real_t rollerR = (real_t)0.05;
        const Quaternion4r rollerQ(Eigen::AngleAxis<real_t>((real_t)(-M_PI_2), Vector3r::UnitY()));

        // Plate spans the grid with margin
        const real_t margin = (real_t)0.05;
        const real_t plateHalfLength = ((real_t)(nx - 1) * spacing) * (real_t)0.5 + margin;
        const real_t plateHalfWidth  = ((real_t)(ny - 1) * spacing) * (real_t)0.5 + margin;
        const real_t rollerHalfL = plateHalfLength + (real_t)0.05; // rollers slightly shorter than plate length
        const real_t rollerOffsetY = plateHalfWidth * (real_t)0.4;

        const real_t styrofoamRho = (real_t)50.0; // approximate density for styrofoam
        sys.addRigidBody(PhysicsSystem::CapsuleShape(rollerR, rollerHalfL), PhysicsSystem::RigidState(Vector3r(0.0, -rollerOffsetY, rollerR), rollerQ), PhysicsSystem::RigidProps(0.015));
        sys.addRigidBody(PhysicsSystem::CapsuleShape(rollerR, rollerHalfL), PhysicsSystem::RigidState(Vector3r(0.0, rollerOffsetY, rollerR), rollerQ), PhysicsSystem::RigidProps(0.015)); 

        // Sheet resting on rollers
        const real_t sheetThickness = (real_t)0.0025;
        const Vector3r sheetHalfExt(plateHalfLength, plateHalfWidth, sheetThickness);
        const real_t sheetCenterZ = 2* rollerR + sheetThickness + (real_t)0.001; // lifted slightly above roller top
        m_sheet = sys.addRigidBody(
            PhysicsSystem::CubeShape(sheetHalfExt),
            PhysicsSystem::RigidState(Vector3r(0, 0, sheetCenterZ)),
            PhysicsSystem::RigidProps::withDensity(styrofoamRho));
        sys.track(m_sheet, "sheet");    

        sys.addConstraint<physics::RotationConstraint>(
            sys.ecs(),
            m_sheet,
            floor,
            physics::JointFrame(m_sheet),
            Vector3r::Constant(std::numeric_limits<real_t>::infinity())
        );

        
        sys.addConstraint<physics::TranslationalConstraint>(
            sys.ecs(),
            m_sheet,
            floor,
            physics::JointFrame(m_sheet),
            Vector3r::Constant(0),
           Vector3r::Constant(0.1)
        );


        

        // Grid placement
        const real_t baseX = -((real_t)(nx - 1) * spacing * (real_t)0.5);
        const real_t baseY = -((real_t)(ny - 1) * spacing * (real_t)0.5);
        const real_t baseZ = sheetCenterZ + sheetThickness; // place bottom center on sheet top

        const Vector3r hingeLocal((real_t)-0.025, (real_t)0.0015, (real_t)0.04);
        const real_t hingeK = (real_t)0.00;
        const real_t hingeD = (real_t)0.0000;

        for (int ix = 0; ix < nx; ++ix) {
            for (int iy = 0; iy < ny; ++iy) {
                const Vector3r pos(baseX + (real_t)ix * spacing, baseY + (real_t)iy * spacing, baseZ + 0.001);

                PhysicsSystem::MeshShape bodyShape("res/meshes/metronome.obj", true);
                auto body = sys.addRigidBody(
                    bodyShape,
                    PhysicsSystem::RigidState(pos),
                    PhysicsSystem::RigidProps((real_t)0.015)); // 100g body

                PhysicsSystem::MeshShape leverShape("res/meshes/metronome_lever.obj", true);
                auto lever = sys.addRigidBody(
                    leverShape,
                    PhysicsSystem::RigidState(pos, Vector3r(0, 0, 0)),
                    PhysicsSystem::RigidProps((real_t)0.05)); // metal lever

                const Vector3r hingeWorld = pos + hingeLocal; // offset from bottom-center in world frame
                physics::JointFrame jf = physics::JointFrame::fromAxis(hingeWorld, Vector3r::UnitX());
                sys.addConstraint<physics::HingeConstraint>(sys.ecs(), body, lever, jf, hingeK, hingeD);
                sys.disableCollisionBetween(body, lever);
                sys.track(lever, "lever_" + std::to_string(ix) + "_" + std::to_string(iy));

                m_levers.push_back(lever);
            }
        }

        // Prepare pseudo-random kick order
        const int total = (int)m_levers.size();
        m_order.resize(total);
        std::iota(m_order.begin(), m_order.end(), 0);
        // std::shuffle(m_order.begin(), m_order.end(), m_rng);
    }

    void updateScene(cardillo::PhysicsSystem& sys, real_t t, real_t dt) override {
        const real_t kickStart = (real_t)0.25;
        const real_t kickEnd = (real_t) 0.5;
        const real_t kickInterval = (kickEnd - kickStart) / (real_t)m_levers.size();
        const real_t kickVelocity = (real_t)0.5; // m/s

        while (m_nextKick < (int)m_order.size() && t >= kickStart + (real_t)m_nextKick * kickInterval) {
            const int idx = m_order[(size_t)m_nextKick];
            if (idx >= 0 && idx < (int)m_levers.size()) {
                entt::entity lever = m_levers[(size_t)idx];
                real_t kickForce = kickVelocity * sys.getMass(lever)(0, 0) / dt;
                if (sys.ecs().valid(lever)) {
                    sys.applyForce(lever, Vector3r((real_t)0.0, kickForce, (real_t)0.0), Vector3r::Zero());
                    sys.applyForce(m_sheet, Vector3r((real_t)0.0, -kickForce, (real_t)0.0), Vector3r::Zero()); // reaction on sheet
                    std::printf("[Metronome] Kick lever %d at t=%.3f with Fy=%.3f N\n", idx, (double)t, (double)kickForce);
                }
            }
            ++m_nextKick;
        }
    }

private:
    std::mt19937 m_rng{1337u};
    entt::entity m_sheet{entt::null};
    std::vector<entt::entity> m_levers;
    std::vector<int> m_order;
    int m_nextKick{0};
};
