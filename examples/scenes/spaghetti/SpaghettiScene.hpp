#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>
#include <vector>
#include <string>

using namespace cardillo;

// Single spaghetti modeled as a straight beam made of capsule segments
class SpaghettiScene : public SceneBase {
public:
    const char* sceneName() const override { return "spaghetti"; }
    SpaghettiScene()  = default;
    ~SpaghettiScene() = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;
        using namespace cardillo::misc;

        // Gravity on (spaghetti should hang, z is up).
        sys.setGravity(Vector3r(0, 0, -9.81));

        // Add a bowl built from 5 flat cubes (floor + 4 walls) aligned with z-up.
        {
            const real_t innerRadius   = (real_t)0.075;   // matches spaghetti layout radius
            const real_t wallThickness = (real_t)0.01;
            const real_t wallHeight    = (real_t)0.06;
            const real_t floorZ        = (real_t)-0.08;
            const real_t floorHalfZ    = wallThickness * (real_t)0.5;
            const real_t wallHalfZ     = wallHeight * (real_t)0.5;
            const real_t wallOffset    = innerRadius + wallThickness * (real_t)0.5;

            const Vector3r floorHalfExtents(innerRadius + wallThickness, innerRadius + wallThickness, floorHalfZ);
            const Vector3r wallHalfExtentsY(innerRadius, wallThickness * (real_t)0.5, wallHalfZ);
            const Vector3r wallHalfExtentsX(wallThickness * (real_t)0.5, innerRadius, wallHalfZ);

            const real_t wallCenterZ = floorZ + floorHalfZ + wallHalfZ;
            
            // Floor
            cardillo::physics::BodyFactory::addStaticBody(sys, PhysicsSystem::CubeShape{floorHalfExtents}, PhysicsSystem::RigidState(Vector3r(0, 0, floorZ)));
            // +Y wall (back)
            cardillo::physics::BodyFactory::addStaticBody(sys, PhysicsSystem::CubeShape{wallHalfExtentsY}, PhysicsSystem::RigidState(Vector3r(0,  wallOffset, wallCenterZ)));
            // -Y wall (front)
            cardillo::physics::BodyFactory::addStaticBody(sys, PhysicsSystem::CubeShape{wallHalfExtentsY}, PhysicsSystem::RigidState(Vector3r(0, -wallOffset, wallCenterZ)));
            // +X wall (right)
            cardillo::physics::BodyFactory::addStaticBody(sys, PhysicsSystem::CubeShape{wallHalfExtentsX}, PhysicsSystem::RigidState(Vector3r( wallOffset, 0, wallCenterZ)));
            // -X wall (left)
            cardillo::physics::BodyFactory::addStaticBody(sys, PhysicsSystem::CubeShape{wallHalfExtentsX}, PhysicsSystem::RigidState(Vector3r(-wallOffset, 0, wallCenterZ)));
        }

        // Cooked spaghetti properties.
        const real_t L  = (real_t)0.22;    // ~22 cm length
        const real_t r  = (real_t)0.0015;  // 1.5 mm radius
        const real_t d  = (real_t)(2 * r);
        const real_t E  = (real_t)1e4;   
        const real_t nu = (real_t)0.30;
        const real_t rho = (real_t)1050.0; // kg/m^3 (water-ish density)

        // Capsule-shaped beam bodies.
        PhysicsSystem::BeamCrossSection section(d, d, PhysicsSystem::BeamBodyType::Capsule);
        auto springs = PhysicsSystem::BeamSpringParams::fromMaterial(E, nu);
        springs.setDampingFromFactor(1.0);

        // Grid layout above the bowl; discard points outside inner radius.
        const real_t bowlInnerRadius = (real_t)0.075; // 7.5 cm inner radius
        const real_t zTop            = (real_t)0.30;  // height of the clamped ends
        const real_t spacing         = d * (real_t)3; // loose packing
        const int    halfCount       = 10;
        const size_t segments = 64;

        PhysicsSystem::RigidProps props = PhysicsSystem::RigidProps::withDensity(rho);

        for (int ix = -halfCount; ix <= halfCount; ++ix) {
            for (int iy = -halfCount; iy <= halfCount; ++iy) {


                Vector3r randomVec = Vector3r::Random();
                randomVec.normalize();

                PhysicsSystem::RigidState stateDefaults(Vector3r::Zero(), randomVec * 0.001, Quaternion4r::Identity());

                const real_t x = (real_t)ix * spacing;
                const real_t y = (real_t)iy * spacing;
                const real_t rXY = std::sqrt(x * x + y * y);
                if (rXY > bowlInnerRadius) {
                    continue; // outside the bowl opening
                }

                const Vector3r pTop(x, y, zTop);
                const Vector3r pBot = pTop + Vector3r(0, 0, -L);
                LinearSpline spline(pTop, pBot);

                auto ends = cardillo::physics::BodyFactory::createBeam(sys, spline, section, springs, stateDefaults, props, segments);
                m_spaghettiEnds.push_back(ends);
            }
        }
    }

    void updateScene(cardillo::PhysicsSystem& /*sys*/, real_t /*t*/, real_t /*dt*/) override {
        // Passive scene: gravity causes the spaghetti to drape over the bowl.
    }

private:
    std::vector<std::pair<entt::entity, entt::entity>> m_spaghettiEnds;
};
