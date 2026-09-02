#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

using namespace cardillo;

class ChuteScene : public SceneBase {
public:
    const char* sceneName() const override { return "chute"; }
    ChuteScene() = default;
    ~ChuteScene() = default;

    void populate(physics::PhysicsEngine& engine) override {
        using namespace cardillo;
        using namespace misc;
        m_elapsedSinceSpawn = 0.0;

        struct ChuteSpec {
            const char* name;
            real_t x, y, z;
            real_t rxDeg, ryDeg, rzDeg;
            real_t sx, sy, sz;
        };

        const std::vector<ChuteSpec> chuteSpecs = {
            {"Cube",     (real_t)0.155224,  (real_t)0.0,       (real_t)-0.372766, (real_t)0.0,        (real_t)0.0,        (real_t)0.0,       (real_t)0.230894,  (real_t)0.230894,  (real_t)0.004592},
            {"Cube.001", (real_t)0.155224,  (real_t)0.229927,  (real_t)-0.14443,  (real_t)90.000003,  (real_t)0.0,        (real_t)0.0,       (real_t)0.230894,  (real_t)0.230894,  (real_t)0.004592},
            {"Cube.002", (real_t)-0.075499, (real_t)0.0,       (real_t)-0.14443,  (real_t)90.000003,  (real_t)0.0,        (real_t)90.000003, (real_t)0.230894,  (real_t)0.230894,  (real_t)0.004592},
            {"Cube.004", (real_t)0.010486,  (real_t)0.0,       (real_t)-0.099402, (real_t)32.253246,  (real_t)0.000002,   (real_t)90.000003, (real_t)0.063601,  (real_t)0.098282,  (real_t)0.001955},
            {"Cube.005", (real_t)0.059472,  (real_t)0.000169,  (real_t)-0.31162,  (real_t)0.0,        (real_t)0.0,        (real_t)0.0,       (real_t)0.067962,  (real_t)0.124295,  (real_t)0.001907},
            {"Cube.006", (real_t)0.203497,  (real_t)0.000169,  (real_t)-0.268026, (real_t)0.0,        (real_t)-30.701526, (real_t)0.0,       (real_t)0.095871,  (real_t)0.124295,  (real_t)0.001907},
            {"Cube.003", (real_t)0.01624,   (real_t)-0.06145,  (real_t)-0.090283, (real_t)32.253243,  (real_t)0.000002,   (real_t)90.000003, (real_t)0.002116,  (real_t)0.098282,  (real_t)0.010967},
            {"Cube.007", (real_t)0.01624,   (real_t)0.064016,  (real_t)-0.090283, (real_t)32.253243,  (real_t)0.000002,   (real_t)90.000003, (real_t)0.002116,  (real_t)0.098282,  (real_t)0.010967},
            {"Cube.008", (real_t)0.059472,  (real_t)-0.123853, (real_t)-0.302663, (real_t)180.000005, (real_t)0.0,        (real_t)0.0,       (real_t)0.067962,  (real_t)0.001216,  (real_t)0.010859},
            {"Cube.009", (real_t)0.059472,  (real_t)0.125495,  (real_t)-0.302663, (real_t)180.000005, (real_t)0.0,        (real_t)0.0,       (real_t)0.067962,  (real_t)0.001216,  (real_t)0.010859},
            {"Cube.010", (real_t)0.197469,  (real_t)-0.123764, (real_t)-0.257874, (real_t)0.0,        (real_t)-30.701526, (real_t)0.0,       (real_t)0.095871,  (real_t)0.001392,  (real_t)0.012183},
            {"Cube.011", (real_t)0.197469,  (real_t)0.125381,  (real_t)-0.257874, (real_t)0.0,        (real_t)-30.701526, (real_t)0.0,       (real_t)0.095871,  (real_t)0.001392,  (real_t)0.012183},
            {"Cube.012", (real_t)-0.008007, (real_t)0.000169,  (real_t)-0.301463, (real_t)180.000005, (real_t)0.0,        (real_t)0.0,       (real_t)-0.002256, (real_t)-0.124295, (real_t)-0.010687},
            {"Cube.013", (real_t)0.281449,  (real_t)0.000169,  (real_t)-0.208213, (real_t)180.000005, (real_t)-30.701526, (real_t)0.0,       (real_t)-0.00253,  (real_t)-0.124295, (real_t)-0.011757},
        };

        auto degToRad = [](real_t deg) { return deg * (real_t)M_PI / (real_t)180.0; };

        real_t maxZ = -std::numeric_limits<real_t>::max();
        real_t maxRadius = (real_t)0.0;

        for (const auto& o : chuteSpecs) {
            const Vector3r position(o.x, o.y, o.z);

            const Quaternion4r qx(Eigen::AngleAxis<real_t>(degToRad(o.rxDeg), Vector3r::UnitX()));
            const Quaternion4r qy(Eigen::AngleAxis<real_t>(degToRad(o.ryDeg), Vector3r::UnitY()));
            const Quaternion4r qz(Eigen::AngleAxis<real_t>(degToRad(o.rzDeg), Vector3r::UnitZ()));
            const Quaternion4r orientation = qz * qy * qx;

            const Vector3r extents(std::abs(o.sx), std::abs(o.sy), std::abs(o.sz));

            physics::CubeShape shape(extents);
            physics::RigidState state{position, orientation};
            engine.addStaticBody(shape, state);

            maxZ = std::max(maxZ, o.z);
            maxRadius = std::max(maxRadius, std::sqrt(o.x * o.x + o.y * o.y));
        }

        spawnRandomShapesBatch(engine);
    }

    void updateScene(physics::PhysicsEngine& engine, real_t t, real_t dt) override {
        m_elapsedSinceSpawn += dt;
        if (m_elapsedSinceSpawn < m_spawnInterval) return;

        m_elapsedSinceSpawn = 0.0;

        if(t < 1.75) spawnRandomShapesBatch(engine);
    }

private:
    struct SpawnCandidate {
        Vector3r pos;
        real_t radius;
    };

    void spawnRandomShapesBatch(physics::PhysicsEngine& engine) {
        using namespace cardillo;

        std::uniform_real_distribution<real_t> angleDist((real_t)0.0, (real_t)(2.0 * M_PI));
        std::uniform_real_distribution<real_t> radiusFrac((real_t)0.0, (real_t)1.0);

        std::uniform_real_distribution<real_t> baseScaleDist((real_t)0.0008, (real_t)0.0018); 
        
        // Variation per dimension: allows shapes to be stretched/squashed around the base size
        // 1.0 means perfect symmetry. E.g., [0.6, 1.4] allows significant independent variation.
        std::uniform_real_distribution<real_t> dimVariation((real_t)0.6, (real_t)1.4);
        
        std::uniform_real_distribution<real_t> fallSpeedDist((real_t)0.2, (real_t)0.6);
        std::uniform_real_distribution<real_t> lateralSpeedDist((real_t)-0.05, (real_t)0.05);
        std::uniform_real_distribution<real_t> angVelDist((real_t)-3.0, (real_t)3.0);
        std::uniform_real_distribution<real_t> axisComponent((real_t)-1.0, (real_t)1.0);
        
        // 0=Cube, 1=Sphere, 2=Cylinder, 3=Cone
        std::uniform_int_distribution<int> shapeDist(0, 3);

        std::vector<SpawnCandidate> spawnedBatch;
        spawnedBatch.reserve(m_objectsPerSpawn);

        const int maxAttemptsPerBatch = m_objectsPerSpawn * 150;
        int attempts = 0;

        while (spawnedBatch.size() < static_cast<size_t>(m_objectsPerSpawn) && attempts < maxAttemptsPerBatch) {
            attempts++;

            const real_t baseScale = baseScaleDist(m_rng);
            
            // Calculate size per dimension independently
            const real_t sx = baseScale * dimVariation(m_rng);
            const real_t sy = baseScale * dimVariation(m_rng);
            const real_t sz = baseScale * dimVariation(m_rng);
            const Vector3r extents(sx, sy, sz);
            
            // Approximate bounding sphere radius for collision avoidance during spawn.
            const real_t approxRadius = extents.maxCoeff() * (real_t)1.8;

            const real_t theta = angleDist(m_rng);
            const real_t r = m_spawnCircleRadius * std::sqrt(radiusFrac(m_rng));
            
            const Vector3r candidatePos(r * std::cos(theta), r * std::sin(theta), m_spawnHeight);

            bool overlaps = false;
            for (const auto& existing : spawnedBatch) {
                const real_t minDistance = approxRadius + existing.radius;
                const real_t dx = candidatePos.x() - existing.pos.x();
                const real_t dy = candidatePos.y() - existing.pos.y();

                if ((dx * dx + dy * dy) < minDistance * minDistance) {
                    overlaps = true;
                    break;
                }
            }

            if (overlaps) continue;

            spawnedBatch.push_back({candidatePos, approxRadius});

            Vector3r axis(axisComponent(m_rng), axisComponent(m_rng), axisComponent(m_rng));
            if (axis.norm() < (real_t)1e-6) axis = Vector3r::UnitZ();
            axis.normalize();
            const Quaternion4r orientation(Eigen::AngleAxis<real_t>(angleDist(m_rng), axis));

            const Vector3r radial(candidatePos.x(), candidatePos.y(), (real_t)0.0);
            Vector3r towardCenter = Vector3r::Zero();
            if (radial.norm() > (real_t)1e-6) towardCenter = -radial.normalized() * fallSpeedDist(m_rng) * (real_t)0.3;
            const Vector3r velocity = Vector3r(lateralSpeedDist(m_rng), lateralSpeedDist(m_rng), -fallSpeedDist(m_rng)) + towardCenter;
            const Vector3r angularVelocity(angVelDist(m_rng), angVelDist(m_rng), angVelDist(m_rng));

            physics::RigidState state(candidatePos, velocity, orientation, angularVelocity);
            physics::RigidProps props = physics::RigidProps::withDensity((real_t)2500.0);

            int shapeType = shapeDist(m_rng);
            
            switch (shapeType) {
                case 0: { // Cube
                    physics::CubeShape shape(extents);
                    engine.addRigidBody(shape, state, props);
                    break;
                }
                case 1: { // Sphere (Averaging dimensions for radius)
                    real_t radius = (sx + sy + sz) / (real_t)3.0;
                    physics::SphereShape shape(radius);
                    engine.addRigidBody(shape, state, props);
                    break;
                }
                case 2: { // Cylinder
                    real_t radius = (sx + sy) / (real_t)2.0; 
                    real_t height = sz; // Removed * 2.0 (engine interprets as half-height)
                    physics::CylinderShape shape(radius, height);
                    engine.addRigidBody(shape, state, props);
                    break;
                }
                case 3: { // Cone
                    real_t radius = (sx + sy) / (real_t)2.0;
                    real_t height = sz; // Removed * 2.0 (engine interprets as half-height)
                    physics::ConeShape shape(radius, height);
                    engine.addRigidBody(shape, state, props);
                    break;
                }
            }
        }

        std::cout << "Spawned " << spawnedBatch.size() << " objects in " << attempts << " attempts.\n";
    }

    std::mt19937 m_rng{1337};
    real_t m_spawnCircleRadius{0.040};
    real_t m_spawnHeight{0.1};
    real_t m_spawnInterval{0.025};
    real_t m_elapsedSinceSpawn{0.0};
    int m_objectsPerSpawn{100};
};