#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>
#include "physics/force_interaction.hpp"

using namespace cardillo;

class SpringTestScene : public SceneBase {
public:
    SpringTestScene() = default;
    ~SpringTestScene() override = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;

        // ground (static obstacle cube)
        PhysicsSystem::Cube ground;
        ground.center = Vector3r(0.0, 0.0, -0.5);
        ground.halfExtents = Vector3r(5.0, 5.0, 0.5);
        ground.q = Quaternion4r::Identity();
        // addObstacleBody returns an integer id; convert to entt::entity to reference in the spring
        index_t ground_id = sys.addObstacleBody(ground);
        entt::entity eGround = entt::entity(static_cast<uint32_t>(ground_id));

        // Flat cube that acts as a "foot" (dynamic)
        PhysicsSystem::Cube foot; 
        foot.halfExtents = Vector3r(0.2, 0.2, 0.02); // wide and flat
        foot.center = Vector3r(0.0, 0.0, 0.3); 
        foot.q = Quaternion4r::Identity();
        real_t footMass = (real_t)0.5;
        entt::entity eFoot = sys.addRigidBody(footMass, foot.center, foot.q, Vector3r::Zero(), Vector3r::Zero(), foot);

        // Rigid sphere above the foot that will be connected by four corner springs
        real_t sphereRadius = (real_t)0.06;
        Vector3r spherePos = Vector3r(0.0, 0.0, 0.7); // dropped from above
        real_t sphereMass = (real_t)1.0;
        entt::entity eSphere = sys.addRigidBodySphere(sphereMass, spherePos, Quaternion4r::Identity(), Vector3r::Zero(), Vector3r::Zero(), sphereRadius);

        // Create four springs from the top corners of the foot to the sphere center
        const real_t k_corner = (real_t)5e2;
        const real_t d_corner = (real_t)0.0;
        // compute local corner offsets on foot (top face)
        Vector3r he = foot.halfExtents;
        std::vector<Vector3r> corners = {
            Vector3r( he.x(),  he.y(),  he.z()),
            Vector3r(-he.x(),  he.y(),  he.z()),
            Vector3r(-he.x(), -he.y(),  he.z()),
            Vector3r( he.x(), -he.y(),  he.z())
        };
        for (const auto& rA_loc : corners) {
            // attach to sphere center (local offset zero)
            cardillo::physics::SpringConstraint sc(sys.ecs(), eFoot, eSphere, rA_loc, Vector3r::Zero());
            sc.addTranslationalSpring(Vector3r::UnitZ(), k_corner, d_corner);
            sc.addTranslationalSpring(Vector3r::UnitX(), k_corner, d_corner);
            sc.addTranslationalSpring(Vector3r::UnitY(), k_corner, d_corner);
            sys.addConstraint(sc);
        }

        // Double pendulum: One obstacle sphere in the air, two dynamic spheres besides it with infinite stiffness springs
        Vector3r obstaclePos = Vector3r(0.5, 0.0, 1.5);
        real_t obstacleRadius = (real_t)0.05;
        PhysicsSystem::Cube cube;
        cube.center = obstaclePos;
        cube.halfExtents = Vector3r(obstacleRadius, obstacleRadius, obstacleRadius);
        cube.q = Quaternion4r::Identity();
        index_t obs_id = sys.addObstacleBody(cube);
        entt::entity eObs = entt::entity(static_cast<uint32_t>(obs_id));  

        // First dynamic sphere
        Vector3r dyn1Pos = obstaclePos + Vector3r(0.0, -0.4, 0.0);
        real_t dyn1Mass = (real_t)0.3;
        entt::entity eDyn1 = sys.addRigidBodySphere(dyn1Mass, dyn1Pos, Quaternion4r::Identity(), Vector3r(0.0,0.0,0.1), Vector3r::Zero(), obstacleRadius * 0.8);

        // Second dynamic sphere
        Vector3r dyn2Pos = dyn1Pos + Vector3r(0.0, -0.4, 0.0);
        real_t dyn2Mass = (real_t)0.3;
        entt::entity eDyn2 = sys.addRigidBodySphere(dyn2Mass, dyn2Pos, Quaternion4r::Identity(), Vector3r(0.0,0.0,1.0), Vector3r::Zero(), obstacleRadius * 0.8);

        // Create infinite stiffness springs to connect them
        const real_t k_inf = std::numeric_limits<real_t>::max();
        {
            // obs to dyn1
            cardillo::physics::SpringConstraint sc(sys.ecs(), eObs, eDyn1, Vector3r::Zero(), Vector3r::Zero());
            sc.addTranslationalSpring(Vector3r::UnitZ(), k_inf, 0);
            sc.addTranslationalSpring(Vector3r::UnitX(), k_inf, 0);
            sc.addTranslationalSpring(Vector3r::UnitY(), k_inf, 0);
            sys.addConstraint(sc);
        }           
        {
            // dyn1 to dyn2
            cardillo::physics::SpringConstraint sc(sys.ecs(), eDyn1, eDyn2, Vector3r::Zero(), Vector3r::Zero());
            sc.addTranslationalSpring(Vector3r::UnitZ(), k_inf, 0);
            sc.addTranslationalSpring(Vector3r::UnitX(), k_inf, 0);
            sc.addTranslationalSpring(Vector3r::UnitY(), k_inf, 0);
            sys.addConstraint(sc);
        }

        // rope between second anker point and lower dynamic sphere
        // anker point
        Vector3r ankerPos = obstaclePos + Vector3r(0.0, -1.2, 0.0);
        cube.center = ankerPos;
        cube.halfExtents = Vector3r(0.001, 0.001, 0.001);
        index_t anker_id = sys.addObstacleBody(cube);
        entt::entity eAnker = entt::entity(static_cast<uint32_t>(anker_id));
        // create rope
        createRope(sys, eAnker, eDyn2, 50, (real_t)8.0, (real_t)0.001, Vector3r::Zero(), Vector3r( 0.0, 0.0, -obstacleRadius * 0.8));
    }

private:
    // Create a rope between two existing bodies (entities) eA and eB.
    // - numSegments: number of internal segments (positive integer). The rope will have numSegments+1 links between endpoints.
    // - slack: a small non-negative value that controls initial sag/extra length (fraction of straight distance).
    // - segmentMass: mass for each intermediate segment body.
    // This implementation creates small spherical segments and connects consecutive bodies with
    // very large stiffness translational springs (approximate infinite stiffness).
    void createRope(cardillo::PhysicsSystem& sys, entt::entity eA, entt::entity eB,
                    int numSegments, real_t slack = (real_t)0.0, real_t segmentMass = (real_t)0.05,
                    const Vector3r& attachA_local = Vector3r::Zero(),
                    const Vector3r& attachB_local = Vector3r::Zero()) {
        using namespace cardillo;
        using namespace cardillo::physics;
        if (numSegments <= 0) return;

        entt::registry& reg = sys.ecs();

        // Read endpoint body centers and orientations (fall back to origin/identity if not present)
        Vector3r centerA = Vector3r::Zero();
        Vector3r centerB = Vector3r::Zero();
        Quaternion4r qA = Quaternion4r::Identity();
        Quaternion4r qB = Quaternion4r::Identity();
        if (reg.all_of<cardillo::PhysicsSystem::C_Position3>(eA)) centerA = reg.get<cardillo::PhysicsSystem::C_Position3>(eA).value;
        if (reg.all_of<cardillo::PhysicsSystem::C_Position3>(eB)) centerB = reg.get<cardillo::PhysicsSystem::C_Position3>(eB).value;
        if (reg.all_of<cardillo::PhysicsSystem::C_Orientation>(eA)) qA = reg.get<cardillo::PhysicsSystem::C_Orientation>(eA).q;
        if (reg.all_of<cardillo::PhysicsSystem::C_Orientation>(eB)) qB = reg.get<cardillo::PhysicsSystem::C_Orientation>(eB).q;

        // Determine local attachment offsets. For spheres, ensure attachment lies on surface.
        Vector3r rA_local = attachA_local;
        Vector3r rB_local = attachB_local;
        const real_t eps = (real_t)1e-9;
        // If attachment not provided and A is a sphere, point toward B; if provided and A is sphere, project to surface
        if (reg.all_of<cardillo::PhysicsSystem::C_RB_Sphere>(eA) && reg.all_of<cardillo::PhysicsSystem::C_Radius>(eA)) {
            real_t radiusA = reg.get<cardillo::PhysicsSystem::C_Radius>(eA).r;
            if (rA_local.squaredNorm() < eps) {
                // default: direction from center to other body center in world frame
                Vector3r worldDir = (centerB - centerA);
                if (worldDir.norm() < eps) worldDir = Vector3r::UnitZ(); else worldDir.normalize();
                // convert to local frame
                Vector3r localDir = qA.conjugate() * worldDir;
                rA_local = localDir.normalized() * radiusA;
            } else {
                rA_local = rA_local.normalized() * radiusA;
            }
        }
        // Same for B
        if (reg.all_of<cardillo::PhysicsSystem::C_RB_Sphere>(eB) && reg.all_of<cardillo::PhysicsSystem::C_Radius>(eB)) {
            real_t radiusB = reg.get<cardillo::PhysicsSystem::C_Radius>(eB).r;
            if (rB_local.squaredNorm() < eps) {
                Vector3r worldDir = (centerA - centerB);
                if (worldDir.norm() < eps) worldDir = Vector3r::UnitZ(); else worldDir.normalize();
                Vector3r localDir = qB.conjugate() * worldDir;
                rB_local = localDir.normalized() * radiusB;
            } else {
                rB_local = rB_local.normalized() * radiusB;
            }
        }

        // Compute world-space attachment points to use as rope endpoints
        Vector3r posA = centerA + (qA * rA_local);
        Vector3r posB = centerB + (qB * rB_local);

        // Straight-line distance and direction
        Vector3r d = posB - posA;
        real_t L = d.norm();
        Vector3r dir;
        if (L > (real_t)0) dir = d / L; else dir = Vector3r::UnitZ();

        // Construct intermediate node world positions with a parabolic sag and
        // resample so nodes are equally spaced along the curve (arc-length)
        // Heuristic: desired extra length = slack * L. We approximate sag height h proportional to slack.
        real_t h = (real_t)0.25 * slack * L; // heuristic factor

        // Parametric curve p(t) = posA + t * d + (0,0, sag(t)) with sag(t) = -4*h*t*(1-t)
        // We'll sample the curve and compute cumulative arc-length, then pick t values
        // that correspond to equally spaced arc-length positions.
        const int samples = std::max<int>(200, (numSegments + 1) * 8);
        std::vector<real_t> ts(samples + 1);
        std::vector<Vector3r> sample_pos(samples + 1);
        ts[0] = 0.0; sample_pos[0] = posA;
        for (int i = 1; i <= samples; ++i) {
            real_t t = (real_t)i / (real_t)samples;
            ts[i] = t;
            Vector3r p = posA + t * d;
            real_t sag = -4.0 * h * t * (1.0 - t);
            p.z() += sag;
            sample_pos[i] = p;
        }

        // cumulative arc-length
        std::vector<real_t> cum(samples + 1);
        cum[0] = 0.0;
        for (int i = 1; i <= samples; ++i) {
            cum[i] = cum[i-1] + (sample_pos[i] - sample_pos[i-1]).norm();
        }
        real_t total_arc = cum.back();

        // Build nodes: endpoints plus numSegments interior nodes spaced equally by arc-length
        std::vector<Vector3r> nodes; nodes.reserve((size_t)numSegments + 2);
        nodes.push_back(posA);
        for (int k = 1; k <= numSegments; ++k) {
            real_t target_s = (real_t)k * (total_arc / (real_t)(numSegments + 1));
            // find interval j such that cum[j] <= target_s <= cum[j+1]
            int j = 0;
            // binary search
            int lo = 0, hi = samples;
            while (lo <= hi) {
                int mid = (lo + hi) / 2;
                if (cum[mid] <= target_s) { j = mid; lo = mid + 1; }
                else { hi = mid - 1; }
            }
            // interpolate t between j and j+1
            real_t t_interp = ts[j];
            if (j < samples) {
                real_t seg_s = cum[j+1] - cum[j];
                if (seg_s > (real_t)0) {
                    real_t alpha = (target_s - cum[j]) / seg_s;
                    t_interp = ts[j] + alpha * (ts[j+1] - ts[j]);
                } else {
                    t_interp = ts[j];
                }
            }
            // evaluate position at t_interp
            Vector3r p = posA + t_interp * d;
            real_t sag = -4.0 * h * t_interp * (1.0 - t_interp);
            p.z() += sag;
            nodes.push_back(p);
        }
        nodes.push_back(posB);

        // Create rigid bodies for intermediate nodes (skip endpoints which are existing bodies)
        std::vector<entt::entity> bodies; bodies.reserve(nodes.size());
        bodies.push_back(eA);
        
        // approximate segment length for radius heuristic
        real_t segmentLength = (real_t)1.0 / (numSegments + 1) * L;

        const real_t nodeRadius = std::min<real_t>((real_t)0.03, segmentLength * (real_t)0.75);
        for (size_t i = 1; i + 1 < nodes.size(); ++i) {
            const Vector3r& p = nodes[i];
            // add small point-mass segment (visual + collidable radius provided)
            index_t seg_id = sys.addPointMass(segmentMass, p, Vector3r::Zero(), nodeRadius);
            entt::entity seg = entt::entity(static_cast<uint32_t>(seg_id));
            bodies.push_back(seg);
        }
        bodies.push_back(eB);

        // Connect consecutive bodies with very large stiffness translational springs
        const real_t k_inf = std::numeric_limits<real_t>::max();
        const real_t d_inf = (real_t)1.0;
        for (size_t i = 0; i + 1 < bodies.size(); ++i) {
            entt::entity a = bodies[i];
            entt::entity b = bodies[i+1];
            // For the first and last constraint, use the provided local attachment offsets for endpoints
            Vector3r rA = Vector3r::Zero();
            Vector3r rB = Vector3r::Zero();
            if (i == 0) {
                // connect endpoint A (a == eA) to first segment: use rA_local for A
                rA = rA_local;
            }
            if (i + 1 == bodies.size() - 1) {
                // connect last segment to endpoint B (b == eB): use rB_local for B
                rB = rB_local;
            }
            SpringConstraint sc(sys.ecs(), a, b, rA, rB);
            // constrain all translational axes to approximate an inextensible link
            sc.addTranslationalSpring(Vector3r::UnitX(), k_inf, 0);
            sc.addTranslationalSpring(Vector3r::UnitY(), k_inf, 0);
            sc.addTranslationalSpring(Vector3r::UnitZ(), k_inf, 0);
            sys.addConstraint(sc);
        }
    }
};
