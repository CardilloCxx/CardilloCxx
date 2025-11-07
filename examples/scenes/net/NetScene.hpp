#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <vector>
#include "physics/constraints.hpp"

using namespace cardillo;

// NetScene: creates a static square frame (4 elongated cubes), a rectangular grid
// of point masses connected by springs, and drops a sphere from above onto the net.
// The net stiffness and geometry are configurable through the main scene config
// (keys: net.k, net.rows, net.cols, net.spacing, net.node_mass, net.sphere_*).
class NetScene : public SceneBase {
public:
    NetScene() = default;
    ~NetScene() override = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;
        using namespace cardillo::physics;

        const auto& cfg = sys.config();
        const real_t k = 10000.0;
        const int rows = 20;
        const int cols = 20;
        // Keep the frame size fixed and derive spacing from rows/cols so the mesh scales to fit.
        const real_t frameSide = (real_t)1.0; // total frame side length (meters)
        const real_t spacing = frameSide / (real_t)(std::max(2, cols) - 1);
        const real_t nodeMass = 0.001;

        // Frame geometry: square centered at origin on the XY-plane (z = 0)
        const real_t side = frameSide;
        const real_t halfSide = side * (real_t)0.5;
        const real_t frame_thickness = std::max<real_t>(spacing * 0.2, (real_t)0.01);
        const real_t frame_height = std::max<real_t>(spacing * 0.1, (real_t)0.01);

        // All frame cubes placed at z = 0 (slightly below the net plane)
        const real_t frame_z = (real_t)0.0;

    PhysicsSystem::CubeShape edgeX{Vector3r(halfSide, frame_thickness * (real_t)0.5, frame_height * (real_t)0.5)};
    PhysicsSystem::CubeShape edgeY{Vector3r(frame_thickness * (real_t)0.5, halfSide, frame_height * (real_t)0.5)};

        // Create four edges as obstacle bodies (top, bottom, left, right)
    PhysicsSystem::RigidProps staticProps; PhysicsSystem::RigidState st; st.orientation = Quaternion4r::Identity();
    st.position = Vector3r(0.0, halfSide, frame_z); entt::entity eTop = sys.addRigidBody(edgeX, st, staticProps);

    st.position = Vector3r(0.0, -halfSide, frame_z); entt::entity eBottom = sys.addRigidBody(edgeX, st, staticProps);

    st.position = Vector3r(halfSide, 0.0, frame_z); entt::entity eRight = sys.addRigidBody(edgeY, st, staticProps);

    st.position = Vector3r(-halfSide, 0.0, frame_z); entt::entity eLeft = sys.addRigidBody(edgeY, st, staticProps);

        std::vector<entt::entity> frameEdges = { eTop, eBottom, eLeft, eRight };

        // No global cross-bracing across the frame; diagonal springs follow the interior pattern

        // Net plane z coordinate slightly above frame (so springs attach to frame faces)
        const real_t net_z = frame_z + frame_height * (real_t)0.5 + spacing * (real_t)0.01;

        // Create grid of point masses centered in the frame
        std::vector<std::vector<entt::entity>> nodes;
        nodes.resize(rows);
        for (int i = 0; i < rows; ++i) nodes[i].resize(cols, entt::entity{entt::null});

        const Vector3r origin = Vector3r(-halfSide, -halfSide, net_z);
        // Use flat cube nodes (wider than before) so springs can attach to corners/edges
        const real_t node_hx = spacing * (real_t)0.15; // half-extent x
        const real_t node_hy = spacing * (real_t)0.15; // half-extent y
        const real_t node_hz = std::max<real_t>(spacing * (real_t)0.06, (real_t)0.001); // thin

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                // Only create physical rigid-cube nodes for interior nodes (not the boundary)
                if (i == 0 || j == 0 || i == rows - 1 || j == cols - 1) {
                    // leave border entries as null; springs will connect interior nodes to frame directly
                    nodes[i][j] = entt::entity{entt::null};
                    continue;
                }
                Vector3r p = origin + Vector3r((real_t)j * spacing, (real_t)i * spacing, 0.0);
                PhysicsSystem::CubeShape nodeShape{Vector3r(node_hx, node_hy, node_hz)};
                PhysicsSystem::RigidState nodeState; nodeState.position = p; nodeState.orientation = Quaternion4r::Identity();
                PhysicsSystem::RigidProps nodeProps; nodeProps.mass = nodeMass;
                entt::entity eNode = sys.addRigidBody(nodeShape, nodeState, nodeProps);
                nodes[i][j] = eNode;
            }
        }

        // Frame top surface z-coordinate for spring attachments
        const real_t frame_surface_z = frame_z + frame_height * (real_t)0.5;

        // small helper to add a corner-to-corner spring between two entities using local corner offsets
        auto addCornerSpring = [&](entt::entity A, entt::entity B, const Vector3r& rA_local, const Vector3r& rB_local, real_t stiffMul) {
            sys.addConstraint<cardillo::physics::LinearDistanceConstraint>(sys.ecs(), A, B, rA_local, rB_local, k * stiffMul, (real_t)0.1);
        };

        // Connect interior nodes using corner-to-corner pattern. When a neighbor lies on the border,
        // connect to the corresponding anchor corner instead so the pattern is preserved at the frame.
        for (int i = 1; i < rows - 1; ++i) {
            for (int j = 1; j < cols - 1; ++j) {
                entt::entity a = nodes[i][j];
                Vector3r a_tl = Vector3r(-node_hx, node_hy, 0.0);
                Vector3r a_tr = Vector3r(node_hx, node_hy, 0.0);
                Vector3r a_bl = Vector3r(-node_hx, -node_hy, 0.0);
                Vector3r a_br = Vector3r(node_hx, -node_hy, 0.0);

                // right neighbor / attach to right frame edge
                if (j + 1 < cols - 1) {
                    entt::entity b = nodes[i][j+1];
                    addCornerSpring(a, b, a_tr, Vector3r(-node_hx, node_hy, 0.0), 1.0);
                    addCornerSpring(a, b, a_br, Vector3r(-node_hx, -node_hy, 0.0), 1.0);
                } else {
                    // attach to right frame edge at positions aligned with a's top-right and bottom-right corners
                    entt::registry& reg = sys.ecs();
                    Vector3r edgeCenter = reg.get<PhysicsSystem::C_Position3>(eRight).value;
                    Quaternion4r qEdge = reg.get<PhysicsSystem::C_Orientation>(eRight).q;
                    Vector3r aCenter = reg.get<PhysicsSystem::C_Position3>(a).value;
                    Vector3r p_tr = Vector3r( halfSide, aCenter.y() + node_hy, frame_surface_z);
                    Vector3r p_br = Vector3r( halfSide, aCenter.y() - node_hy, frame_surface_z);
                    Vector3r r_edge_tr = qEdge.conjugate() * (p_tr - edgeCenter);
                    Vector3r r_edge_br = qEdge.conjugate() * (p_br - edgeCenter);
                    addCornerSpring(a, eRight, a_tr, r_edge_tr, 1.0);
                    addCornerSpring(a, eRight, a_br, r_edge_br, 1.0);
                }

                // up neighbor / anchor
                if (i + 1 < rows - 1) {
                    entt::entity b = nodes[i+1][j];
                    addCornerSpring(a, b, a_tl, Vector3r(-node_hx, -node_hy, 0.0), 1.0);
                    addCornerSpring(a, b, a_tr, Vector3r(node_hx, -node_hy, 0.0), 1.0);
                } else {
                    // attach to top frame edge at positions aligned with a's top-left and top-right corners
                    entt::registry& reg = sys.ecs();
                    Vector3r edgeCenter = reg.get<PhysicsSystem::C_Position3>(eTop).value;
                    Quaternion4r qEdge = reg.get<PhysicsSystem::C_Orientation>(eTop).q;
                    Vector3r aCenter = reg.get<PhysicsSystem::C_Position3>(a).value;
                    Vector3r p_tl = Vector3r( aCenter.x() - node_hx,  halfSide, frame_surface_z);
                    Vector3r p_tr = Vector3r( aCenter.x() + node_hx,  halfSide, frame_surface_z);
                    Vector3r r_edge_tl = qEdge.conjugate() * (p_tl - edgeCenter);
                    Vector3r r_edge_tr = qEdge.conjugate() * (p_tr - edgeCenter);
                    addCornerSpring(a, eTop, a_tl, r_edge_tl, 1.0);
                    addCornerSpring(a, eTop, a_tr, r_edge_tr, 1.0);
                }

                // down-right diagonal (i+1, j+1)
                if (i + 1 < rows - 1 && j + 1 < cols - 1) {
                    entt::entity b = nodes[i+1][j+1];
                    // connect a.top-right -> b.bottom-left
                    addCornerSpring(a, b, a_tr, Vector3r(-node_hx, -node_hy, 0.0), (real_t)0.7);
                } else if (i + 1 == rows - 1 || j + 1 == cols - 1) {
                    // do not attach diagonal springs to the frame
                }

                // down-left diagonal (i+1, j-1)
                if (i + 1 < rows - 1 && j - 1 > 0) {
                    entt::entity b = nodes[i+1][j-1];
                    // connect a.top-left -> b.bottom-right
                    addCornerSpring(a, b, a_tl, Vector3r(node_hx, -node_hy, 0.0), (real_t)0.7);
                } else if (i + 1 == rows - 1 || j - 1 == 0) {
                    // do not attach diagonal springs to the frame
                }

                // left border attachments (symmetric to right)
                if (j - 1 == 0) {
                    entt::registry& reg = sys.ecs();
                    Vector3r edgeCenter = reg.get<PhysicsSystem::C_Position3>(eLeft).value;
                    Quaternion4r qEdge = reg.get<PhysicsSystem::C_Orientation>(eLeft).q;
                    Vector3r aCenter = reg.get<PhysicsSystem::C_Position3>(a).value;
                    Vector3r p_tl = Vector3r( -halfSide, aCenter.y() + node_hy, frame_surface_z);
                    Vector3r p_bl = Vector3r( -halfSide, aCenter.y() - node_hy, frame_surface_z);
                    Vector3r r_edge_tl = qEdge.conjugate() * (p_tl - edgeCenter);
                    Vector3r r_edge_bl = qEdge.conjugate() * (p_bl - edgeCenter);
                    addCornerSpring(a, eLeft, a_tl, r_edge_tl, 1.0);
                    addCornerSpring(a, eLeft, a_bl, r_edge_bl, 1.0);
                }

                // bottom border attachments (symmetric to top)
                if (i - 1 == 0) {
                    entt::registry& reg = sys.ecs();
                    Vector3r edgeCenter = reg.get<PhysicsSystem::C_Position3>(eBottom).value;
                    Quaternion4r qEdge = reg.get<PhysicsSystem::C_Orientation>(eBottom).q;
                    Vector3r aCenter = reg.get<PhysicsSystem::C_Position3>(a).value;
                    Vector3r p_bl = Vector3r( aCenter.x() - node_hx, -halfSide, frame_surface_z);
                    Vector3r p_br = Vector3r( aCenter.x() + node_hx, -halfSide, frame_surface_z);
                    Vector3r r_edge_bl = qEdge.conjugate() * (p_bl - edgeCenter);
                    Vector3r r_edge_br = qEdge.conjugate() * (p_br - edgeCenter);
                    addCornerSpring(a, eBottom, a_bl, r_edge_bl, 1.0);
                    addCornerSpring(a, eBottom, a_br, r_edge_br, 1.0);
                }
            }
        }

    //     // Drop a bolder from above
    //     const std::string bolderPath = "res/meshes/rock.obj";
    //     const real_t bolderMass = (real_t)5.0;
    //     const real_t bolderHeight = (real_t)1.0;
    //     Vector3r bolderPos = Vector3r(0.0, 0.0, bolderHeight);
    //     entt::entity eBolder = sys.addRigidBodyMesh(
    //         bolderMass, 
    //         bolderPos, 
    //         Quaternion4r::Identity(),
    //         Vector3r::Zero(),
    //         Vector3r::Zero(), 
    //         bolderPath, Vector3r(0.1, 0.1, 0.1));
    // }

        // Drop a sphere from above
        const real_t sphereRadius = (real_t)0.05;
        const real_t sphereMass = (real_t)0.5;
        const real_t sphereHeight = (real_t)0.5;
        Vector3r spherePos = Vector3r(0.0, 0.0, sphereHeight);
        {
            PhysicsSystem::SphereShape s{sphereRadius}; PhysicsSystem::RigidState st2; st2.position = spherePos; st2.orientation = Quaternion4r::Identity(); PhysicsSystem::RigidProps pr2; pr2.mass = sphereMass; (void)sys.addRigidBody(s, st2, pr2);
        }
    }
};
