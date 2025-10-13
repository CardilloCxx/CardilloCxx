#pragma once

#include <vector>
#include <optional>
#include <variant>
#include <entt/entt.hpp>
#include "../misc/types.hpp"
#include "../physics/physics_system.hpp"

namespace cardillo::collision {

// Contact information for a pairwise collision
struct Contact {
    entt::entity a;       // entity id of first colliding object
    entt::entity b;       // entity id of second colliding object
    Vector3r point;       // contact point (on surface or midpoint between)
    Vector3r normal;      // normal pointing from a to b
    real_t penetration;   // overlap distance (> 0 means interpenetration)
};

// Collider shapes resolved from ECS
struct SphereCollider {
    entt::entity e;
    Vector3r center;
    real_t radius;
};

struct PlaneCollider {
    entt::entity e;
    Vector3r point;   // plane origin/center
    Vector3r normal;  // unit normal
    Vector3r u;       // unit in-plane axis U (orthonormal to normal)
    Vector3r v;       // unit in-plane axis V (orthonormal to normal and U)
    real_t hx;        // half extent along U
    real_t hy;        // half extent along V
};

struct AabbCollider {
    entt::entity e;
    Vector3r center;
    Vector3r halfExtents;
};

struct ObbCollider {
    entt::entity e;
    Vector3r center;
    Vector3r halfExtents;
    Matrix33r R; // columns are local axes
};

using Collider = std::variant<SphereCollider, PlaneCollider, AabbCollider, ObbCollider>;

// Naive O(N^2) collision detector over all C_Collidable objects.
// Currently supports:
// - Sphere (point mass with C_Radius) vs sphere
// - Sphere vs finite oriented plane (rectangle) using center/normal/up/sizeX/sizeY
// - Sphere vs axis-aligned cube (AABB)
// - Sphere vs oriented cube (OBB) with rotation R
// Returns a collection of detected contacts.
std::vector<Contact> detectAll(const PhysicsSystem& sys);

}
