#pragma once

#include <vector>
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
    MatrixXXr wA;         // contact Force directions for body A     (For point masses, this is just -normal)
    MatrixXXr wB;         // contact Force directions for body B     (For point masses, this is just normal)
    real_t penetration;   // overlap distance (> 0 means interpenetration)
};

// Collider shapes resolved from ECS
struct SphereCollider { entt::entity e; Vector3r center; real_t radius; };

struct PlaneCollider {
    entt::entity e;
    Vector3r point;   // plane origin/center
    Vector3r normal;  // unit normal
    Vector3r u;       // unit in-plane axis U (orthonormal to normal)
    Vector3r v;       // unit in-plane axis V (orthonormal to normal and U)
    real_t hx;        // half extent along U
    real_t hy;        // half extent along V
};

struct AabbCollider { entt::entity e; Vector3r center; Vector3r halfExtents; };

struct ObbCollider { entt::entity e; Vector3r center; Vector3r halfExtents; Matrix33r R; };

}
