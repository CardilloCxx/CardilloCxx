#pragma once

#include <vector>
#include <utility>
#include "types.hpp"
#include "../physics/physics_system.hpp"

namespace cardillo::collision {

struct BroadPhaseData {
    // Typed collider lists collected from the ECS
    std::vector<SphereCollider> spheres;
    std::vector<PlaneCollider> planes;
    std::vector<ObbCollider> obbs;
};

// Unified AABB proxy to abstract all collider kinds for broad-phase
enum class ColliderType : unsigned char { Sphere = 0, Plane = 1, Obb = 2 };

struct AABB { Vector3r min, max; };

struct AabbProxy {
    AABB box;
    ColliderType type;
    int index; // index into BroadPhaseData vector of given type
};

struct Pair { int i, j; }; // indices into proxies vector

class BroadPhase {
public:
    // Collect colliders from ECS into typed arrays
    BroadPhaseData collect(const PhysicsSystem& sys) const;

    // Build proxies from data (one per collider), converting all to AABBs
    std::vector<AabbProxy> buildProxies(const BroadPhaseData& data) const;

    // Compute candidate pairs using 3-axis sweep-and-prune (X primary, filtered by Y,Z)
    std::vector<Pair> makePairs(const std::vector<AabbProxy>& proxies) const;

private:
    // Axis-specific cached endpoints to reuse ordering between frames
    struct Endpoint { real_t value; int proxyIndex; bool start; };
    struct AxisData { std::vector<Endpoint> endpoints; };
    mutable AxisData m_axisX, m_axisY, m_axisZ;

    // Reusable buffers to avoid per-frame allocations during sweep
    mutable std::vector<int> m_active;
    mutable std::vector<int> m_pos;
    // Reuse candidate/filtered pair buffers across calls
    mutable std::vector<Pair> m_candidates;
    mutable std::vector<Pair> m_filtered;

    // aabb converters
    template <typename T> static AABB aabbOf(const T&);
    static AABB aabbOf(const SphereCollider& s);
    static AABB aabbOf(const PlaneCollider& p);
    static AABB aabbOf(const ObbCollider& b);

    // helpers
    static inline bool overlap1D(real_t aMin, real_t aMax, real_t bMin, real_t bMax) {
        return !(aMax < bMin || bMax < aMin);
    }
    static inline bool overlap3D(const AABB& a, const AABB& b) {
        return overlap1D(a.min.x(), a.max.x(), b.min.x(), b.max.x()) &&
               overlap1D(a.min.y(), a.max.y(), b.min.y(), b.max.y()) &&
               overlap1D(a.min.z(), a.max.z(), b.min.z(), b.max.z());
    }

    // Update endpoints from proxies along a given axis, then stable insertion sort
    void updateAndSortAxis(AxisData& axis, const std::vector<AabbProxy>& proxies, int axisIdx) const;

    // Sweep on a single axis filling candidate pairs (by proxy indices) into 'out'
    void sweepAxis(const AxisData& axis, std::vector<Pair>& out) const;
};

}
