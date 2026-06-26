#pragma once

#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>
// Only need public Contact type here
#include "contact_tracker.hpp"
#include "types.hpp"

// Include COAL headers for complete types (we use unique_ptr to these types)
#include <coal/BVH/BVH_model.h>
#include <coal/broadphase/broadphase.h>
#include <coal/collision_object.h>
#include <coal/mesh_loader/loader.h>

// Forward declare World and DynamicsAssembler in the correct namespaces
namespace cardillo {
class World;
}
namespace cardillo {
namespace misc {
class TimingManager;
}
}  // namespace cardillo
namespace cardillo {
namespace physics {
class DynamicsAssembler;
}
}  // namespace cardillo

namespace cardillo::collision {

class CollisionCoal {
   private:
    enum class ColliderKind { Box, Sphere, Halfspace, Mesh, Capsule, Cylinder, Cone };

    // Helpers
    void ensureBroadphaseFromConfig_();
    ColliderKind inferKind_(entt::entity e) const;
    std::shared_ptr<coal::CollisionGeometry> makeGeometryFor_(ColliderKind kind, entt::entity e) const;
    coal::Transform3s makeTransformFromQ_(const VectorXr& q) const;

    // Backrefs
    cardillo::World* m_world = nullptr;                  // not owned
    cardillo::misc::TimingManager* m_timings = nullptr;  // optional timings pointer
    cardillo::config::Config& m_cfg;                     // reference to global config for easy access

   public:
    CollisionCoal(cardillo::World& world, cardillo::misc::TimingManager* timings, cardillo::config::Config& cfg);
    ~CollisionCoal();

    // (Re)build the COAL scene from ECS collidables (creates geometries & objects and registers
    // them in broadphase)
    void rebuild();

    // Update per-object transforms from World state and notify broadphase
    void applyTransforms();

    // Run broadphase + narrowphase and return a reference to the authoritative contact buffer
    std::vector<Contact>& detectAll();

    std::vector<Contact> m_prev_flattened;  // contacts from previous step
    std::vector<Contact> m_flattened;       // contacts from current step (authoritative)

    // Clear all internal caches/objects
    void clear();

    // Disable/enable collisions between a specific pair of entities (order-independent)
    void disablePair(entt::entity a, entt::entity b);
    void enablePair(entt::entity a, entt::entity b);
    bool isPairDisabled(entt::entity a, entt::entity b) const;

    // COAL scene storage
    std::unique_ptr<coal::BroadPhaseCollisionManager> m_broadphase;           // manager chosen from config
    std::vector<std::shared_ptr<coal::CollisionGeometry>> m_geoms;            // one per object
    std::vector<std::unique_ptr<coal::CollisionObject>> m_objects;            // one per object
    std::vector<std::shared_ptr<coal::CollisionGeometry>> m_meshSphereGeoms;  // mesh precheck spheres (optional)
    std::vector<std::unique_ptr<coal::CollisionObject>> m_meshSphereObjects;  // mesh precheck spheres (optional)
    std::vector<entt::entity> m_entities;                                     // index -> entity
    std::unordered_map<uint32_t, std::size_t> m_index_from_entity;            // entity id -> index
    std::vector<ColliderKind> m_kinds;                                        // index -> kind

    // Contact Tracker
    ContactTracker m_contactTracker;

    // Pairs to skip in collision (canonicalized ContactPairKey)
    std::unordered_set<ContactPairKey, ContactPairKeyHash> m_disabledPairs;
    mutable uint64_t m_seen_structure_version = 0;
};

}  // namespace cardillo::collision
