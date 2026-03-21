#pragma once

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <stdexcept>
// Only need public Contact type here
#include "types.hpp"

// Include COAL headers for complete types (we use unique_ptr to these types)
#include <coal/broadphase/broadphase.h>
#include <coal/collision_object.h>
#include <coal/mesh_loader/loader.h>
#include <coal/BVH/BVH_model.h>

// Forward declare World in the correct namespace
namespace cardillo { class World; }

namespace cardillo::collision {

class CollisionCoal {
public:
    CollisionCoal() = default;
    ~CollisionCoal();

    // Attach to a World (not owned)
    void registerSystem(const cardillo::World* sys) { m_sys = sys; }

    // (Re)build the COAL scene from ECS collidables (creates geometries & objects and registers them in broadphase)
    void rebuild();

    // Update per-object transforms from World state and notify broadphase
    void applyTransforms();

    // Run broadphase + narrowphase and return contacts
    std::vector<Contact> detectAll() const;

    // Access the last flattened contact vector produced by the most recent detectAll() call.
    // Returns an empty vector if detectAll() has not been run yet.
    const std::vector<Contact>& lastFlattenedContacts() const { return m_last_flattened; }

    // Clear all internal caches/objects
    void clear();

    // Disable/enable collisions between a specific pair of entities (order-independent)
    void disablePair(entt::entity a, entt::entity b);
    void enablePair(entt::entity a, entt::entity b);
    bool isPairDisabled(entt::entity a, entt::entity b) const;

private:
    enum class ColliderKind { Box, Sphere, Halfspace, Mesh, HeightField, Capsule, Cylinder, Cone };

    // Helpers
    void ensureBroadphaseFromConfig_();
    ColliderKind inferKind_(entt::entity e) const;
    std::shared_ptr<coal::CollisionGeometry> makeGeometryFor_(ColliderKind kind, entt::entity e) const;
    coal::Transform3s makeTransformFromQ_(const VectorXr& q) const;

    // Backrefs
    const cardillo::World* m_sys = nullptr; // not owned

    // COAL scene storage
    std::unique_ptr<coal::BroadPhaseCollisionManager> m_broadphase; // manager chosen from config
    std::vector<std::shared_ptr<coal::CollisionGeometry>> m_geoms;   // one per object
    std::vector<std::unique_ptr<coal::CollisionObject>> m_objects;   // one per object
    std::vector<std::shared_ptr<coal::CollisionGeometry>> m_meshSphereGeoms; // mesh precheck spheres (optional)
    std::vector<std::unique_ptr<coal::CollisionObject>> m_meshSphereObjects; // mesh precheck spheres (optional)
    std::vector<entt::entity> m_entities;                            // index -> entity
    std::unordered_map<uint32_t, std::size_t> m_index_from_entity;   // entity id -> index
    std::vector<ColliderKind> m_kinds;                                // index -> kind

    // Last generation grouped contacts (for potential warmstarting)
    mutable ContactMap m_prevContactMap; // updated at end of detectAll
    mutable std::vector<Contact> m_last_flattened; // flattened contacts from last detectAll()

    // Pairs to skip in collision (canonicalized ContactPairKey)
    std::unordered_set<ContactPairKey, ContactPairKeyHash> m_disabledPairs;
    mutable uint64_t m_seen_structure_version = 0;
};

} // namespace cardillo::collision
