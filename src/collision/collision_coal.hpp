#pragma once

#include <vector>
#include <unordered_map>
#include <memory>
#include <stdexcept>
// Only need public Contact type here
#include "types.hpp"

// Include COAL headers for complete types (we use unique_ptr to these types)
#include <coal/broadphase/broadphase.h>
#include <coal/collision_object.h>
#include <coal/mesh_loader/loader.h>
#include <coal/BVH/BVH_model.h>

// Forward declare PhysicsSystem in the correct namespace
namespace cardillo { class PhysicsSystem; }

namespace cardillo::collision {

class CollisionCoal {
public:
    CollisionCoal() = default;
    ~CollisionCoal();

    // Attach to a PhysicsSystem (not owned)
    void registerSystem(const cardillo::PhysicsSystem* sys) { m_sys = sys; }

    // (Re)build the COAL scene from ECS collidables (creates geometries & objects and registers them in broadphase)
    void rebuild();

    // Update per-object transforms from PhysicsSystem state and notify broadphase
    void applyTransforms();

    // Run broadphase + narrowphase and return contacts
    std::vector<Contact> detectAll() const;

    // Clear all internal caches/objects
    void clear();

private:
    enum class ColliderKind { Box, Sphere, Halfspace, Mesh, HeightField };

    // Helpers
    void ensureBroadphaseFromConfig_();
    ColliderKind inferKind_(entt::entity e) const;
    std::shared_ptr<coal::CollisionGeometry> makeGeometryFor_(ColliderKind kind, entt::entity e) const;
    coal::Transform3s makeTransformFromQ_(const VectorXr& q) const;

    // Backrefs
    const cardillo::PhysicsSystem* m_sys = nullptr; // not owned

    // COAL scene storage
    std::unique_ptr<coal::BroadPhaseCollisionManager> m_broadphase; // manager chosen from config
    std::vector<std::shared_ptr<coal::CollisionGeometry>> m_geoms;   // one per object
    std::vector<std::unique_ptr<coal::CollisionObject>> m_objects;   // one per object
    std::vector<entt::entity> m_entities;                            // index -> entity
    std::unordered_map<uint32_t, std::size_t> m_index_from_entity;   // entity id -> index
    std::vector<ColliderKind> m_kinds;                                // index -> kind

    // Last generation grouped contacts (for potential warmstarting)
    mutable ContactMap m_prevContactMap; // updated at end of detectAll
};

} // namespace cardillo::collision
