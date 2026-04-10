#pragma once

#include <utility>
#include <vector>

#include <entt/entt.hpp>

#include "../world.hpp"

namespace cardillo::physics {

class BeamFactory {
   public:
    static std::pair<entt::entity, entt::entity> createBeam(World& system, const misc::SplinePattern& spline, const BeamCrossSection& section, const BeamSpringParams& springs,
                                                            const RigidState& stateDefaults, const RigidProps& propsDefaults, size_t segments,
                                                            cardillo::collision::CollisionCoal* collision_mgr = nullptr);

    static std::pair<entt::entity, entt::entity> createBeams(World& system, const std::vector<const misc::SplinePattern*>& splines, const BeamCrossSection& section, const BeamSpringParams& springs,
                                                             const RigidState& stateDefaults, const RigidProps& propsDefaults, size_t segments,
                                                             cardillo::collision::CollisionCoal* collision_mgr = nullptr);
};

}  // namespace cardillo::physics
