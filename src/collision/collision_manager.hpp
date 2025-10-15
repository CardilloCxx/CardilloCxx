#pragma once

#include <vector>
#include "types.hpp"
#include "broad_phase.hpp"
#include "narrow_phase.hpp"

namespace cardillo::collision {

class CollisionManager {
public:
    std::vector<Contact> detectAll(const PhysicsSystem& sys) const;
private:
    BroadPhase m_broad;
    NarrowPhase m_narrow;
};

}
