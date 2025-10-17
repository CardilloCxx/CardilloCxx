#pragma once

#include <vector>
#include "types.hpp"

namespace cardillo::collision {

class NarrowPhase {
public:
    // Accumulate contacts from provided typed colliders and index pairs
    void sphereSphere(const std::vector<SphereCollider>& spheres,
                      const std::vector<std::pair<int,int>>& pairs,
                      std::vector<Contact>& out) const;

    void spherePlane(const std::vector<SphereCollider>& spheres,
                     const std::vector<PlaneCollider>& planes,
                     const std::vector<std::pair<int,int>>& pairs,
                     std::vector<Contact>& out) const;

    void sphereObb(const std::vector<SphereCollider>& spheres,
                   const std::vector<ObbCollider>& obbs,
                   const std::vector<std::pair<int,int>>& pairs,
                   std::vector<Contact>& out) const;

    void obbObb(const std::vector<ObbCollider>& obbs,
                const std::vector<std::pair<int,int>>& pairs,
                std::vector<Contact>& out) const;
};

}
