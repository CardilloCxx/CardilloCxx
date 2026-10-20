#pragma once

#include <coal/BV/AABB.h>
#include <coal/collision_object.h>
#include <coal/shape/geometric_shapes.h>


#include <cstddef>
#include <vector>
#include "../misc/types.hpp"
#include "coal_types.hpp"
#include "../physics/world.hpp"

namespace cardillo::collision {

class BeamHullShape : public coal::ShapeBase {
public:
    struct Disc {
        coal::Vec3s c;  coal::Vec3s n;  real_t r;
        bool operator==(const Disc& o) const { return c == o.c && n == o.n && r == o.r; }
    };
    
    BeamHullShape(const C_Collider_BeamHull& collider, const entt::registry& reg);

    void updateRings(const C_Collider_BeamHull& collider, const entt::registry& reg);

    const std::vector<coal::Vec3s>& ringA() const { return ringA_; }
    const std::vector<coal::Vec3s>& ringB() const { return ringB_; }

    // Return the support witness on one endpoint cross-section
    coal::Vec3s endpointSupport(const coal::Vec3s& dir, int endpoint) const;

    coal::NODE_TYPE getNodeType() const override { return coal::GEOM_CUSTOM; }

    void computeShapeSupport(const coal::Vec3s& dir, coal::Vec3s& support, int& hint, coal::details::ShapeSupportData& data) const override;

    // Support set (contact face) in the direction of supportSet.getNormal()
    void computeShapeSupportSet(coal::SupportSet& supportSet, int& hint,
                                coal::details::ShapeSupportData& data,
                                std::size_t numSamples, real_t tol) const override;

    void computeLocalAABB() override;
    CollisionGeometry* clone() const override;
    bool isEqual(const coal::CollisionGeometry& other) const override;

private:
    void recomputeLocalAABB_();

    bool circular_ = false;
    std::vector<coal::Vec3s> ringA_, ringB_;      // for arbitrary cross-sections
    Disc discA_, discB_; 
};

}  // namespace cardillo::collision