#pragma once

#include <coal/BV/AABB.h>
#include <coal/collision_object.h>
#include <coal/shape/geometric_shapes.h>


#include <vector>
#include "../misc/types.hpp"
#include "../physics/world.hpp"

namespace cardillo::collision {

class BeamHullShape : public coal::ShapeBase {
public:
    struct Disc { coal::Vec3s c;  coal::Vec3s n;  coal::CoalScalar r; };
    
    BeamHullShape(const C_Collider_BeamHull& collider, const entt::registry& reg);

    void updateRings(const C_Collider_BeamHull& collider, const entt::registry& reg);

    const std::vector<coal::Vec3s>& ringA() const { return ringA_; }
    const std::vector<coal::Vec3s>& ringB() const { return ringB_; }

    // Return the support witness on one endpoint cross-section. This is separate from the
    // aggregate GJK support because contact Jacobians need endpoint provenance.
    coal::Vec3s endpointSupport(const coal::Vec3s& dir, int endpoint) const;

    coal::NODE_TYPE getNodeType() const override { return coal::GEOM_CUSTOM; }

    void computeShapeSupport(const coal::Vec3s& dir, coal::Vec3s& support, int& hint, coal::details::ShapeSupportData& data) const override;

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


// TODO:
// Maybe use oriented capsule to skip direct support function evals?
// What about round crossections?