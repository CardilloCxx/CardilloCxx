// beam_hull_shape.hpp
//
// A custom coal::ShapeBase representing the convex hull of two identical 2D
// cross-sections, each independently positioned/oriented in 3D -- i.e. a
// generalized (twisted/tapered/bent) beam segment between two arbitrarily
// posed rings of points.
//
// The support function of the convex hull of two point sets A, B is exactly
//     support(dir) = argmax_{p in A u B} dot(dir, p)
// so no hull construction (no qhull call, no triangulation) is needed at
// all -- this is a plain O(|A|+|B|) scan, called once per GJK/EPA iteration,
// same cost class as any other built-in primitive's support function.
//
// Requires coal built from the GEOM_CUSTOM fork (coal-library/coal PR #822,
// rebased onto v3.0.4 -- see the top-level CMakeLists.txt / README for the
// exact fork commit pinned by this project). GEOM_CUSTOM adds
// ShapeBase::computeShapeSupport() as the one virtual escape hatch in
// Coal's otherwise NODE_TYPE-switch-based dispatch. Before building against
// a new fork commit, diff its
//   include/coal/shape/geometric_shapes.h
//   include/coal/narrowphase/support_functions.h
// against what's assumed in beam_hull_shape.cpp, since the exact virtual
// signature is this fork's own design, not a stable upstream API.
#pragma once

#include <coal/BV/AABB.h>
#include <coal/collision_object.h>
#include <coal/shape/geometric_shapes.h>

#include <vector>

namespace cardillo::collision {

class BeamHullShape : public coal::ShapeBase {
public:
    BeamHullShape(std::vector<coal::Vec3s> ringA, std::vector<coal::Vec3s> ringB);

    // Called every frame for a dynamic beam (see CollisionCoal::applyTransforms's Beam branch):
    // overwrites both rings' points in place -- cheap, no allocation once capacity has settled --
    // and refreshes the local AABB. No geometry-pointer swap, no qhull call, ever.
    void updateRings(std::vector<coal::Vec3s> ringA, std::vector<coal::Vec3s> ringB);

    const std::vector<coal::Vec3s>& ringA() const { return ringA_; }
    const std::vector<coal::Vec3s>& ringB() const { return ringB_; }

    coal::NODE_TYPE getNodeType() const override { return coal::GEOM_CUSTOM; }

    void computeShapeSupport(const coal::Vec3s& dir, coal::Vec3s& support, int& hint, coal::details::ShapeSupportData& data) const override;

    void computeLocalAABB() override;
    CollisionGeometry* clone() const override;
    bool isEqual(const coal::CollisionGeometry& other) const override;

private:
    void recomputeLocalAABB_();

    std::vector<coal::Vec3s> ringA_, ringB_;
};

}  // namespace cardillo::collision


// TODO:
// Maybe use oriented capsule to skip direct support function evals?
// What about round crossections?