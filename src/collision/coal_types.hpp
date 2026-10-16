#pragma once

#include <memory>
#include <type_traits>
#include <vector>

#include <Eigen/Geometry>
#include <coal/BVH/BVH_model.h>
#include <coal/data_types.h>

#include "../misc/types.hpp"

namespace cardillo::collision {
// Scalar, quaternion and triangle types of Coal, independent of the Coal version: coal::CoalScalar,
// coal::Quatf and coal::Triangle of Coal 3.0 are deprecated/removed in later versions (coal::Scalar,
// coal::Quats, coal::Triangle32), whereas coal::Vec3s and the BVH's triangle storage exist in all.
using CoalReal = coal::Vec3s::Scalar;
using CoalQuat = Eigen::Quaternion<CoalReal>;
using CoalTriangle = std::decay_t<decltype(std::declval<coal::BVHModelBase&>().tri_indices)>::element_type::value_type;

// CardilloCxx passes its vectors and quaternions to Coal without conversion, so the scalar types
// must agree (Coal uses double unless built with COAL_USE_FLOAT_PRECISION).
static_assert(std::is_same_v<real_t, CoalReal>, "real_t must equal Coal's scalar type (Coal built with COAL_USE_FLOAT_PRECISION?)");
static_assert(std::is_same_v<Vector3r, coal::Vec3s>, "Vector3r must equal coal::Vec3s");
}  // namespace cardillo::collision
