#pragma once

#include <type_traits>

#include <coal/data_types.h>

#include "../misc/types.hpp"

namespace cardillo::collision {
// CardilloCxx passes its scalars, vectors and quaternions to Coal without conversion, so the
// types must agree (Coal uses double unless it is built with COAL_USE_FLOAT_PRECISION).
static_assert(std::is_same_v<real_t, coal::Scalar>, "real_t must equal coal::Scalar (Coal built with COAL_USE_FLOAT_PRECISION?)");
static_assert(std::is_same_v<Vector3r, coal::Vec3s>, "Vector3r must equal coal::Vec3s");
}  // namespace cardillo::collision
