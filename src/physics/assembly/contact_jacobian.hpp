#pragma once

#include "../../misc/types.hpp"

namespace cardillo::physics {

// Builds one dense row (length `dof`, 3 or 6 -- the only two body-DOF shapes in this engine: a
// point mass or a full rigid body) of a contact's Jacobian for one side (body) and one direction
// (normal or tangent). `dir_inertial`/`dir_body` are the same physical direction expressed in the
// inertial basis and in the body-fixed basis respectively; `r_body` is the contact point relative
// to the body's center of mass, in body-local coordinates. `sign` is +1/-1 selecting which side of the
// contact this row is for (bodies A and B pull the constraint in opposite directions).
// Returns a fixed-size (stack-allocated) buffer rather than VectorXr since `dof` is always one of
// exactly two compile-time-known values here -- called up to 6x per contact per timestep (see
// dynamics_assembler.cpp's accumulateDirForSide), so this avoided a heap allocation per call.
// Callers read only the first `dof` entries (row.head(dof)); entries beyond `dof` stay zero.
inline Vector6r buildContactRowByDof(int dof, const Vector3r& dir_inertial, const Vector3r& r_body, const Vector3r& dir_body, real_t sign) {
    Vector6r row = Vector6r::Zero();
    if (dof < 3) return row;
    row.head<3>() = sign * dir_inertial;
    if (dof >= 6) row.tail<3>() = r_body.cross(sign * dir_body);
    return row;
}

}  // namespace cardillo::physics
