#pragma once

#include "../../misc/types.hpp"

namespace cardillo::physics {

// Builds one dense row (length `dof`, 3 or 6 -- the only two body-DOF shapes in this engine: a
// point mass or a full rigid body) of a contact's Jacobian for one side (body) and one direction
// (normal or tangent). `dir_world`/`dir_body` are the same physical direction expressed in the
// world frame and in the body's local frame respectively; `r_body` is the contact point relative
// to the body's center of mass, in body-local coordinates. `s` is +1/-1 selecting which side of the
// contact this row is for (bodies A and B pull the constraint in opposite directions).
// Returns a fixed-size (stack-allocated) buffer rather than VectorXr since `dof` is always one of
// exactly two compile-time-known values here -- called up to 6x per contact per timestep (see
// dynamics_assembler.cpp's accumulateDirForSide), so this avoided a heap allocation per call.
// Callers read only the first `dof` entries (row.head(dof)); entries beyond `dof` stay zero.
inline Vectorr<6> buildContactRowByDof(int dof, const Vector3r& dir_world, const Vector3r& r_body, const Vector3r& dir_body, real_t s) {
    Vectorr<6> row = Vectorr<6>::Zero();
    if (dof < 3) return row;

    row[0] = s * dir_world.x();
    row[1] = s * dir_world.y();
    row[2] = s * dir_world.z();

    if (dof >= 6) {
        const Vector3r t_body = r_body.cross(s * dir_body);
        row[3] = t_body.x();
        row[4] = t_body.y();
        row[5] = t_body.z();
    }
    return row;
}

}  // namespace cardillo::physics
