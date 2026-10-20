# Collision (Coal) and geometry

Code: `src/collision/collision_coal.{hpp,cpp}`, `beam_hull_shape.{hpp,cpp}`, `contact_tracker.*`,
`src/physics/assets/{rigid_body_factory,beam_factory,assets}.cpp`, contact rows in
`DynamicsAssembler::rebuildW_()`.

## Coal version

`src/CMakeLists.txt` fetches the fork `https://github.com/41pha1/coal.git` (branch
`beam-support-set-hook`, commit `a90ca6e`): upstream devel including coal-library/coal#908 (contact
patch fixes by the project owner) plus custom shapes (`GEOM_CUSTOM`, `computeShapeSupportSet` hook)
needed by `BeamHullShape`. Switch back to an upstream release once these are merged there. Coal's
headers are marked as system headers (Coal uses members it has deprecated itself).

`src/collision/coal_types.hpp` only holds `static_assert`s (`real_t == coal::Scalar`,
`Vector3r == coal::Vec3s`); it is included once, in `collision_coal.cpp`. Pass Eigen values to Coal
directly; use `coal::Triangle32` for triangles.

## Contact patches

Coal computes contacts per pair and then contact patches (`collision.use_patch_vertices`,
`collision.max_patches`, `collision.patch_num_samples`). A vertex belongs to a patch if it lies within
the patch tolerance of the separating plane; the engine raises the tolerance per pair to
`max(collision.patch_tolerance, penetration depth)`, otherwise slightly penetrating corners of a
tilted box are dropped. Contact matching for warm starts: `ContactTracker`
(`collision.match_max_dist`).

## Mesh bodies

- **Frames:** a dynamic mesh body (mass > 0) is normalized at creation: its pose is shifted to the
  centre of mass and rotated into the principal axes, so it differs from the requested `RigidState`.
  Everything given in body coordinates (constraint attachment points, body-frame rotations of
  trajectories) refers to this principal frame. Read the actual pose from the ECS
  (`C_Position3`, `C_Orientation`) and convert world-defined points. `C_Mesh::normalized` records the
  choice so the normalized mesh is still used after `makeStatic()`.
- **Triangle meshes (BVH)** yield one contact point per intersecting triangle, i.e. poor support for
  boxes lying on large triangles. For closed convex meshes set `MeshShape::convex = true`: Coal then
  uses a `coal::Convex` built from the mesh's own vertices/triangles and produces full patches.

## Beam colliders

- `BeamColliderMode::RigidBodyPrimitive`: one primitive (box/capsule) per node.
- `BeamColliderMode::InterSegmentHull`: `BeamHullShape` (custom Coal shape) spanning two neighbouring
  node cross sections; contacts get several weighted attachments (both nodes), handled generically by
  the contact Jacobian (`ContactSide::attachments`). Since #77 hull beams place nodes at the spline
  points (segments + 1 per spline, half segment length at the ends); `createBeams()` merges the
  coinciding junction nodes of consecutive splines. End nodes lie exactly at the beam ends, so
  disable collisions with bodies they are attached to (see cable scene).

## Surface velocity

`engine.setSurfaceVelocity(body, v_body)` (`C_SurfaceVelocity`, body coordinates): conveyor belts with
static geometry. Added to the relative contact velocity of each attachment and projected onto the
contact tangent plane, so it only drives friction. Example: `examples/scenes/conveyor_merge`.
