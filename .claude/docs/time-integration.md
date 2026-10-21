# Time integration: Moreau-θ with the position-level force law

Code: `src/physics/integration/moreau.cpp`, `src/physics/assembly/dynamics_assembler.{hpp,cpp}`
(`springBias`, `damperBias`, `*ComplianceScale`, `advancePrescribedMotion`), `src/config/config.hpp`
(header comment). User docs: `docs/chapters/moreau_time_stepping.rst`, `docs/chapters/solvers/*.rst`.
Theory: `LaTeX/theory_doc.tex` if present locally (gitignored, separate git repository).

## Step

```
q_{n+θ} = q_n + (1−θ) h u_n            (explicitPositionUpdate)
advance prescribed rates by (1−θ) h      (ConstraintPattern::advancePrescribedMotion)
assemble at t_{n+θ} = t_n + (1−θ) h      (drivers, contacts, W, Wg, Wγ, g, χ)
solve for u_{n+1} and the impulses
q_{n+1} = q_{n+θ} + θ h u_{n+1};  advance prescribed rates by θ h
```

## Compliant rows (springs, perfect constraints) and dampers

Unknown: the θ-weighted force λ_θ. The position-level law is linearized at q_{n+θ}:

```
C λ_θ + g(t_{n+θ}, q_{n+θ}) + h (θ² Wᵀu_{n+1} − (1−θ)² Wᵀu_n) + (2θ − 1) h χ = 0
```

With the impulses Λ = h λ_θ every solver sees the row `Wᵢᵀu_{n+1} + bᵢ + cᵢ Λᵢ = 0` with

- spring compliance factor `c_g = C / (θ² h²)` (`springComplianceScale`)
- damper compliance factor `c_γ = A / (θ h)` (`damperComplianceScale`)
- `springBias = g / (θ² h) − ((1−θ)/θ)² Wᵀu_n + (2θ−1)/θ² χ`
- `damperBias = ((1−θ)/θ) Wᵀu_n + χ/θ`

Saddle-point assemblers (PJ, cone programs) use the negative bias as right-hand side. No multiplier of
the previous step enters (Λ_n is only a warm start), no stabilization parameter, no drift.

- θ = 1/2 is the trapezoidal rule (energy-conserving for linear springs, `test_position_level`), but
  violations of perfect constraints are not damped and grow → warning; scenes should use θ > 1/2.
- Errors are first order in (θ − 1/2) h (`test_torsional_oscillator` checks the ratio).
- Removed on purpose: Baumgarte bias (`constraint_bias_factor`, `biasFactor`), `moreau.lambda_theta`,
  the velocity-level law, and the variant that evaluated g at t_n.

## Prescribed motion and drivers

- Constraint rates (`setConstraintScalarVelocity`/`Linear`/`Angular`) are integrated into the
  position-level reference by `advancePrescribedMotion()` (distance constraint: L0; joints: g0 and the
  joint frame of B), otherwise the position law would pull back against the prescribed motion.
- Trajectories (`C_StaticTrajectory`, `src/physics/synchronization/trajectory.cpp`) set the pose at
  t_n + (1−θ)h and the velocity to `(pose(t_n+h) − pose(t_n)) / h` (angular velocity in the body
  frame). Only pose functions exist; a velocity v(t) is prescribed via its integral.
- `addTrajectory` makes the body static. For dynamic mesh bodies, read the created pose from the ECS
  first (principal-axes frame, see `collision.md`).

## Joints (`TranslationRotationConstraint`)

- Rows 0..2: translational gap in the joint frame of A.
- Row 3 (x, hinge axis): the unwrapped twist angle −φ, φ = atan2(s, c) with
  s = (e1z·e2y − e1y·e2z)/2, c = (e1y·e2y + e1z·e2z)/2, unwrapped against `m_twist_ref` (rotation per
  evaluation must stay below π). Its Jacobian is the exact derivative −(c n_s − s n_c)/(s² + c²).
- Rows 4, 5: −e1x·e2z and e1x·e2y (sines; stay zero under rotation about e1x, never lose rank).
- `computeJointTriads()` evaluates A_IJ1, A_IJ2, s, c once per `getConstraint()`; Jacobian and
  position error share it.

## Rod elements (`BeamConstraint`)

Discrete Cosserat rod (rigid nodes, six-row compliant elements, strains γ, κ from the unnormalized
mean quaternion). Rest state: `BeamSpringParams::relaxedInInitialConfiguration()` (default) or
`relaxedInStraightConfiguration()`, optional `withAxialPretension(strain)`. Element stiffness scales
with 1/l₀, so zero-length elements are singular (see the junction merge in `BeamFactory::createBeams`).
Large axial stiffness × pretension strain quickly yields enormous forces (3DPrinter belt: EA = 5e9 N,
pretension kept at 0.17 %).

## Contact impulses

`Contact::last_impulse = (pn, pt1, pt2)`: percussion on body B = pn·normal + pt1·tangent1 +
pt2·tangent2, normal from A to B, pn ≥ 0. All solvers work with it directly (no conversion). Bilateral
percussions are π = h·λ (stored in `Lambda_g`/`Lambda_gamma`); every block acts as +W·π,
u = u_free + M⁻¹Wπ, with W unchanged (theory manual, eq. (3.4), (3.12), Section 4). VTK contact
manifolds export `pn`, `pt1`, `pt2`, `pt_mag`, `percussion`, `percussion_normal`, `percussion_tangent`.
