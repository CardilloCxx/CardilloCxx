Beams and Flexible Rods
=======================

Cardillo models elastic rods as a chain of rigid segments connected by
**beam constraints**, a Cosserat-rod spring model with six scalar rows
(three for stretch/shear, three for torsion/bending). The segments are
ordinary ECS entities, so all body factory methods, trajectories, and
post-processing tools work on them directly.

.. contents:: On this page
   :local:
   :depth: 2

Cross-section
---------------------------------

:cpp:struct:`BeamCrossSection <cardillo::physics::BeamCrossSection>` defines the
gometry of each beam segment and provides the derived section properties
used to compute stiffness:

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Field
     - Default
     - Meaning
   * - ``width``
     - ``0``
     - Cross-section width in metres. For circular sections (``Capsule`` or
       ``Cylinder`` type) this is treated as the **diameter** and the radius
       is ``min(width, height) / 2``.
   * - ``height``
     - ``0``
     - Cross-section height in metres. For rectangular sections this is the
       dimension along the beam's y-axis.
   * - ``type``
     - ``Cube``
     - Collision/visual shape used for each segment. One of ``BeamBodyType::Cube``,
       ``BeamBodyType::Capsule``, or ``BeamBodyType::Cylinder``.

Derived section properties computed from these fields:

- :cpp:member:`area() <cardillo::physics::BeamCrossSection::area>` -- cross-sectional area A (m²)
- :cpp:member:`Iy() <cardillo::physics::BeamCrossSection::Iy>` -- second moment of area about the y-axis (m⁴)
- :cpp:member:`Iz() <cardillo::physics::BeamCrossSection::Iz>` -- second moment of area about the z-axis (m⁴)
- :cpp:member:`Jp() <cardillo::physics::BeamCrossSection::Jp>` -- polar moment of inertia ``Iy + Iz`` (m⁴)
- :cpp:member:`sectionModulus() <cardillo::physics::BeamCrossSection::sectionModulus>` -- minimum section modulus (m³); used for yield checks

.. code-block:: cpp

   using namespace cardillo::physics;

   // 2 cm × 2 cm square cross-section
   BeamCrossSection squareSection(0.02, 0.02, BeamBodyType::Cube);

   // 1 cm diameter circular cross-section
   BeamCrossSection circSection(0.01, 0.01, BeamBodyType::Capsule);

Material stiffness
--------------------------------------

:cpp:struct:`BeamSpringParams <cardillo::physics::BeamSpringParams>` describes the
elastic and damping parameters of the beam.

Fields
~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 22 15 63

   * - Field
     - Default
     - Meaning
   * - ``E``
     - ``0``
     - Young's modulus (Pa). Used together with ``nu`` to derive ``Ke`` and
       ``Kf`` from the cross-section geometry. Set to 0 if you override with
       ``Ke_direct`` / ``Kf_direct``.
   * - ``nu``
     - ``0``
     - Poisson's ratio (dimensionless). Used to compute shear modulus
       ``G = E / (2*(1 + nu))``.
   * - ``scaleKe``
     - ``(1, 1, 1)``
     - Per-component scale factors applied to the computed extensional/shear
       stiffness vector ``Ke = (E*A/L, G*A/L, G*A/L)``. Component order:
       ``[axial, shear-y, shear-z]``.
   * - ``scaleKf``
     - ``(1, 1, 1)``
     - Per-component scale factors applied to the computed torsion/bending
       stiffness vector ``Kf = (G*Jp/L, E*Iy/L, E*Iz/L)``. Component order:
       ``[torsion, bend-y, bend-z]``.
   * - ``Ke_direct``
     - ``nullopt``
     - Optional direct override for ``Ke`` (N/m, beam-local frame).
       When set, ``E``, ``nu``, and ``scaleKe`` are ignored for extensional/shear.
   * - ``Kf_direct``
     - ``nullopt``
     - Optional direct override for ``Kf`` (Nm/rad). When set, ``E``, ``nu``,
       and ``scaleKf`` are ignored for torsion/bending.
   * - ``gamma0``
     - ``nullopt``
     - Optional rest-state translational strain ``γ₀``. When unset the
       factory initialises it from the relative pose of the two connected
       segments at creation time (natural configuration).
   * - ``kappa0``
     - ``nullopt``
     - Optional rest-state curvature ``κ₀``. When unset the factory
       initialises it from the initial segment orientation.
   * - ``dampingFactor``
     - ``0``
     - Rayleigh-type damping factor ``d`` applied to all stiffness terms.
       The effective damping stiffness is ``K * d`` (units: s).
       A value of 0.01 provides light numerical damping.

Named constructors
~~~~~~~~~~~~~~~~~~

**From a material**:

.. code-block:: cpp

   BeamSpringParams springs = BeamSpringParams::fromMaterial(
       real_t E_in,              // Young's modulus (Pa)
       real_t nu_in,             // Poisson's ratio
       real_t axialScale  = 1,   // scale on axial stiffness Ke[0]
       real_t shearScale  = 1,   // scale on shear stiffnesses Ke[1], Ke[2]
       real_t torsionScale = 1,  // scale on torsional stiffness Kf[0]
       real_t bendYScale  = 1,   // scale on bending stiffness Kf[1]
       real_t bendZScale  = 1,   // scale on bending stiffness Kf[2]
       real_t dampingFactor = 0  // Rayleigh damping factor
   );

   // Steel rod, slight damping
   BeamCrossSection section(0.01, 0.01, BeamBodyType::Capsule);
   BeamSpringParams springs = BeamSpringParams::fromMaterial(
       200e9,   // E = 200 GPa (steel)
       0.3,     // nu
       1.0, 1.0, 1.0, 1.0, 1.0, // all scales = 1
       0.005    // 0.5% Rayleigh damping
   );

**From direct stiffness values** (bypasses material formulas):

.. code-block:: cpp

   BeamSpringParams springs(
       Vector3r(1e4, 5e3, 5e3),   // Ke: axial, shear-y, shear-z (N/m)
       Vector3r(2e3, 8e2, 8e2),   // Kf: torsion, bend-y, bend-z (Nm/rad)
       0.01                        // damping factor
   );

How Ke and Kf are computed
~~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``Ke_direct`` and ``Kf_direct`` are not set, the stiffnesses are
derived per segment from the material constants, cross-section, and segment
length ``L``:

.. math::

   K_e = \begin{pmatrix} EA/L \\ GA/L \\ GA/L \end{pmatrix}
   \odot \text{scaleKe}
   \qquad
   K_f = \begin{pmatrix} GJ_p/L \\ EI_y/L \\ EI_z/L \end{pmatrix}
   \odot \text{scaleKf}

where :math:`G = E / (2(1+\nu))`.

Creating beams
---------------

:cpp:func:`createBeam <cardillo::physics::PhysicsEngine::createBeam>`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Samples a spline into ``segments`` rigid bodies and connects adjacent ones
with beam constraints. Returns the root and tip entity handles.

.. code-block:: cpp

   std::pair<entt::entity, entt::entity> [root, tip] = engine.createBeam(
       const misc::SplinePattern& spline,  // geometry of the rod centreline
       const BeamCrossSection& section,    // cross-section geometry + type
       const BeamSpringParams& springs,    // material / stiffness params
       const RigidState& stateDefaults,    // default state applied to all segments
       const RigidProps& propsDefaults,    // default props (mass/density, flags)
       size_t segments                     // number of rigid segments
   );

   auto [root, tip] = engine.createBeam(spline, section, springs,
                                         defaultState, defaultProps, 48);

The factory:

1. Samples the spline at ``segments`` evenly-spaced parameter values.
2. Creates a rigid body at each sample point (shape from ``section.type``).
3. Connects every adjacent pair with :cpp:func:`PhysicsEngine::addBeamConstraint <cardillo::physics::PhysicsEngine::addBeamConstraint>`.

:cpp:func:`createBeams <cardillo::physics::PhysicsEngine::createBeams>`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Chains several splines into a single beam, with the tip of each spline
connected to the root of the next. Returns the overall root and tip.

.. code-block:: cpp

   std::pair<entt::entity, entt::entity> [root, tip] = engine.createBeams(
       const std::vector<const misc::SplinePattern*>& splines,
       const BeamCrossSection& section,
       const BeamSpringParams& springs,
       const RigidState& stateDefaults,
       const RigidProps& propsDefaults,
       size_t segments   // total segments distributed across all splines
   );

Complete example
~~~~~~~~~~~~~~~~

.. code-block:: cpp

   using namespace cardillo;
   using namespace cardillo::physics;

   // 1-metre straight rod along the z-axis
   misc::LinearSpline spline(Vector3r(0, 0, 0), Vector3r(0, 0, 1));

   // 5 mm diameter circular cross-section
   BeamCrossSection section(0.005, 0.005, BeamBodyType::Capsule);

   // Nylon-like material: E ≈ 3 GPa, nu = 0.4
   BeamSpringParams springs = BeamSpringParams::fromMaterial(
       3e9, 0.4,
       1.0, 1.0, 1.0, 1.0, 1.0,
       0.01   // light Rayleigh damping
   );

   RigidState defaultState;                  // segments positioned by spline
   RigidProps defaultProps = RigidProps::withDensity(1150.0); // nylon density

   auto [beamRoot, beamTip] = engine.createBeam(spline, section, springs,
                                                 defaultState, defaultProps, 32);

   // Pin the root to the world
   engine.addHingeConstraint(beamRoot, entt::null,
       JointFrame::fromAxis(Vector3r::Zero(), Vector3r::UnitZ()),
       0.0, 0.0);

   // Track tip for output
   engine.track(beamTip, "beam_tip");

.. tip::
   Because beam segments are ordinary ECS entities, you can attach a
  trajectory to the root segment to drive the beam's base, or call
  :cpp:func:`applyForce() <cardillo::physics::PhysicsEngine::applyForce>` to impose point loads on any segment.
