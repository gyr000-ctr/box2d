# Soft Bodies

Box2D is primarily a rigid body engine. Rigid bodies are fast and stable but cannot deform. Flags,
ropes, jelly, and soft characters cannot be represented by a single rigid body; the traditional
workaround is to chain many rigid bodies with joints, which is expensive and hard to tune.

Box2D therefore includes an *optional* soft body solver built on top of the rigid body world. A
soft body is a particle mesh connected by distance and bending constraints and simulated with an
iterative position-based dynamics (PBD) solver.

The feature is strictly opt-in. Soft bodies do not exist unless you explicitly call
`b2CreateSoftBody()`. When no soft body is created, the rigid simulation, its API, and its
deterministic results are unchanged.

## Mesh construction

You provide a triangle and/or quadrilateral mesh in world space.

```c
b2SoftBodyDef def = b2DefaultSoftBodyDef();
def.vertices = vertices;
def.vertexCount = vertexCount;
def.triangles = triangles;
def.triangleCount = triangleCount;
b2SoftBodyId softId = b2CreateSoftBody( worldId, &def );
```

From the mesh Box2D derives:

- **Distance constraints** for every unique structural edge.
- **Bending constraints** connecting the opposite vertices of pairs of triangles that share an
  edge. Boundary edges (a single neighboring triangle) produce no bending constraint.
- **Particles**, one per vertex, each with the configured mass and collision radius.

Quads are triangulated internally along a diagonal and both child triangles must be non-degenerate.

## Validation and error handling

Soft body creation validates its inputs and returns `b2_nullSoftBodyId` (after logging a specific
error) rather than producing an unstable simulation. The following are rejected up front:

- Fewer than three vertices or no triangle/quad elements.
- Duplicate or non-finite vertices.
- Out-of-range or repeated indices inside an element.
- Degenerate (zero-area) triangles and concave or self-intersecting quads.
- Non-positive or non-finite mass and particle radius.
- Stiffness, damping, or friction outside `[0, 1]`.
- Tearing enabled with a threshold outside `(1, 100]`.
- Anchors that reference an invalid particle or body.

At run time the solver guards against extreme stretch, non-finite positions and velocities, and
insufficient sub-steps. Instead of exploding, it clamps the offending particle and reports the
condition through `b2SoftBodyStats.numericalFailure`.

## Constraints and stiffness

Each PBD iteration projects every active distance and bending constraint back toward its rest
length. The positional correction is scaled by a per-constraint stiffness in `[0, 1]`, so a value
of one is nearly inextensible while lower values allow stretching. More iterations improve
constraint satisfaction at a linear cost and are set with `b2SoftBodyDef.iterations`.

Tearing is optional. When enabled, a distance constraint whose length exceeds
`tearMaxStretch` times its rest length is permanently removed; bending constraints never tear.
The cumulative number of removed constraints is reported as `tornCount`.

## Gravity, damping, and anchors

Soft bodies use the world gravity unless `enableGravity` is false. A per-sub-step velocity damping
factor removes energy and helps the mesh settle.

Anchors pin a particle either to a fixed world point or to a rigid body:

- A null `bodyId` creates a world-space anchor.
- A valid `bodyId` creates a body anchor whose target follows the body pose every sub-step using
  the stored local point.

Anchored particles are given infinite mass and always win the constraint projection. If the body
backing an anchor is destroyed, the anchor is detached automatically and keeps pinning its particle
at the last known world point rather than dereferencing a stale body. Anchors can also be added and
removed at run time; removing the last anchor of a particle restores its finite mass.

## Collision and filtering

Particles collide with rigid shapes as inflated circles of radius `particleRadius`. Contact
detection uses the world broad phase and the standard category/mask/group filtering, including the
world custom filter callback. The broad-phase query uses a swept box spanning the particle's
sub-step motion, which prevents fast particles from tunneling through thin shapes. Deep overlaps
are resolved with a surface-based fallback normal.

Particle-versus-particle self-collision is off by default and enabled with
`enableSelfCollision`. It uses a spatial hash grid sized once at creation; self-collision keeps
particles at least two radii apart.

Two-way coupling applies an equal and opposite impulse to any dynamic rigid body a particle
pushes against, including an angular impulse about the body center of mass. Rigid bodies and soft
bodies share the same `b2World_Step` call and the same sub-step length.

## Memory and performance

All soft body arrays — particles, constraints, anchors, the self-collision grid, and the rigid
contact scratch — are allocated once at creation. Stepping performs **no heap allocation**, and the
per-body coupling impulse arrays grow only when a new rigid body is created, never during a step.

The solver cost is approximately linear in the number of particles and constraints. Measured
per-step cost on a release build (GCC 15, `-O2`, 8 iterations, 4 sub-steps, single worker) for a
hanging cloth grid is:

| Particles | Self-collision off | Self-collision on |
|-----------|--------------------|-------------------|
| 256       | ~0.49 ms (~1.9 µs/particle) | ~0.58 ms (~2.3 µs/particle) |
| 961       | ~2.0 ms (~2.1 µs/particle)  | ~2.4 ms (~2.5 µs/particle)  |
| 3721      | ~8.2 ms (~2.2 µs/particle)  | ~13.4 ms (~3.6 µs/particle) |

Without self-collision the per-particle cost is essentially constant (~2 microseconds per particle
per step) across these sizes, confirming linear scaling. Self-collision uses a spatial hash built
each sub-step with an intrusive linked list (also linear to build). Its constant factor grows with
mesh density because a folded cloth presents more genuinely overlapping neighbor pairs, which is
inherent to the number of contacts rather than an algorithmic super-linearity.

## Stability and energy diagnostics

`b2SoftBody_GetStats` returns per-step diagnostics:

- particle and active constraint counts,
- rigid and self contact counts,
- number of torn constraints,
- kinetic, strain, and gravitational potential energy in joules,
- the maximum constraint stretch observed, and
- a `numericalFailure` flag that is set whenever a guard had to clamp the solution.

These let an application monitor stability and energy rather than discovering a blow-up after the
fact. A hanging cloth, for example, should show finite and bounded energy with `numericalFailure`
remaining false.
