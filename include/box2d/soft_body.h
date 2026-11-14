// SPDX-FileCopyrightText: 2026 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "base.h"
#include "id.h"
#include "math_types.h"
#include "types.h"

#include <stdbool.h>

/**
 * @defgroup softBody Soft Body
 * @brief Optional position-based soft bodies (cloth, ropes, jelly, soft characters).
 *
 * Soft bodies are an optional feature built on top of the rigid body world. A soft body is a
 * particle mesh (triangles and/or quadrilaterals) connected by distance constraints and bending
 * constraints. It is solved with an iterative position-based dynamics (PBD) solver that shares the
 * world time step and the world contact filtering.
 *
 * Soft bodies do not exist unless explicitly created with b2CreateSoftBody(). Creating no soft
 * bodies leaves the rigid body simulation, its API and its deterministic results unchanged.
 *
 * All soft body arrays are allocated up front at creation time. Stepping a soft body does not
 * perform any heap allocation and the solver cost is approximately linear in the number of
 * particles and constraints.
 * @{
 */

/// A triangle element in a soft body mesh. Indices reference the vertex array of the soft body
/// definition. Winding order does not matter.
typedef struct b2SoftBodyTriangle
{
	int vertex1, vertex2, vertex3;
} b2SoftBodyTriangle;

/// A quadrilateral element in a soft body mesh. Indices reference the vertex array of the soft
/// body definition. The quad is stored as given (it is triangulated internally) and should be
/// convex.
typedef struct b2SoftBodyQuad
{
	int vertex1, vertex2, vertex3, vertex4;
} b2SoftBodyQuad;

/// An anchor attaches a soft body particle to a point in the world or to a rigid body. Anchors with
/// a null body id are world anchors: the particle is pinned at the anchor world position. Anchors
/// with a valid body id follow the body: the local point is transformed with the body pose every
/// sub-step.
typedef struct b2SoftBodyAnchor
{
	/// Particle index within the soft body
	int particleIndex;

	/// Body to attach to, or b2_nullBodyId for a static world-space anchor
	b2BodyId bodyId;

	/// Anchor point. For a world anchor this is the world position. For a body anchor this is the
	/// anchor point in body local coordinates.
	b2Vec2 localPoint;
} b2SoftBodyAnchor;

/// Soft body definition. Initialize with b2DefaultSoftBodyDef(), then fill in the mesh. At least
/// one vertex and one constraint are required.
typedef struct b2SoftBodyDef
{
	/// The mesh vertices in world space. Duplicate vertices (coincident points) are rejected.
	const b2Vec2* vertices;
	int vertexCount;

	/// Optional triangle elements. May be NULL when triangleCount is zero.
	const b2SoftBodyTriangle* triangles;
	int triangleCount;

	/// Optional quadrilateral elements. May be NULL when quadCount is zero.
	const b2SoftBodyQuad* quads;
	int quadCount;

	/// Per-particle mass in kilograms. Must be greater than zero.
	float particleMass;

	/// Particle collision radius in meters. Must be greater than zero. Particles collide with
	/// rigid shapes and, when enabled, with each other as inflated circles of this radius.
	float particleRadius;

	/// Contact filter used for particle versus rigid shape collision. Reuses the standard Box2D
	/// category/mask/group filtering rules, including the custom filter callback of the world.
	b2Filter filter;

	/// Stiffness of the distance constraints, in [0, 1]. A value of 1 makes edges nearly
	/// inextensible; lower values allow stretching.
	float stiffness;

	/// Stiffness of the bending constraints generated across triangle and quad diagonals, in
	/// [0, 1]. Set to zero to disable bending resistance.
	float bendStiffness;

	/// Damping applied to particle velocities each sub-step, in [0, 1]. Zero is no damping.
	float damping;

	/// Particle friction coefficient used when sliding along rigid shapes, in [0, 1].
	float friction;

	/// Number of PBD solver iterations per sub-step. More iterations improve constraint
	/// satisfaction at a linear cost. Must be at least one.
	int iterations;

	/// Enable particle versus particle self-collision.
	bool enableSelfCollision;

	/// Enable tearing. When enabled a distance constraint whose length exceeds
	/// tearMaxStretch times its rest length is permanently removed.
	bool enableTearing;

	/// Tearing threshold as a multiple of the rest length. Must be greater than one.
	float tearMaxStretch;

	/// Optional anchors. May be NULL when anchorCount is zero.
	const b2SoftBodyAnchor* anchors;
	int anchorCount;

	/// Use the world gravity. When false the soft body is weightless (other forces still apply).
	bool enableGravity;

	/// User data pointer stored on the soft body.
	void* userData;

	/// Used internally to detect a valid definition. DO NOT SET.
	int internalValue;
} b2SoftBodyDef;

/// Use this to initialize your soft body definition with safe defaults.
B2_API b2SoftBodyDef b2DefaultSoftBodyDef( void );

/// Per-step stability and energy diagnostics for a soft body. All energies are in joules.
typedef struct b2SoftBodyStats
{
	/// Number of particles currently simulated
	int particleCount;

	/// Number of distance and bending constraints still active (torn constraints are removed)
	int constraintCount;

	/// Number of constraints removed by tearing so far
	int tornCount;

	/// Number of particle versus rigid shape contacts applied during the last step
	int rigidContactCount;

	/// Number of self-collision contacts applied during the last step
	int selfContactCount;

	/// Total particle kinetic energy: 0.5 * m * |v|^2
	float kineticEnergy;

	/// Potential energy stored in the distance and bending constraints
	float strainEnergy;

	/// Gravitational potential energy relative to the lowest particle
	float potentialEnergy;

	/// Maximum constraint stretch observed during the last step as a multiple of rest length
	float maxStretch;

	/// True when the last step detected a numerical failure (NaN/Inf or runaway stretch). The
	/// solver clamps the failure instead of exploding; this flag explains the clamp.
	bool numericalFailure;
} b2SoftBodyStats;

/// Create a soft body in the given world from a triangle/quad mesh. Returns b2_nullSoftBodyId and
/// logs an error when the definition is invalid (degenerate mesh, duplicate vertices, out of range
/// indices, non-finite data, bad parameters).
B2_API b2SoftBodyId b2CreateSoftBody( b2WorldId worldId, const b2SoftBodyDef* def );

/// Destroy a soft body and release all its memory. Safe to call with a null id.
B2_API void b2DestroySoftBody( b2SoftBodyId softBodyId );

/// True when the soft body id is valid and refers to a live soft body.
B2_API bool b2SoftBody_IsValid( b2SoftBodyId softBodyId );

/// Get the number of particles in a soft body.
B2_API int b2SoftBody_GetParticleCount( b2SoftBodyId softBodyId );

/// Get the current world position of a particle.
B2_API b2Vec2 b2SoftBody_GetParticlePosition( b2SoftBodyId softBodyId, int particleIndex );

/// Get the current velocity of a particle in meters per second.
B2_API b2Vec2 b2SoftBody_GetParticleVelocity( b2SoftBodyId softBodyId, int particleIndex );

/// Override the position of a particle immediately. The corresponding velocity is left unchanged.
B2_API void b2SoftBody_SetParticlePosition( b2SoftBodyId softBodyId, int particleIndex, b2Vec2 position );

/// Apply an impulse (kg * m/s) to a particle. Applied at the beginning of the next sub-step.
B2_API void b2SoftBody_ApplyParticleImpulse( b2SoftBodyId softBodyId, int particleIndex, b2Vec2 impulse );

/// Get the number of anchors on a soft body.
B2_API int b2SoftBody_GetAnchorCount( b2SoftBodyId softBodyId );

/// Remove an anchor by index. A particle whose anchors are all removed becomes free. Removing the
/// last anchor of a particle does not destroy the particle.
B2_API void b2SoftBody_DestroyAnchor( b2SoftBodyId softBodyId, int anchorIndex );

/// Add an anchor to an existing soft body. Returns false and logs an error when the particle
/// index or body id is invalid.
B2_API bool b2SoftBody_AddAnchor( b2SoftBodyId softBodyId, const b2SoftBodyAnchor* anchor );

/// Enable or disable self-collision on an existing soft body.
B2_API void b2SoftBody_EnableSelfCollision( b2SoftBodyId softBodyId, bool enable );

/// Enable or disable tearing and update the tear threshold.
B2_API void b2SoftBody_EnableTearing( b2SoftBodyId softBodyId, bool enable, float maxStretch );

/// Update the solver iteration count (at least one).
B2_API void b2SoftBody_SetIterations( b2SoftBodyId softBodyId, int iterations );

/// Update the user data pointer.
B2_API void b2SoftBody_SetUserData( b2SoftBodyId softBodyId, void* userData );

/// Get the user data pointer.
B2_API void* b2SoftBody_GetUserData( b2SoftBodyId softBodyId );

/// Get stability and energy diagnostics from the most recent step.
B2_API b2SoftBodyStats b2SoftBody_GetStats( b2SoftBodyId softBodyId );

/**@}*/
