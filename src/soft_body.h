// SPDX-FileCopyrightText: 2026 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "container.h"
#include "id_pool.h"

#include "box2d/soft_body.h"
#include "box2d/types.h"

typedef struct b2World b2World;

b2DeclareArray( b2Vec2 );
b2DeclareArray( b2SoftBodyTriangle );

// A soft body particle. Position is stored as a b2Pos so large world coordinates keep their
// precision; every per-sub-step delta is float.
typedef struct b2SoftParticle
{
	b2Pos position;
	b2Pos prevPosition;
	b2Vec2 velocity;
	b2Vec2 accumulatedImpulse;
	float invMass;
} b2SoftParticle;
b2DeclareArray( b2SoftParticle );

// Distance constraints are structural edges. Bending constraints resist folding between pairs of
// particles that share an element but are not structural edges.
typedef enum b2SoftConstraintType
{
	b2_softDistanceConstraint = 0,
	b2_softBendingConstraint = 1,
} b2SoftConstraintType;

typedef struct b2SoftConstraint
{
	int particle1;
	int particle2;
	float restLength;

	// Per-iteration PBD stiffness factor in [0,1] (1 = inextensible, 0 = no resistance).
	float stiffnessFactor;

	// Bending constraints are never torn.
	uint8_t type;
	bool active;
} b2SoftConstraint;
b2DeclareArray( b2SoftConstraint );

typedef struct b2SoftAnchor
{
	int particleIndex;

	// B2_NULL_INDEX for a world anchor, otherwise the internal body id.
	int bodyId;
	b2Vec2 localPoint;

	// Set when the anchored body has been destroyed. The anchor then pins at lastWorldPoint.
	bool detached;
	b2Pos lastWorldPoint;
} b2SoftAnchor;
b2DeclareArray( b2SoftAnchor );

// Spatial hash grid for self-collision, built each sub-step with a counting sort. All arrays are
// sized once at creation and never grow during a step, so a rebuild is strictly O(particles).
typedef struct b2SoftGrid
{
	// Number of hash buckets (power of two, at least twice the particle count).
	int bucketCount;

	// bucketHead[b] is the first particle in bucket b, or -1 when empty.
	b2Array( int ) bucketHead;

	// nextLink[i] is the next particle sharing i's bucket, forming an intrusive linked list.
	b2Array( int ) nextLink;

	// Scratch cell coordinates and hash bucket per particle, recomputed each rebuild.
	b2Array( int ) cellX;
	b2Array( int ) cellY;
	b2Array( int ) particleBucket;

	float cellSize;
	int bucketMask;
} b2SoftGrid;

// The full soft body. All arrays are allocated once at creation and never grow during a step.
typedef struct b2SoftBody
{
	int id;
	uint16_t generation;
	bool inUse;

	b2Array( b2SoftParticle ) particles;
	b2Array( b2SoftConstraint ) constraints;
	b2Array( b2SoftAnchor ) anchors;

	// Element triangulation, kept for debug drawing and surface queries
	b2Array( b2SoftBodyTriangle ) triangles;

	b2SoftGrid grid;

	b2Filter filter;
	float particleRadius;
	float particleMass;
	float stiffness;
	float bendStiffness;
	float damping;
	float friction;
	int iterations;

	bool enableSelfCollision;
	bool enableTearing;
	float tearMaxStretch;
	bool enableGravity;

	// Reusable scratch data for rigid collision, sized at creation
	b2Array( int ) shapeCandidates;

	// Diagnostics from the most recent step
	b2SoftBodyStats stats;

	// Bounding box of the particles, refreshed each step
	b2AABB bounds;

	void* userData;
} b2SoftBody;
b2DeclareArray( b2SoftBody );

// World lifecycle integration
void b2CreateSoftBodyPool( b2World* world );
void b2DestroySoftBodyPool( b2World* world );

// Advance every soft body in the world by one full step. Called after the rigid body solver so
// particles see the final rigid body poses of this step. The soft body uses the same h/dt.
void b2StepSoftBodies( b2World* world, float dt, float invDt, float h, float invH, int subStepCount );

// Lookup with validation. Returns NULL for a stale or null id.
b2SoftBody* b2GetSoftBodyFullId( b2World* world, b2SoftBodyId id );
b2SoftBodyId b2MakeSoftBodyId( b2World* world, int id );

// Ensure the world's soft body impulse staging arrays cover the given internal body id. This may
// allocate, so it is only called outside of the step loop.
void b2EnsureSoftImpulseCapacity( b2World* world, int bodyId );

// Mark every anchor referencing the given internal body id as detached. The anchor keeps pinning
// its particle at the last known world point instead of dereferencing the destroyed body.
void b2DetachSoftBodyAnchors( b2World* world, int bodyId );

// Draw every soft body in the world using the debug draw interface.
struct b2DebugDraw;
void b2DrawSoftBodies( b2World* world, struct b2DebugDraw* draw );
