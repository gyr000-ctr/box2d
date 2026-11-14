// SPDX-FileCopyrightText: 2026 Erin Catto
// SPDX-License-Identifier: MIT

#include "soft_body.h"

#include "aabb.h"
#include "body.h"
#include "broad_phase.h"
#include "core.h"
#include "dynamic_tree.h"
#include "physics_world.h"
#include "shape.h"

#include "box2d/box2d.h"
#include "box2d/constants.h"

#include <float.h>
#include <math.h>
#include <string.h>

// Soft body limits. These keep integer indexing and grid hashing safe while staying far above any
// practical mesh size.
#define B2_SOFT_MAX_PARTICLES ( 1 << 22 )
#define B2_SOFT_MAX_CONSTRAINTS ( 1 << 23 )
#define B2_SOFT_MAX_ELEMENTS ( 1 << 22 )

// A particle that would move more than this many radii in one sub-step is considered runaway. The
// solver clamps it and reports a numerical failure instead of producing NaN positions.
#define B2_SOFT_RUNAWAY_RADII 100.0f

// Stretch beyond this multiple of rest length is treated as numerical failure even without
// tearing enabled. Tearing thresholds above this are rejected at creation.
#define B2_SOFT_MAX_SANE_STRETCH 100.0f

typedef struct b2SoftEdge
{
	int particle1, particle2;
} b2SoftEdge;

// -------------------------------------------------------------------------------------------------
// Small helpers
// -------------------------------------------------------------------------------------------------

static uint32_t b2SoftHashPair( int a, int b )
{
	// Deterministic, order independent hash. The multiplier constants are arbitrary primes.
	uint32_t x = (uint32_t)a, y = (uint32_t)b;
	x = x * 73856093u;
	y = y * 19349663u;
	return x ^ y;
}

static bool b2SoftEdgeEquals( b2SoftEdge e1, b2SoftEdge e2 )
{
	return ( e1.particle1 == e2.particle1 && e1.particle2 == e2.particle2 ) ||
		   ( e1.particle1 == e2.particle2 && e1.particle2 == e2.particle1 );
}

// Per-iteration PBD stiffness in [0,1]. The position correction scales linearly by this factor,
// which is unconditionally stable for any positive time step. The constraint stores this value
// directly (1 => inextensible, 0 => no resistance).
static float b2SoftCompliance( float stiffness )
{
	return b2ClampFloat( stiffness, 0.0f, 1.0f );
}

// -------------------------------------------------------------------------------------------------
// Id and pool
// -------------------------------------------------------------------------------------------------

b2SoftBodyId b2MakeSoftBodyId( b2World* world, int id )
{
	B2_ASSERT( 0 <= id && id < world->softBodies.count );
	b2SoftBody* softBody = world->softBodies.data + id;
	b2SoftBodyId result = { id + 1, world->worldId, softBody->generation };
	return result;
}

b2SoftBody* b2GetSoftBodyFullId( b2World* world, b2SoftBodyId id )
{
	if ( B2_IS_NULL( id ) || id.world0 != world->worldId || id.generation == 0 )
	{
		return NULL;
	}

	int index = id.index1 - 1;
	if ( index < 0 || index >= world->softBodies.count )
	{
		return NULL;
	}

	b2SoftBody* softBody = world->softBodies.data + index;
	if ( softBody->inUse == false || softBody->generation != id.generation )
	{
		return NULL;
	}
	return softBody;
}

static b2SoftBody* b2GetSoftBody( b2SoftBodyId id )
{
	if ( B2_IS_NULL( id ) )
	{
		b2Log( "invalid input: null soft body id in %s\n", __func__ );
		return NULL;
	}

	b2World* world = b2GetWorld( id.world0 );
	if ( world == NULL )
	{
		b2Log( "invalid input: bad world index in %s\n", __func__ );
		return NULL;
	}

	b2SoftBody* softBody = b2GetSoftBodyFullId( world, id );
	if ( softBody == NULL )
	{
		b2Log( "invalid input: stale soft body id in %s\n", __func__ );
	}
	return softBody;
}

void b2CreateSoftBodyPool( b2World* world )
{
	b2Array_Create( world->softBodies );
	world->softBodyIdPool = b2CreateIdPool();
	b2Array_Create( world->softLinearImpulses );
	b2Array_Create( world->softAngularImpulses );
	b2Array_Create( world->softImpulseBodyIds );
}

static void b2DestroyGrid( b2SoftGrid* grid )
{
	b2Array_Destroy( grid->bucketHead );
	b2Array_Destroy( grid->nextLink );
	b2Array_Destroy( grid->cellX );
	b2Array_Destroy( grid->cellY );
	b2Array_Destroy( grid->particleBucket );
}

static void b2DestroySoftBodyInternal( b2SoftBody* softBody )
{
	if ( softBody->inUse == false )
	{
		return;
	}
	b2Array_Destroy( softBody->particles );
	b2Array_Destroy( softBody->constraints );
	b2Array_Destroy( softBody->anchors );
	b2Array_Destroy( softBody->triangles );
	b2Array_Destroy( softBody->shapeCandidates );
	b2DestroyGrid( &softBody->grid );
	memset( softBody, 0, sizeof( *softBody ) );
}

void b2DestroySoftBodyPool( b2World* world )
{
	for ( int i = 0; i < world->softBodies.count; ++i )
	{
		b2DestroySoftBodyInternal( world->softBodies.data + i );
	}
	b2Array_Destroy( world->softBodies );
	b2DestroyIdPool( &world->softBodyIdPool );
	b2Array_Destroy( world->softLinearImpulses );
	b2Array_Destroy( world->softAngularImpulses );
	b2Array_Destroy( world->softImpulseBodyIds );
}

// -------------------------------------------------------------------------------------------------
// Mesh validation and construction
// -------------------------------------------------------------------------------------------------

static bool b2SoftFindEdge( const b2SoftEdge* edges, int count, b2SoftEdge edge )
{
	for ( int i = 0; i < count; ++i )
	{
		if ( b2SoftEdgeEquals( edges[i], edge ) )
		{
			return true;
		}
	}
	return false;
}

static bool b2SoftPushEdge( b2SoftEdge* edges, int* count, int capacity, int v1, int v2 )
{
	b2SoftEdge edge = { v1, v2 };
	if ( b2SoftFindEdge( edges, *count, edge ) )
	{
		return true;
	}
	if ( *count >= capacity )
	{
		return false;
	}
	edges[ ( *count )++ ] = edge;
	return true;
}

static float b2SoftTriangleArea( const b2Vec2* vertices, const b2SoftBodyTriangle* t )
{
	b2Vec2 a = vertices[t->vertex1], b = vertices[t->vertex2], c = vertices[t->vertex3];
	return 0.5f * b2Cross( b2Sub( b, a ), b2Sub( c, a ) );
}

static bool b2SoftValidateTriangle( const b2SoftBodyTriangle* t, int vertexCount )
{
	return t->vertex1 >= 0 && t->vertex2 >= 0 && t->vertex3 >= 0 && t->vertex1 < vertexCount &&
		   t->vertex2 < vertexCount && t->vertex3 < vertexCount &&
		   t->vertex1 != t->vertex2 && t->vertex2 != t->vertex3 && t->vertex3 != t->vertex1;
}

static bool b2SoftValidateQuad( const b2SoftBodyQuad* q, int vertexCount )
{
	return q->vertex1 >= 0 && q->vertex2 >= 0 && q->vertex3 >= 0 && q->vertex4 >= 0 &&
		   q->vertex1 < vertexCount && q->vertex2 < vertexCount && q->vertex3 < vertexCount &&
		   q->vertex4 < vertexCount && q->vertex1 != q->vertex2 && q->vertex2 != q->vertex3 &&
		   q->vertex3 != q->vertex4 && q->vertex4 != q->vertex1 && q->vertex1 != q->vertex3 &&
		   q->vertex2 != q->vertex4;
}

b2SoftBodyDef b2DefaultSoftBodyDef( void )
{
	b2SoftBodyDef def = { 0 };
	def.particleMass = 0.1f;
	def.particleRadius = 0.05f;
	def.filter = b2DefaultFilter();
	def.stiffness = 0.9f;
	def.bendStiffness = 0.3f;
	def.damping = 0.02f;
	def.friction = 0.3f;
	def.iterations = 8;
	def.enableSelfCollision = false;
	def.enableTearing = false;
	def.tearMaxStretch = 2.0f;
	def.enableGravity = true;
	def.internalValue = B2_SECRET_COOKIE;
	return def;
}

b2SoftBodyId b2CreateSoftBody( b2WorldId worldId, const b2SoftBodyDef* def )
{
	b2SoftBodyId nullId = b2_nullSoftBodyId;

	B2_CHECK_INPUT_RETURN( def != NULL, nullId );
	B2_CHECK_INPUT_RETURN( def->internalValue == B2_SECRET_COOKIE, nullId );

	b2World* world = b2GetWorldFromId( worldId );
	B2_CHECK_INPUT_RETURN( world != NULL, nullId );
	B2_CHECK_INPUT_RETURN( world->locked == false, nullId );

	// --- Parameter validation ---------------------------------------------------------------
	if ( def->vertices == NULL || def->vertexCount <= 2 || def->vertexCount > B2_SOFT_MAX_PARTICLES )
	{
		b2Log( "soft body error: vertex count must be in [3, %d]\n", B2_SOFT_MAX_PARTICLES );
		B2_ASSERT( false );
		return nullId;
	}
	if ( def->triangleCount < 0 || def->quadCount < 0 ||
		 ( def->triangleCount == 0 && def->quadCount == 0 ) ||
		 def->triangleCount > B2_SOFT_MAX_ELEMENTS || def->quadCount > B2_SOFT_MAX_ELEMENTS )
	{
		b2Log( "soft body error: at least one triangle or quad element is required\n" );
		B2_ASSERT( false );
		return nullId;
	}
	if ( ( def->triangleCount > 0 && def->triangles == NULL ) ||
		 ( def->quadCount > 0 && def->quads == NULL ) )
	{
		b2Log( "soft body error: element count > 0 but element array is NULL\n" );
		B2_ASSERT( false );
		return nullId;
	}
	if ( def->particleMass <= 0.0f || b2IsValidFloat( def->particleMass ) == false )
	{
		b2Log( "soft body error: particle mass must be positive and finite\n" );
		B2_ASSERT( false );
		return nullId;
	}
	if ( def->particleRadius <= 0.0f || b2IsValidFloat( def->particleRadius ) == false )
	{
		b2Log( "soft body error: particle radius must be positive and finite\n" );
		B2_ASSERT( false );
		return nullId;
	}
	if ( def->stiffness < 0.0f || def->stiffness > 1.0f || b2IsValidFloat( def->stiffness ) == false ||
		 def->bendStiffness < 0.0f || def->bendStiffness > 1.0f ||
		 b2IsValidFloat( def->bendStiffness ) == false ||
		 def->damping < 0.0f || def->damping >= 1.0f || def->friction < 0.0f || def->friction > 1.0f )
	{
		b2Log( "soft body error: stiffness/damping/friction parameters out of range\n" );
		B2_ASSERT( false );
		return nullId;
	}
	if ( def->iterations < 1 || def->iterations > 256 )
	{
		b2Log( "soft body error: iterations must be in [1, 256]\n" );
		B2_ASSERT( false );
		return nullId;
	}
	if ( def->enableTearing &&
		 ( def->tearMaxStretch <= 1.0f || def->tearMaxStretch > B2_SOFT_MAX_SANE_STRETCH ||
		   b2IsValidFloat( def->tearMaxStretch ) == false ) )
	{
		b2Log( "soft body error: tearMaxStretch must be in (1, %g]\n", B2_SOFT_MAX_SANE_STRETCH );
		B2_ASSERT( false );
		return nullId;
	}
	if ( def->anchorCount < 0 || ( def->anchorCount > 0 && def->anchors == NULL ) ||
		 def->anchorCount > def->vertexCount )
	{
		b2Log( "soft body error: invalid anchor array\n" );
		B2_ASSERT( false );
		return nullId;
	}

	int vertexCount = def->vertexCount;
	int elementCount = def->triangleCount + 2 * def->quadCount;

	// --- Vertex validation: finite data and no duplicates -----------------------------------
	b2Vec2 minVertex = def->vertices[0];
	b2Vec2 maxVertex = def->vertices[0];
	// Duplicate tolerance: relative to the particle radius, much smaller than any valid edge.
	float duplicateTolerance = 0.25f * def->particleRadius;
	float duplicateToleranceSq = duplicateTolerance * duplicateTolerance;

	for ( int i = 0; i < vertexCount; ++i )
	{
		b2Vec2 v = def->vertices[i];
		if ( b2IsValidVec2( v ) == false )
		{
			b2Log( "soft body error: vertex %d is not finite\n", i );
			B2_ASSERT( false );
			return nullId;
		}
		minVertex.x = b2MinFloat( minVertex.x, v.x );
		minVertex.y = b2MinFloat( minVertex.y, v.y );
		maxVertex.x = b2MaxFloat( maxVertex.x, v.x );
		maxVertex.y = b2MaxFloat( maxVertex.y, v.y );
		for ( int j = 0; j < i; ++j )
		{
			if ( b2LengthSquared( b2Sub( v, def->vertices[j] ) ) <= duplicateToleranceSq )
			{
				b2Log( "soft body error: duplicate vertices %d and %d\n", j, i );
				B2_ASSERT( false );
				return nullId;
			}
		}
	}

	// --- Element validation -----------------------------------------------------------------
	b2SoftBodyTriangle* elements = B2_ALLOC_ARRAY( elementCount, b2SoftBodyTriangle );
	if ( elements == NULL )
	{
		b2Log( "soft body error: out of memory\n" );
		return nullId;
	}

	int elementIndex = 0;
	float minEdgeLength = FLT_MAX;
	for ( int i = 0; i < def->triangleCount; ++i )
	{
		const b2SoftBodyTriangle* t = def->triangles + i;
		if ( b2SoftValidateTriangle( t, vertexCount ) == false )
		{
			b2Log( "soft body error: triangle %d has invalid vertex indices\n", i );
			B2_FREE_ARRAY( elements, elementCount, b2SoftBodyTriangle );
			B2_ASSERT( false );
			return nullId;
		}
		float area = b2SoftTriangleArea( def->vertices, t );
		if ( b2IsValidFloat( area ) == false || b2AbsFloat( area ) <= duplicateToleranceSq )
		{
			b2Log( "soft body error: degenerate triangle %d (zero or invalid area)\n", i );
			B2_FREE_ARRAY( elements, elementCount, b2SoftBodyTriangle );
			B2_ASSERT( false );
			return nullId;
		}
		elements[elementIndex++] = *t;
	}

	for ( int i = 0; i < def->quadCount; ++i )
	{
		const b2SoftBodyQuad* q = def->quads + i;
		if ( b2SoftValidateQuad( q, vertexCount ) == false )
		{
			b2Log( "soft body error: quad %d has invalid vertex indices\n", i );
			B2_FREE_ARRAY( elements, elementCount, b2SoftBodyTriangle );
			B2_ASSERT( false );
			return nullId;
		}
		// Split each quad along the 1-3 diagonal. Both child triangles must be non-degenerate,
		// which also rejects concave and self-intersecting quads.
		b2SoftBodyTriangle t1 = { q->vertex1, q->vertex2, q->vertex3 };
		b2SoftBodyTriangle t2 = { q->vertex1, q->vertex3, q->vertex4 };
		float area1 = b2SoftTriangleArea( def->vertices, &t1 );
		float area2 = b2SoftTriangleArea( def->vertices, &t2 );
		if ( b2IsValidFloat( area1 ) == false || b2IsValidFloat( area2 ) == false ||
			 b2AbsFloat( area1 ) <= duplicateToleranceSq || b2AbsFloat( area2 ) <= duplicateToleranceSq ||
			 ( area1 > 0.0f ) != ( area2 > 0.0f ) )
		{
			b2Log( "soft body error: degenerate or non-convex quad %d\n", i );
			B2_FREE_ARRAY( elements, elementCount, b2SoftBodyTriangle );
			B2_ASSERT( false );
			return nullId;
		}
		elements[elementIndex++] = t1;
		elements[elementIndex++] = t2;
	}

	// --- Collect structural edges (unique) and remember element adjacency for bending -------
	int edgeCapacity = 3 * elementCount;
	b2SoftEdge* edges = B2_ALLOC_ARRAY( edgeCapacity, b2SoftEdge );
	if ( edges == NULL )
	{
		B2_FREE_ARRAY( elements, elementCount, b2SoftBodyTriangle );
		b2Log( "soft body error: out of memory\n" );
		return nullId;
	}
	int edgeCount = 0;

	for ( int i = 0; i < elementCount; ++i )
	{
		const b2SoftBodyTriangle* t = elements + i;
		int triEdges[3][2] = {
			{ t->vertex1, t->vertex2 }, { t->vertex2, t->vertex3 }, { t->vertex3, t->vertex1 }
		};
		for ( int e = 0; e < 3; ++e )
		{
			int v1 = triEdges[e][0], v2 = triEdges[e][1];
			float length = b2Length( b2Sub( def->vertices[v1], def->vertices[v2] ) );
			minEdgeLength = b2MinFloat( minEdgeLength, length );
			if ( b2SoftPushEdge( edges, &edgeCount, edgeCapacity, v1, v2 ) == false )
			{
				B2_FREE_ARRAY( edges, edgeCapacity, b2SoftEdge );
				B2_FREE_ARRAY( elements, elementCount, b2SoftBodyTriangle );
				b2Log( "soft body error: constraint count exceeded internal limit\n" );
				return nullId;
			}
		}
	}

	if ( edgeCount == 0 || minEdgeLength <= duplicateTolerance )
	{
		B2_FREE_ARRAY( edges, edgeCapacity, b2SoftEdge );
		B2_FREE_ARRAY( elements, elementCount, b2SoftBodyTriangle );
		b2Log( "soft body error: mesh has no valid edges\n" );
		B2_ASSERT( false );
		return nullId;
	}

	// Bending pairs: the unconnected vertex pairs inside each triangle. Shared edges are already
	// structural, so adding all three missing pairs of every triangle covers every diagonal.
	int bendCapacity = 3 * elementCount;
	b2SoftEdge* bendEdges = B2_ALLOC_ARRAY( bendCapacity, b2SoftEdge );
	if ( bendEdges == NULL )
	{
		B2_FREE_ARRAY( edges, edgeCapacity, b2SoftEdge );
		B2_FREE_ARRAY( elements, elementCount, b2SoftBodyTriangle );
		b2Log( "soft body error: out of memory\n" );
		return nullId;
	}
	int bendCount = 0;

	// Bending constraints connect the two opposite vertices of every pair of triangles sharing a
	// structural edge. Each edge is scanned against all triangles, which is O(edges * elements)
	// and therefore linear in the mesh size. Boundary edges have a single neighbor and produce no
	// bending constraint.
	for ( int e = 0; e < edgeCount; ++e )
	{
		int shared1 = edges[e].particle1;
		int shared2 = edges[e].particle2;
		int opposite[2] = { B2_NULL_INDEX, B2_NULL_INDEX };
		int oppositeCount = 0;
		for ( int i = 0; i < elementCount && oppositeCount < 2; ++i )
		{
			const b2SoftBodyTriangle* t = elements + i;
			int tv[3] = { t->vertex1, t->vertex2, t->vertex3 };
			bool has1 = false, has2 = false;
			int other = B2_NULL_INDEX;
			for ( int k = 0; k < 3; ++k )
			{
				if ( tv[k] == shared1 )
				{
					has1 = true;
				}
				else if ( tv[k] == shared2 )
				{
					has2 = true;
				}
				else
				{
					other = tv[k];
				}
			}
			if ( has1 && has2 && other != B2_NULL_INDEX )
			{
				// Skip duplicates when more than two triangles share an edge (non-manifold mesh).
				if ( oppositeCount == 0 || opposite[0] != other )
				{
					opposite[oppositeCount++] = other;
				}
			}
		}
		if ( oppositeCount == 2 && opposite[0] != opposite[1] )
		{
			b2SoftEdge bend = { opposite[0], opposite[1] };
			if ( b2SoftFindEdge( edges, edgeCount, bend ) == false &&
				 b2SoftFindEdge( bendEdges, bendCount, bend ) == false )
			{
				if ( bendCount < bendCapacity )
				{
					bendEdges[bendCount++] = bend;
				}
			}
		}
	}

	// --- Allocate the soft body --------------------------------------------------------------
	int softId = b2AllocId( &world->softBodyIdPool );
	if ( softId >= world->softBodies.count )
	{
		b2Array_Resize( world->softBodies, softId + 1 );
	}
	b2SoftBody* softBody = world->softBodies.data + softId;
	memset( softBody, 0, sizeof( *softBody ) );
	softBody->id = softId;
	softBody->generation = (uint16_t)( softBody->generation + ( softBody->generation == 0 ? 1 : 0 ) );
	if ( softBody->generation == 0 )
	{
		softBody->generation = 1;
	}
	softBody->inUse = true;

	// Typical sub-step at 60Hz with 4 sub-steps.

	b2Array_CreateN( softBody->particles, vertexCount );
	b2Array_Resize( softBody->particles, vertexCount );
	float invMass = 1.0f / def->particleMass;
	for ( int i = 0; i < vertexCount; ++i )
	{
		b2SoftParticle* p = softBody->particles.data + i;
		p->position = (b2Pos){ def->vertices[i].x, def->vertices[i].y };
		p->prevPosition = p->position;
		p->velocity = (b2Vec2){ 0.0f, 0.0f };
		p->accumulatedImpulse = (b2Vec2){ 0.0f, 0.0f };
		p->invMass = invMass;
	}

	int totalConstraints = edgeCount + bendCount;
	if ( totalConstraints > B2_SOFT_MAX_CONSTRAINTS )
	{
		b2FreeId( &world->softBodyIdPool, softId );
		softBody->inUse = false;
		B2_FREE_ARRAY( bendEdges, bendCapacity, b2SoftEdge );
		B2_FREE_ARRAY( edges, edgeCapacity, b2SoftEdge );
		B2_FREE_ARRAY( elements, elementCount, b2SoftBodyTriangle );
		b2Log( "soft body error: constraint count %d exceeds limit %d\n", totalConstraints,
			   B2_SOFT_MAX_CONSTRAINTS );
		return nullId;
	}

	b2Array_CreateN( softBody->constraints, totalConstraints );
	b2Array_Resize( softBody->constraints, totalConstraints );
	int constraintIndex = 0;
	for ( int i = 0; i < edgeCount; ++i )
	{
		b2SoftConstraint* c = softBody->constraints.data + constraintIndex++;
		c->particle1 = edges[i].particle1;
		c->particle2 = edges[i].particle2;
		c->restLength = b2Length( b2Sub( def->vertices[c->particle1], def->vertices[c->particle2] ) );
		c->stiffnessFactor = b2SoftCompliance( def->stiffness );
		c->type = b2_softDistanceConstraint;
		c->active = true;
	}
	for ( int i = 0; i < bendCount; ++i )
	{
		b2SoftConstraint* c = softBody->constraints.data + constraintIndex++;
		c->particle1 = bendEdges[i].particle1;
		c->particle2 = bendEdges[i].particle2;
		c->restLength = b2Length( b2Sub( def->vertices[c->particle1], def->vertices[c->particle2] ) );
		c->stiffnessFactor = b2SoftCompliance( def->bendStiffness );
		c->type = b2_softBendingConstraint;
		c->active = def->bendStiffness > 0.0f;
	}

	b2Array_CreateN( softBody->triangles, elementCount );
	b2Array_Resize( softBody->triangles, elementCount );
	memcpy( softBody->triangles.data, elements, (size_t)elementCount * sizeof( b2SoftBodyTriangle ) );

	// --- Anchors -----------------------------------------------------------------------------
	if ( def->anchorCount > 0 )
	{
		b2Array_CreateN( softBody->anchors, def->anchorCount );
		b2Array_Resize( softBody->anchors, def->anchorCount );
	}
	else
	{
		b2Array_Create( softBody->anchors );
	}
	for ( int i = 0; i < def->anchorCount; ++i )
	{
		const b2SoftBodyAnchor* src = def->anchors + i;
		b2SoftAnchor* dst = softBody->anchors.data + i;
		if ( src->particleIndex < 0 || src->particleIndex >= vertexCount ||
			 b2IsValidVec2( src->localPoint ) == false )
		{
			// Roll back the whole creation rather than silently dropping an anchor.
			b2DestroySoftBodyInternal( softBody );
			b2FreeId( &world->softBodyIdPool, softId );
			B2_FREE_ARRAY( bendEdges, bendCapacity, b2SoftEdge );
			B2_FREE_ARRAY( edges, edgeCapacity, b2SoftEdge );
			B2_FREE_ARRAY( elements, elementCount, b2SoftBodyTriangle );
			b2Log( "soft body error: anchor %d references an invalid particle\n", i );
			B2_ASSERT( false );
			return nullId;
		}
		dst->particleIndex = src->particleIndex;
		dst->detached = false;
		if ( B2_IS_NULL( src->bodyId ) )
		{
			dst->bodyId = B2_NULL_INDEX;
			dst->localPoint = src->localPoint;
			dst->lastWorldPoint = (b2Pos){ src->localPoint.x, src->localPoint.y };
		}
		else
		{
			b2Body* body = b2GetBodyFullId( world, src->bodyId );
			if ( body == NULL || body->id == B2_NULL_INDEX )
			{
				b2DestroySoftBodyInternal( softBody );
				b2FreeId( &world->softBodyIdPool, softId );
				B2_FREE_ARRAY( bendEdges, bendCapacity, b2SoftEdge );
				B2_FREE_ARRAY( edges, edgeCapacity, b2SoftEdge );
				B2_FREE_ARRAY( elements, elementCount, b2SoftBodyTriangle );
				b2Log( "soft body error: anchor %d references an invalid body\n", i );
				B2_ASSERT( false );
				return nullId;
			}
			dst->bodyId = body->id;
			dst->localPoint = src->localPoint;
			b2WorldTransform transform = b2GetBodyTransformQuick( world, body );
			dst->lastWorldPoint = b2TransformWorldPoint( transform, src->localPoint );
		}
	}

	// --- Remaining scalar state --------------------------------------------------------------
	softBody->filter = def->filter;
	softBody->particleRadius = def->particleRadius;
	softBody->stiffness = def->stiffness;
	softBody->bendStiffness = def->bendStiffness;
	softBody->damping = def->damping;
	softBody->friction = def->friction;
	softBody->iterations = def->iterations;
	softBody->enableSelfCollision = def->enableSelfCollision;
	softBody->enableTearing = def->enableTearing;
	softBody->tearMaxStretch = def->tearMaxStretch;
	softBody->enableGravity = def->enableGravity;
	softBody->userData = def->userData;
	softBody->bounds.lowerBound = minVertex;
	softBody->bounds.upperBound = maxVertex;
	softBody->stats = (b2SoftBodyStats){ 0 };
	softBody->stats.particleCount = vertexCount;
	softBody->stats.constraintCount = totalConstraints;

	// Pinned particles get infinite mass.
	for ( int i = 0; i < softBody->anchors.count; ++i )
	{
		softBody->particles.data[softBody->anchors.data[i].particleIndex].invMass = 0.0f;
	}

	// --- Self-collision grid, sized once -----------------------------------------------------
	// A fixed power-of-two bucket table at a low load factor plus an intrusive per-particle link.
	// Both stay resident; rebuilding each sub-step only rewrites their contents and never
	// allocates.
	int bucketCount = 8;
	while ( bucketCount < 2 * vertexCount )
	{
		bucketCount *= 2;
	}
	b2SoftGrid* grid = &softBody->grid;
	b2Array_CreateN( grid->bucketHead, bucketCount );
	b2Array_Resize( grid->bucketHead, bucketCount );
	b2Array_CreateN( grid->nextLink, vertexCount );
	b2Array_Resize( grid->nextLink, vertexCount );
	b2Array_CreateN( grid->cellX, vertexCount );
	b2Array_Resize( grid->cellX, vertexCount );
	b2Array_CreateN( grid->cellY, vertexCount );
	b2Array_Resize( grid->cellY, vertexCount );
	b2Array_CreateN( grid->particleBucket, vertexCount );
	b2Array_Resize( grid->particleBucket, vertexCount );
	grid->bucketCount = bucketCount;
	grid->bucketMask = bucketCount - 1;
	grid->cellSize = b2MaxFloat( 2.0f * def->particleRadius, minEdgeLength );

	// Rigid collision candidate scratch. Preallocated generously so a broad-phase query never
	// allocates while stepping. b2SoftTreeCallback stops appending once the capacity is reached.
	int candidateCapacity = b2MaxInt( 32, 4 * vertexCount );
	b2Array_CreateN( softBody->shapeCandidates, candidateCapacity );
	softBody->particleMass = def->particleMass;

	// Size the per-body coupling impulse arrays to cover every body currently in the world so the
	// step loop never has to grow them.
	b2EnsureSoftImpulseCapacity( world, b2MaxInt( 0, world->bodies.capacity - 1 ) );

	B2_FREE_ARRAY( bendEdges, bendCapacity, b2SoftEdge );
	B2_FREE_ARRAY( edges, edgeCapacity, b2SoftEdge );
	B2_FREE_ARRAY( elements, elementCount, b2SoftBodyTriangle );

	return b2MakeSoftBodyId( world, softId );
}

// -------------------------------------------------------------------------------------------------
// Solver
// -------------------------------------------------------------------------------------------------

typedef struct b2SoftQueryContext
{
	b2World* world;
	b2SoftBody* softBody;
} b2SoftQueryContext;

static bool b2SoftTreeCallback( int proxyId, uint64_t userData, void* context )
{
	B2_UNUSED( proxyId );
	b2SoftQueryContext* queryContext = context;
	b2World* world = queryContext->world;
	int shapeId = (int)userData;
	b2Shape* shape = world->shapes.data + shapeId;

	// Standard category/mask/group filtering, matching rigid contacts.
	if ( b2ShouldShapesCollide( queryContext->softBody->filter, shape->filter ) == false )
	{
		return true;
	}

	// World custom filter callback. The soft body has no shape id of its own, so the rigid shape
	// id is passed as both arguments; applications distinguish soft body queries by user data.
	if ( shape->enableCustomFiltering && world->customFilterFcn != NULL )
	{
		b2ShapeId publicId = { shapeId + 1, world->worldId, shape->generation };
		if ( world->customFilterFcn( publicId, publicId, world->customFilterContext ) == false )
		{
			return true;
		}
	}

	// Bounded push: the scratch array is preallocated at creation, so stepping never allocates.
	// An overflow is recorded as a numerical failure instead of dropping silently or exploding.
	if ( queryContext->softBody->shapeCandidates.count < queryContext->softBody->shapeCandidates.capacity )
	{
		b2Array_Push( queryContext->softBody->shapeCandidates, shapeId );
	}
	else
	{
		queryContext->softBody->stats.numericalFailure = true;
	}
	return true;
}

// Resolve particle versus one rigid shape using GJK distance. On contact, fills the world-space
// outward normal (pointing from shape to particle) and the penetration depth.
static bool b2SoftCollideShape( b2World* world, int shapeId, b2Vec2 particleCenter, float particleRadius,
								b2Vec2* normalOut, float* penetrationOut )
{
	b2Shape* shape = world->shapes.data + shapeId;
	b2Body* body = world->bodies.data + shape->bodyId;
	b2WorldTransform shapeTransform = b2GetBodyTransformQuick( world, body );

	b2ShapeProxy shapeProxy = b2MakeShapeDistanceProxy( shape );

	// Distance runs in the shape's local frame, so the particle center is expressed relative to it.
	b2WorldTransform particleFrame = { .p = (b2Pos){ particleCenter.x, particleCenter.y },
									   .q = { 1.0f, 0.0f } };
	b2Transform relative = b2InvMulWorldTransforms( shapeTransform, particleFrame );

	// First query with radii so the inflated particle detects grazing contact.
	b2ShapeProxy particleProxy = { 0 };
	particleProxy.count = 1;
	particleProxy.radius = particleRadius;
	particleProxy.points[0] = (b2Vec2){ 0.0f, 0.0f };

	b2DistanceInput input = { 0 };
	input.proxyA = shapeProxy;
	input.proxyB = particleProxy;
	input.transform = relative;
	input.useRadii = true;

	b2SimplexCache cache = { 0 };
	b2DistanceOutput output = b2ComputeShapeDistance( &input, &cache, NULL, 0 );
	if ( b2IsValidFloat( output.distance ) == false )
	{
		return false;
	}

	if ( output.distance > 0.0f )
	{
		// Inflated proxies are separated: no contact.
		return false;
	}

	b2Vec2 localNormal;
	float penetration;

	if ( b2LengthSquared( output.normal ) > FLT_EPSILON )
	{
		// Shallow overlap: GJK returned a usable direction. The inflated distance is zero, so the
		// penetration follows from how far the bare particle center sits inside the shape surface.
		localNormal = output.normal;
		b2ShapeProxy bareParticle = particleProxy;
		bareParticle.radius = 0.0f;
		b2DistanceInput bareInput = { 0 };
		bareInput.proxyA = shapeProxy;
		bareInput.proxyB = bareParticle;
		bareInput.transform = relative;
		bareInput.useRadii = true;
		b2SimplexCache bareCache = { 0 };
		b2DistanceOutput bare = b2ComputeShapeDistance( &bareInput, &bareCache, NULL, 0 );
		// Bare center distance d means the inflated particle overlaps by (particleRadius - d).
		penetration = particleRadius - b2MaxFloat( bare.distance, 0.0f );
	}
	else
	{
		// Deep overlap: the radius-based GJK cannot produce a direction. Re-run without the shape
		// radius to get the closest point on the actual shape surface, then derive the normal from
		// that surface point toward the particle center.
		b2ShapeProxy bareShape = shapeProxy;
		bareShape.radius = 0.0f;
		b2ShapeProxy bareParticle = particleProxy;
		bareParticle.radius = 0.0f;
		b2DistanceInput deepInput = { 0 };
		deepInput.proxyA = bareShape;
		deepInput.proxyB = bareParticle;
		deepInput.transform = relative;
		deepInput.useRadii = false;
		b2SimplexCache deepCache = { 0 };
		b2DistanceOutput deep = b2ComputeShapeDistance( &deepInput, &deepCache, NULL, 0 );
		if ( b2IsValidFloat( deep.distance ) == false )
		{
			return false;
		}

		if ( b2LengthSquared( deep.normal ) > FLT_EPSILON )
		{
			localNormal = deep.normal;
			penetration = particleRadius - b2MaxFloat( deep.distance, 0.0f );
		}
		else
		{
			// The particle center is inside the shape core. Fall back to the vector from the
			// closest shape surface point to the particle center.
			b2Vec2 toCenter = b2Sub( relative.p, deep.pointA );
			localNormal = b2Normalize( toCenter );
			if ( b2LengthSquared( localNormal ) < FLT_EPSILON )
			{
				// Center coincides with the closest point; use the body frame up axis.
				localNormal = (b2Vec2){ 0.0f, 1.0f };
			}
			penetration = particleRadius + b2Length( toCenter );
		}
	}

	// Rotate the local-frame normal into world space.
	b2Vec2 normal = b2RotateVector( shapeTransform.q, localNormal );
	if ( b2LengthSquared( normal ) < FLT_EPSILON )
	{
		normal = (b2Vec2){ 0.0f, 1.0f };
	}

	*normalOut = normal;
	*penetrationOut = b2MaxFloat( penetration, 0.0f );
	return true;
}

void b2EnsureSoftImpulseCapacity( b2World* world, int bodyId )
{
	int needed = bodyId + 1;
	if ( world->softLinearImpulses.capacity >= needed )
	{
		return;
	}
	int oldCapacity = world->softLinearImpulses.capacity;
	int newCapacity = oldCapacity == 0 ? 64 : oldCapacity;
	while ( newCapacity < needed )
	{
		newCapacity *= 2;
	}
	b2Array_Reserve( world->softLinearImpulses, newCapacity );
	b2Array_Reserve( world->softAngularImpulses, newCapacity );
	b2Array_Reserve( world->softImpulseBodyIds, newCapacity );
	b2Vec2 zeroVec = { 0.0f, 0.0f };
	for ( int i = oldCapacity; i < newCapacity; ++i )
	{
		b2Array_Push( world->softLinearImpulses, zeroVec );
		b2Array_Push( world->softAngularImpulses, 0.0f );
		b2Array_Push( world->softImpulseBodyIds, 0 );
	}
}

// Rebuild the spatial hash from current particle positions using an intrusive linked list per
// bucket. O(particles): one pass to clear the heads and one pass to link every particle. No
// allocation and no probing, so the cost stays independent of mesh density.
static void b2SoftBuildGrid( b2SoftBody* softBody )
{
	b2SoftGrid* grid = &softBody->grid;
	int particleCount = softBody->particles.count;

	for ( int b = 0; b < grid->bucketCount; ++b )
	{
		grid->bucketHead.data[b] = -1;
	}

	for ( int i = 0; i < particleCount; ++i )
	{
		b2Vec2 p = b2SubPos( softBody->particles.data[i].position, (b2Pos){ 0 } );
		int cx = (int)floorf( p.x / grid->cellSize );
		int cy = (int)floorf( p.y / grid->cellSize );
		int bucket = (int)( b2SoftHashPair( cx, cy ) & (uint32_t)grid->bucketMask );

		grid->cellX.data[i] = cx;
		grid->cellY.data[i] = cy;
		grid->particleBucket.data[i] = bucket;
		grid->nextLink.data[i] = grid->bucketHead.data[bucket];
		grid->bucketHead.data[bucket] = i;
	}
}

// PBD distance constraint solve. Returns the absolute stretch ratio for diagnostics and tearing.
static float b2SoftSolveConstraint( b2SoftBody* softBody, b2SoftConstraint* c )
{
	b2SoftParticle* p1 = softBody->particles.data + c->particle1;
	b2SoftParticle* p2 = softBody->particles.data + c->particle2;

	b2Vec2 delta = b2Sub( b2SubPos( p2->position, (b2Pos){ 0 } ), b2SubPos( p1->position, (b2Pos){ 0 } ) );
	if ( b2IsValidVec2( delta ) == false )
	{
		softBody->stats.numericalFailure = true;
		return 1.0f;
	}
	float length = b2Length( delta );
	if ( length < FLT_EPSILON || c->restLength < FLT_EPSILON )
	{
		// Coincident particles cannot be separated along a direction; skip the constraint rather
		// than dividing by zero.
		return 1.0f;
	}

	float stretch = length / c->restLength;
	float C = length - c->restLength;

	float w1 = p1->invMass, w2 = p2->invMass;
	float wSum = w1 + w2;
	if ( wSum <= 0.0f )
	{
		return stretch;
	}

	// Classic PBD position correction. n points from p1 to p2. When stretched (C > 0) the
	// particles are pulled together: p1 moves along +n, p2 along -n, each weighted by its
	// inverse-mass share. Scaled by the per-iteration stiffness factor.
	float scale = C * c->stiffnessFactor / wSum;
	b2Vec2 n = b2MulSV( 1.0f / length, delta );

	if ( w1 > 0.0f )
	{
		p1->position = b2OffsetPos( p1->position, b2MulSV( w1 * scale, n ) );
	}
	if ( w2 > 0.0f )
	{
		p2->position = b2OffsetPos( p2->position, b2MulSV( -w2 * scale, n ) );
	}
	return stretch;
}

// Rigid collision for one particle. The broad-phase query uses the swept AABB spanning the
// particle's motion since the previous sub-step (via prevPosition), so a fast particle cannot
// tunnel through thin shapes. The contact itself is an XPBD-style position correction with the
// rigid body treated as having finite inverse mass. Two-way impulses are staged on the world.
static void b2SoftSolveRigidContacts( b2World* world, b2SoftBody* softBody, int particleIndex, float h )
{
	b2SoftParticle* particle = softBody->particles.data + particleIndex;
	if ( particle->invMass <= 0.0f )
	{
		return;
	}

	b2Vec2 center = b2SubPos( particle->position, (b2Pos){ 0 } );
	b2Vec2 prevCenter = b2SubPos( particle->prevPosition, (b2Pos){ 0 } );
	float radius = softBody->particleRadius;

	// Swept box covering the whole sub-step motion, inflated by the particle radius.
	b2Vec2 boxLower = {
		b2MinFloat( center.x, prevCenter.x ) - radius,
		b2MinFloat( center.y, prevCenter.y ) - radius,
	};
	b2Vec2 boxUpper = {
		b2MaxFloat( center.x, prevCenter.x ) + radius,
		b2MaxFloat( center.y, prevCenter.y ) + radius,
	};
	b2AABB sweptBox;
	sweptBox.lowerBound = boxLower;
	sweptBox.upperBound = boxUpper;

	// Reuse the query machinery with a swept box instead of a point-sized box.
	b2Array_Clear( softBody->shapeCandidates );
	b2SoftQueryContext context = { world, softBody };
	for ( int i = 0; i < b2_bodyTypeCount; ++i )
	{
		b2DynamicTree_OverlapAABB( world->broadPhase.trees + i, sweptBox, softBody->filter.maskBits,
								   b2SoftTreeCallback, &context );
	}

	for ( int k = 0; k < softBody->shapeCandidates.count; ++k )
	{
		int shapeId = softBody->shapeCandidates.data[k];
		b2Vec2 normal;
		float penetration;
		if ( b2SoftCollideShape( world, shapeId, center, radius, &normal, &penetration ) == false )
		{
			// The swept box overlaps the shape but the current particle center does not: the
			// particle crossed the surface this sub-step (tunneling). Recover the surface plane
			// from a query at the previous (outside) position and clamp onto it.
			b2Vec2 prevNormal;
			float prevPen;
			if ( b2SoftCollideShape( world, shapeId, prevCenter, radius, &prevNormal, &prevPen ) )
			{
				normal = prevNormal;
				penetration = 0.0f;
				b2Vec2 clamped = b2Add( prevCenter, b2MulSV( radius, prevNormal ) );
				particle->position = (b2Pos){ clamped.x, clamped.y };
			}
			continue;
		}

		softBody->stats.rigidContactCount += 1;

		b2Shape* shape = world->shapes.data + shapeId;
		b2Body* body = world->bodies.data + shape->bodyId;
		float bodyInvMass = 0.0f;
		if ( body->type == b2_dynamicBody && body->setIndex != b2_disabledSet )
		{
			b2SolverSet* set = b2Array_Get( world->solverSets, body->setIndex );
			if ( body->setIndex == b2_awakeSet )
			{
				b2BodySim* bodySim = b2Array_Get( set->bodySims, body->localIndex );
				bodyInvMass = bodySim->invMass;
			}
		}

		float wSum = particle->invMass + bodyInvMass;
		if ( wSum <= 0.0f )
		{
			// Static/kinematic wall: the particle takes the whole correction.
			wSum = particle->invMass;
		}

		float correctionScale = penetration / wSum;
		b2Vec2 particleCorrection = b2MulSV( particle->invMass * correctionScale, normal );
		particle->position = b2OffsetPos( particle->position, particleCorrection );

		if ( bodyInvMass > 0.0f )
		{
			// Stage the equal-and-opposite impulse on the body. Impulse = m * dx / h. The impulse
			// arrays are sized to the world body array at soft body creation, so this never
			// allocates during a step. Bodies created later grow the array outside of stepping.
			b2Vec2 bodyCorrection = b2MulSV( -bodyInvMass * correctionScale, normal );
			b2Vec2 impulse = b2MulSV( 1.0f / h, b2MulSV( 1.0f / bodyInvMass, bodyCorrection ) );

			B2_ASSERT( body->id < world->softLinearImpulses.capacity );
			if ( body->id < world->softLinearImpulses.capacity )
			{
				world->softLinearImpulses.data[body->id] =
					b2Add( world->softLinearImpulses.data[body->id], impulse );

				// Angular impulse about the body center of mass.
				b2SolverSet* set = b2Array_Get( world->solverSets, body->setIndex );
				b2BodySim* bodySim = b2Array_Get( set->bodySims, body->localIndex );
				b2Vec2 r = b2Sub( center, b2SubPos( bodySim->center, (b2Pos){ 0 } ) );
				world->softAngularImpulses.data[body->id] += b2Cross( r, impulse );
				world->softImpulseBodyIds.data[body->id] = 1;
			}
			else
			{
				softBody->stats.numericalFailure = true;
			}
		}

		// Friction: project the tangential motion onto the surface using a Coulomb clamp.
		b2Vec2 tangent = b2CrossSV( 1.0f, normal );
		b2Vec2 relVelocity = particle->velocity;
		if ( bodyInvMass > 0.0f )
		{
			b2SolverSet* set = b2Array_Get( world->solverSets, body->setIndex );
			b2BodyState* state = b2Array_Get( set->bodyStates, body->localIndex );
			relVelocity = b2Sub( relVelocity, state->linearVelocity );
		}
		float vn = b2Dot( relVelocity, normal );
		if ( vn < 0.0f )
		{
			float vt = b2Dot( relVelocity, tangent );
			float maxFriction = softBody->friction * b2AbsFloat( vn );
			float frictionImpulse = b2ClampFloat( -vt, -maxFriction, maxFriction );
			particle->velocity = b2MulAdd( particle->velocity, particle->invMass,
										   b2MulSV( frictionImpulse, tangent ) );
		}
	}
}

// Particle-particle contact using the spatial hash. For each particle the nine neighboring cells
// are examined by walking their bucket lists; hash collisions are rejected by comparing cell
// coordinates. Each pair is handled once (only higher-indexed neighbors), keeping the result
// deterministic and O(particles).
static void b2SoftSolveSelfContacts( b2SoftBody* softBody )
{
	b2SoftGrid* grid = &softBody->grid;
	float minDist = 2.0f * softBody->particleRadius;
	float minDistSq = minDist * minDist;
	int particleCount = softBody->particles.count;

	for ( int i = 0; i < particleCount; ++i )
	{
		b2SoftParticle* p1 = softBody->particles.data + i;
		if ( p1->invMass <= 0.0f )
		{
			continue;
		}
		b2Vec2 pos1 = b2SubPos( p1->position, (b2Pos){ 0 } );
		int cx = grid->cellX.data[i];
		int cy = grid->cellY.data[i];

		for ( int ox = -1; ox <= 1; ++ox )
		{
			for ( int oy = -1; oy <= 1; ++oy )
			{
				int nx = cx + ox;
				int ny = cy + oy;
				int bucket = (int)( b2SoftHashPair( nx, ny ) & (uint32_t)grid->bucketMask );

				int j = grid->bucketHead.data[bucket];
				while ( j != -1 )
				{
					// Reject particles that merely share a hash bucket but live in another cell.
					if ( grid->cellX.data[j] == nx && grid->cellY.data[j] == ny && j > i )
					{
						b2SoftParticle* p2 = softBody->particles.data + j;
						b2Vec2 pos2 = b2SubPos( p2->position, (b2Pos){ 0 } );
						b2Vec2 d = b2Sub( pos2, pos1 );
						float distSq = b2LengthSquared( d );
						if ( distSq < minDistSq && distSq >= FLT_EPSILON )
						{
							softBody->stats.selfContactCount += 1;
							float dist = sqrtf( distSq );
							b2Vec2 n = b2MulSV( 1.0f / dist, d );
							float overlap = minDist - dist;
							float wSum = p1->invMass + p2->invMass;
							if ( wSum > 0.0f )
							{
								float scale = overlap / wSum;
								p1->position = b2OffsetPos( p1->position, b2MulSV( -p1->invMass * scale, n ) );
								p2->position = b2OffsetPos( p2->position, b2MulSV( p2->invMass * scale, n ) );
							}
						}
					}
					j = grid->nextLink.data[j];
				}
			}
		}
	}
}

// -------------------------------------------------------------------------------------------------
// Per-step integration
// -------------------------------------------------------------------------------------------------

// Apply anchor targets before the constraint solve. World anchors pin at their stored point, body
// anchors follow the current body pose. Anchors whose body was destroyed fall back to the last
// known world point instead of dereferencing a stale body id.
static void b2SoftApplyAnchors( b2World* world, b2SoftBody* softBody )
{
	for ( int i = 0; i < softBody->anchors.count; ++i )
	{
		b2SoftAnchor* anchor = softBody->anchors.data + i;
		b2SoftParticle* particle = softBody->particles.data + anchor->particleIndex;

		b2Pos target = anchor->lastWorldPoint;
		if ( anchor->detached == false && anchor->bodyId != B2_NULL_INDEX )
		{
			b2Body* body = world->bodies.data + anchor->bodyId;
			if ( body->id == B2_NULL_INDEX )
			{
				// The body was destroyed while the soft body lived on. Pin at the last point and
				// mark the anchor failed so no stale dereference can happen later.
				anchor->detached = true;
			}
			else
			{
				b2WorldTransform transform = b2GetBodyTransformQuick( world, body );
				target = b2TransformWorldPoint( transform, anchor->localPoint );
				anchor->lastWorldPoint = target;
			}
		}

		particle->position = target;
	}
}

// One PBD sub-step.
static void b2SoftStep( b2World* world, b2SoftBody* softBody, float h, float invH )
{
	int particleCount = softBody->particles.count;
	int constraintCount = softBody->constraints.count;
	float dampingFactor = 1.0f - softBody->damping;

	// 1. Integrate: apply gravity, external impulses, and damping; predict positions.
	for ( int i = 0; i < particleCount; ++i )
	{
		b2SoftParticle* p = softBody->particles.data + i;
		p->prevPosition = p->position;

		if ( p->invMass <= 0.0f )
		{
			continue;
		}

		if ( softBody->enableGravity )
		{
			p->velocity.y += world->gravity.y * h;
			p->velocity.x += world->gravity.x * h;
		}

		p->velocity = b2MulSV( dampingFactor, p->velocity );
		p->velocity = b2Add( p->velocity, b2MulSV( p->invMass, p->accumulatedImpulse ) );
		p->accumulatedImpulse = (b2Vec2){ 0.0f, 0.0f };

		// Runaway guard before the position is committed.
		float moveSq = b2LengthSquared( b2MulSV( h, p->velocity ) );
		float maxMove = B2_SOFT_RUNAWAY_RADII * softBody->particleRadius;
		if ( b2IsValidVec2( p->velocity ) == false || moveSq > maxMove * maxMove )
		{
			softBody->stats.numericalFailure = true;
			if ( b2IsValidVec2( p->velocity ) == false )
			{
				p->velocity = (b2Vec2){ 0.0f, 0.0f };
			}
			else
			{
				float move = sqrtf( moveSq );
				p->velocity = b2MulSV( maxMove / move, p->velocity );
			}
		}

		p->position = b2OffsetPos( p->position, b2MulSV( h, p->velocity ) );
		if ( b2IsValidPos( p->position ) == false )
		{
			softBody->stats.numericalFailure = true;
			p->position = p->prevPosition;
			p->velocity = (b2Vec2){ 0.0f, 0.0f };
		}
	}

	// 2. Pin anchored particles to their current targets.
	b2SoftApplyAnchors( world, softBody );

	// Refresh the per-constraint stiffness factors for this sub-step.
	float distanceStiffness = b2SoftCompliance( softBody->stiffness );
	float bendStiffness = b2SoftCompliance( softBody->bendStiffness );
	for ( int i = 0; i < constraintCount; ++i )
	{
		b2SoftConstraint* c = softBody->constraints.data + i;
		c->stiffnessFactor = ( c->type == b2_softDistanceConstraint ) ? distanceStiffness : bendStiffness;
	}

	// 3+4. Constraint and collision iterations.
	float stepMaxStretch = 1.0f;
	if ( softBody->enableSelfCollision )
	{
		b2SoftBuildGrid( softBody );
	}

	for ( int iteration = 0; iteration < softBody->iterations; ++iteration )
	{
		for ( int i = 0; i < constraintCount; ++i )
		{
			b2SoftConstraint* c = softBody->constraints.data + i;
			if ( c->active == false )
			{
				continue;
			}

			float stretch = b2SoftSolveConstraint( softBody, c );
			if ( c->type == b2_softDistanceConstraint )
			{
				stepMaxStretch = b2MaxFloat( stepMaxStretch, stretch );

				// Tearing: remove the constraint permanently once it crosses the threshold.
				if ( softBody->enableTearing && stretch > softBody->tearMaxStretch )
				{
					c->active = false;
					softBody->stats.tornCount += 1;
				}
			}
		}

		// Re-pin after the constraint pass so anchors always win.
		b2SoftApplyAnchors( world, softBody );

		// Rigid collisions. Each particle queries the broad phase against its current predicted
		// position; the candidate array is preallocated so the query never allocates.
		for ( int i = 0; i < particleCount; ++i )
		{
			b2SoftSolveRigidContacts( world, softBody, i, h );
		}
		b2SoftApplyAnchors( world, softBody );

		if ( softBody->enableSelfCollision )
		{
			b2SoftSolveSelfContacts( softBody );
			b2SoftApplyAnchors( world, softBody );
		}
	}

	// 5. Derive velocities from the position change, then clamp against runaway feedback. The
	// PBD position solve can produce large corrections on stiff meshes; without this clamp the
	// derived velocity feeds back into the next sub-step and diverges.
	float maxSubStepSpeed = B2_SOFT_RUNAWAY_RADII * softBody->particleRadius * invH;
	for ( int i = 0; i < particleCount; ++i )
	{
		b2SoftParticle* p = softBody->particles.data + i;
		b2Vec2 delta = b2SubPos( p->position, p->prevPosition );
		p->velocity = b2MulSV( invH, delta );

		if ( b2IsValidVec2( p->velocity ) == false )
		{
			softBody->stats.numericalFailure = true;
			p->velocity = (b2Vec2){ 0.0f, 0.0f };
			continue;
		}

		float speedSq = b2LengthSquared( p->velocity );
		if ( speedSq > maxSubStepSpeed * maxSubStepSpeed )
		{
			softBody->stats.numericalFailure = true;
			float speed = sqrtf( speedSq );
			p->velocity = b2MulSV( maxSubStepSpeed / speed, p->velocity );
		}
	}

	softBody->stats.maxStretch = stepMaxStretch;
}

// Update the cached bounds used by debug drawing and queries.
static void b2SoftUpdateBounds( b2SoftBody* softBody )
{
	b2Vec2 lower = { FLT_MAX, FLT_MAX };
	b2Vec2 upper = { -FLT_MAX, -FLT_MAX };
	float r = softBody->particleRadius;
	for ( int i = 0; i < softBody->particles.count; ++i )
	{
		b2Vec2 p = b2SubPos( softBody->particles.data[i].position, (b2Pos){ 0 } );
		lower.x = b2MinFloat( lower.x, p.x - r );
		lower.y = b2MinFloat( lower.y, p.y - r );
		upper.x = b2MaxFloat( upper.x, p.x + r );
		upper.y = b2MaxFloat( upper.y, p.y + r );
	}
	softBody->bounds.lowerBound = lower;
	softBody->bounds.upperBound = upper;
}

// Compute the diagnostic energies from the current state. Strain energy uses a per-constraint
// spring model 0.5 * k * (L - L0)^2 with k derived as 1/(m * compliance) when compliance is finite.
static void b2SoftComputeEnergies( b2World* world, b2SoftBody* softBody, float h )
{
	b2SoftBodyStats* stats = &softBody->stats;
	float kinetic = 0.0f;
	float potential = 0.0f;
	float lowestY = FLT_MAX;

	for ( int i = 0; i < softBody->particles.count; ++i )
	{
		b2SoftParticle* p = softBody->particles.data + i;
		b2Vec2 pos = b2SubPos( p->position, (b2Pos){ 0 } );
		lowestY = b2MinFloat( lowestY, pos.y );
		kinetic += 0.5f * softBody->particleMass * b2LengthSquared( p->velocity );
	}

	float strain = 0.0f;
	float h2 = h * h;
	for ( int i = 0; i < softBody->constraints.count; ++i )
	{
		b2SoftConstraint* c = softBody->constraints.data + i;
		if ( c->active == false )
		{
			continue;
		}
		b2SoftParticle* p1 = softBody->particles.data + c->particle1;
		b2SoftParticle* p2 = softBody->particles.data + c->particle2;
		float length = b2Length( b2Sub( b2SubPos( p2->position, (b2Pos){ 0 } ),
										b2SubPos( p1->position, (b2Pos){ 0 } ) ) );
		float extension = length - c->restLength;

		// PBD has no explicit spring constant; an equivalent energy uses the positional
		// correction stiffness k and the sub-step length: k_eff = k * m / h^2. This stays finite
		// for compliant constraints and vanishes as stiffness goes to zero.
		if ( c->stiffnessFactor > 0.0f )
		{
			float effectiveStiffness = c->stiffnessFactor * softBody->particleMass / h2;
			strain += 0.5f * effectiveStiffness * extension * extension;
		}
	}

	if ( softBody->enableGravity )
	{
		for ( int i = 0; i < softBody->particles.count; ++i )
		{
			b2SoftParticle* p = softBody->particles.data + i;
			b2Vec2 pos = b2SubPos( p->position, (b2Pos){ 0 } );
			potential += softBody->particleMass * ( -world->gravity.y ) * ( pos.y - lowestY );
		}
	}

	stats->kineticEnergy = kinetic;
	stats->strainEnergy = strain;
	stats->potentialEnergy = potential;
}

void b2StepSoftBodies( b2World* world, float dt, float invDt, float h, float invH, int subStepCount )
{
	B2_UNUSED( invDt );

	if ( world->softBodies.count == 0 )
	{
		return;
	}

	uint64_t ticks = b2GetTicks();

	// Clear per-step coupling accumulators. Only bodies marked in softImpulseBodyIds are touched.
	for ( int i = 0; i < world->softImpulseBodyIds.count; ++i )
	{
		world->softImpulseBodyIds.data[i] = 0;
		world->softLinearImpulses.data[i] = (b2Vec2){ 0.0f, 0.0f };
		world->softAngularImpulses.data[i] = 0.0f;
	}

	for ( int s = 0; s < world->softBodies.count; ++s )
	{
		b2SoftBody* softBody = world->softBodies.data + s;
		if ( softBody->inUse == false )
		{
			continue;
		}

		softBody->stats = (b2SoftBodyStats){ 0 };
		softBody->stats.particleCount = softBody->particles.count;

		for ( int sub = 0; sub < subStepCount; ++sub )
		{
			b2SoftStep( world, softBody, h, invH );
		}

		softBody->stats.constraintCount = 0;
		for ( int i = 0; i < softBody->constraints.count; ++i )
		{
			if ( softBody->constraints.data[i].active )
			{
				softBody->stats.constraintCount += 1;
			}
		}

		b2SoftUpdateBounds( softBody );
		b2SoftComputeEnergies( world, softBody, dt > 0.0f ? dt : h );
	}

	// Apply two-way coupling impulses to rigid bodies. The rigid solver already ran this step, so
	// the impulses take effect next step; this keeps the rigid body state itself deterministic and
	// avoids mutating solver arrays mid-solve.
	for ( int bodyId = 0; bodyId < world->softImpulseBodyIds.count; ++bodyId )
	{
		if ( world->softImpulseBodyIds.data[bodyId] == 0 )
		{
			continue;
		}
		b2Body* body = world->bodies.data + bodyId;
		if ( body->id == B2_NULL_INDEX || body->type != b2_dynamicBody ||
			 body->setIndex != b2_awakeSet )
		{
			continue;
		}

		b2SolverSet* set = b2Array_Get( world->solverSets, b2_awakeSet );
		b2BodyState* state = b2Array_Get( set->bodyStates, body->localIndex );
		b2BodySim* bodySim = b2Array_Get( set->bodySims, body->localIndex );

		b2Vec2 linearImpulse = world->softLinearImpulses.data[bodyId];
		state->linearVelocity = b2MulAdd( state->linearVelocity, bodySim->invMass, linearImpulse );
		state->angularVelocity += bodySim->invInertia * world->softAngularImpulses.data[bodyId];
		b2LimitVelocity( state, world->maxLinearSpeed );
	}

	world->profile.softBodies = b2GetMilliseconds( ticks );
}

// -------------------------------------------------------------------------------------------------
// Public API
// -------------------------------------------------------------------------------------------------

void b2DestroySoftBody( b2SoftBodyId softBodyId )
{
	if ( B2_IS_NULL( softBodyId ) )
	{
		return;
	}

	b2World* world = b2GetWorldLocked( softBodyId.world0 );
	if ( world == NULL )
	{
		return;
	}

	b2SoftBody* softBody = b2GetSoftBodyFullId( world, softBodyId );
	if ( softBody == NULL )
	{
		b2Log( "invalid input: stale soft body id in %s\n", __func__ );
		B2_ASSERT( false );
		return;
	}

	int id = softBody->id;
	b2DestroySoftBodyInternal( softBody );
	b2FreeId( &world->softBodyIdPool, id );
}

bool b2SoftBody_IsValid( b2SoftBodyId softBodyId )
{
	if ( B2_IS_NULL( softBodyId ) )
	{
		return false;
	}
	b2World* world = b2GetWorld( softBodyId.world0 );
	if ( world == NULL )
	{
		return false;
	}
	return b2GetSoftBodyFullId( world, softBodyId ) != NULL;
}

int b2SoftBody_GetParticleCount( b2SoftBodyId softBodyId )
{
	b2SoftBody* softBody = b2GetSoftBody( softBodyId );
	if ( softBody == NULL )
	{
		return 0;
	}
	return softBody->particles.count;
}

b2Vec2 b2SoftBody_GetParticlePosition( b2SoftBodyId softBodyId, int particleIndex )
{
	b2Vec2 zero = { 0.0f, 0.0f };
	b2SoftBody* softBody = b2GetSoftBody( softBodyId );
	if ( softBody == NULL || particleIndex < 0 || particleIndex >= softBody->particles.count )
	{
		b2Log( "invalid input: particle index out of range in %s\n", __func__ );
		return zero;
	}
	return b2SubPos( softBody->particles.data[particleIndex].position, (b2Pos){ 0 } );
}

b2Vec2 b2SoftBody_GetParticleVelocity( b2SoftBodyId softBodyId, int particleIndex )
{
	b2Vec2 zero = { 0.0f, 0.0f };
	b2SoftBody* softBody = b2GetSoftBody( softBodyId );
	if ( softBody == NULL || particleIndex < 0 || particleIndex >= softBody->particles.count )
	{
		b2Log( "invalid input: particle index out of range in %s\n", __func__ );
		return zero;
	}
	return softBody->particles.data[particleIndex].velocity;
}

void b2SoftBody_SetParticlePosition( b2SoftBodyId softBodyId, int particleIndex, b2Vec2 position )
{
	B2_CHECK_INPUT( b2IsValidVec2( position ) );
	b2SoftBody* softBody = b2GetSoftBody( softBodyId );
	if ( softBody == NULL || particleIndex < 0 || particleIndex >= softBody->particles.count )
	{
		b2Log( "invalid input: particle index out of range in %s\n", __func__ );
		return;
	}
	b2SoftParticle* p = softBody->particles.data + particleIndex;
	p->position = (b2Pos){ position.x, position.y };
	p->prevPosition = p->position;
}

void b2SoftBody_ApplyParticleImpulse( b2SoftBodyId softBodyId, int particleIndex, b2Vec2 impulse )
{
	B2_CHECK_INPUT( b2IsValidVec2( impulse ) );
	b2SoftBody* softBody = b2GetSoftBody( softBodyId );
	if ( softBody == NULL || particleIndex < 0 || particleIndex >= softBody->particles.count )
	{
		b2Log( "invalid input: particle index out of range in %s\n", __func__ );
		return;
	}
	b2SoftParticle* p = softBody->particles.data + particleIndex;
	if ( p->invMass <= 0.0f )
	{
		b2Log( "soft body warning: impulse on anchored (infinite mass) particle ignored in %s\n", __func__ );
		return;
	}
	p->accumulatedImpulse = b2Add( p->accumulatedImpulse, impulse );
}

int b2SoftBody_GetAnchorCount( b2SoftBodyId softBodyId )
{
	b2SoftBody* softBody = b2GetSoftBody( softBodyId );
	if ( softBody == NULL )
	{
		return 0;
	}
	return softBody->anchors.count;
}

void b2SoftBody_DestroyAnchor( b2SoftBodyId softBodyId, int anchorIndex )
{
	b2SoftBody* softBody = b2GetSoftBody( softBodyId );
	if ( softBody == NULL || anchorIndex < 0 || anchorIndex >= softBody->anchors.count )
	{
		b2Log( "invalid input: anchor index out of range in %s\n", __func__ );
		return;
	}

	int particleIndex = softBody->anchors.data[anchorIndex].particleIndex;
	b2Array_RemoveSwap( softBody->anchors, anchorIndex );

	// Restore finite mass if no remaining anchor pins this particle.
	for ( int i = 0; i < softBody->anchors.count; ++i )
	{
		if ( softBody->anchors.data[i].particleIndex == particleIndex )
		{
			return;
		}
	}
	softBody->particles.data[particleIndex].invMass = 1.0f / softBody->particleMass;
}

bool b2SoftBody_AddAnchor( b2SoftBodyId softBodyId, const b2SoftBodyAnchor* anchor )
{
	b2SoftBody* softBody = b2GetSoftBody( softBodyId );
	if ( softBody == NULL || anchor == NULL )
	{
		b2Log( "invalid input: null anchor in %s\n", __func__ );
		return false;
	}
	if ( anchor->particleIndex < 0 || anchor->particleIndex >= softBody->particles.count ||
		 b2IsValidVec2( anchor->localPoint ) == false )
	{
		b2Log( "invalid input: anchor references an invalid particle in %s\n", __func__ );
		return false;
	}

	b2World* world = b2GetWorld( softBodyId.world0 );
	b2SoftAnchor internal = { 0 };
	internal.particleIndex = anchor->particleIndex;
	internal.detached = false;

	if ( B2_IS_NULL( anchor->bodyId ) )
	{
		internal.bodyId = B2_NULL_INDEX;
		internal.localPoint = anchor->localPoint;
		internal.lastWorldPoint = (b2Pos){ anchor->localPoint.x, anchor->localPoint.y };
	}
	else
	{
		b2Body* body = b2GetBodyFullId( world, anchor->bodyId );
		if ( body == NULL || body->id == B2_NULL_INDEX )
		{
			b2Log( "invalid input: anchor references an invalid body in %s\n", __func__ );
			return false;
		}
		internal.bodyId = body->id;
		internal.localPoint = anchor->localPoint;
		b2WorldTransform transform = b2GetBodyTransformQuick( world, body );
		internal.lastWorldPoint = b2TransformWorldPoint( transform, anchor->localPoint );
	}

	b2Array_Push( softBody->anchors, internal );
	softBody->particles.data[anchor->particleIndex].invMass = 0.0f;
	return true;
}

void b2SoftBody_EnableSelfCollision( b2SoftBodyId softBodyId, bool enable )
{
	b2SoftBody* softBody = b2GetSoftBody( softBodyId );
	if ( softBody == NULL )
	{
		return;
	}
	softBody->enableSelfCollision = enable;
}

void b2SoftBody_EnableTearing( b2SoftBodyId softBodyId, bool enable, float maxStretch )
{
	b2SoftBody* softBody = b2GetSoftBody( softBodyId );
	if ( softBody == NULL )
	{
		return;
	}
	if ( enable && ( maxStretch <= 1.0f || maxStretch > B2_SOFT_MAX_SANE_STRETCH ||
					 b2IsValidFloat( maxStretch ) == false ) )
	{
		b2Log( "invalid input: tear threshold must be in (1, %g] in %s\n", B2_SOFT_MAX_SANE_STRETCH,
			   __func__ );
		return;
	}
	softBody->enableTearing = enable;
	softBody->tearMaxStretch = maxStretch;
}

void b2SoftBody_SetIterations( b2SoftBodyId softBodyId, int iterations )
{
	b2SoftBody* softBody = b2GetSoftBody( softBodyId );
	if ( softBody == NULL )
	{
		return;
	}
	if ( iterations < 1 || iterations > 256 )
	{
		b2Log( "invalid input: iterations must be in [1, 256] in %s\n", __func__ );
		return;
	}
	softBody->iterations = iterations;
}

void b2SoftBody_SetUserData( b2SoftBodyId softBodyId, void* userData )
{
	b2SoftBody* softBody = b2GetSoftBody( softBodyId );
	if ( softBody == NULL )
	{
		return;
	}
	softBody->userData = userData;
}

void* b2SoftBody_GetUserData( b2SoftBodyId softBodyId )
{
	b2SoftBody* softBody = b2GetSoftBody( softBodyId );
	if ( softBody == NULL )
	{
		return NULL;
	}
	return softBody->userData;
}

b2SoftBodyStats b2SoftBody_GetStats( b2SoftBodyId softBodyId )
{
	b2SoftBodyStats empty = { 0 };
	b2SoftBody* softBody = b2GetSoftBody( softBodyId );
	if ( softBody == NULL )
	{
		return empty;
	}
	return softBody->stats;
}

void b2DetachSoftBodyAnchors( b2World* world, int bodyId )
{
	for ( int s = 0; s < world->softBodies.count; ++s )
	{
		b2SoftBody* softBody = world->softBodies.data + s;
		if ( softBody->inUse == false )
		{
			continue;
		}
		for ( int i = 0; i < softBody->anchors.count; ++i )
		{
			b2SoftAnchor* anchor = softBody->anchors.data + i;
			if ( anchor->bodyId == bodyId )
			{
				anchor->detached = true;
				anchor->bodyId = B2_NULL_INDEX;
			}
		}
	}
}

// -------------------------------------------------------------------------------------------------
// Debug drawing
// -------------------------------------------------------------------------------------------------

void b2DrawSoftBodies( b2World* world, b2DebugDraw* draw )
{
	for ( int s = 0; s < world->softBodies.count; ++s )
	{
		b2SoftBody* softBody = world->softBodies.data + s;
		if ( softBody->inUse == false )
		{
			continue;
		}

		// Active constraints as lines. Distance edges are pale blue, bending constraints gold.
		for ( int i = 0; i < softBody->constraints.count; ++i )
		{
			const b2SoftConstraint* c = softBody->constraints.data + i;
			if ( c->active == false )
			{
				continue;
			}
			b2Pos p1 = softBody->particles.data[c->particle1].position;
			b2Pos p2 = softBody->particles.data[c->particle2].position;
			b2HexColor color =
				( c->type == b2_softDistanceConstraint ) ? b2_colorLightBlue : b2_colorGold;
			draw->DrawLineFcn( p1, p2, color, draw->context );
		}

		// Particles as points, anchored ones highlighted.
		for ( int i = 0; i < softBody->particles.count; ++i )
		{
			b2Pos p = softBody->particles.data[i].position;
			bool anchored = softBody->particles.data[i].invMass <= 0.0f;
			draw->DrawPointFcn( p, anchored ? 7.0f : 4.0f,
								anchored ? b2_colorRed : b2_colorSkyBlue, draw->context );
		}
	}
}
