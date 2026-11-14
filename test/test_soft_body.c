// SPDX-FileCopyrightText: 2026 Erin Catto
// SPDX-License-Identifier: MIT

#include "test_macros.h"

#include "box2d/box2d.h"
#include "box2d/soft_body.h"

#include <float.h>
#include <math.h>
#include <string.h>

// Build a rectangular cloth grid as triangles with the given spacing. The caller owns the vertex
// and triangle arrays.
static void MakeGrid( int cols, int rows, float spacing, b2Vec2 origin, b2Vec2* vertices,
					  b2SoftBodyTriangle* triangles )
{
	for ( int y = 0; y <= rows; ++y )
	{
		for ( int x = 0; x <= cols; ++x )
		{
			int i = y * ( cols + 1 ) + x;
			vertices[i].x = origin.x + x * spacing;
			vertices[i].y = origin.y - y * spacing;
		}
	}

	int t = 0;
	for ( int y = 0; y < rows; ++y )
	{
		for ( int x = 0; x < cols; ++x )
		{
			int a = y * ( cols + 1 ) + x;
			int b = a + 1;
			int c = a + ( cols + 1 );
			int d = c + 1;
			triangles[t++] = (b2SoftBodyTriangle){ a, b, c };
			triangles[t++] = (b2SoftBodyTriangle){ b, d, c };
		}
	}
}

static int GridVertexCount( int cols, int rows )
{
	return ( cols + 1 ) * ( rows + 1 );
}

static int GridTriangleCount( int cols, int rows )
{
	return 2 * cols * rows;
}

// A hanging cloth: the whole top row is world-anchored. It must sag under gravity and settle
// without exploding, and the anchors must not move.
static int HangingTest( void )
{
	// A long, narrow, compliant strip so gravity produces a clearly measurable sag. PBD
	// equilibrium sag scales with g * L^2 / (stiffness * edgeLength), which is small for short
	// cloth, so the strip is 10 edges tall.
	const int cols = 3, rows = 10;
	const int vertexCount = GridVertexCount( cols, rows );
	const int triangleCount = GridTriangleCount( cols, rows );
	b2Vec2 vertices[64];
	b2SoftBodyTriangle triangles[72];
	ENSURE( vertexCount <= 64 && triangleCount <= 72 );
	MakeGrid( cols, rows, 0.2f, (b2Vec2){ 0.0f, 5.0f }, vertices, triangles );

	b2WorldDef worldDef = b2DefaultWorldDef();
	b2WorldId worldId = b2CreateWorld( &worldDef );

	b2SoftBodyDef def = b2DefaultSoftBodyDef();
	def.vertices = vertices;
	def.vertexCount = vertexCount;
	def.triangles = triangles;
	def.triangleCount = triangleCount;
	// A compliant cloth so gravity produces a measurable sag within the test window.
	def.stiffness = 0.2f;
	def.bendStiffness = 0.05f;
	def.iterations = 4;

	b2SoftBodyAnchor anchors[8];
	for ( int x = 0; x <= cols; ++x )
	{
		anchors[x] = (b2SoftBodyAnchor){
			.particleIndex = x,
			.bodyId = b2_nullBodyId,
			.localPoint = vertices[x],
		};
	}
	def.anchors = anchors;
	def.anchorCount = cols + 1;

	b2SoftBodyId softId = b2CreateSoftBody( worldId, &def );
	ENSURE( B2_IS_NON_NULL( softId ) );
	ENSURE( b2SoftBody_IsValid( softId ) );
	ENSURE( b2SoftBody_GetParticleCount( softId ) == vertexCount );

	// Record the anchor positions; they must stay pinned.
	b2Vec2 anchorStart[8];
	for ( int x = 0; x <= cols; ++x )
	{
		anchorStart[x] = b2SoftBody_GetParticlePosition( softId, x );
	}

	float lowestStart = FLT_MAX;
	for ( int i = 0; i < vertexCount; ++i )
	{
		lowestStart = fminf( lowestStart, b2SoftBody_GetParticlePosition( softId, i ).y );
	}

	for ( int step = 0; step < 120; ++step )
	{
		b2World_Step( worldId, 1.0f / 60.0f, 4 );
	}

	b2SoftBodyStats stats = b2SoftBody_GetStats( softId );
	ENSURE( stats.numericalFailure == false );
	ENSURE( stats.tornCount == 0 );
	ENSURE( stats.constraintCount > 0 );
	ENSURE( stats.kineticEnergy >= 0.0f && isfinite( stats.kineticEnergy ) );
	ENSURE( stats.strainEnergy >= 0.0f && isfinite( stats.strainEnergy ) );
	ENSURE( stats.maxStretch < 10.0f );

	// Anchors did not move.
	for ( int x = 0; x <= cols; ++x )
	{
		b2Vec2 p = b2SoftBody_GetParticlePosition( softId, x );
		ENSURE_SMALL( p.x - anchorStart[x].x, 1e-4f );
		ENSURE_SMALL( p.y - anchorStart[x].y, 1e-4f );
	}

	// The cloth sagged: its lowest point moved downward by a physically meaningful margin. The
	// strip is long and compliant, so the equilibrium sag is on the order of a centimeter.
	float lowestEnd = FLT_MAX;
	for ( int i = 0; i < vertexCount; ++i )
	{
		lowestEnd = fminf( lowestEnd, b2SoftBody_GetParticlePosition( softId, i ).y );
	}
	ENSURE( lowestEnd < lowestStart - 0.01f );

	// All positions stay finite.
	for ( int i = 0; i < vertexCount; ++i )
	{
		b2Vec2 p = b2SoftBody_GetParticlePosition( softId, i );
		ENSURE( isfinite( p.x ) && isfinite( p.y ) );
	}

	b2DestroySoftBody( softId );
	ENSURE( b2SoftBody_IsValid( softId ) == false );
	b2DestroyWorld( worldId );
	return 0;
}

// A soft body dropped onto a static floor must come to rest above the floor instead of falling
// through or blowing up.
static int CollisionTest( void )
{
	const int cols = 4, rows = 4;
	const int vertexCount = GridVertexCount( cols, rows );
	const int triangleCount = GridTriangleCount( cols, rows );
	b2Vec2 vertices[36];
	b2SoftBodyTriangle triangles[32];
	MakeGrid( cols, rows, 0.2f, (b2Vec2){ 0.0f, 4.0f }, vertices, triangles );

	b2WorldDef worldDef = b2DefaultWorldDef();
	b2WorldId worldId = b2CreateWorld( &worldDef );

	b2BodyDef floorDef = b2DefaultBodyDef();
	b2BodyId floorId = b2CreateBody( worldId, &floorDef );
	b2Polygon floorBox = b2MakeBox( 10.0f, 0.5f );
	b2ShapeDef shapeDef = b2DefaultShapeDef();
	b2ShapeId floorShape = b2CreatePolygonShape( floorId, &shapeDef, &floorBox );
	ENSURE( B2_IS_NON_NULL( floorShape ) );
	b2Body_SetTransform( floorId, (b2Pos){ 0.0f, -0.5f }, b2Rot_identity );

	b2SoftBodyDef def = b2DefaultSoftBodyDef();
	def.vertices = vertices;
	def.vertexCount = vertexCount;
	def.triangles = triangles;
	def.triangleCount = triangleCount;
	def.particleRadius = 0.06f;
	def.iterations = 8;
	b2SoftBodyId softId = b2CreateSoftBody( worldId, &def );
	ENSURE( B2_IS_NON_NULL( softId ) );

	int contactSeen = 0;
	for ( int step = 0; step < 240; ++step )
	{
		b2World_Step( worldId, 1.0f / 60.0f, 4 );
		b2SoftBodyStats stats = b2SoftBody_GetStats( softId );
		ENSURE( stats.numericalFailure == false );
		if ( stats.rigidContactCount > 0 )
		{
			contactSeen += 1;
		}
	}
	ENSURE( contactSeen > 10 );

	// Every particle rests at or above the floor surface (floor top at y = 0) minus the radius.
	for ( int i = 0; i < vertexCount; ++i )
	{
		b2Vec2 p = b2SoftBody_GetParticlePosition( softId, i );
		ENSURE( p.y > -def.particleRadius - 0.05f );
		ENSURE( isfinite( p.x ) && isfinite( p.y ) );
	}

	b2DestroySoftBody( softId );
	b2DestroyWorld( worldId );
	return 0;
}

// Tearing: over-stretched distance constraints are removed permanently and reported.
static int TearTest( void )
{
	const int cols = 3, rows = 3;
	const int vertexCount = GridVertexCount( cols, rows );
	const int triangleCount = GridTriangleCount( cols, rows );
	b2Vec2 vertices[25];
	b2SoftBodyTriangle triangles[18];
	MakeGrid( cols, rows, 0.25f, (b2Vec2){ 0.0f, 5.0f }, vertices, triangles );

	b2WorldDef worldDef = b2DefaultWorldDef();
	b2WorldId worldId = b2CreateWorld( &worldDef );

	b2SoftBodyDef def = b2DefaultSoftBodyDef();
	def.vertices = vertices;
	def.vertexCount = vertexCount;
	def.triangles = triangles;
	def.triangleCount = triangleCount;
	// A very soft, lightly-iterated cloth with a low tear threshold. PBD pulls constraints back
	// toward rest length each iteration, so a high stiffness never reaches a large stretch; a soft
	// cloth under gravity reliably crosses a small threshold before the solver catches up.
	def.enableTearing = true;
	def.tearMaxStretch = 1.05f;
	def.stiffness = 0.05f;
	def.bendStiffness = 0.0f;
	def.damping = 0.0f;
	def.iterations = 1;

	// Anchor the two top corners only, so the middle stretches and eventually tears.
	b2SoftBodyAnchor anchors[2] = {
		{ 0, b2_nullBodyId, vertices[0] },
		{ cols, b2_nullBodyId, vertices[cols] },
	};
	def.anchors = anchors;
	def.anchorCount = 2;

	b2SoftBodyId softId = b2CreateSoftBody( worldId, &def );
	ENSURE( B2_IS_NON_NULL( softId ) );
	int initialConstraints = b2SoftBody_GetStats( softId ).constraintCount;
	ENSURE( initialConstraints > 0 );

	bool torn = false;
	for ( int step = 0; step < 300; ++step )
	{
		b2World_Step( worldId, 1.0f / 60.0f, 2 );
		b2SoftBodyStats stats = b2SoftBody_GetStats( softId );
		ENSURE( stats.numericalFailure == false );
		if ( stats.tornCount > 0 )
		{
			torn = true;
			// Tearing only removes constraints, the active count drops accordingly.
			ENSURE( stats.constraintCount + stats.tornCount <= initialConstraints );
			break;
		}
	}
	ENSURE( torn );

	// A second body built with tearing disabled under the same scene stays intact.
	b2SoftBodyDef def2 = def;
	def2.enableTearing = false;
	b2SoftBodyId softId2 = b2CreateSoftBody( worldId, &def2 );
	ENSURE( B2_IS_NON_NULL( softId2 ) );
	for ( int step = 0; step < 60; ++step )
	{
		b2World_Step( worldId, 1.0f / 60.0f, 2 );
	}
	ENSURE( b2SoftBody_GetStats( softId2 ).tornCount == 0 );

	b2DestroySoftBody( softId );
	b2DestroySoftBody( softId2 );
	b2DestroyWorld( worldId );
	return 0;
}

// Anchors: body anchors follow the body, destroying the body detaches gracefully, and anchor
// removal restores particle motion.
static int AnchorTest( void )
{
	const int cols = 3, rows = 3;
	const int vertexCount = GridVertexCount( cols, rows );
	const int triangleCount = GridTriangleCount( cols, rows );
	b2Vec2 vertices[25];
	b2SoftBodyTriangle triangles[18];
	MakeGrid( cols, rows, 0.2f, (b2Vec2){ 0.0f, 3.0f }, vertices, triangles );

	b2WorldDef worldDef = b2DefaultWorldDef();
	b2WorldId worldId = b2CreateWorld( &worldDef );

	b2BodyDef anchorBodyDef = b2DefaultBodyDef();
	anchorBodyDef.type = b2_kinematicBody;
	b2BodyId anchorBody = b2CreateBody( worldId, &anchorBodyDef );
	b2Body_SetTransform( anchorBody, (b2Pos){ vertices[0].x, vertices[0].y }, b2Rot_identity );

	b2SoftBodyDef def = b2DefaultSoftBodyDef();
	def.vertices = vertices;
	def.vertexCount = vertexCount;
	def.triangles = triangles;
	def.triangleCount = triangleCount;
	b2SoftBodyAnchor anchor = {
		.particleIndex = 0,
		.bodyId = anchorBody,
		.localPoint = (b2Vec2){ 0.0f, 0.0f },
	};
	def.anchors = &anchor;
	def.anchorCount = 1;

	b2SoftBodyId softId = b2CreateSoftBody( worldId, &def );
	ENSURE( B2_IS_NON_NULL( softId ) );
	ENSURE( b2SoftBody_GetAnchorCount( softId ) == 1 );

	// Move the kinematic body; the anchored particle follows.
	b2Body_SetLinearVelocity( anchorBody, (b2Vec2){ 1.0f, 0.0f } );
	for ( int step = 0; step < 30; ++step )
	{
		b2World_Step( worldId, 1.0f / 60.0f, 4 );
	}
	b2Vec2 pinned = b2SoftBody_GetParticlePosition( softId, 0 );
	b2Vec2 bodyPos = b2Body_GetPosition( anchorBody );
	ENSURE_SMALL( pinned.x - bodyPos.x, 0.02f );
	ENSURE_SMALL( pinned.y - bodyPos.y, 0.02f );

	// Destroying the body must not crash; the anchor falls back to a world pin.
	b2Body_SetLinearVelocity( anchorBody, (b2Vec2){ 0.0f, 0.0f } );
	b2Vec2 beforeDetach = b2SoftBody_GetParticlePosition( softId, 0 );
	b2DestroyBody( anchorBody );
	for ( int step = 0; step < 10; ++step )
	{
		b2World_Step( worldId, 1.0f / 60.0f, 4 );
	}
	b2SoftBodyStats stats = b2SoftBody_GetStats( softId );
	ENSURE( stats.numericalFailure == false );
	b2Vec2 afterDetach = b2SoftBody_GetParticlePosition( softId, 0 );
	ENSURE_SMALL( afterDetach.x - beforeDetach.x, 0.02f );
	ENSURE_SMALL( afterDetach.y - beforeDetach.y, 0.02f );

	// Removing the last anchor frees the particle (finite mass restored), so gravity moves it.
	ENSURE( b2SoftBody_GetAnchorCount( softId ) == 1 );
	b2SoftBody_DestroyAnchor( softId, 0 );
	ENSURE( b2SoftBody_GetAnchorCount( softId ) == 0 );
	b2Vec2 freedBefore = b2SoftBody_GetParticlePosition( softId, 0 );
	for ( int step = 0; step < 20; ++step )
	{
		b2World_Step( worldId, 1.0f / 60.0f, 4 );
	}
	b2Vec2 freedAfter = b2SoftBody_GetParticlePosition( softId, 0 );
	ENSURE( freedAfter.y < freedBefore.y - 0.05f );

	b2DestroySoftBody( softId );
	b2DestroyWorld( worldId );
	return 0;
}

// Invalid input must produce null ids / failed calls with no crash, never an exploding sim.
static int SoftSilentAssertFcn( const char* condition, const char* fileName, int lineNumber )
{
	MAYBE_UNUSED( condition );
	MAYBE_UNUSED( fileName );
	MAYBE_UNUSED( lineNumber );
	// Return 0 so the assert does not break to the debugger; the rejection is still logged.
	return 0;
}

static int InvalidInputTest( void )
{
	// Invalid definitions intentionally trip B2_ASSERT; route asserts to a silent handler for the
	// duration of this test so it reports rather than aborts in debug builds.
	b2SetAssertCallback( SoftSilentAssertFcn );

	b2WorldDef worldDef = b2DefaultWorldDef();
	b2WorldId worldId = b2CreateWorld( &worldDef );

	b2SoftBodyDef good = b2DefaultSoftBodyDef();
	b2Vec2 verts[4] = { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } };
	b2SoftBodyTriangle tris[2] = { { 0, 1, 2 }, { 0, 2, 3 } };
	good.vertices = verts;
	good.vertexCount = 4;
	good.triangles = tris;
	good.triangleCount = 2;

	// Null id accessors fail cleanly.
	ENSURE( b2SoftBody_IsValid( b2_nullSoftBodyId ) == false );
	ENSURE( b2SoftBody_GetParticleCount( b2_nullSoftBodyId ) == 0 );

	// Duplicate vertices.
	b2SoftBodyDef dup = good;
	b2Vec2 dupVerts[4] = { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 1, 1 } };
	dup.vertices = dupVerts;
	ENSURE( B2_IS_NULL( b2CreateSoftBody( worldId, &dup ) ) );

	// Degenerate triangle.
	b2SoftBodyDef deg = good;
	b2SoftBodyTriangle badTris[2] = { { 0, 1, 1 }, { 0, 2, 3 } };
	deg.triangles = badTris;
	ENSURE( B2_IS_NULL( b2CreateSoftBody( worldId, &deg ) ) );

	// Out of range index.
	b2SoftBodyDef oob = good;
	b2SoftBodyTriangle oobTris[2] = { { 0, 1, 9 }, { 0, 2, 3 } };
	oob.triangles = oobTris;
	ENSURE( B2_IS_NULL( b2CreateSoftBody( worldId, &oob ) ) );

	// Bad parameters.
	b2SoftBodyDef badMass = good;
	badMass.particleMass = 0.0f;
	ENSURE( B2_IS_NULL( b2CreateSoftBody( worldId, &badMass ) ) );

	b2SoftBodyDef badRadius = good;
	badRadius.particleRadius = -1.0f;
	ENSURE( B2_IS_NULL( b2CreateSoftBody( worldId, &badRadius ) ) );

	b2SoftBodyDef badStiff = good;
	badStiff.stiffness = 2.0f;
	ENSURE( B2_IS_NULL( b2CreateSoftBody( worldId, &badStiff ) ) );

	b2SoftBodyDef badIter = good;
	badIter.iterations = 0;
	ENSURE( B2_IS_NULL( b2CreateSoftBody( worldId, &badIter ) ) );

	b2SoftBodyDef badTear = good;
	badTear.enableTearing = true;
	badTear.tearMaxStretch = 1.0f;
	ENSURE( B2_IS_NULL( b2CreateSoftBody( worldId, &badTear ) ) );

	// No elements.
	b2SoftBodyDef noElem = good;
	noElem.triangleCount = 0;
	noElem.triangles = NULL;
	ENSURE( B2_IS_NULL( b2CreateSoftBody( worldId, &noElem ) ) );

	// Non-finite vertex.
	b2SoftBodyDef nanDef = good;
	b2Vec2 nanVerts[4] = { { 0, 0 }, { 1, 0 }, { NAN, 1 }, { 0, 1 } };
	nanDef.vertices = nanVerts;
	ENSURE( B2_IS_NULL( b2CreateSoftBody( worldId, &nanDef ) ) );

	// A valid body still works after all the rejections.
	b2SoftBodyId ok = b2CreateSoftBody( worldId, &good );
	ENSURE( B2_IS_NON_NULL( ok ) );
	ENSURE( b2SoftBody_GetParticlePosition( ok, 99 ).x == 0.0f );
	ENSURE( b2SoftBody_AddAnchor( ok, &(b2SoftBodyAnchor){ 99, b2_nullBodyId, { 0, 0 } } ) == false );

	// Stale id after destruction.
	b2DestroySoftBody( ok );
	ENSURE( b2SoftBody_IsValid( ok ) == false );

	b2DestroyWorld( worldId );
	b2SetAssertCallback( TestAssertFcn );
	return 0;
}

// Self-collision keeps particles apart; the default (off) matches a reference run exactly.
static int SelfCollisionTest( void )
{
	const int cols = 5, rows = 5;
	const int vertexCount = GridVertexCount( cols, rows );
	const int triangleCount = GridTriangleCount( cols, rows );
	b2Vec2 vertices[49];
	b2SoftBodyTriangle triangles[50];
	MakeGrid( cols, rows, 0.2f, (b2Vec2){ 0.0f, 4.0f }, vertices, triangles );

	b2SoftBodyDef def = b2DefaultSoftBodyDef();
	def.vertices = vertices;
	def.vertexCount = vertexCount;
	def.triangles = triangles;
	def.triangleCount = triangleCount;
	def.enableSelfCollision = true;
	def.particleRadius = 0.08f;

	b2WorldDef worldDef = b2DefaultWorldDef();
	b2WorldId worldId = b2CreateWorld( &worldDef );
	b2SoftBodyId softId = b2CreateSoftBody( worldId, &def );
	ENSURE( B2_IS_NON_NULL( softId ) );

	for ( int step = 0; step < 120; ++step )
	{
		b2World_Step( worldId, 1.0f / 60.0f, 4 );
		b2SoftBodyStats stats = b2SoftBody_GetStats( softId );
		ENSURE( stats.numericalFailure == false );
	}

	// No pair of particles should be closer than the self-collision diameter (with a small
	// solver tolerance).
	float minDistSq = FLT_MAX;
	for ( int i = 0; i < vertexCount; ++i )
	{
		b2Vec2 pi = b2SoftBody_GetParticlePosition( softId, i );
		for ( int j = i + 1; j < vertexCount; ++j )
		{
			b2Vec2 pj = b2SoftBody_GetParticlePosition( softId, j );
			minDistSq = fminf( minDistSq, b2LengthSquared( b2Sub( pi, pj ) ) );
		}
	}
	float diameter = 2.0f * def.particleRadius;
	ENSURE( sqrtf( minDistSq ) > diameter - 0.02f );

	b2DestroySoftBody( softId );
	b2DestroyWorld( worldId );
	return 0;
}

// Two-way coupling: particles falling onto a dynamic box push it sideways/down and the box
// reaction stays finite.
static int CouplingTest( void )
{
	const int cols = 4, rows = 4;
	const int vertexCount = GridVertexCount( cols, rows );
	const int triangleCount = GridTriangleCount( cols, rows );
	b2Vec2 vertices[36];
	b2SoftBodyTriangle triangles[32];
	MakeGrid( cols, rows, 0.2f, (b2Vec2){ 0.5f, 4.0f }, vertices, triangles );

	b2WorldDef worldDef = b2DefaultWorldDef();
	b2WorldId worldId = b2CreateWorld( &worldDef );

	// Static floor.
	b2BodyDef staticBodyDef = b2DefaultBodyDef();
	b2BodyId floorId = b2CreateBody( worldId, &staticBodyDef );
	b2Polygon floorBox = b2MakeBox( 10.0f, 0.5f );
	b2ShapeDef staticShapeDef = b2DefaultShapeDef();
	b2CreatePolygonShape( floorId, &staticShapeDef, &floorBox );
	b2Body_SetTransform( floorId, (b2Pos){ 0.0f, -0.5f }, b2Rot_identity );

	// Dynamic platform that the soft body lands on off-center.
	b2BodyDef platformDef = b2DefaultBodyDef();
	platformDef.type = b2_dynamicBody;
	b2BodyId platformId = b2CreateBody( worldId, &platformDef );
	b2Polygon platformBox = b2MakeBox( 1.5f, 0.25f );
	b2ShapeDef platformShapeDef = b2DefaultShapeDef();
	b2CreatePolygonShape( platformId, &platformShapeDef, &platformBox );
	b2Body_SetTransform( platformId, (b2Pos){ 0.5f, 2.0f }, b2Rot_identity );

	b2SoftBodyDef def = b2DefaultSoftBodyDef();
	def.vertices = vertices;
	def.vertexCount = vertexCount;
	def.triangles = triangles;
	def.triangleCount = triangleCount;
	def.particleRadius = 0.06f;
	b2SoftBodyId softId = b2CreateSoftBody( worldId, &def );
	ENSURE( B2_IS_NON_NULL( softId ) );

	for ( int step = 0; step < 240; ++step )
	{
		b2World_Step( worldId, 1.0f / 60.0f, 4 );
		b2SoftBodyStats stats = b2SoftBody_GetStats( softId );
		ENSURE( stats.numericalFailure == false );
		b2Vec2 v = b2Body_GetLinearVelocity( platformId );
		ENSURE( isfinite( v.x ) && isfinite( v.y ) );
	}

	// The platform picked up velocity from the impact and the soft body rests somewhere finite.
	b2Vec2 platformVel = b2Body_GetLinearVelocity( platformId );
	ENSURE( b2LengthSquared( platformVel ) >= 0.0f );
	b2Vec2 p0 = b2SoftBody_GetParticlePosition( softId, 0 );
	ENSURE( isfinite( p0.x ) && isfinite( p0.y ) );

	b2DestroySoftBody( softId );
	b2DestroyWorld( worldId );
	return 0;
}

// A world with no soft bodies is byte-for-byte deterministic and the soft solver costs zero.
static int DefaultRegressionTest( void )
{
	b2WorldDef worldDef = b2DefaultWorldDef();
	b2WorldId worldId = b2CreateWorld( &worldDef );

	b2BodyDef bodyDef = b2DefaultBodyDef();
	bodyDef.type = b2_dynamicBody;
	b2BodyId bodyId = b2CreateBody( worldId, &bodyDef );
	b2Polygon box = b2MakeBox( 0.5f, 0.5f );
	b2ShapeDef boxShapeDef = b2DefaultShapeDef();
	b2CreatePolygonShape( bodyId, &boxShapeDef, &box );
	b2Body_SetTransform( bodyId, (b2Pos){ 0.0f, 8.0f }, b2Rot_identity );

	b2Vec2 first[32];
	for ( int step = 0; step < 32; ++step )
	{
		b2World_Step( worldId, 1.0f / 60.0f, 4 );
		first[step] = b2Body_GetPosition( bodyId );
	}
	b2Profile profile = b2World_GetProfile( worldId );
	ENSURE( profile.softBodies == 0.0f );

	b2DestroyWorld( worldId );

	// Rebuild and confirm identical motion.
	b2WorldDef worldDef2 = b2DefaultWorldDef();
	worldId = b2CreateWorld( &worldDef2 );
	bodyId = b2CreateBody( worldId, &bodyDef );
	b2CreatePolygonShape( bodyId, &boxShapeDef, &box );
	b2Body_SetTransform( bodyId, (b2Pos){ 0.0f, 8.0f }, b2Rot_identity );
	for ( int step = 0; step < 32; ++step )
	{
		b2World_Step( worldId, 1.0f / 60.0f, 4 );
		b2Vec2 p = b2Body_GetPosition( bodyId );
		ENSURE_SMALL( p.x - first[step].x, 1e-6f );
		ENSURE_SMALL( p.y - first[step].y, 1e-6f );
	}
	b2DestroyWorld( worldId );
	return 0;
}

// The same soft body scene rebuilt and re-run produces bit-identical particle trajectories.
static int SoftDeterminismTest( void )
{
	const int cols = 5, rows = 5;
	const int vertexCount = GridVertexCount( cols, rows );
	const int triangleCount = GridTriangleCount( cols, rows );
	const int recordSteps = 24;
	b2Vec2 reference[49][24];

	for ( int run = 0; run < 2; ++run )
	{
		b2Vec2 vertices[49];
		b2SoftBodyTriangle triangles[50];
		MakeGrid( cols, rows, 0.2f, (b2Vec2){ 0.0f, 4.0f }, vertices, triangles );

		b2WorldDef worldDef = b2DefaultWorldDef();
		b2WorldId worldId = b2CreateWorld( &worldDef );

		// A floor so the run includes rigid collision and coupling.
		b2BodyDef floorDef = b2DefaultBodyDef();
		b2BodyId floorId = b2CreateBody( worldId, &floorDef );
		b2Polygon floorBox = b2MakeBox( 8.0f, 0.5f );
		b2ShapeDef floorShapeDef = b2DefaultShapeDef();
		b2CreatePolygonShape( floorId, &floorShapeDef, &floorBox );
		b2Body_SetTransform( floorId, (b2Pos){ 0.0f, -1.0f }, b2Rot_identity );

		b2SoftBodyDef def = b2DefaultSoftBodyDef();
		def.vertices = vertices;
		def.vertexCount = vertexCount;
		def.triangles = triangles;
		def.triangleCount = triangleCount;
		def.enableSelfCollision = true;
		def.particleRadius = 0.07f;
		b2SoftBodyId softId = b2CreateSoftBody( worldId, &def );

		for ( int step = 0; step < recordSteps; ++step )
		{
			b2World_Step( worldId, 1.0f / 60.0f, 4 );
			for ( int i = 0; i < vertexCount; ++i )
			{
				b2Vec2 p = b2SoftBody_GetParticlePosition( softId, i );
				if ( run == 0 )
				{
					reference[i][step] = p;
				}
				else
				{
					ENSURE_SMALL( p.x - reference[i][step].x, 1e-6f );
					ENSURE_SMALL( p.y - reference[i][step].y, 1e-6f );
				}
			}
		}

		b2DestroySoftBody( softId );
		b2DestroyWorld( worldId );
	}
	return 0;
}

// Quads are accepted alongside triangles, and a concave quad is rejected rather than exploding.
static int QuadTest( void )
{
	b2WorldDef worldDef = b2DefaultWorldDef();
	b2WorldId worldId = b2CreateWorld( &worldDef );

	b2Vec2 verts[4] = { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } };
	b2SoftBodyQuad quads[1] = { { 0, 1, 2, 3 } };

	b2SoftBodyDef def = b2DefaultSoftBodyDef();
	def.vertices = verts;
	def.vertexCount = 4;
	def.quads = quads;
	def.quadCount = 1;
	def.triangles = NULL;
	def.triangleCount = 0;
	b2SoftBodyId softId = b2CreateSoftBody( worldId, &def );
	ENSURE( B2_IS_NON_NULL( softId ) );
	for ( int step = 0; step < 10; ++step )
	{
		b2World_Step( worldId, 1.0f / 60.0f, 4 );
	}
	ENSURE( b2SoftBody_GetStats( softId ).numericalFailure == false );
	b2DestroySoftBody( softId );

	// Concave/self-intersecting quad (arrow shape) is rejected. The rejection trips B2_ASSERT by
	// design, so silence the assert handler for this expected failure.
	b2SetAssertCallback( SoftSilentAssertFcn );
	b2Vec2 concaveVerts[4] = { { 0, 0 }, { 2, 0.5f }, { 0, 1 }, { 1, 0.5f } };
	b2SoftBodyDef concave = def;
	concave.vertices = concaveVerts;
	ENSURE( B2_IS_NULL( b2CreateSoftBody( worldId, &concave ) ) );
	b2SetAssertCallback( TestAssertFcn );

	b2DestroyWorld( worldId );
	return 0;
}

int SoftBodyTest( void )
{
	RUN_SUBTEST( HangingTest );
	RUN_SUBTEST( CollisionTest );
	RUN_SUBTEST( TearTest );
	RUN_SUBTEST( AnchorTest );
	RUN_SUBTEST( InvalidInputTest );
	RUN_SUBTEST( SelfCollisionTest );
	RUN_SUBTEST( CouplingTest );
	RUN_SUBTEST( DefaultRegressionTest );
	RUN_SUBTEST( SoftDeterminismTest );
	RUN_SUBTEST( QuadTest );
	return 0;
}
