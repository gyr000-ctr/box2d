// SPDX-FileCopyrightText: 2026 Erin Catto
// SPDX-License-Identifier: MIT

#include "draw.h"
#include "sample.h"
#include "utils.h"

#include "box2d/box2d.h"
#include "box2d/soft_body.h"

#include <GLFW/glfw3.h>

#include <stdio.h>
#include <vector>

// A soft cloth anchored at its two top corners and draped over a dynamic rigid platform.
// Demonstrates the optional soft body feature, self-collision, tearing, and two-way coupling
// between soft particles and rigid bodies.
class SoftCloth : public Sample
{
public:
	explicit SoftCloth( SampleContext* context )
		: Sample( context )
	{
		if ( m_context->restart == false )
		{
			m_context->camera.center = { 0.0f, 1.0f };
			m_context->camera.zoom = 14.0f;
		}

		m_softBodyId = b2_nullSoftBodyId;
		m_selfCollision = true;
		m_tearing = false;
		CreateScene();
	}

	~SoftCloth() override
	{
		if ( B2_IS_NON_NULL( m_softBodyId ) )
		{
			b2DestroySoftBody( m_softBodyId );
		}
	}

	static Sample* Create( SampleContext* context )
	{
		return new SoftCloth( context );
	}

	void CreateScene()
	{
		if ( B2_IS_NON_NULL( m_softBodyId ) )
		{
			b2DestroySoftBody( m_softBodyId );
			m_softBodyId = b2_nullSoftBodyId;
		}

		// Static floor.
		b2BodyDef floorDef = b2DefaultBodyDef();
		b2BodyId floorId = b2CreateBody( m_worldId, &floorDef );
		b2Polygon floorBox = b2MakeBox( 12.0f, 0.5f );
		b2ShapeDef shapeDef = b2DefaultShapeDef();
		b2CreatePolygonShape( floorId, &shapeDef, &floorBox );
		b2Body_SetTransform( floorId, (b2Pos){ 0.0f, -6.0f }, b2Rot_identity );

		// A dynamic platform the cloth lands on, free to be pushed and rotated.
		b2BodyDef platformDef = b2DefaultBodyDef();
		platformDef.type = b2_dynamicBody;
		b2BodyId platformId = b2CreateBody( m_worldId, &platformDef );
		b2Polygon platformBox = b2MakeBox( 1.5f, 0.2f );
		b2CreatePolygonShape( platformId, &shapeDef, &platformBox );
		b2Body_SetTransform( platformId, (b2Pos){ 0.0f, 1.0f }, b2Rot_identity );

		constexpr int cols = 14;
		constexpr int rows = 12;
		constexpr int vertexCount = ( cols + 1 ) * ( rows + 1 );
		constexpr int triangleCount = 2 * cols * rows;
		constexpr float spacing = 0.25f;

		std::vector<b2Vec2> vertices( vertexCount );
		std::vector<b2SoftBodyTriangle> triangles( triangleCount );

		for ( int y = 0; y <= rows; ++y )
		{
			for ( int x = 0; x <= cols; ++x )
			{
				int i = y * ( cols + 1 ) + x;
				vertices[i].x = ( x - 0.5f * cols ) * spacing;
				vertices[i].y = 5.0f - y * spacing;
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

		b2SoftBodyDef softDef = b2DefaultSoftBodyDef();
		softDef.vertices = vertices.data();
		softDef.vertexCount = vertexCount;
		softDef.triangles = triangles.data();
		softDef.triangleCount = triangleCount;
		softDef.enableSelfCollision = m_selfCollision;
		softDef.enableTearing = m_tearing;
		softDef.tearMaxStretch = 1.6f;
		softDef.iterations = 8;

		// Anchor the two top corners to fixed world points.
		b2SoftBodyAnchor anchors[2] = {
			{ 0, b2_nullBodyId, vertices[0] },
			{ cols, b2_nullBodyId, vertices[cols] },
		};
		softDef.anchors = anchors;
		softDef.anchorCount = 2;

		m_softBodyId = b2CreateSoftBody( m_worldId, &softDef );
	}

	void Step() override
	{
		Sample::Step();

		if ( B2_IS_NULL( m_softBodyId ) )
		{
			return;
		}

		int particleCount = b2SoftBody_GetParticleCount( m_softBodyId );
		for ( int i = 0; i < particleCount; ++i )
		{
			b2Vec2 p = b2SoftBody_GetParticlePosition( m_softBodyId, i );
			DrawPoint( m_draw, (b2Pos){ p.x, p.y }, 4.0f, b2_colorSkyBlue );
		}

		b2SoftBodyStats stats = b2SoftBody_GetStats( m_softBodyId );
		DrawScreenTextLine(
			"particles=%d constraints=%d torn=%d rigidContacts=%d selfContacts=%d fail=%d",
			stats.particleCount, stats.constraintCount, stats.tornCount, stats.rigidContactCount,
			stats.selfContactCount, stats.numericalFailure );
	}

	bool DrawControls() override
	{
		ImGui::SeparatorText( "Soft Cloth" );
		ImGui::Checkbox( "Self Collision", &m_selfCollision );
		ImGui::Checkbox( "Tearing", &m_tearing );
		if ( ImGui::Button( "Reset" ) )
		{
			CreateScene();
		}
		return true;
	}

private:
	b2SoftBodyId m_softBodyId;
	bool m_selfCollision;
	bool m_tearing;
};

static int softClothSample = RegisterSample( "Soft Body", "Cloth Coupling", SoftCloth::Create );
