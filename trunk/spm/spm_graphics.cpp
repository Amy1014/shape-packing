#include "spm_graphics.h"

namespace Geex
{
	const GLfloat SPM_Graphics::c1 = 0.35f;
	const GLfloat SPM_Graphics::c2 = 0.5f;
	const GLfloat SPM_Graphics::c3 = 1.0f;
	GLfloat SPM_Graphics::color_table[24][3] = 
	{
		{c3, c2, c2}, {c2, c3, c2},	{c2, c2, c3}, {c2, c3, c3},
		{c3, c2, c3}, {c3, c3, c2},	{c1, c2, c2}, {c2, c1, c2},
		{c2, c2, c1}, {c2, c1, c1},	{c1, c2, c1}, {c1, c1, c2},
		{c1, c3, c3}, {c3, c1, c3}, {c3, c3, c1}, {c3, c1, c1},
		{c1, c3, c1}, {c1, c1, c3},	{c1, c2, c3}, {c1, c3, c2}, 
		{c2, c1, c3}, {c2, c3, c1}, {c3, c1, c2}, {c3, c2, c1}
	};
	GLfloat SPM_Graphics::surf_diff[4] = {0.95f, 0.95f, 0.0f, 1.0f};
	GLfloat SPM_Graphics::surf_spec[4] = {0.8f, 0.8f, 0.8f, 1.0f} ;
	GLfloat SPM_Graphics::surf_edge[4] = {0.1f, 0.1f, 0.1f, 1.0f};
	//GLfloat SPM_Graphics::pgn_edge[4] = {0.1f, 0.1f, 0.1f, 1.0f};
	GLfloat SPM_Graphics::shininess[] = {75.0f};

	SPM_Graphics::SPM_Graphics(Packer *_packer): packer(_packer)
	{
		how_to_draw_polygons = OUTLINE_DRAW;
	}

	inline void SPM_Graphics::gl_table_color(int index)
	{
		int i = index % (sizeof(color_table)/sizeof(color_table[0]));
		glColor3f(color_table[i][0], color_table[i][1], color_table[i][2]);
	}

	void SPM_Graphics::draw()
	{
		glUseProgramObjectARB(0);
		if (show_mesh_)
			draw_mesh();
		if (show_polygons_)
			draw_polygons();
		//packer->draw_RDT();
	}

	void SPM_Graphics::draw_mesh()
	{
		glCullFace(GL_FRONT) ;
		glEnable(GL_LIGHTING);
		glPushMatrix();
		glMaterialfv(GL_FRONT, GL_DIFFUSE, surf_diff);
		glMaterialfv(GL_FRONT, GL_SPECULAR, surf_spec);
		glMaterialfv(GL_FRONT, GL_SHININESS, shininess);
		const TriMesh& m = packer->mesh_domain();
		for(unsigned int i=0; i<m.size(); i++) 
		{
			vec3 n = m[i].normal();
			glBegin(GL_TRIANGLES);
			glNormal(n);
			glVertex(m[i].vertex[0]) ;
			glVertex(m[i].vertex[1]) ;
			glVertex(m[i].vertex[2]) ;
			glEnd();
		}
		glPopMatrix();
		glDisable(GL_LIGHTING);
	}

	void SPM_Graphics::draw_polygons()
	{
		switch (how_to_draw_polygons)
		{
		case OUTLINE_DRAW:
		case FILL_DRAW:
		case TEXTURE_DRAW:
			outline_draw_polygons();
		}
	}

	void SPM_Graphics::outline_draw_polygons()
	{
		static GLfloat pgn_edge[4] = {0.1f, 0.1f, 0.1f, 1.0f};
		const vector<Packing_object>& tiles = packer->get_tiles();
		glDisable(GL_LIGHTING);
		glColor4fv(pgn_edge);
		glLineWidth(1.5f);
		for (unsigned int i = 0; i < tiles.size(); i++)
		{
			glBegin(GL_LINE_LOOP);
			for (unsigned int j = 0; j < tiles[i].size(); j++)
				glPoint_3(tiles[i].vertex(j));
			glEnd();
		}
		glEnable(GL_LIGHTING);
	}
}