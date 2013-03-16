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
		show_voronoi_cell_ = false;
		show_triangulation_ = false;
		highlighted_group = -1;
	}

	SPM_Graphics::~SPM_Graphics()
	{
		glDeleteLists(triangulation_displist, 1);
		glDeleteLists(vc_displist, 1);
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
		if (show_triangulation_)
			draw_triangulation();
		if (show_voronoi_cell_)
			draw_voronoi_cell();
		//draw_all_vertices();
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
		//glColor4fv(pgn_edge);
		glLineWidth(1.5f);
		for (unsigned int i = 0; i < tiles.size(); i++)
		{
			if (i != highlighted_group)
				glColor4fv(pgn_edge);
			else
				glColor4f(0.0f, 0.0f, 1.0f, 1.0f);
			glBegin(GL_LINE_LOOP);
			for (unsigned int j = 0; j < tiles[i].size(); j++)
				glPoint_3(tiles[i].vertex(j));
			glEnd();
		}
		glEnable(GL_LIGHTING);
	}

	void SPM_Graphics::draw_all_vertices()
	{
		const RDT_data_structure& rdt = packer->rpvd.get_rdt();
		
		glColor3f(0.0f, 1.0f, 1.0f);
		glPointSize(4.0f);
		glBegin(GL_POINTS);
		for (RDT_data_structure::Vertex_iterator vi = rdt.vertices_begin(); vi != rdt.vertices_end(); ++vi)
		{
			Point_3 p = vi->point_3();
			glPoint_3(p);
		}
		glEnd();
	}

	void SPM_Graphics::draw_triangulation()
	{
		if (!glIsList(triangulation_displist))
		{
			triangulation_displist = glGenLists(1);
			glNewList(triangulation_displist, GL_COMPILE);
			glLineWidth(1.0f);
			glColor3f(0.6f, 0.0f, 0.0f);
			glDisable(GL_LIGHTING);
			typedef RestrictedPolygonVoronoiDiagram RPVD;
			const RPVD& rpvd = packer->get_rpvd();
			glBegin(GL_LINES);
			for (RPVD::Edge_iterator eit = rpvd.edges_begin(); eit != rpvd.edges_end(); ++eit)
			{
				RPVD::Face_iterator f = eit->first;
				int i = eit->second;
				RPVD::Vertex_iterator vh0 = f->vertex(f->cw(i)), vh1 = f->vertex(f->ccw(i));
				glPoint_3(vh0->point_3());
				glPoint_3(vh1->point_3());
			}
			glEnd();
			glEnable(GL_LIGHTING);
			glEndList();
		}
		glCallList(triangulation_displist);
	}

	void SPM_Graphics::draw_voronoi_cell()
	{
		typedef RestrictedPolygonVoronoiDiagram RPVD;
		if (!glIsList(vc_displist))
		{
			vc_displist = glGenLists(1);
			glNewList(vc_displist);
			glDisable(GL_LIGHTING);
			const RPVD& rpvd = packer->get_rpvd();
			for (unsigned int i = 0; i < rpvd.number_of_groups(); i++)
			{
				const RPVD::VertGroup& vg = rpvd.sample_points_group(i);
				gl_table_color(i);
				glBegin(GL_TRIANGLES);
				glEnd();
			}
			glEnable(GL_LIGHTING);
			glEndList();
		}
	}
}