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
		show_mesh_ = true;
		show_polygons_ = true;
		show_voronoi_cell_ = true;
		show_triangulation_ = false;
		show_vertices_ = false;
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
		if (show_vertices_)
			draw_all_vertices();
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
		RDT_data_structure& rdt = packer->rpvd.get_rdt();
		
		glDisable(GL_LIGHTING);
		glColor3f(0.0f, 0.0f, 0.0f);
		glPointSize(2.0f);
		glBegin(GL_POINTS);
		for (RDT_data_structure::Vertex_iterator vi = rdt.vertices_begin(); vi != rdt.vertices_end(); ++vi)
		{
			//Point_3 p = vi->point_3();
			Point_3 p = vi->mp;
			glPoint_3(p);
		}
		glEnd();
		glEnable(GL_LIGHTING);
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
			RPVD& rpvd = packer->get_rpvd();
			glBegin(GL_LINES);
			for (RPVD::Halfedge_iterator eit = rpvd.edges_begin(); eit != rpvd.edges_end(); ++eit)
			{
				if (!eit->is_border()&& !eit->is_border() && rpvd.is_delaunay_edge(eit))
					glColor3f(0.6f, 0.0f, 0.0f);
				else
					glColor3f(0.0f, 0.0f, 0.4f);
				RDT_data_structure::Vertex_handle vh0 = eit->vertex(), vh1 = eit->opposite()->vertex();
				glPoint_3(vh0->mp);
				glPoint_3(vh1->mp);
				
			}
			glEnd();
			glEnable(GL_LIGHTING);
			glEndList();
		}
		glCallList(triangulation_displist);
	}

	void SPM_Graphics::draw_voronoi_cell()
	{
		//std::cout<<"Drawing Voronoi cells...\n";
		typedef RestrictedPolygonVoronoiDiagram RPVD;
		//if (!glIsList(vc_displist))
		//{
			//vc_displist = glGenLists(1);
			//glNewList(vc_displist, GL_COMPILE);
			glDisable(GL_LIGHTING);
			const RPVD& rpvd = packer->get_rpvd();
			const std::vector<Packing_object>& po = packer->get_tiles();
			for (unsigned int i = 0; i < rpvd.number_of_groups(); i++)
			{
				//if (i == highlighted_group) {
				const RPVD::VertGroup& vg = rpvd.sample_points_group(i);
				gl_table_color(i);
				//glBegin(GL_TRIANGLES);
				glLineWidth(1.5f);
				glBegin(GL_LINES);
				Point_3 c = po[i].centroid();
				for (unsigned int j = 0; j < vg.size(); j++)
				{
					const std::vector<Point_3>& verts = vg[j]->vd_vertices;
					unsigned int sz = verts.size();
					if (sz > 1)
						for (unsigned int k = 0; k < sz-1; k++)
						{
							//glPoint_3(c);
							glPoint_3(verts[k]);
							//glPoint_3(verts[(k+1)%sz]);
							glPoint_3(verts[(k+1)]);
						}
				}
				glEnd();
				//}
			}
			glEnable(GL_LIGHTING);
			//glEndList();
		//}
		//glCallList(vc_displist);
	}
}