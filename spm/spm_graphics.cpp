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
	GLfloat SPM_Graphics::surf_diff[4] = {0.543f, 0.543f, 1.0f, 1.0f};
	GLfloat SPM_Graphics::surf_spec[4] = {0.9f, 0.9f, 0.9f, 1.0f} ;
	GLfloat SPM_Graphics::surf_edge[4] = {0.1f, 0.1f, 0.1f, 1.0f};
	//GLfloat SPM_Graphics::pgn_edge[4] = {0.1f, 0.1f, 0.1f, 1.0f};
	GLfloat SPM_Graphics::shininess[] = {150.0f};

	SPM_Graphics::SPM_Graphics(Packer *_packer): packer(_packer)
	{
		how_to_draw_polygons = OUTLINE_DRAW;
		show_mesh_ = true;
		show_polygons_ = true;
		show_voronoi_cell_ = true;
		show_triangulation_ = false;
		show_vertices_ = false;
		show_hole_triangles_ = false;
		show_local_frame_ = false;
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
		if (show_hole_triangles_)
			draw_hole_triangles();
		if (show_holes_)
			draw_holes();
		if (show_local_frame_)
			draw_local_frames();
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
			outline_draw_polygons();
			break;
		case FILL_DRAW:
			fill_draw_polygons();
			break;
		case TEXTURE_DRAW:
			outline_draw_polygons();
			break;
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

	void SPM_Graphics::fill_draw_polygons()
	{
		static GLfloat amb_mat[] = {0.19225f, 0.19225f, 0.19225f, 1.0f};
		static GLfloat diff_mat[] = {0.50754f, 0.50754f, 0.50754f, 1.0f};
		static GLfloat spec_mat[] = {0.608273f, 0.608273f, 0.608273f, 1.0f};
		static GLfloat pgn_shininess = 0.4;

		const vector<Packing_object>& tiles = packer->get_tiles();

		glCullFace(GL_FRONT);
		glEnable(GL_LIGHTING);
		glPushMatrix();
		//glMaterialfv(GL_FRONT, GL_AMBIENT, amb_mat);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, diff_mat);
		glMaterialfv(GL_FRONT, GL_SPECULAR, spec_mat);
		glMaterialf(GL_FRONT, GL_SHININESS, 50.0f);
		for (unsigned int i = 0; i < tiles.size(); i++)
		{
			glBegin(GL_TRIANGLES);
			Point_3 c = tiles[i].centroid();
			Vector_3 n = tiles[i].norm();
			glNormal3d(n.x(), n.y(), n.z());
			for (unsigned int j = 0; j < tiles[i].size(); j++)
			{
				glPoint_3(c);
				glPoint_3(tiles[i].vertex(j));
				glPoint_3(tiles[i].vertex((j+1)%tiles[i].size()));
			}
			glEnd();
		}
		glPopMatrix();
		glDisable(GL_LIGHTING);
	}
	void SPM_Graphics::draw_all_vertices()
	{
		RDT_data_structure& rdt = packer->rpvd.get_rdt();
		
		glDisable(GL_LIGHTING);
		glColor3f(0.0f, 0.0f, 0.0f);
		glPointSize(2.0f);
		glBegin(GL_POINTS);
		for (Vertex_iterator vi = rdt.vertices_begin(); vi != rdt.vertices_end(); ++vi)
		{
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
			RPVD& rpvd = packer->get_rpvd();
			glBegin(GL_LINES);
			for (Halfedge_iterator eit = rpvd.edges_begin(); eit != rpvd.edges_end(); ++eit)
			{
				if (!eit->is_border()&& !eit->opposite()->is_border() && rpvd.is_delaunay_edge(eit))
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
		if (!glIsList(vc_displist))
		{
			vc_displist = glGenLists(1);
			glNewList(vc_displist, GL_COMPILE);
			glDisable(GL_LIGHTING);
			const RPVD& rpvd = packer->get_rpvd();
			const std::vector<Packing_object>& po = packer->get_tiles();
			for (unsigned int i = 0; i < rpvd.number_of_groups(); i++)
			{
				if (i == highlighted_group || highlighted_group == unsigned int(-1)) {
				const RPVD::VertGroup& vg = rpvd.sample_points_group(i);
				glColor3f(0.0f, 0.9f, 0.0f);
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
				}
			}
			glEnable(GL_LIGHTING);
			glEndList();
		}
		glCallList(vc_displist);
	}

	void SPM_Graphics::draw_hole_triangles()
	{
		RPVD& rpvd = packer->get_rpvd();
		glColor3f(0.3f, 0.3f, 0.3f);
		glDisable(GL_LIGHTING);
		for (Facet_iterator fit = rpvd.faces_begin(); fit != rpvd.faces_end(); ++fit)
		{
			if (fit->vacant)
			{
				RPVD::Halfedge_handle eh = fit->halfedge();
				Vertex_handle v0 = eh->vertex();
				eh = eh->next();
				Vertex_handle v1 = eh->vertex();
				eh = eh->next();
				Vertex_handle v2 = eh->vertex();
				glBegin(GL_TRIANGLES);
				glPoint_3(v0->mp);
				glPoint_3(v1->mp);
				glPoint_3(v2->mp);
				glEnd();
			}
		}
		glEnable(GL_LIGHTING);
	}

	void SPM_Graphics::draw_holes()
	{
		std::vector<Packer::Hole>& holes = packer->get_holes();
		glDisable(GL_LIGHTING);
		//glColor3f(0.0f, 0.0f, 1.0f);
		
		for (unsigned int i = 0; i < holes.size(); i++)
		{
			Packer::Hole h = holes[i];
			gl_table_color(i);
			glBegin(GL_LINES);
			for (unsigned int j = 0; j < h.size(); j++)
			{
				Halfedge_handle e = h[j];
				Halfedge_handle oe = e->opposite();
				glPoint_3(e->vertex()->mp);
				glPoint_3(oe->vertex()->mp);
			}
			glEnd();
		}
		
		glEnable(GL_LIGHTING);
	}

	void SPM_Graphics::draw_local_frames()
	{
		const std::vector<Packing_object>& tiles = packer->get_tiles();
		glDisable(GL_LIGHTING);
		glColor3f(0.0f, 0.0f, 1.0f);
		for (unsigned int i = 0; i < tiles.size(); i++)
		{
			Point_3 o = tiles[i].centroid();
			Vector_3 u(o, tiles[i].vertex(0));
			u = u / CGAL::sqrt(u.squared_length());
			Vector_3 w = tiles[i].norm();
			w = w / CGAL::sqrt(w.squared_length());
			Vector_3 v = CGAL::cross_product(w, u);
			glBegin(GL_LINES);
			glPoint_3(o);
			glPoint_3(o + 0.1 * u);
			glPoint_3(o);
			glPoint_3(o + 0.1 * v);
			glPoint_3(o);
			glPoint_3(o + 0.1 * w);
			glEnd();
		}
		glEnable(GL_LIGHTING);
	}
}