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
	GLfloat SPM_Graphics::surf_diff[4] = {0.3922f, 0.5843f, 0.9294f, 1.0f};
	GLfloat SPM_Graphics::surf_spec[4] = {0.9f, 0.9f, 0.9f, 1.0f} ;
	GLfloat SPM_Graphics::surf_edge[4] = {0.1f, 0.1f, 0.1f, 1.0f};
	//GLfloat SPM_Graphics::pgn_edge[4] = {0.1f, 0.1f, 0.1f, 1.0f};
	GLfloat SPM_Graphics::shininess[] = {150.0f};

	SPM_Graphics::SPM_Graphics(Packer *_packer): packer(_packer)
	{
		how_to_draw_polygons = OUTLINE_DRAW;
		show_mesh_ = true;
		show_polygons_ = true;
		show_voronoi_cell_ = false;
		show_triangulation_ = false;
		show_vertices_ = false;
		show_hole_triangles_ = false;
		show_local_frame_ = false;
		show_curvatures_ = false;
		textured = false;
		show_holes_ = false;
		show_cdt_ = false;
		highlighted_group = -1;
	}

	SPM_Graphics::~SPM_Graphics()
	{
		glDeleteLists(triangulation_displist, 1);
		glDeleteLists(vc_displist, 1);
		glDeleteLists(cur_color_displist, 1);
		for (unsigned int i = 0; i < texture_lib.size(); i++)
			glDeleteTextures(1, &texture_lib[i]);
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
		if (show_smoothed_voronoi_cell_)
			draw_smoothed_voronoi_cell();
		if (show_vertices_)
			draw_all_vertices();
		if (show_hole_triangles_)
			draw_hole_triangles();
		if (show_holes_)
			draw_holes();
		if (show_local_frame_)
			draw_local_frames();
		if (show_cdt_)
			draw_cdt();
		if (show_curvatures_)
			glCallList(cur_color_displist);

		if (show_multi_submeshes_)
			draw_multi_submeshes();
		if (show_multi_tiles_)
			draw_multi_tiles();
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

	void SPM_Graphics::draw_multi_submeshes()
	{
		show_mesh_ = false;

		glCullFace(GL_FRONT) ;
		glEnable(GL_LIGHTING);
		glPushMatrix();
		glMaterialfv(GL_FRONT, GL_DIFFUSE, surf_diff);
		glMaterialfv(GL_FRONT, GL_SPECULAR, surf_spec);
		glMaterialfv(GL_FRONT, GL_SHININESS, shininess);
		const std::vector<TriMesh>& submeshes = packer->get_multigroup_submeshes();
		for (unsigned int i = 0; i < submeshes.size(); i++)
		{
			const TriMesh& m = submeshes[i];
			glBegin(GL_TRIANGLES);
			for (unsigned int j = 0; j < m.size(); j++)
			{
				vec3 n = m[j].normal();			
				glNormal(n);
				glVertex(m[j].vertex[0]) ;
				glVertex(m[j].vertex[1]) ;
				glVertex(m[j].vertex[2]) ;
			}
			glEnd();
		}
		glPopMatrix();
		glDisable(GL_LIGHTING);

	}

	void SPM_Graphics::draw_multi_tiles()
	{
		const std::vector<std::vector<Packing_object>>& multi_tiles = packer->res_pack_objects;

		static GLfloat amb_mat[][4] = {{0.19225f, 0.19225f, 0.19225f, 1.0f}, {0.19225f, 0.19225f, 0.19225f, 1.0f}, {0.19225f, 0.19225f, 0.19225f, 1.0f}};
		static GLfloat diff_mat[][4] = {{0.2f, 0.2f, 0.2f, 1.0f}, {1.0f, 1.0f, 1.0f, 1.0f}, {0.50754f, 0.50754f, 0.50754f, 1.0f}};
		static GLfloat spec_mat[][4] = {{0.608273f, 0.608273f, 0.608273f, 1.0f}, {0.608273f, 0.608273f, 0.608273f, 1.0f}, {0.608273f, 0.608273f, 0.608273f, 1.0f}};
		static GLfloat pgn_shininess = 0.4;

		for (unsigned int k = 0; k < multi_tiles.size(); k++)
		{
			glCullFace(GL_FRONT);
			glEnable(GL_LIGHTING);
			const std::vector<Packing_object>& tiles = multi_tiles[k];
			glPushMatrix();
			//glMaterialfv(GL_FRONT, GL_AMBIENT, amb_mat);
			glMaterialfv(GL_FRONT, GL_DIFFUSE, diff_mat[k%3]);
			glMaterialfv(GL_FRONT, GL_SPECULAR, spec_mat[k%3]);
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
			if (textured)
				texture_draw_polygons();
			else
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

	void SPM_Graphics::texture_draw_polygons()
	{
		static GLfloat mat[] = {0.9f, 0.9f, 0.9f, 1.0f};
		const std::vector<Packing_object>& tiles = packer->get_tiles();
		glEnable(GL_TEXTURE_2D);
		GLboolean old_cull_face_config = glIsEnabled(GL_CULL_FACE);
		glEnable(GL_CULL_FACE);
		glCullFace(GL_BACK);
		//glEnable(GL_LIGHTING);
		glDisable(GL_LIGHTING);
		glPushMatrix();
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat);
		for (unsigned int i = 0; i < tiles.size(); i++)
		{
			Point_3 c = tiles[i].centroid();
			const std::vector<Point_2>& texture_coord = tiles[i].texture_coord;
			Point_2 tex_c = CGAL::centroid(texture_coord.begin(), texture_coord.end(), CGAL::Dimension_tag<0>());
			Vector_3 n = tiles[i].norm();
			glBindTexture(GL_TEXTURE_2D, texture_lib[tiles[i].texture_id]);
			glBegin(GL_TRIANGLES);
			unsigned int sz = tiles[i].size();
			for (unsigned int j = 0; j < sz; j++)
			{
				glNormal3d(n.x(), n.y(), n.z());
				glTexCoord2d(texture_coord[j].x(), texture_coord[j].y());
				glPoint_3(tiles[i].vertex(j));
				//glNormal3d(-n.x(), -n.y(), -n.z());
				glTexCoord2d(texture_coord[(j+1)%sz].x(), texture_coord[(j+1)%sz].y());
				glPoint_3(tiles[i].vertex((j+1)%sz));
				//glNormal3d(-n.x(), -n.y(), -n.z());
				glTexCoord2d(tex_c.x(), tex_c.y());
				glPoint_3(c);
			}
			glEnd();
		}
		glPopMatrix();
		//glDisable(GL_LIGHTING);
		if (!old_cull_face_config)
			glDisable(GL_CULL_FACE);
		glDisable(GL_TEXTURE_2D);
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
			GLboolean old_cull_face_config = glIsEnabled(GL_CULL_FACE);
			glEnable(GL_CULL_FACE);
			glCullFace(GL_BACK);
			glLineWidth(1.0f);
			glColor3f(0.6f, 0.0f, 0.0f);
			glDisable(GL_LIGHTING);
			RPVD& rpvd = packer->get_rpvd();
			glBegin(GL_LINES);
			for (Halfedge_iterator eit = rpvd.edges_begin(); eit != rpvd.edges_end(); ++eit)
			{
				//if (eit->is_border() || eit->opposite()->is_border())
				//	continue;
				//if (!eit->is_border()&& !eit->opposite()->is_border() && rpvd.is_delaunay_edge(eit))
					glColor3f(0.6f, 0.0f, 0.0f);
				//else
				//	glColor3f(0.0f, 0.0f, 0.4f);
				RDT_data_structure::Vertex_handle vh0 = eit->vertex(), vh1 = eit->opposite()->vertex();
				glPoint_3(vh0->mp);
				glPoint_3(vh1->mp);
				
			}
			glEnd();
			glEnable(GL_LIGHTING);
			if (!old_cull_face_config)
				glDisable(GL_CULL_FACE);
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
		
		for (unsigned int i = holes.size()-1; i < holes.size(); i++)
		{
			Packer::Hole h = holes[i];
			//gl_table_color(i);
			glColor3f(1.0, 0.0f, 0.0f);
			glLineWidth(1.5f);
			glBegin(GL_LINES);
			for (unsigned int j = 0; j < h.size(); j++)
			{
				//Halfedge_handle e = h[j];
				//Halfedge_handle oe = e->opposite();
				glPoint_3(h[j].first->mp);
				glPoint_3(h[j].second->mp);
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

	void SPM_Graphics::cur_color_map(double cur, double min_cur, double max_cur, GLfloat& r, GLfloat& g, GLfloat& b)
	{
		if ( std::fabs(max_cur - min_cur) < 1.0 || min_cur > max_cur)
		{
			r = 0.0f;
			g = 0.0f;
			b = 1.0f;
		}
		else
		{
			double mean = (min_cur + max_cur)/2.0;
			if (cur <= mean)
			{
				r = 0.0f ;
				b = ( mean - cur ) / ( mean - min_cur );
				g = ( cur - min_cur ) / ( mean - min_cur );
			}
			else
			{
				b = 0.0f;
				g = ( max_cur - cur ) / ( max_cur - mean );
				r = ( cur - mean ) / ( max_cur - mean );
			}
		}
	}
	void SPM_Graphics::build_curv_color_list()
	{
	
		const TriMesh& mesh = packer->mesh_domain();
		double max_curvature = std::numeric_limits<double>::min();
		double min_curvature = std::numeric_limits<double>::max();
		for (unsigned int i = 0; i < mesh.nb_vertices(); i++)
		{
			const MeshVertex& v = mesh.vertex(i);
			max_curvature = std::max(v.curvature, max_curvature);
			min_curvature = std::min(v.curvature, min_curvature);
		}

		double xmin, xmax, ymin, ymax, zmin, zmax;
		packer->get_bbox(xmin, ymin, zmin, xmax, ymax, zmax);

		cur_color_displist = glGenLists(1);
		glNewList(cur_color_displist, GL_COMPILE);
		GLint old_shade_model;
		glGetIntegerv(GL_SHADE_MODEL, &old_shade_model);
		glDisable(GL_LIGHTING);
		glShadeModel(GL_SMOOTH);
		GLfloat r, g, b;

		// draw color bar
		glMatrixMode(GL_PROJECTION);
		glPushMatrix();	
		glLoadIdentity();
		int glut_viewer_W = 800;
		int glut_viewer_H = 800;
		glut_viewer_get_screen_size(&glut_viewer_W, &glut_viewer_H);
		glOrtho(0, glut_viewer_W, glut_viewer_H, 0, 0, 1);
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();
		glBegin(GL_QUADS);	
		glColor3f(1.0f, 0.0f, 0.0f);
		glVertex2i(glut_viewer_W - 50, 30); 
		glVertex2i(glut_viewer_W - 30, 30);

		glColor3f(0.0f, 1.0f, 0.0f);
		glVertex2i(glut_viewer_W - 30, (glut_viewer_H)/2);
		glVertex2i(glut_viewer_W - 50, (glut_viewer_H)/2); 

		glColor3f(0.0f, 1.0f, 0.0f);
		glVertex2i(glut_viewer_W - 50, (glut_viewer_H)/2); 
		glVertex2i(glut_viewer_W - 30, (glut_viewer_H)/2);

		glColor3f(0.0f, 0.0f, 1.0f);
		glVertex2i(glut_viewer_W - 30, glut_viewer_H - 30); 
		glVertex2i(glut_viewer_W - 50, glut_viewer_H - 30); 
		glEnd();
		glMatrixMode(GL_MODELVIEW);
		glPopMatrix();
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();

		glBegin(GL_TRIANGLES);
		for (unsigned int i = 0; i < mesh.size(); i++)	
		{
			const Facet& f = mesh[i];
			for (int j = 0; j < 3; j++)
			{
				double v_cur = mesh.curvature_at_vertex(f.vertex_index[j]);
				cur_color_map(v_cur, min_curvature, max_curvature, r, g, b);
				glColor3f(r, g, b);
				glVertex(f.vertex[j]);
			}
		}
		glEnd();
		glShadeModel(old_shade_model);
		glEndList();
		glEnable(GL_LIGHTING);
	}

	void SPM_Graphics::build_texture_lib()
	{
		if (!packer->get_project_ioer().texture_specified())
			return;
		std::vector<std::string> texture_files;
		packer->get_project_ioer().read_texture_files(texture_files);
		if (texture_files.empty())
			return;
		texture_lib.resize(texture_files.size());
		for (unsigned int i = 0; i < texture_files.size(); i++)
		{
			cv::Mat m = cv::imread(texture_files[i]);
			if (!m.data)
			{
				std::cout<<"Error: cannot read texture file \""<<texture_files[i]<<"\"\n";
				return;
			}
			glGenTextures(1, &texture_lib[i]);
			glBindTexture(GL_TEXTURE_2D, texture_lib[i]);
			int r, c;
			int e = 0;
			while ((1<<e) < m.rows)
				e++;
			r = (1<<e);
			e = 0;
			while ((1<<e) < m.cols)
				e++;
			c = (1<<e);
			GLubyte *scaleTexels = (GLubyte*)malloc(r*c*3*sizeof(GLubyte));
			GLubyte *pixels = (GLubyte*)malloc(m.rows*m.cols*3*sizeof(GLubyte));
			GLubyte *ptr = pixels;
			// load content in m to 1D array pixels
			for (int y = 0; y < m.rows; y++)
				for (int x = 0; x < m.cols; x++)
				{
					cv::Vec3b pixel = m.at<Vec3b>(y, x);
					// to rgb model
					*ptr++ = pixel.val[2];
					*ptr++ = pixel.val[1];
					*ptr++ = pixel.val[0];
				}
			if (gluScaleImage(GL_RGB, m.cols, m.rows, GL_UNSIGNED_BYTE, pixels, c, r, GL_UNSIGNED_BYTE, scaleTexels))
			{
				std::cout<<"Error: Scale image failed.\n";
				return;
			}
			glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
			glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
			glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, c, r, 0, GL_RGB, GL_UNSIGNED_BYTE, scaleTexels);
			free(scaleTexels);
			free(pixels);
		}
		textured = true;
	}

	void SPM_Graphics::draw_cdt()
	{
		glColor3f(1.0f, 0.0f, 0.0f);
		glLineWidth(1.5f);
		glBegin(GL_LINES);
		CDT& cdt = packer->get_cdt();
		//for (CDT::All_edges_iterator eit = cdt.all_edges_begin(); eit != cdt.all_edges_end(); ++eit)
		//{
		//	if (cdt.is_infinite(eit))
		//		continue;
		//	CDT::Face_handle f = eit->first;
		//	int vi = eit->second;
		//	CDT::Vertex_handle v0 = f->vertex(f->cw(vi)), v1 = f->vertex(f->ccw(vi));
		//	if (!cdt.is_constrained(*eit) && v0->already_exist_in_rdt && v1->already_exist_in_rdt)
		//		continue;
		//	Point_3 p0 = v0->geo_info.prj_pnt, p1 = v1->geo_info.prj_pnt;
		//	glPoint_3(p0);
		//	glPoint_3(p1);
		//}
		for (CDT::All_faces_iterator fit = cdt.all_faces_begin(); fit != cdt.all_faces_end(); ++fit)
		{
			if (!fit->inside_hole)
				continue;
			CDT::Vertex_handle v[] = {fit->vertex(0), fit->vertex(1), fit->vertex(2)};
			
			for (int i = 0; i < 3; i++)
			{
				glPoint_3(v[i]->geo_info.prj_pnt);
				glPoint_3(v[(i+1)%3]->geo_info.prj_pnt);
			}
		}
		glEnd();
	}

	void SPM_Graphics::draw_smoothed_voronoi_cell()
	{
		const RestrictedPolygonVoronoiDiagram& rpvd = packer->get_rpvd();
		const std::vector<std::vector<Point_3>>& smoothed_regions = rpvd.get_smoothed_voronoi_regions();
		glDisable(GL_LIGHTING);
		glColor3f(0.0f, 0.0f, 1.0f);
		glLineWidth(1.5f);
		for (unsigned int i = 0; i < smoothed_regions.size(); i++)
		{
			const std::vector<Point_3>& smoothed_region = smoothed_regions[i];
			glBegin(GL_LINE_LOOP);
			for (unsigned int j = 0; j < smoothed_region.size(); j++)
				glPoint_3(smoothed_region[j]);
			glEnd();
		}
		glEnable(GL_LIGHTING);
	}
}