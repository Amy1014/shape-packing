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
	//GLfloat SPM_Graphics::surf_diff[4] = {0.3922f, 0.5843f, 0.9294f, 0.5f}; // shadow blue
	//GLfloat SPM_Graphics::surf_diff[4] = {0.94117647f, 0.8705882f, 0.725490196f, 1.0f}; // shadow brown
	GLfloat SPM_Graphics::surf_diff[4] = {0.912745f, 0.912745f, 0.912745f, 1.0f}; // shadow gray
	GLfloat SPM_Graphics::surf_spec[4] = {0.5f, 0.5f, 0.5f, 0.5f} ;
	GLfloat SPM_Graphics::surf_edge[4] = {0.1f, 0.1f, 0.1f, 1.0f};
	//GLfloat SPM_Graphics::pgn_edge[4] = {0.1f, 0.1f, 0.1f, 1.0f};
	GLfloat SPM_Graphics::shininess[] = {150.0f};

	SPM_Graphics::SPM_Graphics(Packer *_packer): packer(_packer)
	{
		how_to_draw_polygons = OUTLINE_DRAW;
		show_mesh_ = true;
		show_feature_lines_ = false;
		show_polygons_ = true;
		show_vector_field_ = false;
		show_voronoi_cell_ = false;
		show_smoothed_voronoi_cell_ = false;
		show_midpoint_cell_ = true;
		show_triangulation_ = false;
		show_vertices_ = false;
		show_hole_triangles_ = false;
		show_local_frame_ = false;
		show_curvatures_ = false;
		textured = false;
		show_holes_ = false;
		//show_cdt_ = false;
		highlighted_group = -1;
		show_multi_tiles_ = false;
		show_multi_submeshes_ = false;
		show_inactive_ = true;
		alpha_texture = false;
		sub_pack_id = 0;

		tunable_red = surf_diff[0];
		tunable_green = surf_diff[1];
		tunable_blue = surf_diff[2];

		show_vicinity_ = false;
		show_swapped_polygons_ = false;
	}

	SPM_Graphics::~SPM_Graphics()
	{
		glDeleteLists(triangulation_displist, 1);
		glDeleteLists(vc_displist, 1);
		glDeleteLists(cur_color_displist, 1);
		clear_all_textures();
	}

	void SPM_Graphics::draw_feature_line()
	{
		const TriMesh& mesh = packer->mesh_domain();
		const vector<pair<int, int>>& feature_edges = mesh.getFeatureEdges();
		glDisable(GL_LIGHTING);
		glLineWidth(1.5f);
		glColor3f(1.0f, 0.0f, 0.0f);
		glBegin(GL_LINES);
		for (unsigned int i = 0; i < feature_edges.size(); i++)
		{
			const vec3& v0 = mesh.vertex(feature_edges[i].first).pos_;
			const vec3& v1 = mesh.vertex(feature_edges[i].second).pos_;
			glVertex(v0);
			glVertex(v1);
		}
		glEnd();
		glEnable(GL_LIGHTING);
	}

	void SPM_Graphics::clear_all_textures()
	{
		if (multi_texture_libs.empty())
		{
			for (unsigned int i = 0; i < texture_lib.size(); i++)
				glDeleteTextures(1, &texture_lib[i]);		
		}
		else
		{
			for (unsigned int i = 0; i < multi_texture_libs.size(); i++)
				for (unsigned int j = 0; j < multi_texture_libs[i].size(); j++)
					glDeleteTextures(1, &multi_texture_libs[i][j]);
			std::for_each(multi_texture_libs.begin(), multi_texture_libs.end(), std::mem_fun_ref(&std::vector<GLuint>::clear));
			multi_texture_libs.clear();
		}
		texture_lib.clear();
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
		if (show_feature_lines_)
			draw_feature_line();

		if (show_swapped_polygons_)
			draw_swapped_polygons();

		if (show_polygons_)
			draw_polygons();
		if (show_triangulation_)
			draw_triangulation();
		if (show_voronoi_cell_)
			draw_voronoi_cell();
		if (show_smoothed_voronoi_cell_)
			draw_smoothed_voronoi_cell();
		if (show_midpoint_cell_)
			draw_midpoint_cell();

		if (show_vertices_)
			draw_all_vertices();
		if (show_hole_triangles_)
			draw_hole_triangles();
		if (show_holes_)
			draw_holes();
		if (show_local_frame_)
			draw_local_frames();
		//if (show_cdt_)
		//	draw_cdt();
		if (show_curvatures_)
			glCallList(cur_color_displist);

		if (show_multi_submeshes_)
			draw_multi_submeshes();
		if (show_multi_tiles_)
			draw_multi_tiles();
		if (show_vicinity_)
			draw_vicinity();

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
		vec3 mc = m.mesh_center();
		for(unsigned int i=0; i<m.size(); i++) 
		{
			vec3 n = m[i].normal();
			glBegin(GL_TRIANGLES);
			glNormal(n);
//#ifdef _DEMO_
			glVertex((m[i].vertex[0]-mc)*(1-3.0e-3) + mc) ;
			glVertex((m[i].vertex[1]-mc)*(1-3.0e-3) + mc) ;
			glVertex((m[i].vertex[2]-mc)*(1-3.0e-3) + mc) ;
// #else
// 			glVertex(m[i].vertex[0]) ;
// 			glVertex(m[i].vertex[1]) ;
// 			glVertex(m[i].vertex[2]) ;
// #endif
			glEnd();
			
			//glDisable(GL_LIGHTING);
			//glBegin(GL_LINE_LOOP);
			//glVertex(m[i].vertex[0]) ;
			//glVertex(m[i].vertex[1]) ;
			//glVertex(m[i].vertex[2]) ;
			//glEnd();
			//glEnable(GL_LIGHTING);
		}
		glPopMatrix();
		glDisable(GL_LIGHTING);	
// #ifdef _DEMO_
// 		const TriMesh::BoundaryEdgeSet& bd_edges = m.getBoundary();
// 		glColor3f(0.0f, 0.0f, 0.0f);
// 		glBegin(GL_LINES);
// 		for (TriMesh::BoundaryEdgeSet::const_iterator it = bd_edges.begin(); it != bd_edges.end(); ++it)
// 		{
// 			int v0 = it->first, v1 = it->second;
// 			glVertex((m.vertex(v0).pos_ - mc) * (1-1.0e-4) + mc);
// 			glVertex((m.vertex(v1).pos_ - mc) * (1-1.0e-4) + mc);
// 		}
// 		glEnd();
// #endif

		if (m.with_vector_field() && show_vector_field_)
		{
			glColor3f(0.0f, 0.0f, 1.0f);
			glBegin(GL_LINES);
			for (unsigned int i = 0; i < m.nb_vertices(); i++)
			{
				glVertex(m.vertex(i).pos_);
				glVertex(m.vertex(i).pos_ + m.vector_field_at(i));
			}
			glEnd();
		}
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
		const std::vector< std::vector<Packing_object> >& multi_tiles = packer->get_multigroup_tiles();

		static GLfloat amb_mat[][4] = {{0.19225f, 0.19225f, 0.19225f, 1.0f}, {0.19225f, 0.19225f, 0.19225f, 1.0f}, {0.19225f, 0.19225f, 0.19225f, 1.0f}};
		static GLfloat diff_mat[][4] = {{0.2f, 0.2f, 0.2f, 1.0f}, {1.0f, 1.0f, 1.0f, 1.0f}, {0.50754f, 0.50754f, 0.50754f, 1.0f}};
		static GLfloat spec_mat[][4] = {{0.608273f, 0.608273f, 0.608273f, 1.0f}, {0.608273f, 0.608273f, 0.608273f, 1.0f}, {0.608273f, 0.608273f, 0.608273f, 1.0f}};
		static GLfloat pgn_shininess = 0.4;

		for (unsigned int k = 0; k < multi_tiles.size(); k++)
		{
			switch (how_to_draw_polygons)
			{
			case OUTLINE_DRAW:
				outline_draw_polygons(multi_tiles[k]);
				break;
			case FILL_DRAW:
				fill_draw_polygons(multi_tiles[k]);
				break;
			case TEXTURE_DRAW:
				if (!textured)
					outline_draw_polygons(multi_tiles[k]);
				else
				{
					texture_lib.clear();
					texture_lib = multi_texture_libs[k];
					texture_draw_polygons(multi_tiles[k]);
				}
				break;
			}
		}	
	}
	void SPM_Graphics::draw_polygons()
	{
		const vector<Packing_object>& tiles = packer->get_tiles();
		switch (how_to_draw_polygons)
		{
		case OUTLINE_DRAW:
			outline_draw_polygons(tiles);
			break;
		case FILL_DRAW:
			fill_draw_polygons(tiles);
			break;
		case TEXTURE_DRAW:
			if (textured)
				texture_draw_polygons(tiles);
			else
				outline_draw_polygons(tiles);
			break;
		//case DISCRETE_SCALE_DRAWE:
		//	discrete_draw_polygons(tiles);
		//	break;
		}
	}

	void SPM_Graphics::outline_draw_polygons(const std::vector<Packing_object>& tiles)
	{
		static GLfloat pgn_edge[4] = {0.1f, 0.1f, 0.1f, 1.0f};
		//const vector<Packing_object>& tiles = packer->get_tiles();
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

	void SPM_Graphics::fill_draw_polygons(const std::vector<Packing_object>& tiles)
	{
		//static GLfloat amb_mat[] = {0.19225f, 0.19225f, 0.19225f, 1.0f};
		//static GLfloat diff_mat[] = {0.50754f, 0.50754f, 0.50754f, 1.0f};
		static GLfloat diff_mat[] = {0.8f, 0.5078431f, 0.31568627f, 1.0f};
		static GLfloat spec_mat[] = {0.608273f, 0.608273f, 0.608273f, 1.0f};
		//static GLfloat spec_mat[] = {0.256777f,	0.137622f, 0.086014f, 1.0f};
		static GLfloat pgn_shininess = 1500.0f;

		//const vector<Packing_object>& tiles = packer->get_tiles();

		glCullFace(GL_FRONT);
		glEnable(GL_LIGHTING);
		glPushMatrix();
		//glMaterialfv(GL_FRONT, GL_AMBIENT, amb_mat);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, diff_mat);
		glMaterialfv(GL_FRONT, GL_SPECULAR, spec_mat);
		glMaterialf(GL_FRONT, GL_SHININESS, 50.0f);
		for (unsigned int i = 0; i < tiles.size(); i++)
		{
			Vector_3 n = tiles[i].norm();
			glNormal3d(n.x(), n.y(), n.z());
//			glBegin(GL_TRIANGLES);
//			Point_3 c = tiles[i].centroid();
//			Vector_3 n = tiles[i].norm();
//			for (unsigned int j = 0; j < tiles[i].size(); j++)
//			{
//#ifdef _DEMO_
//				glPoint_3(c + 1.0e-2*n);
//				glPoint_3(tiles[i].vertex(j) + 1.0e-2*n);
// 				glPoint_3(tiles[i].vertex((j+1)%tiles[i].size()) + 1.0e-2*n);
//#else
//				glPoint_3(c);
//				glPoint_3(tiles[i].vertex(j));
//				glPoint_3(tiles[i].vertex((j+1)%tiles[i].size()));
//#endif
//			}
//			glEnd();
			typedef CGAL::Partition_traits_2< CGAL::Exact_predicates_inexact_constructions_kernel >::Polygon_2 Partition_polygon;

			std::list<Partition_polygon> convex_partitions;
			Polygon_2 pgn2d;
			Local_frame lf = tiles[i].local_frame();
			for (unsigned int j = 0; j < tiles[i].size(); j++)
				pgn2d.push_back(lf.to_uv(tiles[i].vertex(j)));
			CGAL::approx_convex_partition_2( pgn2d.vertices_begin(), pgn2d.vertices_end(), std::back_inserter(convex_partitions) );

			for (std::list<Partition_polygon>::iterator it = convex_partitions.begin(); it != convex_partitions.end(); ++it)
			{
				glBegin(GL_POLYGON);
				for (Partition_polygon::Vertex_const_iterator vit = it->vertices_begin(); vit != it->vertices_end(); ++vit)
					glPoint_3(lf.to_xy(*vit) + 1.0e-2*n );
				glEnd();
			}
		}

		glPopMatrix();
		glDisable(GL_LIGHTING);
		glLineWidth(0.6f);
		for (unsigned int i = 0; i < tiles.size(); i++)
		{
			glColor3f(0.0f, 0.0f, 0.7f);
			glBegin(GL_LINE_LOOP);
			for (unsigned int j = 0;j < tiles[i].size(); j++)
			{
				glPoint_3(tiles[i].vertex(j));
			}
			glEnd();
		}
		glLineWidth(1.0f);

	}

	void SPM_Graphics::texture_draw_polygons(const std::vector<Packing_object>& tiles)
	{
		static GLfloat mat[] = {0.9f, 0.9f, 0.9f, 1.0f};
		//const std::vector<Packing_object>& tiles = packer->get_tiles();
		glEnable(GL_TEXTURE_2D);
		if (alpha_texture)
		{
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);          
			glEnable(GL_BLEND); 
			glEnable(GL_ALPHA_TEST);
			glAlphaFunc(GL_GREATER, 0.5f);  
		}
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
		if (!old_cull_face_config)
			glDisable(GL_CULL_FACE);
		glDisable(GL_TEXTURE_2D);
		if (alpha_texture)
		{
			glDisable(GL_BLEND);
			glDisable(GL_ALPHA_TEST);
		}
	}

	//void SPM_Graphics::discrete_draw_polygons(const std::vector<Packing_object>& tiles)
	//{
	//	// draw color bar
	//	GLint old_shade_model;
	//	glGetIntegerv(GL_SHADE_MODEL, &old_shade_model);
	//	glDisable(GL_LIGHTING);
	//	glShadeModel(GL_SMOOTH);
	//	GLfloat r, g, b;

	//	// draw color bar
	//	glMatrixMode(GL_PROJECTION);
	//	glPushMatrix();	
	//	glLoadIdentity();
	//	int glut_viewer_W = 800;
	//	int glut_viewer_H = 800;
	//	glut_viewer_get_screen_size(&glut_viewer_W, &glut_viewer_H);
	//	glOrtho(0, glut_viewer_W, glut_viewer_H, 0, 0, 1);
	//	glMatrixMode(GL_MODELVIEW);
	//	glPushMatrix();
	//	glLoadIdentity();
	//	glBegin(GL_QUADS);	
	//	glColor3f(1.0f, 0.0f, 0.0f);
	//	glVertex2i(glut_viewer_W - 50, 30); 
	//	glVertex2i(glut_viewer_W - 30, 30);

	//	glColor3f(0.0f, 1.0f, 0.0f);
	//	glVertex2i(glut_viewer_W - 30, (glut_viewer_H)/2);
	//	glVertex2i(glut_viewer_W - 50, (glut_viewer_H)/2); 

	//	glColor3f(0.0f, 1.0f, 0.0f);
	//	glVertex2i(glut_viewer_W - 50, (glut_viewer_H)/2); 
	//	glVertex2i(glut_viewer_W - 30, (glut_viewer_H)/2);

	//	glColor3f(0.0f, 0.0f, 1.0f);
	//	glVertex2i(glut_viewer_W - 30, glut_viewer_H - 30); 
	//	glVertex2i(glut_viewer_W - 50, glut_viewer_H - 30); 
	//	glEnd();
	//	glMatrixMode(GL_MODELVIEW);
	//	glPopMatrix();
	//	glMatrixMode(GL_PROJECTION);
	//	glPopMatrix();

	//	double min_scale = packer->min_scale_factor(), max_scale = packer->max_scale_factor();
	//	glBegin(GL_TRIANGLES);
	//	for (unsigned int i = 0; i < tiles.size(); i++)
	//	{
	//		double scale = tiles[i].factor;
	//		if (show_inactive_ && !tiles[i].active)
	//		{
	//			r = 0.6f;
	//			g = 0.6f;
	//			b = 0.6f;
	//		}
	//		else
	//			cur_color_map(scale, min_scale, max_scale, r, g, b);

	//		glColor3f(r, g, b);
	//		unsigned int nb_verts = tiles[i].size();
	//		Point_3 c = tiles[i].centroid();
	//		for (unsigned int j = 0; j < nb_verts; j++)
	//		{
	//			glPoint_3(c);
	//			glPoint_3(tiles[i].vertex(j));
	//			glPoint_3(tiles[i].vertex((j+1)%nb_verts));
	//		}
	//	}
	//	glEnd();
	//	glShadeModel(old_shade_model);
	//	glEnable(GL_LIGHTING);
	//}
	void SPM_Graphics::draw_all_vertices()
	{
		const RDT_data_structure& rdt = packer->get_rpvd().get_rdt();
		
		glDisable(GL_LIGHTING);
		glColor3f(0.0f, 0.0f, 0.0f);
		glPointSize(2.0f);
		glBegin(GL_POINTS);
		for (Vertex_const_iterator vi = rdt.vertices_begin(); vi != rdt.vertices_end(); ++vi)
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
			const RPVD& rpvd = packer->get_rpvd();
			glBegin(GL_LINES);
			for (Halfedge_const_iterator eit = rpvd.edges_begin(); eit != rpvd.edges_end(); ++eit)
			{
				//if (eit->is_border() || eit->opposite()->is_border())
				//	continue;
				//if (!eit->is_border()&& !eit->opposite()->is_border() && rpvd.is_delaunay_edge(eit))
					glColor3f(0.6f, 0.0f, 0.0f);
				//else
				//	glColor3f(0.0f, 0.0f, 0.4f);
				RDT_data_structure::Vertex_const_handle vh0 = eit->vertex();
				RDT_data_structure::Vertex_const_handle vh1 = eit->opposite()->vertex();
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
		const RPVD& rpvd = packer->get_rpvd();
		const std::vector<Packing_object>& po = packer->get_tiles();

//#ifndef _DEMO_
//#pragma message ("Non-demo code compiled!")
		////////////////////////////////  Outline draw ////////////////////////////////////////
		glDisable(GL_LIGHTING);
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
		//////////////////////////////////////////////////////////////////////////
// #else
// 		///////////////////////////////////  Demo draw ///////////////////////////
// 		GLfloat mat_diff[] = {tunable_red, tunable_green, tunable_blue, 1.0};
// 		glCullFace(GL_FRONT) ;
// 		glShadeModel(GL_SMOOTH);
// 		glPushMatrix();
// 		//glMaterialfv(GL_FRONT, GL_DIFFUSE, surf_diff);
// 		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diff);
// 		glMaterialfv(GL_FRONT, GL_SPECULAR, surf_spec);
// 		GLfloat shn[] = { 1500.0f };
// 		glMaterialfv(GL_FRONT, GL_SHININESS, shn);
// 		/////////////////////////////////////////////////////////////////////////		
// 		for (unsigned int i = 0; i < rpvd.number_of_groups(); i++)
// 		{
// 			const RestrictedPolygonVoronoiDiagram::VertGroup& vg = rpvd.sample_points_group(i);
// 			Point_3 c = po[i].centroid();
// 			Vector_3 n = po[i].norm();			
// 			for (unsigned int j = 0; j < vg.size(); j++)
// 			{
// 				const std::vector<Point_3>& vd_verts = vg[j]->no_prj_vd_vertices;				
// 				if (vd_verts.size() > 1)
// 				{
// 					glEnable(GL_LIGHTING);
// 					glBegin(GL_TRIANGLES);
// 					glNormal3d(n.x(), n.y(), n.z());
// 					for (unsigned int k = 0; k < vd_verts.size()-1; k++)
// 					{
// 						glPoint_3(c);
// 						glPoint_3(vd_verts[k]);
// 						glPoint_3(vd_verts[k+1]);
// 					}
// 					glEnd();
// 					glDisable(GL_LIGHTING);
// 					glColor3f(0.0f, 0.0f, 0.0f);
// 					glLineWidth(1.5f);
// 					glBegin(GL_LINES);
// 					for (unsigned int k = 0; k < vd_verts.size()-1; k++)
// 					{
// 						glPoint_3(vd_verts[k]);
// 						glPoint_3(vd_verts[k+1]);
// 					}
// 					glEnd();
// 				}
// 			}			
// 		}
// #endif
	}

	void SPM_Graphics::draw_hole_triangles()
	{
		const RPVD& rpvd = packer->get_rpvd();
		glColor3f(0.3f, 0.3f, 0.3f);
		glDisable(GL_LIGHTING);
		for (Facet_const_iterator fit = rpvd.faces_begin(); fit != rpvd.faces_end(); ++fit)
		{
			if (fit->vacant)
			{
				RPVD::Halfedge_const_handle eh = fit->halfedge();
				Vertex_const_handle v0 = eh->vertex();
				eh = eh->next();
				Vertex_const_handle v1 = eh->vertex();
				eh = eh->next();
				Vertex_const_handle v2 = eh->vertex();
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
		const std::vector<Packer::Hole>& holes = packer->get_holes();
		glDisable(GL_LIGHTING);
		//glColor3f(0.0f, 0.0f, 1.0f);
		
		for (unsigned int i = 0; i < holes.size(); i++)
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
		if ( std::fabs(max_cur - min_cur) < 0.001 || min_cur > max_cur)
		{
			r = 0.0f;
			g = 0.0f;
			b = 1.0f;
		}
		else
		{
			double mean = (min_cur + max_cur)/2.0;
			if (cur < min_cur)
			{
				r = 0.0f;
				g = 0.0f;
				b = 1.0f;
			}
			else if ( cur > max_cur )
			{
				r = 1.0f;
				g = 0.0f;
				b = 0.0f;
			}
			else if (cur <= mean)
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

	bool SPM_Graphics::build_texture_from_files(std::vector<GLuint>& text_indices, const std::vector<std::string>& text_files, int cv_flag)
	{
		text_indices.resize(text_files.size());
		for (unsigned int i = 0; i < text_files.size(); i++)
		{
			cv::Mat m = cv::imread(text_files[i], cv_flag);
			if (!m.data)
			{
				std::cout<<"Error: cannot read texture file \""<<text_files[i]<<"\"\n";
				return false;
			}
			glGenTextures(1, &text_indices[i]);
			glBindTexture(GL_TEXTURE_2D, text_indices[i]);
			int r, c;
			int e = 0;
			while ((1<<e) < m.rows)
				e++;
			r = (1<<e);
			e = 0;
			while ((1<<e) < m.cols)
				e++;
			c = (1<<e);
			int channel_sz = 3;
			if (cv_flag < 0)
				channel_sz = 4;
			GLubyte *scaleTexels = (GLubyte*)malloc(r*c*channel_sz*sizeof(GLubyte));
			GLubyte *pixels = (GLubyte*)malloc(m.rows*m.cols*channel_sz*sizeof(GLubyte));
			GLubyte *ptr = pixels;
			// load content in m to 1D array pixels
			if (cv_flag > 0)
				for (int y = 0; y < m.rows; y++)
					for (int x = 0; x < m.cols; x++)
					{
						cv::Vec3b pixel = m.at<Vec3b>(y, x);
						// to rgb model
						*ptr++ = pixel.val[2];
						*ptr++ = pixel.val[1];
						*ptr++ = pixel.val[0];
					}
			else
				for (int y = 0; y < m.rows; y++)
					for (int x = 0; x < m.cols; x++)
					{
						cv::Vec4b pixel = m.at<Vec4b>(y, x);
						// to rgb model
						*ptr++ = pixel.val[2];
						*ptr++ = pixel.val[1];
						*ptr++ = pixel.val[0];
						*ptr++ = pixel.val[3];
					}
			GLuint rgb_mode = GL_RGB;
			if (cv_flag < 0)
				rgb_mode = GL_RGBA;
			if (gluScaleImage(rgb_mode, m.cols, m.rows, GL_UNSIGNED_BYTE, pixels, c, r, GL_UNSIGNED_BYTE, scaleTexels))
			{
				std::cout<<"Error: Scale image failed.\n";
				return false;
			}

			glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
			glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
			glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
			glTexImage2D(GL_TEXTURE_2D, 0, rgb_mode, c, r, 0, rgb_mode, GL_UNSIGNED_BYTE, scaleTexels);
			free(scaleTexels);
			free(pixels);
		}
		return true;
	}
	void SPM_Graphics::build_texture_lib()
	{
		int cv_flag = 1;
		if (!packer->get_project_ioer().texture_specified())
			return;
		if (!packer->get_project_ioer().mesh_polygon_coupled())
		{
			std::vector<std::string> texture_files;
			packer->get_project_ioer().read_texture_files(texture_files);
			if (texture_files.empty())
				return;
			
			if (Geex::FileSystem::extension(texture_files[0]) == "png" || Geex::FileSystem::extension(texture_files[0]) == "PNG")
			{
				cv_flag = -1;
				std::cout<<"Loading texture images with alpha mode...\n";
			}
			if (! build_texture_from_files(texture_lib, texture_files, cv_flag) )
			{
				clear_all_textures();
				return;
			}
		}
		else
		{
			std::vector<std::vector<std::string>> multi_file_sets;
			packer->get_project_ioer().read_texture_files(multi_file_sets);
			if (multi_file_sets.empty())
				return;
			multi_texture_libs.resize(multi_file_sets.size());
			for (unsigned int i = 0; i < multi_file_sets.size(); i++)
			{
				if (Geex::FileSystem::extension(multi_file_sets[i][0]) == "png" || Geex::FileSystem::extension(multi_file_sets[i][0]) == "PNG")
					cv_flag = -1;
				if (! build_texture_from_files(multi_texture_libs[i], multi_file_sets[i], cv_flag) )
				{
					clear_all_textures();
					return;
				}
			}
		}

		textured = true;

		if (cv_flag < 0)
		{
			alpha_texture = true;  
		}
	}

	void SPM_Graphics::save_materials()
	{
		if (!packer->get_project_ioer().texture_specified())
			return;
		std::vector<std::string> texture_files;
		packer->get_project_ioer().read_texture_files(texture_files);
		if (texture_files.empty())
			return;
		// pick out only the file name
		// save material file
		std::string output_dir = packer->get_project_ioer().attribute_value("OutputDir");
		std::string mat_fn = output_dir + "/mat.mtl";
		std::ofstream of(mat_fn.c_str());
		for (unsigned int i = 0; i < texture_files.size(); i++)
		{
			unsigned int last_slash_pos = texture_files[i].rfind('/');
			if (last_slash_pos == std::string::npos)
				last_slash_pos = texture_files[i].rfind('\\');
			if (last_slash_pos == std::string::npos)
				last_slash_pos = 0;
			std::string fn(texture_files[i], last_slash_pos+1);
			unsigned int dot_pos = fn.find('.');
			std::string mat_name(fn, 0, dot_pos);
			of << "newmtl " << mat_name << std::endl;
			of << "Ka " << "0.8 0.8 0.8 \n";
			of << "Kd " << "0.8 0.8 0.8 \n";
			of << "illum 2 \n";
			of << "map_Ka " << fn << std::endl;
			of << "map_Kd " << fn << std::endl;
			of << std::endl;
		}
	}
	void SPM_Graphics::load_next_texture_lib()
	{
		texture_lib.clear();
		texture_lib = multi_texture_libs[sub_pack_id];
		sub_pack_id++;
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

	void SPM_Graphics::draw_midpoint_cell()
	{
		const RestrictedPolygonVoronoiDiagram& rpvd = packer->get_rpvd();
		const vector<Packing_object>& tiles = packer->get_tiles();
		glLineWidth(1.5f);

//#ifndef _DEMO_
		/////////////		Line Draw		//////////////////
 		glDisable(GL_LIGHTING);
 		glColor3f(0.1f, 0.1f, 0.1f);
 		for (unsigned int i = 0; i < rpvd.number_of_groups(); i++)
 		{
 			const RestrictedPolygonVoronoiDiagram::VertGroup& vg = rpvd.sample_points_group(i);
 			//Point_3 c = tiles[i].centroid();
 			//Vector_3 n = tiles[i].norm();
 		
 			glBegin(GL_LINES);
 			for (unsigned int j = 0; j < vg.size(); j++)
 			{
 				const std::vector<Segment_3>& med_segs = vg[j]->med_segs;
 				//const std::vector<Segment_3>& med_segs = vg[j]->no_prj_med_segs;
 				for (unsigned int k = 0; k < med_segs.size(); k++)
 				{
 					glPoint_3(med_segs[k].source());
 					glPoint_3(med_segs[k].target());
 				}
 			}
 			glEnd();
 			
 		}
 		glEnable(GL_LIGHTING);
		//////////////////////////////////////////////////////
//#else
//		/////////////////////////////// For demo /////////////////////////////////////
		//GLfloat mat_diff[] = {tunable_red, tunable_green, tunable_blue, 1.0};
 	//	glCullFace(GL_FRONT) ;
 	//	glShadeModel(GL_SMOOTH);
 	//	glPushMatrix();
 	//	//glMaterialfv(GL_FRONT, GL_DIFFUSE, surf_diff);
		//glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diff);
 	//	glMaterialfv(GL_FRONT, GL_SPECULAR, surf_spec);
 	//	GLfloat shn[] = { 1500.0f };
 	//	glMaterialfv(GL_FRONT, GL_SHININESS, shn);
 	//	/////////////////////////////////////////////////////////////////////////		
 	//	for (unsigned int i = 0; i < rpvd.number_of_groups(); i++)
 	//	{
 	//		const RestrictedPolygonVoronoiDiagram::VertGroup& vg = rpvd.sample_points_group(i);
 	//		Point_3 c = tiles[i].centroid();
 	//		Vector_3 n = tiles[i].norm();			
 	//		for (unsigned int j = 0; j < vg.size(); j++)
 	//		{
 	//			const std::vector<Segment_3>& med_segs = vg[j]->no_prj_med_segs;
 	//			//gl_table_color(i);
 	//			glEnable(GL_LIGHTING);
 	//			glBegin(GL_TRIANGLES);
 	//			for (unsigned int k = 0; k < med_segs.size(); k++)
 	//			{
 	//				glNormal3d(n.x(),n.y(), n.z());
 	//				glPoint_3(c);
 	//				glPoint_3(med_segs[k].source());
 	//				glPoint_3(med_segs[k].target());
 	//			}
 	//			glEnd();
 	//			glDisable(GL_LIGHTING);
 	//			glColor3f(0.0f, 0.0f, 0.0f);
 	//			glBegin(GL_LINES);
 	//			for (unsigned int k = 0; k < med_segs.size(); k++)
 	//			{
 	//				glPoint_3(med_segs[k].source());
 	//				glPoint_3(med_segs[k].target());
 	//			}
 	//			glEnd();
 	//		}			
 	//	}
//#endif
 		glLineWidth(1.0f);
	}

	void SPM_Graphics::draw_vicinity()
	{
		glDisable(GL_LIGHTING);
		glColor3f(0.0f, 0.8f, 0.0f);
		glLineWidth(4.0f);
		const std::list< std::list<Segment_3> >& vicinity_holes = packer->vicinity_holes;
		for (std::list< std::list<Segment_3> >::const_iterator it = vicinity_holes.begin(); it != vicinity_holes.end(); ++it)
		{
			glBegin(GL_LINES);
			for (std::list<Segment_3>::const_iterator sit = it->begin(); sit != it->end(); ++sit)
			{
				glVertex3d(sit->source().x(), sit->source().y(), sit->source().z());
				glVertex3d(sit->target().x(), sit->target().y(), sit->target().z());
			}
			glEnd();
		}
		glEnable(GL_LIGHTING);
	}

	void SPM_Graphics::draw_swapped_polygons()
	{
		//static GLfloat pgn_edge[4] = {0.1f, 0.1f, 0.1f, 1.0f};
		//glDisable(GL_LIGHTING);
		//glLineWidth(1.5f);
		const std::vector<Packing_object>& tiles = packer->pack_objects;
		//glColor4f(0.0f, 0.0f, 1.0f, 1.0f);
		//for (std::list<unsigned int>::const_iterator it = packer->swapped_this_time.begin(); it != packer->swapped_this_time.end(); ++it)
		//{			
		//	glBegin(GL_LINE_LOOP);
		//	for (unsigned int j = 0; j < tiles[*it].size(); j++)
		//		glPoint_3(tiles[*it].vertex(j));
		//	glEnd();
		//}
		//glEnable(GL_LIGHTING);

		static GLfloat diff_mat[] = {0.0f, 0.0f, 0.3f, 1.0f};
		static GLfloat spec_mat[] = {0.608273f, 0.608273f, 0.608273f, 1.0f};
		static GLfloat pgn_shininess = 1500.0f;

		//const vector<Packing_object>& tiles = packer->get_tiles();

		glCullFace(GL_FRONT);
		glEnable(GL_LIGHTING);
		glPushMatrix();
		glMaterialfv(GL_FRONT, GL_DIFFUSE, diff_mat);
		glMaterialfv(GL_FRONT, GL_SPECULAR, spec_mat);
		glMaterialf(GL_FRONT, GL_SHININESS, 50.0f);
		for (std::list<unsigned int>::const_iterator it = packer->swapped_this_time.begin(); it != packer->swapped_this_time.end(); ++it)
		{
			unsigned int i = *it;
			//glBegin(GL_TRIANGLES);
			//Point_3 c = tiles[i].centroid();
			Vector_3 n = tiles[i].norm();
			glNormal3d(n.x(), n.y(), n.z());
//			for (unsigned int j = 0; j < tiles[i].size(); j++)
//			{
//#ifdef _DEMO_
//				glPoint_3(c + 1.0e-2*n);
//				glPoint_3(tiles[i].vertex(j) + 1.0e-2*n);
//				glPoint_3(tiles[i].vertex((j+1)%tiles[i].size()) + 1.0e-2*n);
//#else
//				glPoint_3(c);
//				glPoint_3(tiles[i].vertex(j));
//				glPoint_3(tiles[i].vertex((j+1)%tiles[i].size()));
//#endif
//			}
			//glEnd();
			typedef CGAL::Partition_traits_2< CGAL::Exact_predicates_inexact_constructions_kernel >::Polygon_2 Partition_polygon;

			std::list<Partition_polygon> convex_partitions;
			Polygon_2 pgn2d;
			Local_frame lf = tiles[i].local_frame();
			for (unsigned int j = 0; j < tiles[i].size(); j++)
				pgn2d.push_back(lf.to_uv(tiles[i].vertex(j)));
			CGAL::approx_convex_partition_2( pgn2d.vertices_begin(), pgn2d.vertices_end(), std::back_inserter(convex_partitions) );

			for (std::list<Partition_polygon>::iterator it = convex_partitions.begin(); it != convex_partitions.end(); ++it)
			{
				glBegin(GL_POLYGON);
				for (Partition_polygon::Vertex_const_iterator vit = it->vertices_begin(); vit != it->vertices_end(); ++vit)
					glPoint_3(lf.to_xy(*vit));
				glEnd();
			}
		}
	}
}