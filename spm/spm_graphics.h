#ifndef _SPM_GRAPHICS_H_
#define  _SPM_GRAPHICS_H_

//#undef _DEMO_

#ifdef WIN32
#include <windows.h>
#endif

#include <GL/GL.h>
#include <GL/glut.h>
#include <glut_viewer/glut_viewer.h>
#include <Geex/graphics/opengl.h>
#include <cv.h>
#include <CGAL/enum.h>
#include <CGAL/partition_2.h>
#include <highgui.h>
#include "packer.h"
#include "polygon_diagram.h"


namespace Geex
{
	class SPM_Graphics
	{

		typedef RestrictedPolygonVoronoiDiagram RPVD;
		typedef RPVD::Vertex_handle Vertex_handle;
		typedef RPVD::Vertex_const_handle Vertex_const_handle;
		typedef RPVD::Halfedge_handle Halfedge_handle;
		typedef RPVD::Halfedge_iterator	Halfedge_iterator;
		typedef RPVD::Vertex_iterator Vertex_iterator;
		typedef RPVD::Facet_iterator Facet_iterator;
		typedef RPVD::Vertex_const_iterator Vertex_const_iterator;
		typedef RPVD::Halfedge_const_iterator	Halfedge_const_iterator;
		typedef RPVD::Facet_const_iterator	Facet_const_iterator;

	public:

		typedef enum {OUTLINE_DRAW, FILL_DRAW, TEXTURE_DRAW, DISCRETE_SCALE_DRAWE} Polygon_draw_type;

		SPM_Graphics(Packer *_packer);

		~SPM_Graphics();

		void draw();

		GLboolean& show_mesh()	{ return show_mesh_; }

		GLboolean& show_multi_submeshes() { return show_multi_submeshes_; }

		GLboolean& show_feature_lines() { return show_feature_lines_; }

		GLboolean& show_vector_field() { return show_vector_field_ ; }

		GLboolean& show_multi_tiles() { return show_multi_tiles_; }

		GLboolean& show_polygons()	{ return show_polygons_; }
		GLboolean& show_inactive() { return show_inactive_; }

		GLboolean& show_voronoi_cell() { return show_voronoi_cell_; }
		GLboolean& show_midpoint_cell() { return show_midpoint_cell_; }

		GLboolean& show_smoothed_voronoi_cell() { return show_smoothed_voronoi_cell_; }

		void redraw_voronoi_cell() { if (glIsList(vc_displist))	glDeleteLists(vc_displist, 1); }

		GLboolean& show_triangulation() { return show_triangulation_; }
		void redraw_triangulation() { if (glIsList(triangulation_displist))	glDeleteLists(triangulation_displist, 1); }

		GLboolean& show_vertices() { return show_vertices_; }

		GLboolean& show_tiles() { return show_polygons_; }

		GLboolean& show_hole_triangles() { return show_hole_triangles_; }
		GLboolean& show_holes() { return show_holes_; }
		GLboolean& show_vicinity() { return show_vicinity_; }
		GLboolean& show_swapped_polygons() { return show_swapped_polygons_; }
		GLboolean& show_curvatures() { return show_curvatures_; }

		Polygon_draw_type& get_polygon_draw_type() { return how_to_draw_polygons; }

		int& highlighted_group_id() { return highlighted_group; }

		GLboolean& show_local_frame() { return show_local_frame_; }
		//GLboolean& show_cdt() { return show_cdt_; }

		void load_graphics() 
		{
			build_curv_color_list();
			build_texture_lib();
			//if (textured)
			//	glEnable(GL_TEXTURE_2D);
		}

		// file io
		void save_materials();

	protected:
		// for multiple polygon sets display
		void load_next_texture_lib();

	private:
		inline void gl_table_color(int index);

		void outline_draw_polygons(const std::vector<Packing_object>& tiles);
		void fill_draw_polygons(const std::vector<Packing_object>& tiles);
		void texture_draw_polygons(const std::vector<Packing_object>& tiles);
		//void discrete_draw_polygons(const std::vector<Packing_object>& tiles);

		void draw_mesh();
		void draw_multi_submeshes(); 
		void draw_feature_line();

		void draw_multi_tiles();

		void draw_polygons();
		void draw_triangulation();

		void draw_voronoi_cell();
		void draw_midpoint_cell();

		void draw_smoothed_voronoi_cell();

		void draw_hole_triangles();
		void draw_holes();
		void draw_vicinity();
		void draw_swapped_polygons();
		void draw_curvature();
		void build_curv_color_list(); 
		void build_texture_lib();

		bool build_texture_from_files(std::vector<GLuint>& text_indices, const std::vector<std::string>& text_files, int cv_flag = 1);

		void clear_all_textures();

		void cur_color_map(double cur, double min_cur, double max_cur, GLfloat& r, GLfloat& g, GLfloat& b);
		// debug
		void draw_all_vertices();
		//void draw_cdt(); 

		void draw_local_frames();

		/** call opengl **/
		inline void glPoint_3(const Point_3& p) { glVertex3d(p.x(), p.y(), p.z()); }

	public:
		GLfloat tunable_red;
		GLfloat tunable_green;
		GLfloat tunable_blue;

	private:
		Packer* packer;
		Polygon_draw_type how_to_draw_polygons;
		GLboolean show_mesh_;
		GLboolean show_multi_submeshes_;
		GLboolean show_multi_tiles_;
		GLboolean show_vector_field_;

		GLboolean show_feature_lines_;

		GLboolean show_polygons_;
		GLboolean show_voronoi_cell_;
		GLboolean show_midpoint_cell_;
		GLboolean show_smoothed_voronoi_cell_;

		GLboolean show_triangulation_;
		GLboolean show_vertices_;
		GLboolean show_hole_triangles_;
		GLboolean show_holes_;
		GLboolean show_curvatures_;
		GLboolean show_vicinity_;

		GLboolean show_inactive_;
		GLboolean show_swapped_polygons_;
		GLboolean textured;
		GLboolean alpha_texture;
		std::vector<GLuint> texture_lib;

		// for multiple sets of polygons
		std::vector<std::vector<GLuint>> multi_texture_libs;
		unsigned int sub_pack_id;

		/** display list **/
		GLuint triangulation_displist;
		GLuint vc_displist;//voronoi cell display list
		GLuint cur_color_displist;

		// debug
		GLboolean show_local_frame_;
		//GLboolean show_cdt_;
		int highlighted_group;
		// colors
		int random_color_index;
		static const GLfloat c1, c2, c3;
		static GLfloat color_table[24][3];

		static GLfloat surf_diff[4];
		static GLfloat surf_spec[4];
		static GLfloat surf_edge[4];
		//static GLfloat pgn_edge[4];
		static GLfloat shininess[];
	};
}

#endif
