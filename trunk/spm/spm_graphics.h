#ifndef _SPM_GRAPHICS_H_
#define  _SPM_GRAPHICS_H_

#ifdef WIN32
#include <windows.h>
#endif

#include <GL/GL.h>
#include <GL/glut.h>
#include <glut_viewer/glut_viewer.h>
#include <Geex/graphics/opengl.h>
#include "packer.h"
#include "polygon_diagram.h"

namespace Geex
{
	class SPM_Graphics
	{

		typedef RestrictedPolygonVoronoiDiagram RPVD;
		typedef RPVD::Vertex_handle Vertex_handle;
		typedef RPVD::Halfedge_handle Halfedge_handle;
		typedef RPVD::Halfedge_iterator	Halfedge_iterator;
		typedef RPVD::Vertex_iterator Vertex_iterator;
		typedef RPVD::Facet_iterator Facet_iterator;

	public:

		typedef enum {OUTLINE_DRAW, FILL_DRAW, TEXTURE_DRAW} Polygon_draw_type;

		SPM_Graphics(Packer *_packer);

		~SPM_Graphics();

		void draw();

		GLboolean& show_mesh()	{ return show_mesh_; }

		GLboolean& show_polygons()	{ return show_polygons_; }

		GLboolean& show_voronoi_cell() { return show_voronoi_cell_; }
		void redraw_voronoi_cell() { glDeleteLists(vc_displist, 1); }

		GLboolean& show_triangulation() { return show_triangulation_; }
		void redraw_triangulation() { glDeleteLists(triangulation_displist, 1); }

		GLboolean& show_vertices() { return show_vertices_; }

		GLboolean& show_tiles() { return show_polygons_; }

		GLboolean& show_hole_triangles() { return show_hole_triangles_; }
		GLboolean& show_holes() { return show_holes_; }

		Polygon_draw_type& get_polygon_draw_type() { return how_to_draw_polygons; }

		int& highlighted_group_id() { return highlighted_group; }

		GLboolean& show_local_frame() { return show_local_frame_; }

	private:
		inline void gl_table_color(int index);
		void outline_draw_polygons();
		void fill_draw_polygons();
		void draw_mesh();
		void draw_polygons();
		void draw_triangulation();
		void draw_voronoi_cell();
		void draw_hole_triangles();
		void draw_holes();
		// debug
		void draw_all_vertices();

		void draw_local_frames();
		/** call opengl **/
		inline void glPoint_3(const Point_3& p) { glVertex3d(p.x(), p.y(), p.z()); }

	private:
		Packer* packer;
		Polygon_draw_type how_to_draw_polygons;
		GLboolean show_mesh_;
		GLboolean show_polygons_;
		GLboolean show_voronoi_cell_;
		GLboolean show_triangulation_;
		GLboolean show_vertices_;
		GLboolean show_hole_triangles_;
		GLboolean show_holes_;

		GLboolean show_local_frame_;

		/** display list **/
		GLuint triangulation_displist;
		GLuint vc_displist;//voronoi cell display list

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
