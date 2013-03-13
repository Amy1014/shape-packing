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

	public:
		SPM_Graphics(Packer *_packer);
		void draw();
		GLboolean& show_mesh()	{ return show_mesh_; }
		void draw_mesh();
		GLboolean& show_polygons()	{ return show_polygons_; }
		void draw_polygons();
		GLboolean& show_voronoi_cell() { return show_voronoi_cell_; }
		GLboolean& show_triangulation() { return show_triangulation_; }

		int& highlighted_group_id() { return highlighted_group; }

		// debug
		void draw_all_vertices();

		typedef enum {OUTLINE_DRAW, FILL_DRAW, TEXTURE_DRAW} Polygon_draw_type;

	private:
		inline void gl_table_color(int index);
		void outline_draw_polygons();

		/** call opengl **/
		inline void glPoint_3(const Point_3& p) { glVertex3d(p.x(), p.y(), p.z()); }

	private:
		Packer* packer;
		Polygon_draw_type how_to_draw_polygons;
		GLboolean show_mesh_;
		GLboolean show_polygons_;
		GLboolean show_voronoi_cell_;
		GLboolean show_triangulation_;

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
