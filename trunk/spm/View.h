#pragma  once

#include <vector>
#include <Geex/graphics/opengl.h>
#include <glut_viewer/glut_viewer.h>
#include <GL/glut.h>
#include "Doc.h"
#include "spm_cgal.h"
#include "power_diagram.h"
#include "SphericalBounder.h"

using std::vector;
using CGAL::object_cast;

namespace Geex
{
	//class Doc;

	class View
	{
	public:
		View(Doc* doc) ;
		void print_screen();
		void draw() ;
		GLboolean& show_domain_mesh() { return show_domain_mesh_ ; }
		void draw_domain_mesh();
		GLboolean& show_mesh_skeleton() { return show_mesh_skeleton_; }
		void draw_mesh_skeleton();
		GLboolean& show_polygons() {return show_polygons_; }
		void draw_polygons();
		GLboolean& show_bounding_spheres() { return show_bounding_spheres_; }
		void draw_bounding_spheres();
		GLboolean& show_bisectors()	{ return show_bisectors_; }
		void draw_bisectors();
		GLboolean& show_regular_triangulation() {return show_regular_triangulation_;}
		void draw_triangulation();
	private:
		//void check_edge_flag(const Facet& f, int e);
		void draw_trimesh(const TriMesh& m);
		bool check_edge_flag(const Facet& F, int e);
		void draw_mesh_skeleton(const TriMesh&);
		inline void draw_segment(const Segment_3 *s);
		inline void draw_ray(const Ray_3 *r);
		inline void glPoint_3 (const Point_3& p) { glVertex3d(p.x(), p.y(), p.z());}
		//debug
		void draw_selected_triangles();
	private:
		Doc* pDoc_ ;
		GLboolean show_domain_mesh_;
		GLboolean show_mesh_skeleton_;
		GLboolean show_polygons_;
		GLboolean show_bounding_spheres_;
		GLboolean show_bisectors_;
		GLboolean show_regular_triangulation_;
	};

}//namespace Geex