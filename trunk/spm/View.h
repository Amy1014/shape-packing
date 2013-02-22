#pragma  once

#include <vector>
#include <Geex/graphics/opengl.h>
#include <glut_viewer/glut_viewer.h>
#include <GL/glut.h>
#include <cv.h>
#include <highgui.h>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include "Doc.h"
//#include "spm_cgal.h"
#include "power_diagram.h"
//#include "SphericalBounder.h"
//debug
//#include "dbscan/clusters.h"
//#include "dbscan/db_scan.h"

using std::vector;
using CGAL::object_cast;
using namespace cv;

namespace Geex
{
	//class Doc;
	class View
	{
	public:
		View(Doc* doc) ;
		~View();
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
		GLboolean& show_holes() { return show_holes_; }
		void draw_holes();
		int& setHighlightedGroup()	{return highlighted;}
		void draw_power_cell();
		GLboolean& show_domain_boundary() { return show_boundary_; }
		void draw_domain_boundary();
		GLboolean& show_features() { return show_features_; }
		void draw_features();
		GLboolean& show_high_curvature() { return show_high_curvature_; }
		void draw_high_curvature();
		GLboolean& fill_polygon() { return fill_polygon_; }
		void fill_draw_polygons();
		GLboolean& apply_texture() { return apply_texture_; }
		void polygon_draw_with_texture();
		//rendering
		void draw_fractal_polygons();
		void thicken_polygons();
		//debug
		void draw_project_facets();
		void draw_undetected_holes();
		//void draw_centroids();

		void read_texture(const vector<string>& filenames);
		void read_backgroud_texture(const string& bg_tex_filename);
		void read_backgroud_texture(const vector<string>& bg_tex_filenames);
		void draw_fitting_points();

	private:
		//void check_edge_flag(const Facet& f, int e);
		void draw_trimesh(const TriMesh& m);
		bool check_edge_flag(const Facet& F, int e);
		void draw_mesh_skeleton(const TriMesh&);
		inline void draw_segment(const Segment_3 *s);
		inline void draw_ray(const Ray_3 *r);
		inline void glPoint_3 (const Point_3& p) { glVertex3d(p.x(), p.y(), p.z());}
		inline void triangulate_draw(const vector<Segment_3>& p, const vec3& n);//triangulate a polygon and draw
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
		GLboolean show_holes_;
		GLboolean show_boundary_;
		GLboolean show_features_;
		GLboolean show_high_curvature_;
		GLboolean fill_polygon_;
		GLboolean apply_texture_;
		int highlighted;

		//vector<Mat> texLib;//texture library
		vector<GLuint> texNames;
		vector<GLuint> bgTexNames;
		GLuint bgTexName;//back ground texture index
		vector<vector<pair<GLfloat, GLfloat>>> texCoord;

		// for fractal 
		typedef boost::variate_generator<boost::mt19937&, boost::random::uniform_real_distribution<>> Real_generator;
		//Real_generator rg;
		vector<double> random_offset;
	};

}//namespace Geex
