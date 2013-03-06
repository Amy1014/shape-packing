#ifndef _POLYGON_DIAGRAM_H_
#define _POLYGON_DIAGRAM_H_

/************************************************************************/
/*  Approximate restricted Voronoi Diagram between polygons on a mesh   */
/*                   													*/
/************************************************************************/

#ifdef WIN32
#include <windows.h>
#endif
#include <GL/GL.h>
#include <GL/GLU.h>

#include <vector>
//#include <LpCVT/combinatorics/delaunay.h>
#include <Geex/cvt/delaunay.h>
#include <Geex/cvt/RVD.h>
//#include <LpCVT/combinatorics/RVD.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include "spm_cgal.h"
 
namespace Geex
{

	/** representation of restricted Delaunay triangulation **/
	struct Embedded_point_2 : public K
	{
		typedef K::Point_3 Point_2;
	};

	template < class Gt, class Vb = CGAL::Triangulation_vertex_base_2<Gt> >
	class Embedded_vertex_2 : public Vb // vertex of 2d triangulation embedding in 3d space
	{
	public:
		typedef typename Vb::Vertex_handle Vertex_handle;
		typedef typename Vb::Face_handle Face_handle;
		typedef typename Vb::Point		Point;

		template < typename TDS2 > 
		struct Rebind_TDS
		{
			typedef typename Vb::template Rebind_TDS<TDS2>::Other Vb2;
			typedef Embedded_vertex_2<Gt, Vb2>	Other;
		};
	public:
		Embedded_vertex_2() : Vb() {}
		Embedded_vertex_2(const Point& p) : Vb(p) {}
		Embedded_vertex_2(const Point& p, Face_handle f) : Vb(f, p) {}
		Embedded_vertex_2(Face_handle f) : Vb(f) {}
	};

	template < class Gt, class Fb = CGAL::Triangulation_face_base_2<Gt> >
	class Embedded_face_2 : public Fb
	{					};

	typedef CGAL::Triangulation_data_structure_2<Embedded_vertex_2<Embedded_point_2>, Embedded_face_2<K>> RDT_data_structure; // representation of RDT


class RestrictedPolygonVoronoiDiagram
{
	typedef RestrictedVoronoiDiagram_poly RestrictedVoronoiDiagram;
	

public:
	RestrictedPolygonVoronoiDiagram();
	~RestrictedPolygonVoronoiDiagram();

	/** set functions **/
	void set_mesh(const std::string& mesh_file_name);

	void begin_insert() ;
	// insert one polygon
	void insert_polygons(const Polygon_3& polygon, unsigned int group_id, unsigned int samp_nb = 24); 
	// insert a group of polygons
	template<class InputIterator> void insert_polygons(InputIterator first, InputIterator last, unsigned int samp_nb = 24);

	void end_insert() ;

	void iDT_update(); // edge flip to preserve intrinsic Delaunay structure

	/** access functions **/

	/** debug **/
	void draw_DT();

//private:
	//void uniform_sample(const Polygon_3& p, unsigned int samp_nb); // sample polygons

private:
	bool open;
	double *pnt_coord; // point coordinates in array form
	std::vector<unsigned int> which_group; // belong to which polygon. element is polygon id

	// handles to sample points
	//template <class Gt, class Vb> class Embedded_vertex_2;
	typedef RDT_data_structure::Vertex_handle Vertex_handle;
	std::vector<Vertex_handle> samp_pnts;
	
	// std::vector<Point_3> samp_pnts;

	/** geometry **/
	TopoPolyMesh *mesh;
	Delaunay *del;
	RestrictedVoronoiDiagram *rvd;

	

	RDT_data_structure rdt_ds;

	/** construct of RDT **/
	struct Construct_RDT_structure
	{
		Construct_RDT_structure(RDT_data_structure& _rdt_ds, std::vector<Vertex_handle>& _samp_pnts) : rdt_ds( _rdt_ds ), samp_pnts(_samp_pnts) { }
		
		void operator()(unsigned int i, unsigned int j, unsigned int k)
		{
			rdt_ds.create_face(samp_pnts[i], samp_pnts[j], samp_pnts[k]);
		}
		RDT_data_structure& rdt_ds;
		std::vector<Vertex_handle>& samp_pnts;
	};

	

	/** debug for DT draw **/
	GLuint RDT_disp_list;
	struct draw_primal_triangles
	{
		draw_primal_triangles(double *coordinates) : coord(coordinates) {}

		void operator()(unsigned int i, unsigned int j, unsigned int k) const
		{
			//std::cout<<"("<<i<<", "<<j<<", "<<k<<")\n";
			glBegin(GL_LINE_LOOP);
			glVertex3d(coord[i*3], coord[i*3+1], coord[i*3+2]);
			glVertex3d(coord[j*3], coord[j*3+1], coord[j*3+2]);
			glVertex3d(coord[k*3], coord[k*3+1], coord[k*3+2]);
			glEnd();
		}
		double *coord;
	};
};

template <class InputIterator>
void RestrictedPolygonVoronoiDiagram::insert_polygons(InputIterator first, InputIterator last, unsigned int samp_nb)
{
	assert(open);
	unsigned int group_id = 0;
	while (first != last)
	{
		insert_polygons(*first, group_id, samp_nb);
		first++;
		group_id++;
	}
}

}

#endif