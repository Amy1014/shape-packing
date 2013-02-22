#ifndef _SPM_CGAL_
#define _SPM_CGAL_

#include <vector>

#include <CGAL/Polygon_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAl/Regular_triangulation_filtered_traits_3.h>
#include <CGAL/Triangulation_hierarchy_vertex_base_3.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Extended_cartesian.h>
#include <CGAL/Object.h>
#include <CGAL/Triangle_3.h>
/** for exact computation */
#include <CGAL/MP_Float.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <Geex/CVT/geometry.h>


namespace Geex {
	/* basic kernels */
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef K::Point_2		Point_2;
	typedef K::Vector_2		Vector_2;
	typedef K::Line_2		Line_2;
	typedef K::Segment_2	Segment_2;
	typedef CGAL::Polygon_2<K>	Polygon_2;
	typedef CGAL::Aff_transformation_2<K>	Transformation_2;
	typedef CGAL::Triangle_2<K>				Triangle_2;
	typedef K::Point_3    Point_3;
	typedef K::Vector_3   Vector_3;
	typedef K::Segment_3	Segment_3;
	typedef K::Ray_3		Ray_3;
	typedef K::Direction_3	Direction_3;
	typedef CGAL::Plane_3<K>	Plane_3;
	typedef CGAL::Object	Object;
	typedef CGAL::Aff_transformation_3<K>	Transformation_3;
	typedef CGAL::Iso_rectangle_2<K>	Bounding_box_2;
	/** For exact computation **/
	typedef CGAL::MP_Float	MP_Float;
	typedef CGAL::Lazy_exact_nt<CGAL::Quotient<MP_Float>>	NT;
	typedef CGAL::Cartesian<NT> Ext_kernel;
	typedef Ext_kernel::Point_3	Ext_point_3;
	typedef CGAL::Vector_3<Ext_kernel> Ext_vector_3;
	typedef CGAL::Segment_3<Ext_kernel> Ext_segment_3;
	typedef CGAL::Ray_3<Ext_kernel>	Ext_ray_3;
	typedef CGAL::Plane_3<Ext_kernel> Ext_plane_3;
	typedef CGAL::Triangle_3<Ext_kernel> Ext_triangle_3;
	//using std::std::vector;

	// ------------- My stuff to extend CGAL --------------------------
	//typedef CGAL::Nef_polyhedron_3<CGAL::Extended_cartesian<CGAL::Lazy_exact_nt<CGAL::Gmpq>>> Nef_polyhedron_3;

	//class Vertex;

	template <class Vertex_base>
	class VertexDecoration {	
	public:
		VertexDecoration() {
			index = -1 ;
			radius = 0.0;
			group_id = -1;
		}
	public:
		double radius;  // the radius of sphere, weight() = radius^2
		int index ;
		//bool dual_intersects_boundary ;
		int group_id;   
		//Nef_polyhedron_3 vCell;//associated voronoi cell
		std::vector<std::vector<Object>> vCell;
		std::vector<typename Vertex_base::Vertex_handle> dualVertHandle;
		std::vector<int> dualVert;
		std::vector<Segment_3> tangentIntersection;//intersection with the tangent plane of the faces in the vCell, correspondingly
		//std::vector<int> tanIntersecDualVert;//correspond to the above tangentIntersection, recording whether the intersection segment is formed with which vertex
		bool included;//for hole detection with spheres
	} ;

	class CellDecoration {
	public:
		std::vector<bool> visited;//modified for openmp parallelization
	} ;


	// ------------- CGAL stuff ---------------------------------

	// ----------------------- A CGAL::Vertex with decoration ------------------
	template < class Gt, class Vb = CGAL::Triangulation_vertex_base_3<Gt> >
	class Vertex : public  Vb, public VertexDecoration<Vb>  {
		typedef Vb superclass;
	public:
		typedef typename Vb::Vertex_handle      Vertex_handle;
		typedef typename Vb::Cell_handle		Cell_handle;
		typedef typename Vb::Point				Point;

		template < typename TDS2 >
		struct Rebind_TDS {
			typedef typename Vb::template Rebind_TDS<TDS2>::Other Vb2;
			typedef Vertex<Gt,Vb2> Other;
		} ;
	public:
		Vertex() : superclass() {}
		Vertex(const Point & p) : superclass(p) {}
		Vertex(const Point & p, Cell_handle f) : superclass(f,p) {}
		Vertex(Cell_handle f) : superclass(f) {}
	} ;


	// ----------------------- A CGAL::Cell with decoration ------------------
	template < class Gt, class Cb = CGAL::Triangulation_cell_base_3<Gt> >
	class Cell : public Cb, public CellDecoration {
		typedef Cb superclass;
	public:
		typedef typename Cb::Vertex_handle      Vertex_handle;
		typedef typename Cb::Cell_handle        Cell_handle;
		template < typename TDS2 >
		struct Rebind_TDS {
			typedef typename Cb::template Rebind_TDS<TDS2>::Other Cb2;
			typedef Cell<Gt,Cb2> Other;
		};


		Cell() : superclass() {  }
		Cell(
			Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3
			) : superclass(v0, v1, v2, v3) { }

		Cell(
			Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3,
			Cell_handle n0, Cell_handle n1, Cell_handle n2 , Cell_handle n3
			) : superclass(v0,v1,v2,v3,n0,n1,n2,n3) { }

	} ;



	typedef CGAL::Regular_triangulation_filtered_traits_3<K>  Traits;
	typedef Vertex<Traits> Vb ;
	typedef CGAL::Triangulation_hierarchy_vertex_base_3<Vb> Vbh;
	typedef Cell<Traits>   Cb ;
	typedef CGAL::Triangulation_data_structure_3<Vbh,Cb> TDS;
	typedef CGAL::Regular_triangulation_3<Traits, TDS> Regular_triangulation_3;
	typedef CGAL::Delaunay_triangulation_2<K> DT;

	typedef Regular_triangulation_3::Weighted_point Weighted_point;

	inline vec3 to_geex(const Point_3& p) { return vec3(p.x(), p.y(), p.z()) ; }
	inline Point_3 to_cgal(const vec3& p) { return Point_3(p.x, p.y, p.z) ; }
	inline Point_2 to_uv(const Vector_3& uvx, const Vector_3& uvy, const Point_3& o, const Point_3& p)
	{
		Vector_3 v(o, p);
		return Point_2(v*uvx, v*uvy);
	}
	inline Point_3 to_xy(const Vector_3& uvx, const Vector_3& uvy, const Point_3& o, const Point_2& p)
	{
		return o + p.x()*uvx + p.y()*uvy;
	}

	class RegularBase_3 : public Regular_triangulation_3  {
	public:
	} ;

//______________________________________________________________________________________

#define FOR_EACH_VERTEX(TCLASS, INSTANCE, IT)                       \
	for(                                                            \
	TCLASS::Finite_vertices_iterator  IT = (INSTANCE)->finite_vertices_begin() ; \
	IT != (INSTANCE)->finite_vertices_end(); IT++                      \
	)

#define FOR_EACH_FACE(TCLASS, INSTANCE, IT)                      \
	for(                                                         \
	TCLASS::Face_iterator IT = (INSTANCE)->faces_begin() ;   \
	IT != (INSTANCE)->faces_end(); IT++                      \
	)

}

#endif