#pragma  once

#include <CGAL/Simple_cartesian.h>
#include <CGAL/HalfedgeDS_vertex_base.h>
#include <CGAL/HalfedgeDS_halfedge_base.h>
#include <CGAL/HalfedgeDS_face_base.h>
#include <CGAL/Polyhedron_items_3.h>
#include <CGAL/Polyhedron_3.h>

#include "spm_cgal.h"


namespace Geex
{
	/** representation of restricted Delaunay triangulation **/
	template <class Refs, class Point>
	struct MyVertex
		: public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point>
	{
		//CGAL::Color color;
		MyVertex() {} // repeat the required constructors
		MyVertex( const Point& p)
			: CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point>(p) {}
		MyVertex (const Point& p, int id) 
			: CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point>(p), group_id(id) {}

		int group_id;
		Point_3 mp; // projection point on mesh
		std::vector<Point_3> vd_vertices;
		bool contain_non_delaunay_facet;
		//std::vector<double> weights; // weights for optimization
		int idx;
	};

	template <class Refs>
	struct MyHalfEdge : public CGAL::HalfedgeDS_halfedge_base<Refs, CGAL::Tag_true, CGAL::Tag_true, CGAL::Tag_true> 
	{
	//public:
	//	bool visited;
	};

	template <class Refs>
	class MyFace : public CGAL::HalfedgeDS_face_base<Refs>
	{
	public:
		typedef CGAL::HalfedgeDS_face_base<Refs> baseclass;
	public:
		Vector_3 n;
		bool vacant;
		bool is_delaunay;
	};

	class MyItems : public CGAL::Polyhedron_items_3
	{
	public:
		template <class Refs, class Traits>
		struct Vertex_wrapper 
		{
			typedef typename Traits::Point_3 Point_3;
			typedef MyVertex<Refs, Point_3> Vertex;
		};
		template <class Refs, class Traits>
		struct Halfedge_wrapper
		{
			//typedef MyHalfEdge<Refs> HalfEdge;
			typedef CGAL::HalfedgeDS_halfedge_base<Refs,CGAL::Tag_true, CGAL::Tag_true, CGAL::Tag_true>     Halfedge;
		};
		template <class Refs, class Traits>
		struct Face_wrapper 
		{
			typedef MyFace<Refs> Face;
		};
	};


	class RDT_data_structure : public CGAL::Polyhedron_3<K, MyItems>
	{
	
	};
}