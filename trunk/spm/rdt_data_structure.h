#pragma  once

#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include "spm_cgal.h"


namespace Geex
{
	/** representation of restricted Delaunay triangulation **/
	struct Embedded_point_2
	{
		//typedef K::Point_3 Point_2;
		Embedded_point_2() {}
		Embedded_point_2(const Point_3& _p) : p(_p) {}
		Point_3 p;
		Point_3 mp; // projection point on mesh
	};


	template < class Gt, class Vb = CGAL::Triangulation_vertex_base_2<Gt> >
	class Embedded_vertex_2 : public Vb, public Embedded_point_2 // vertex of 2d triangulation embedding in 3d space
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
		Embedded_vertex_2() : Vb(), group_id(-1) { }
		Embedded_vertex_2(const Point_3& _p) : Embedded_point_2(_p), group_id(-1), n(0.0, 0.0, 0.0) {}
		Embedded_vertex_2(const Point_3& _p, int _group_id) : Embedded_point_2(_p), group_id(_group_id), n(0.0, 0.0, 0.0) {}
		inline Point_3 point_3() const { return Embedded_point_2::p; }

	public:
		int group_id;
		Vector_3 n;
		std::vector<Point_3> vd_vertices; // vertices of Voronoi diagram associated with this vertex
	};

	template < class Gt, class Fb = CGAL::Triangulation_face_base_2<Gt> >
	class Embedded_face_2 : public Fb
	{
	public:
		typedef typename Fb::Vertex_handle Vertex_handle;
		typedef typename Fb::Face_handle Face_handle;

		template < typename TDS2 >
		struct Rebind_TDS {
			typedef typename Fb::template Rebind_TDS<TDS2>::Other Fb2;
			typedef Embedded_face_2<Gt,Fb2> Other;
		};
		Embedded_face_2() : Fb() {  }
		Embedded_face_2(
			Vertex_handle v0, Vertex_handle v1, Vertex_handle v2
			) : Fb(v0, v1, v2) { }

		Embedded_face_2(
			Vertex_handle v0, Vertex_handle v1, Vertex_handle v2,
			Face_handle n0, Face_handle n1, Face_handle n2 
			) : Fb(v0,v1,v2,n0,n1,n2) { }

	public:
		bool visited[3]; // edge flipping, whether an edge has been visited
		Vector_3 n;
	};

	typedef CGAL::Triangulation_data_structure_2<Embedded_vertex_2<K>, Embedded_face_2<K>> RDT_data_structure_base;
	class RDT_data_structure: public RDT_data_structure_base// representation of RDT
	{
		typedef RDT_data_structure_base baseclass;
	public:
		typedef baseclass::Vertex_handle Vertex_handle;
		typedef baseclass::Face_handle Face_handle;
		typedef baseclass::Vertex	Vertex;
		typedef baseclass::Face		Face;
		typedef baseclass::Edge		Edge;
		typedef baseclass::Vertex_iterator	Vertex_iterator;
		typedef baseclass::Edge_iterator	Edge_iterator;
		typedef baseclass::Face_iterator	Face_iterator;
	public:
		inline bool is_infinite(Vertex_handle v) const
		{
			return v == _infinite_vertex;
		}
		inline bool is_infinite(Face_handle f) const
		{
			return is_infinite(f->vertex(0)) || is_infinite(f->vertex(1)) || is_infinite(f->vertex(2));
		}
		inline bool is_infinite(const Edge& e) const
		{
			Vertex_handle v0 = e.first->vertex(e.first->ccw(e.second));
			Vertex_handle v1 = e.first->vertex(e.first->cw(e.second));
			return is_infinite(v0) || is_infinite(v1);
		}
		inline size_t number_of_vertices() const { return baseclass::number_of_vertices()-1;}

		inline Vertex_handle infinite_vertex() const { return _infinite_vertex; }

		inline void clear()
		{
			baseclass::clear();
			_infinite_vertex = baseclass::create_vertex(Vertex());
		}

		RDT_data_structure() : baseclass()
		{
			_infinite_vertex = create_vertex(Vertex());
		}
	private:
		Vertex_handle _infinite_vertex;
	};
}