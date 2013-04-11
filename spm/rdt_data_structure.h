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
		MyVertex() {} // repeat the required constructors
		MyVertex( const Point& p)
			: CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point>(p) {}
		MyVertex (const Point& p, int id) 
			: CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point>(p), group_id(id) {}

		int group_id;
		Point_3 mp; // projection point on mesh
		std::vector<Point_3> vd_vertices;
		bool contain_non_delaunay_facet;
		bool penetration;
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

	/************************************************************************/
	/*                        Helper classes                                */
	/************************************************************************/
	struct MyPoint
	{
		Point_3 p;
		int group_id;
		Point_3 prj_pnt;
		MyPoint(const Point_3& _p, const Point_3& _prj, int _group_id) : p(_p), prj_pnt(_prj), group_id(_group_id) {}
	};
	class Builder : public CGAL::Modifier_base<RDT_data_structure::HalfedgeDS>
	{
		typedef RDT_data_structure::HalfedgeDS HDS;
		typedef CGAL::Polyhedron_incremental_builder_3<HDS> Internal_Builder;
		typedef RestrictedVoronoiDiagram_poly RestrictedVoronoiDiagram;
		typedef HDS::Vertex_handle Vertex_handle;
		typedef HDS::Vertex Vertex;

	public:
		Builder(RestrictedVoronoiDiagram& _rvd, std::vector<MyPoint>& verts, TriMesh& _m)
			:rvd(_rvd), points(verts), m(_m), nb_invalid_facets(0) {}
		void operator()(HDS& hds)
		{
			Internal_Builder builder(hds, true);
			builder.begin_surface(points.size(), points.size()/3, CGAL::Polyhedron_incremental_builder_3<HDS>::ABSOLUTE_INDEXING);

			for (unsigned int i = 0; i < points.size(); i++)
			{
				Vertex_handle vh = builder.add_vertex(points[i].p);
				vh->group_id = points[i].group_id;
				vh->mp = points[i].prj_pnt;
				vh->idx = i;
				//of<<"v "<<points[i].prj_pnt.x()<<' '<<points[i].prj_pnt.y()<<' '<<points[i].prj_pnt.z()<<std::endl;
			}
			rvd.for_each_primal_triangle(construct_faces(builder, points, m, nb_invalid_facets));
			builder.end_surface();

			std::cout<<"Number of invalid facets: "<<nb_invalid_facets<<std::endl;
		}
	private:
		RestrictedVoronoiDiagram& rvd;
		std::vector<MyPoint>& points;
		TriMesh& m;
		int nb_invalid_facets;

		struct construct_faces
		{
			Internal_Builder& b;
			std::vector<MyPoint>& points;
			TriMesh& m;
			int& nb_invalid_facets;

			construct_faces(Internal_Builder& _b, std::vector<MyPoint>& _points, TriMesh& _m, int& _nb_invalid_facets)
				: b(_b), points(_points), m(_m), nb_invalid_facets(_nb_invalid_facets) {}
			void operator()(unsigned int i, unsigned int j, unsigned int k) const
			{
				vec3 dv;
				Vector_3 avgp = (points[i].prj_pnt - CGAL::ORIGIN)+(points[j].prj_pnt - CGAL::ORIGIN)+(points[k].prj_pnt - CGAL::ORIGIN); 
				avgp = avgp / 3.0;
				m.project_to_mesh(to_geex_pnt(CGAL::ORIGIN + avgp), dv);
				Vector_3 n = to_cgal_vec(dv);
				n = n / CGAL::sqrt(n.squared_length());
				Vector_3 nij(points[i].prj_pnt, points[j].prj_pnt);
				Vector_3 njk(points[j].prj_pnt, points[k].prj_pnt);
				cgal_vec_normalize(nij);
				cgal_vec_normalize(njk);
				Vector_3 det = CGAL::cross_product(nij, njk);
				cgal_vec_normalize(det);
				if ( n*det < 0.0 )
				{
					unsigned int indices[] = {k, j, i};
					if (b.test_facet(indices, indices+3))
					{
						HDS::Face_handle f = b.begin_facet();
						b.add_vertex_to_facet(k);
						b.add_vertex_to_facet(j);
						b.add_vertex_to_facet(i);
						b.end_facet();
						f->is_delaunay = true;
					}
					else
						nb_invalid_facets++;

				}
				else
				{	
					unsigned int indices[] = {i, j, k};
					if (b.test_facet(indices, indices+3))
					{
						HDS::Face_handle f = b.begin_facet();
						b.add_vertex_to_facet(i);
						b.add_vertex_to_facet(j);
						b.add_vertex_to_facet(k);
						b.end_facet();
						f->is_delaunay = true;
					}
					else
						nb_invalid_facets++;
				}
				//of<<"f "<<i+1<<' '<<j+1<<' '<<k+1<<std::endl;
			}
		};

	};

	class Append : public CGAL::Modifier_base<RDT_data_structure::HalfedgeDS>
	{
		typedef RDT_data_structure::HalfedgeDS HDS;
		typedef CGAL::Polyhedron_incremental_builder_3<HDS> Internal_Builder;
		typedef RestrictedVoronoiDiagram_poly RestrictedVoronoiDiagram;
		typedef HDS::Vertex_handle Vertex_handle;
		typedef HDS::Vertex Vertex;

	public:
		Append(std::vector<MyPoint>& pnts, std::vector<CGAL::Triple<int, int, int>>& adjacency, TriMesh& trimesh) 
			: appended_points(pnts), adjacency_info(adjacency), mesh(trimesh), nb_invalid_facets(0) {}

		void operator()(HDS& hds)
		{
			Internal_Builder builder(hds, true);
			builder.begin_surface(appended_points.size(), appended_points.size()/3, CGAL::Polyhedron_incremental_builder_3<HDS>::RELATIVE_INDEXING);
			for (unsigned int i = 0; i < appended_points.size(); i++)
			{
				Vertex_handle vh = builder.add_vertex(appended_points[i].p);
				vh->group_id = appended_points[i].group_id;
				vh->mp = appended_points[i].prj_pnt;
			}

			for (unsigned int i = 0; i < adjacency_info.size(); i++)
			{
				vec3 dv;
				const Point_3& p0 = appended_points[adjacency_info[i].first].prj_pnt;
				const Point_3& p1 = appended_points[adjacency_info[i].second].prj_pnt;
				const Point_3& p2 = appended_points[adjacency_info[i].third].prj_pnt;
				Point_3 cent = CGAL::centroid(p0, p1, p2);
				mesh.project_to_mesh(to_geex_pnt(cent), dv);
				Vector_3 n = to_cgal_vec(dv);
				cgal_vec_normalize(n);
				Vector_3 v01(p0, p1), v12(p1, p2);
				cgal_vec_normalize(v01);
				cgal_vec_normalize(v12);
				Vector_3 det = CGAL::cross_product(v01, v12);
				unsigned int indices[3];
				if (n*det < 0.0)
				{
					indices[0] = adjacency_info[i].third;
					indices[1] = adjacency_info[i].second;
					indices[2] = adjacency_info[i].first;
				}
				else
				{
					indices[0] = adjacency_info[i].first;;
					indices[1] = adjacency_info[i].second;
					indices[2] = adjacency_info[i].third;
				}
				if (builder.test_facet(indices, indices+3))
				{
					HDS::Face_handle f = builder.begin_facet();
					builder.add_vertex_to_facet(indices[0]);
					builder.add_vertex_to_facet(indices[1]);
					builder.add_vertex_to_facet(indices[2]);
					builder.end_facet();
					f->is_delaunay = true;
				}
				else
					nb_invalid_facets++;
			}
			builder.end_surface();
			std::cout<<"number of invalid facets: "<<nb_invalid_facets<<std::endl;
		}

	private:
		std::vector<MyPoint>& appended_points;
		std::vector<CGAL::Triple<int, int, int>>& adjacency_info;
		TriMesh& mesh;
		int nb_invalid_facets;
	};
}