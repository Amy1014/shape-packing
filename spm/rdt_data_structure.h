#pragma  once

#include <CGAL/Simple_cartesian.h>
#include <CGAL/HalfedgeDS_vertex_base.h>
#include <CGAL/HalfedgeDS_halfedge_base.h>
#include <CGAL/HalfedgeDS_face_base.h>
#include <CGAL/Polyhedron_items_3.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/HalfedgeDS_items_decorator.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>

#include "spm_cgal.h"


namespace Geex
{
	struct MyPoint
	{
		Point_3 p;
		int group_id;
		Point_3 prj_pnt;
		MyPoint() {}
		MyPoint(const Point_3& _p, const Point_3& _prj, int _group_id) : p(_p), prj_pnt(_prj), group_id(_group_id) {}
	};

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

	//////////////////////////////////////////////////////////////////////////
	//					Constrained Delaunay Triangulation					//
	//////////////////////////////////////////////////////////////////////////
	template <class Gt, class Vb = CGAL::Triangulation_vertex_base_2<Gt>>
	class CDTVertex : public Vb
	{
	public:
		typedef typename Vb::Vertex_handle	Vertex_handle;
		typedef typename Vb::Face_handle		Face_handle;
		typedef typename Vb::Point			Point;

		template < typename TDS2 >
		struct Rebind_TDS {
			typedef typename Vb::template Rebind_TDS<TDS2>::Other Vb2;
			typedef CDTVertex<Gt,Vb2> Other;
		} ;
	public:

		CDTVertex() : Vb() {}

		CDTVertex(const Point& p) : Vb(p) {}
		
		//CDTVertex(const Point& p2d, const MyPoint& _geo_info, bool _already_exist_in_triangulation) 
		//	: Vb(p2d), geo_info(_geo_info), already_exist_in_triangulation(_already_exist_in_triangulation) {}

	public:
		MyPoint geo_info;
		bool already_exist_in_rdt;
		RDT_data_structure::Vertex_handle rdt_handle;
	};

	template <class Gt, class Fb = CGAL::Constrained_triangulation_face_base_2<Gt>>
	class CDTFace : public Fb
	{

	};

	typedef CGAL::Triangulation_data_structure_2<CDTVertex<K>, CDTFace<K>> CDTDS;

	typedef CGAL::Exact_predicates_tag Itag;

	class CDT : public CGAL::Constrained_Delaunay_triangulation_2<K, CDTDS, Itag> 
	{

	};

	/************************************************************************/
	/*                        Helper classes                                */
	/************************************************************************/
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

	class CDTtoRDT : public CGAL::Modifier_base<RDT_data_structure::HalfedgeDS>
	{
		typedef RDT_data_structure::HalfedgeDS HDS;
		typedef CGAL::Polyhedron_incremental_builder_3<HDS> Internal_Builder;
		typedef RestrictedVoronoiDiagram_poly RestrictedVoronoiDiagram;
		typedef RDT_data_structure::Vertex_handle Vertex_handle;
		typedef RDT_data_structure::Halfedge_handle Halfedge_handle;
		typedef RDT_data_structure::Facet_handle	Facet_handle;
		typedef RDT_data_structure::Vertex		Vertex;
		typedef RDT_data_structure::Halfedge	Halfedge;
		typedef RDT_data_structure::Facet	Facet;
		typedef RDT_data_structure::Halfedge_around_vertex_circulator Edge_circulator;

	public:
		CDTtoRDT(CDT& _cdt, TriMesh& _trimesh) : cdt(_cdt), trimesh(_trimesh), nb_invalid_facets(0) {}

		void operator()(HDS& hds)
		{
			for (CDT::All_vertices_iterator vit = cdt.all_vertices_begin(); vit != cdt.all_vertices_end(); ++vit)
			{
				if (cdt.is_infinite(vit))
					continue;
				if (!vit->already_exist_in_rdt)
				{
					Vertex_handle rdt_v = hds.vertices_push_back(Vertex(vit->geo_info.p));
					rdt_v->mp = vit->geo_info.prj_pnt;
					vit->already_exist_in_rdt = true;
					vit->rdt_handle = rdt_v;
				}
			}
			// register already halfedges
			for (CDT::All_edges_iterator eit = cdt.all_edges_begin(); eit != cdt.all_edges_end(); ++eit)
			{
				if (cdt.is_infinite(eit))
					continue;
				CDT::Face_handle f = eit->first;
				int vi = eit->second;
				CDT::Vertex_handle v0 = f->vertex(f->cw(vi)), v1 = f->vertex(f->ccw(vi));
				Vertex_handle rdt_v0 = v0->rdt_handle, rdt_v1 = v1->rdt_handle;
				if (cdt.is_constrained(*eit)) // this edge already existed in rdt
				{
					Edge_circulator start_edge = rdt_v0->vertex_begin();
					Edge_circulator current_edge = start_edge;
					bool found = false, border_break = false;
					do 
					{
						if (current_edge->opposite()->vertex() == rdt_v1)
						{
							found = true;
							break;
						}
						if (current_edge->is_border() || current_edge->opposite()->is_border())
						{
							border_break = true;
							break;
						}
						++current_edge;
					} while (current_edge != start_edge);
					if (!found)
					{
						start_edge = current_edge;
						--current_edge;
						do 
						{
							if ( current_edge->opposite()->vertex() == rdt_v1 )
							{
								found = true;
								break;	
							}
							if (current_edge->is_border() || current_edge->opposite()->is_border())
								break;
							--current_edge;
						} while (current_edge != start_edge);
					}
					if (!found)
					{
						std::cout<<"!!! Caution: rdt data structure assertion failed!\n";
						system("pause");
					}
					existing_halfedges[std::make_pair(rdt_v1, rdt_v0)] = current_edge;
					existing_halfedges[std::make_pair(rdt_v0, rdt_v1)] = current_edge->opposite();				
				}
				else
				{
					Halfedge he;
					Halfedge ohe;
					//he.set_vertex(rdt_v0);
					//decorator.set_vertex_halfedge(rdt_v0, he);
					//ohe.set_vertex(rdt_v1);
					//decorator.set_vertex_halfedge(rdt_v1, ohe);
					Halfedge_handle h = hds.edges_push_back(he, ohe);
					existing_halfedges[std::make_pair(rdt_v1, rdt_v0)] = h;
					existing_halfedges[std::make_pair(rdt_v0, rdt_v1)] = h->opposite();
					decorator.set_vertex_halfedge(rdt_v0, h);
					decorator.set_vertex_halfedge(rdt_v1, h->opposite());
				}

			}
			for (CDT::All_faces_iterator fit = cdt.all_faces_begin(); fit != cdt.all_faces_end(); ++fit)
			{
				if (cdt.is_infinite(fit))
					continue;
				CDT::Vertex_handle v[] = { fit->vertex(0), fit->vertex(1), fit->vertex(2) };
				// check orientation
				Point_3 c = CGAL::centroid(v[0]->geo_info.prj_pnt, v[1]->geo_info.prj_pnt, v[2]->geo_info.prj_pnt);
				vec3 dv;
				trimesh.project_to_mesh(to_geex_pnt(c), dv);
				Vector_3 ndv = to_cgal_vec(dv);
				cgal_vec_normalize(ndv);
				Vector_3 v01(v[0]->geo_info.prj_pnt, v[1]->geo_info.prj_pnt), v12(v[1]->geo_info.prj_pnt, v[2]->geo_info.prj_pnt);
				cgal_vec_normalize(v01);
				cgal_vec_normalize(v12);
				Vector_3 det = CGAL::cross_product(v01, v12);
				cgal_vec_normalize(det);
				if (ndv * det >= 0.0)
				{
					Halfedge_handle he01 = existing_halfedges[std::make_pair(v[0]->rdt_handle, v[1]->rdt_handle)];
					Halfedge_handle he12 = existing_halfedges[std::make_pair(v[1]->rdt_handle, v[2]->rdt_handle)];
					Halfedge_handle he20 = existing_halfedges[std::make_pair(v[2]->rdt_handle, v[0]->rdt_handle)];
					//he01->set_prev(he20);
					//he01->set_next(he12);
					//he12->set_prev(he01);
					//he12->set_next(he20);
					//he20->set_prev(he12);
					//he20->set_next(he01);
					decorator.set_prev(he01, he20);
					decorator.set_prev(he12, he01);
					decorator.set_prev(he20, he12);
					Facet_handle f = hds.faces_push_back(Facet());
					//f->set_halfedge(he01);
					//he01->set_face(f);
					//he12->set_face(f);
					//he20->set_face(f);
					decorator.set_face_halfedge(f, he01);
					decorator.set_face(he01, f);
					decorator.set_face(he12, f);
					decorator.set_face(he20, f);
				}
				else
				{
					Halfedge_handle he21 = existing_halfedges[std::make_pair(v[2]->rdt_handle, v[1]->rdt_handle)];
					Halfedge_handle he10 = existing_halfedges[std::make_pair(v[1]->rdt_handle, v[0]->rdt_handle)];
					Halfedge_handle he02 = existing_halfedges[std::make_pair(v[0]->rdt_handle, v[2]->rdt_handle)];
					//he21->set_prev(he02);
					//he21->set_next(he10);
					//he10->set_prev(he21);
					//he10->set_next(he02);
					//he02->set_prev(he10);
					//he02->set_next(he21);
					decorator.set_prev(he21, he02);
					decorator.set_prev(he10, he21);
					decorator.set_prev(he02, he10);
					Facet_handle f = hds.faces_push_back(Facet());
					//f->set_halfedge(he21);
					//he21->set_face(f);
					//he10->set_face(f);
					//he02->set_face(f);
					decorator.set_face_halfedge(f, he21);
					decorator.set_face(he21, f);
					decorator.set_face(he10, f);
					decorator.set_face(he02, f);
				}
			}
			//std::cout<<"number of invalid facets: "<<nb_invalid_facets<<std::endl;
		}

	private:
		CDT& cdt;
		CGAL::HalfedgeDS_items_decorator<HDS> decorator;
		TriMesh& trimesh;
		std::map<std::pair<Vertex_handle, Vertex_handle>, Halfedge_handle> existing_halfedges;
		int nb_invalid_facets;
	};

}