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
#include <map>
#include <set>
#include <queue>
#include <stack>
#include <algorithm>
#include <numeric>
#include <fstream>

#include <Geex/cvt/delaunay.h>
#include <Geex/cvt/RVD.h>

#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Object.h>
#include "spm_cgal.h"
#include "meshes.h"
#include "rdt_data_structure.h"
 
namespace Geex
{

	extern int nb_invalid_edges;
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
						:rvd(_rvd), points(verts), m(_m) {}
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
			}
			rvd.for_each_primal_triangle(construct_faces(builder, points, m));
			builder.end_surface();

			std::cout<<"Number of invalid facets: "<<nb_invalid_edges<<std::endl;
		}
	private:
		RestrictedVoronoiDiagram& rvd;
		std::vector<MyPoint>& points;
		TriMesh& m;
		//std::set<std::pair<unsigned int, unsigned int>> ordering;

		struct construct_faces
		{
			Internal_Builder& b;
			std::vector<MyPoint>& points;
			TriMesh& m;
			//std::set<std::pair<unsigned int, unsigned int>>& ordering;
			construct_faces(Internal_Builder& _b, std::vector<MyPoint>& _points, TriMesh& _m)
				: b(_b), points(_points), m(_m) {}
			void operator()(unsigned int i, unsigned int j, unsigned int k) const
			{
				vec3 dv;
				Vector_3 avgp = (points[i].prj_pnt - CGAL::ORIGIN)+(points[j].prj_pnt - CGAL::ORIGIN)+(points[k].prj_pnt - CGAL::ORIGIN); 
				avgp = avgp / 3.0;
				m.project_to_mesh(to_geex_pnt(CGAL::ORIGIN + avgp), dv);
				Vector_3 n = to_cgal_vec(dv);
				Vector_3 nij(points[i].prj_pnt, points[j].prj_pnt);
				Vector_3 njk(points[j].prj_pnt, points[k].prj_pnt);
				Vector_3 det = CGAL::cross_product(nij, njk);
				if ( n*det < 0.0 )
				{
					unsigned int indices[] = {k, j, i};
					if(!b.test_facet(indices, indices+3))
						nb_invalid_edges++;
					else
					{
						b.begin_facet();
						b.add_vertex_to_facet(k);
						b.add_vertex_to_facet(j);
						b.add_vertex_to_facet(i);
						b.end_facet();
					}

				}
				else
				{	
					unsigned int indices[] = {i, j, k};
					if(!b.test_facet(indices, indices+3))
						nb_invalid_edges++;
					else
					{
						b.begin_facet();
						b.add_vertex_to_facet(i);
						b.add_vertex_to_facet(j);
						b.add_vertex_to_facet(k);
						b.end_facet();
					}

				}
			}
		};
	};
class RestrictedPolygonVoronoiDiagram
{

public:
	typedef RestrictedVoronoiDiagram_poly RestrictedVoronoiDiagram;
	typedef RDT_data_structure::HalfedgeDS::Vertex_handle Vertex_handle;
	//typedef RDT_data_structure::Vertex_handle Vertex_handle;
	typedef RDT_data_structure::Halfedge_handle Halfedge_handle;
	typedef RDT_data_structure::Facet_handle Facet_handle;
	//typedef std::pair<Vertex_handle, Vertex_handle> Vertex_pair;
	typedef std::vector<Vertex_handle> VertGroup; // vertices belonging to the same group
	typedef RDT_data_structure::Halfedge::Vertex Vertex;
	//typedef RDT_data_structure::Vertex Vertex;
	typedef RDT_data_structure::Halfedge Halfedge;
	typedef RDT_data_structure::Face Face;
	//typedef RDT_data_structure::Edge Edge;
	typedef RDT_data_structure::Halfedge_iterator Halfedge_iterator;
	typedef RDT_data_structure::Edge_iterator Edge_iterator;
	typedef RDT_data_structure::Vertex_iterator Vertex_iterator;
	typedef RDT_data_structure::Facet_iterator Facet_iterator;

public:
	RestrictedPolygonVoronoiDiagram();
	~RestrictedPolygonVoronoiDiagram();

	/** set functions **/
	void set_mesh(const std::string& mesh_file_name);
	void set_trimesh(TriMesh *_trimesh) { assert(!trimesh); trimesh = _trimesh; }

	void begin_insert() ;

	// insert a group of polygons
	template<class InputIterator> void insert_polygons(InputIterator first, InputIterator last, unsigned int samp_nb = 24);

	// insert feature and bounding points, samp_nb is number of sampling points on each edge
	void insert_bounding_points(unsigned int samp_nb = 5);

	void end_insert() ;

	/** geometry **/

	void compute_clipped_VD(std::vector<Plane_3>& clipping_planes);

	void iDT_update(); // edge flip to preserve intrinsic Delaunay structure

	/** access functions **/
	const VertGroup& sample_points_group(unsigned int polygon_id) const { return samp_pnts[polygon_id]; }
	VertGroup& sample_points_group(unsigned int polygon_id) { return samp_pnts[polygon_id]; }
	
	Vertex_iterator vertices_begin() { return rdt_ds.vertices_begin(); }
	Vertex_iterator vertices_end() { return rdt_ds.vertices_end(); }
	Edge_iterator edges_begin() { return rdt_ds.edges_begin(); }
	Edge_iterator edges_end() { return rdt_ds.edges_end(); }
	Facet_iterator faces_begin()  { return rdt_ds.facets_begin(); }
	Facet_iterator faces_end()  { return rdt_ds.facets_end(); }
	Halfedge_iterator halfedges_begin() { return rdt_ds.halfedges_begin(); }
	Halfedge_iterator halfedges_end() { return rdt_ds.halfedges_end(); }
	unsigned int number_of_groups() const { return nb_groups; }

	// for debug
	const RDT_data_structure& get_rdt() const { return rdt_ds; }
	void save_triangulation(const std::string fn);

	bool is_delaunay_edge(Halfedge_handle e);

private: // private functions

	/** geometry **/
	
	//inline bool is_interior_edge(const Edge& e);
	//inline void add_quadrilateral_edge(Halfedge_handle e, std::queue<Halfedge_handle>& q, std::set<Halfedge_handle>& s);
	inline void add_quadrilateral_edge(Halfedge_handle e, std::stack<Halfedge_handle>& q, std::set<Halfedge_handle>& s);
	void add_quadrilateral_edge(std::pair<Vertex_handle, Vertex_handle>, 
						std::stack<std::pair<Vertex_handle, Vertex_handle>>&, 
						std::set<std::pair<Vertex_handle, Vertex_handle>>& s);

	// insert one polygon, return how many sample points are actually inserted
	unsigned int insert_polygons(const Polygon_3& polygon, unsigned int group_id, unsigned int samp_nb = 24);

	template <class EdgeInputIterator>
	unsigned int insert_bounding_edges(EdgeInputIterator first, EdgeInputIterator last, unsigned int samp_nb);

private:
	bool open;
	// points && vertices
	std::vector<VertGroup> samp_pnts;
	std::vector<Vertex_handle> bounding_pnts;
	std::vector<MyPoint> all_points;
	/** geometry **/
	TopoPolyMesh *mesh;
	TriMesh *trimesh;//the same as mesh, but with some other necessary functions
	unsigned int nb_groups;

	RDT_data_structure rdt_ds;
};

template <class InputIterator>
void RestrictedPolygonVoronoiDiagram::insert_polygons(InputIterator first, InputIterator last, unsigned int samp_nb)
{
	assert(open);
	unsigned int group_id = 0, nb_pnts = 0;
	samp_pnts.resize(last - first);
	//unsigned int approx_nb_vertices = (last-first)*samp_nb;
	//cents.reserve(last - first);
	while (first != last)
	{
		nb_pnts += insert_polygons(*first, group_id, samp_nb);
		++first;
		group_id++;
	}
	nb_groups = group_id;
	//assert(nb_pnts == rdt_ds.size_of_vertices());
}

template <class EdgeInputIterator>
unsigned int RestrictedPolygonVoronoiDiagram::insert_bounding_edges(EdgeInputIterator first, EdgeInputIterator last, unsigned int samp_nb)
{
	assert(open);
	std::set<int> vert_idx; 
	//bounding_pnts.reserve(last - first);
	unsigned int cnt = 0;
	while (first != last)
	{
		int vi0 = first->first, vi1 = first->second;
		const vec3& v0 = trimesh->vertex(vi0).point();
		const vec3& v1 = trimesh->vertex(vi1).point();
		if ( vert_idx.find(vi0) == vert_idx.end() )
		{
			Point_3 p = to_cgal_pnt(v0);
			all_points.push_back(MyPoint(p, p, -1));
			vert_idx.insert(vi0);
			cnt++;
		}
		if ( vert_idx.find(vi1) == vert_idx.end() )
		{
			Point_3 p = to_cgal_pnt(v1);
			all_points.push_back(MyPoint(p, p, -1));
			vert_idx.insert(vi1);
			cnt++;
		}
		for ( unsigned int i = 1; i < samp_nb; i++ )
		{
			vec3 sv = ( (samp_nb - i)*v0 + i*v1 ) / samp_nb;
			Point_3 p = to_cgal_pnt(sv);
			all_points.push_back(MyPoint(p, p, -1));
			cnt++;
		}	
		++first;
	}
	return cnt;
}

}

#endif