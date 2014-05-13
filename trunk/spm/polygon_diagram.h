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
#include <CGAL/ch_graham_andrew.h>

#include "spm_cgal.h"
#include "meshes.h"
#include "rdt_data_structure.h"
 
namespace Geex
{
	//extern std::ofstream of;


class RestrictedPolygonVoronoiDiagram
{

public:
	typedef RestrictedVoronoiDiagram_poly RestrictedVoronoiDiagram;
	typedef RDT_data_structure::HalfedgeDS::Vertex_handle Vertex_handle;
	typedef RDT_data_structure::Vertex_const_handle	Vertex_const_handle;
	typedef RDT_data_structure::Halfedge_handle Halfedge_handle;
	typedef RDT_data_structure::Halfedge_const_iterator Halfedge_const_handle;
	typedef RDT_data_structure::Facet_handle Facet_handle;
	typedef RDT_data_structure::Halfedge::Vertex Vertex;
	typedef RDT_data_structure::Halfedge Halfedge;
	typedef RDT_data_structure::Face Face;
	typedef RDT_data_structure::Halfedge_iterator Halfedge_iterator;
	typedef RDT_data_structure::Edge_iterator Edge_iterator;
	typedef RDT_data_structure::Vertex_iterator Vertex_iterator;
	typedef RDT_data_structure::Facet_iterator Facet_iterator;
	typedef RDT_data_structure::Vertex_const_iterator Vertex_const_iterator;
	typedef RDT_data_structure::Halfedge_const_iterator Halfedge_const_iterator;
	typedef RDT_data_structure::Facet_const_iterator	Facet_const_iterator;
	typedef RDT_data_structure::Edge_const_iterator		Edge_const_iterator;

	typedef std::vector<Vertex_handle> VertGroup; // vertices belonging to the same group
public:
	RestrictedPolygonVoronoiDiagram();
	~RestrictedPolygonVoronoiDiagram();

	/** set functions **/
	void set_mesh(const std::string& mesh_file_name);
	void set_trimesh(TriMesh *_trimesh) { trimesh = _trimesh; }

	void begin_insert() ;

	// insert a group of polygons
	template<class InputIterator> void insert_polygons(InputIterator first, InputIterator last, unsigned int samp_nb);

	// only insert polygon vertices
	template<class InputIterator> void insert_polygons(InputIterator first, InputIterator last);

	// modify rdt through adding a cdt
	void delegate(CDTtoRDT& cdt2rdt) { rdt_ds.delegate(cdt2rdt); }

	// insert feature and bounding points, samp_nb is number of sampling points on each edge
	void insert_bounding_points(unsigned int samp_nb = 5);

	void end_insert() ;

	/** geometry **/

	void compute_clipped_VD(std::vector<Plane_3>& clipping_planes, std::vector<Point_3>& ref_pnts);

	void compute_approx_VD(std::vector<Plane_3>& clipping_planes, std::vector<Point_3>& ref_pnts);

	void compute_mixed_VD(std::vector<Plane_3>& clipping_planes, std::vector<Point_3>& ref_pnts, std::vector<bool>& use_approx);

	void compute_midpoint_VD(std::vector<Plane_3>& clipping_planes, std::vector<Point_3>& ref_pnts);

	void smooth_VD_region();

	void iDT_update(); // edge flip to preserve intrinsic Delaunay structure

	void erase_facet(Halfedge_handle h) { rdt_ds.erase_facet(h); }

	/** access functions **/
	const VertGroup& sample_points_group(unsigned int polygon_id) const { return samp_pnts[polygon_id]; }
	VertGroup& sample_points_group(unsigned int polygon_id) { return samp_pnts[polygon_id]; }
	void delete_point_group(unsigned int group_id) { samp_pnts[group_id].clear();}
	const std::vector< std::vector<Point_3> >& get_smoothed_voronoi_regions() const { return smoothed_VD_regions; }
	
	Vertex_iterator vertices_begin() { return rdt_ds.vertices_begin(); }
	Vertex_iterator vertices_end() { return rdt_ds.vertices_end(); }
	Vertex_const_iterator vertices_begin() const { return rdt_ds.vertices_begin(); }
	Vertex_const_iterator vertices_end() const { return rdt_ds.vertices_end(); }
	Edge_iterator edges_begin() { return rdt_ds.edges_begin(); }
	Edge_iterator edges_end() { return rdt_ds.edges_end(); }
	Facet_iterator faces_begin()  { return rdt_ds.facets_begin(); }
	Facet_iterator faces_end()  { return rdt_ds.facets_end(); }
	Halfedge_iterator halfedges_begin() { return rdt_ds.halfedges_begin(); }
	Halfedge_iterator halfedges_end() { return rdt_ds.halfedges_end(); }
	Halfedge_const_iterator halfedges_begin() const { return rdt_ds.halfedges_begin(); }
	Halfedge_const_iterator halfedges_end() const { return rdt_ds.halfedges_end(); }
	Facet_const_iterator faces_begin() const { return rdt_ds.facets_begin(); }
	Facet_const_iterator faces_end() const { return rdt_ds.facets_end(); }
	Edge_const_iterator edges_begin() const { return rdt_ds.edges_begin(); }
	Edge_const_iterator edges_end() const { return rdt_ds.edges_end(); }

	unsigned int number_of_groups() const { return nb_groups; }

	void get_neighbor_indices(size_t idx, std::set<size_t>& neighbors);

	// for debug
	const RDT_data_structure& get_rdt() const { return rdt_ds; }
	void save_triangulation(const std::string fn) const;

	bool is_delaunay_edge(Halfedge_handle e);

private: // private functions

	/** geometry **/
	
	inline void add_quadrilateral_edge(Halfedge_handle e, std::stack<Halfedge_handle>& q, std::set<Halfedge_handle>& s);

	// insert one polygon, return how many sample points are actually inserted
	unsigned int insert_polygons(const Polygon_3& polygon, unsigned int group_id, unsigned int samp_nb);

	// insert the vertices of a polygon only
	unsigned int insert_polygons(const Polygon_3& polygon, unsigned int group_id);

	template <class EdgeInputIterator>
	unsigned int insert_bounding_edges(EdgeInputIterator first, EdgeInputIterator last, unsigned int samp_nb);

private:
	bool open;
	// points && vertices
	std::vector<VertGroup> samp_pnts;
	std::vector<Vertex_handle> bounding_pnts;
	std::vector<MyPoint> all_points;
	std::vector< std::vector<Point_3> > smoothed_VD_regions;
	/** geometry **/
	TopoPolyMesh *mesh;
	TriMesh *trimesh;//the same as mesh, but with some other necessary functions
	unsigned int nb_groups;

	RDT_data_structure rdt_ds;

	static double pi;
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

template <class InputIterator>
void RestrictedPolygonVoronoiDiagram::insert_polygons(InputIterator first, InputIterator last)
{
	assert(open);
	unsigned int group_id = 0, nb_pnts = 0;
	samp_pnts.resize(last - first);
	while (first != last)
	{
		nb_pnts += insert_polygons(*first, group_id);
		++first;
		++group_id;
	}
	nb_groups = group_id;
}

template <class EdgeInputIterator>
unsigned int RestrictedPolygonVoronoiDiagram::insert_bounding_edges(EdgeInputIterator first, EdgeInputIterator last, unsigned int samp_nb)
{
	assert(open);
	std::set<int> vert_idx; 
	unsigned int cnt = 0;
	while (first != last)
	{
		int vi0 = first->first, vi1 = first->second;
		const vec3& v0 = trimesh->vertex(vi0).point();
		const vec3& v1 = trimesh->vertex(vi1).point();
// 		vec3 v0 = 0.98*trimesh->vertex(vi0).point();
// 		vec3 v1 = 0.98*trimesh->vertex(vi1).point();
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