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
#include <algorithm>
#include <numeric>

#include <Geex/cvt/delaunay.h>
#include <Geex/cvt/RVD.h>

#include <CGAL/Object.h>
#include "spm_cgal.h"
#include "meshes.h"
#include "rdt_data_structure.h"
 
namespace Geex
{

class RestrictedPolygonVoronoiDiagram
{

public:
	typedef RestrictedVoronoiDiagram_poly RestrictedVoronoiDiagram;
	typedef RDT_data_structure::Vertex_handle Vertex_handle;
	typedef RDT_data_structure::Face_handle Face_handle;
	typedef std::pair<Vertex_handle, Vertex_handle> Vertex_pair;
	typedef std::vector<Vertex_handle> VertGroup; // vertices belonging to the same group
	typedef RDT_data_structure::Vertex Vertex;
	typedef RDT_data_structure::Edge Edge;
	typedef RDT_data_structure::Edge_iterator Edge_iterator;
	typedef RDT_data_structure::Vertex_iterator Vertex_iterator;
	typedef RDT_data_structure::Face_iterator Face_iterator;

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
	
	Vertex_iterator vertices_begin() const { return rdt_ds.vertices_begin(); }
	Vertex_iterator vertices_end() const { return rdt_ds.vertices_end(); }
	Edge_iterator edges_begin() const { return rdt_ds.edges_begin(); }
	Edge_iterator edges_end() const { return rdt_ds.edges_end(); }
	Face_iterator faces_begin() const { return rdt_ds.faces_begin(); }
	Face_iterator faces_end() const { return rdt_ds.faces_end(); }
	unsigned int number_of_groups() const { return nb_groups; }

	// for debug
	const RDT_data_structure& get_rdt() const { return rdt_ds; }


private: // private functions

	/** geometry **/
	bool is_delaunay_edge(const Edge& e);
	inline bool is_interior_edge(const Edge& e);
	void add_quadrilateral_edge(Face_handle f, int i, std::queue<Edge>& q);
	void add_quadrilateral_edge(Vertex_handle v0, Vertex_handle v1, std::queue<Vertex_pair>&, std::set<Vertex_pair>& s);
	/** computation **/
	inline Plane_3 get_bisector(Vertex_handle v0, Vertex_handle v1)
	{
		//return CGAL::bisector(v0->point_3(), v1->point_3());
		return CGAL::bisector(v0->mp, v1->mp);
	}

	// insert one polygon, return how many sample points are actually inserted
	unsigned int insert_polygons(const Polygon_3& polygon, unsigned int group_id, unsigned int samp_nb = 24);

	template <class EdgeInputIterator>
	unsigned int insert_bounding_edges(EdgeInputIterator first, EdgeInputIterator last, unsigned int samp_nb);

private:
	bool open;
	//double *pnt_coord; // point coordinates in array form

	// handles to sample points
	std::vector<VertGroup> samp_pnts;
	std::vector<Vertex_handle> bounding_pnts;

	/** geometry **/
	TopoPolyMesh *mesh;
	TriMesh *trimesh;//the same as mesh, but with some other necessary functions
	//Delaunay *del;
	//RestrictedVoronoiDiagram *rvd;
	//std::vector<std::vector<Point_3>> clipped_VD;//clipped Voronoi diagram for each polygon, expressed as a set of lines
	unsigned int nb_groups;

	RDT_data_structure rdt_ds;
	std::map<Vertex_pair, Face_handle> edge_face_adjacency;
	std::map<Vertex_pair, int> edge_valance;

	/** helper classes **/
	// construct of RDT
	struct Construct_RDT_structure
	{
		Construct_RDT_structure(RDT_data_structure& _rdt_ds, std::vector<Vertex_handle>& _samp_pnts,
						std::map<Vertex_pair, Face_handle>& ef_adj, std::map<Vertex_pair, int>& _edge_valance) 
							: rdt_ds( _rdt_ds ), samp_pnts(_samp_pnts), edge_face_adjacency(ef_adj), edge_valance(_edge_valance) { }
		
		void operator()(unsigned int i, unsigned int j, unsigned int k) const ;

		void find_adjacency(const Vertex_handle& v0, const Vertex_handle& v1, int v2, const Face_handle& f) const;

		RDT_data_structure& rdt_ds;
		std::vector<Vertex_handle>& samp_pnts;
		std::map<Vertex_pair, Face_handle>& edge_face_adjacency;
		std::map<Vertex_pair, int>& edge_valance;
	};
};

template <class InputIterator>
void RestrictedPolygonVoronoiDiagram::insert_polygons(InputIterator first, InputIterator last, unsigned int samp_nb)
{
	assert(open);
	unsigned int group_id = 0, nb_pnts = 0;
	samp_pnts.reserve(last - first);
	//cents.reserve(last - first);
	while (first != last)
	{
		nb_pnts += insert_polygons(*first, group_id, samp_nb);
		++first;
		group_id++;
	}
	nb_groups = group_id;
	assert(nb_pnts == rdt_ds.number_of_vertices());
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
			Vertex_handle vh = rdt_ds.create_vertex(RDT_data_structure::Vertex(to_cgal_pnt(v0), -1));
			bounding_pnts.push_back(vh);
			vh->mp = to_cgal_pnt(v0);
			vert_idx.insert(vi0);
		}
		if ( vert_idx.find(vi1) == vert_idx.end() )
		{
			Vertex_handle vh = rdt_ds.create_vertex(RDT_data_structure::Vertex(to_cgal_pnt(v1), -1));
			bounding_pnts.push_back(vh);
			vh->mp = to_cgal_pnt(v1);
			vert_idx.insert(vi1);
		}
		for ( unsigned int i = 1; i < samp_nb; i++ )
		{
			vec3 sv = ( (samp_nb - i)*v0 + i*v1 ) / samp_nb;
			Vertex_handle vh = rdt_ds.create_vertex(RDT_data_structure::Vertex(to_cgal_pnt(sv), -1));
			bounding_pnts.push_back(vh);
			vh->mp = to_cgal_pnt(sv);
			cnt++;
		}
		
		++first;
	}
	return cnt;
}

}

#endif