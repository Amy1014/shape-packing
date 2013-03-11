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
#include <algorithm>
#include <numeric>
//#include <LpCVT/combinatorics/delaunay.h>
#include <Geex/cvt/delaunay.h>
#include <Geex/cvt/RVD.h>
//#include <LpCVT/combinatorics/RVD.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Object.h>
#include "spm_cgal.h"
#include "meshes.h"
 
namespace Geex
{

	/** representation of restricted Delaunay triangulation **/
	struct Embedded_point_2
	{
		//typedef K::Point_3 Point_2;
		Embedded_point_2() {}
		Embedded_point_2(const Point_3& _p) : p(_p) {}
		Point_3 p;
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
		Embedded_vertex_2(const Point_3& _p) : Embedded_point_2(_p), group_id(-1) {}
		Embedded_vertex_2(const Point_3& _p, int _group_id) : Embedded_point_2(_p), group_id(_group_id) {}
		inline Point_3 point_3() const { return Embedded_point_2::p; }
	public:
		int group_id;
		std::vector<Point_3> vd_vertices; // vertices of Voronoi diagram associated with this vertex
	};

	template < class Gt, class Fb = CGAL::Triangulation_face_base_2<Gt> >
	class Embedded_face_2 : public Fb
	{					};

	typedef CGAL::Triangulation_data_structure_2<Embedded_vertex_2<K>, Embedded_face_2<K>> RDT_data_structure; // representation of RDT


class RestrictedPolygonVoronoiDiagram
{
	typedef RestrictedVoronoiDiagram_poly RestrictedVoronoiDiagram;
	typedef RDT_data_structure::Vertex_handle Vertex_handle;
	typedef RDT_data_structure::Face_handle Face_handle;
	typedef std::pair<Vertex_handle, Vertex_handle> Vertex_pair;
	typedef std::vector<Vertex_handle> VertGroup; // vertices belonging to the same group

public:
	RestrictedPolygonVoronoiDiagram();
	~RestrictedPolygonVoronoiDiagram();

	/** set functions **/
	void set_mesh(const std::string& mesh_file_name);
	void set_trimesh(TriMesh *_trimesh) { assert(!trimesh); trimesh = _trimesh; }

	void begin_insert() ;

	// insert a group of polygons
	template<class InputIterator> void insert_polygons(InputIterator first, InputIterator last, unsigned int samp_nb = 24);

	void end_insert() ;

	void compute_clipped_VD(std::vector<Plane_3>& clipping_planes);

	void iDT_update(); // edge flip to preserve intrinsic Delaunay structure

	/** access functions **/


//#if _DEBUG
	void draw_DT();
	void draw_clipped_VD();
	inline void setGroupColor(unsigned int i, unsigned int nb_groups)
	{
		i = rand();
		GLfloat r = GLfloat(i)/nb_groups,
				g = GLfloat((i+nb_groups/3)%nb_groups)/nb_groups,
				b = GLfloat((i+2*nb_groups/3)%nb_groups)/nb_groups;
		glColor3f(r, g, b);
	}
//#endif


private: // private functions

	/** computation **/
	inline Plane_3 get_bisector(Vertex_handle v0, Vertex_handle v1)
	{
		return CGAL::bisector(v0->point_3(), v1->point_3());
	}

	// insert one polygon, return how many sample points are actually inserted
	unsigned int insert_polygons(const Polygon_3& polygon, unsigned int group_id, unsigned int samp_nb = 24); 

private:
	bool open;
	//double *pnt_coord; // point coordinates in array form

	// handles to sample points
	std::vector<VertGroup> samp_pnts;

	/** geometry **/
	TopoPolyMesh *mesh;
	TriMesh *trimesh;//the same as mesh, but with some other necessary functions
	//Delaunay *del;
	//RestrictedVoronoiDiagram *rvd;
	std::vector<std::vector<Point_3>> clipped_VD;//clipped Voronoi diagram for each polygon, expressed as a set of lines
	unsigned int nb_groups;

	RDT_data_structure rdt_ds;
	std::map<Vertex_pair, Face_handle> edge_face_adjacency;

	/** construct of RDT, helper class **/
	struct Construct_RDT_structure
	{
		Construct_RDT_structure(RDT_data_structure& _rdt_ds, std::vector<Vertex_handle>& _samp_pnts, std::map<Vertex_pair, Face_handle>& ef_adj) 
							: rdt_ds( _rdt_ds ), samp_pnts(_samp_pnts), edge_face_adjacency(ef_adj) { }
		
		void operator()(unsigned int i, unsigned int j, unsigned int k) const ;

		void find_adjacency(const Vertex_handle& v0, const Vertex_handle& v1, int v2, const Face_handle& f) const;

		RDT_data_structure& rdt_ds;
		std::vector<Vertex_handle>& samp_pnts;
		std::map<Vertex_pair, Face_handle>& edge_face_adjacency;
	};

//#if _DEBUG	
	/** debug for DT draw **/
	GLuint RDT_disp_list;
	GLuint clipped_VD_disp_list;
	std::vector<Point_3> cents;
//#endif
};

template <class InputIterator>
void RestrictedPolygonVoronoiDiagram::insert_polygons(InputIterator first, InputIterator last, unsigned int samp_nb)
{
	assert(open);
	unsigned int group_id = 0, nb_pnts = 0;
	samp_pnts.reserve(last - first);
	cents.reserve(last - first);
	while (first != last)
	{
		nb_pnts += insert_polygons(*first, group_id, samp_nb);
		first++;
		group_id++;
	}
	nb_groups = group_id;
	assert(nb_pnts == rdt_ds.number_of_vertices());
}



}

#endif