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
#include "spm_cgal.h"
 
namespace Geex
{

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

private:
	void uniform_sample(const Polygon_3& p, unsigned int samp_nb); // sample polygons

private:
	bool open;
	double *pnt_coord; // point coordinates in array form
	std::vector<unsigned int> which_group; // belong to which polygon. element is polygon id
	std::vector<Point_3> samp_pnts;
	// geometry
	TopoPolyMesh *mesh;
	Delaunay *del;
	RestrictedVoronoiDiagram *rvd;
	//std::vector<Polygon_3> *polygon_sites;

	/** debug for DT draw **/
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