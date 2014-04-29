#ifndef _INIT_POINT_GENERATOR_
#define _INIT_POINT_GENERATOR_

#include <cstdlib>
#include <set>
#include <vector>
#include <string>
#include <CGAL/centroid.h>
#include <CGAL/Dimension.h>
#include <Geex/cvt/delaunay.h>
#include <Geex/cvt/RVD.h>
#include "spm_cgal.h"
#include "meshes.h"

namespace Geex
{

class Init_point_generator
{

public:

	Init_point_generator(const TriMesh& trimesh, unsigned int nb_init_pnts) : mesh(trimesh), nb_pnts(nb_init_pnts) { nb_pnt_fetched = 0; }

	virtual Point_3 next_point() = 0;

	virtual unsigned int nb_points() = 0;

	virtual ~Init_point_generator() { } 
	
	inline double density(unsigned int facet_idx) 
	{
		double cur = mesh.curvature_at_face(facet_idx);
		return std::pow(cur, 2.0);
	}

	inline double weight(unsigned int facet_idx);

protected:
	const TriMesh& mesh;
	unsigned int nb_pnt_fetched;
	unsigned int nb_pnts;
};

class Uniform_point_generator : public Init_point_generator
{

public:
	Uniform_point_generator(const TriMesh& trimesh, unsigned int nb_init_pnts) : Init_point_generator(trimesh, nb_init_pnts) { std::srand(time(NULL)); }

	Point_3 next_point();

	unsigned int nb_points() { return nb_pnts; }

private:
	std::set<int> facets_used;
};

class Density_point_generator : public Init_point_generator
{
public:
	Density_point_generator(const TriMesh& trimesh, unsigned int nb_init_pnts);

	Point_3 next_point();

	unsigned int nb_points() { return pnts.size(); }

private:
	std::vector<Point_3> pnts;

};

class CVT_point_generator : public Init_point_generator
{
	typedef RestrictedVoronoiDiagram_poly RVD;

public:
	CVT_point_generator(const TriMesh& trimesh, unsigned int nb_init_pnts, const std::string& mesh_filename);

	Point_3 next_point();
	
	unsigned int nb_points() ;
	
	~CVT_point_generator() { delete topo_mesh; }

private:
	std::vector<vec3> pnts;
	TopoPolyMesh *topo_mesh;
};

// used by CVT_point_generator, to extract information in order for computation of centroids
class Compute_centroids
{
public:
	Compute_centroids(const TriMesh& _mesh, const std::vector<double>& _facets_density, unsigned int nb_seeds);
	void operator()(unsigned int vi, TopoPolyMesh *clipped_region) const;
	void get_centroids(std::vector<vec3>& centroids) const;
private:
	mutable std::vector<vec3> numerator_sum;
	mutable std::vector<double> denominator_sum;
	// just for density extraction
	const std::vector<double>& facets_density;
	const TriMesh& mesh;
};

}
#endif