#pragma once

// define the exact number type
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
// define the kernels
#include <CGAL/Simple_cartesian.h>

// typedefs for the traits and the algorithm
#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/aff_transformation_tags.h>
#include <list>
#include <CGAL/intersections.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/enum.h>
#include <algorithm>

#include "spm_cgal.h"

using CGAL::bisector;
using CGAL::intersection;
using CGAL::do_intersect;
using CGAL::assign;
using CGAL::make_object;
using CGAL::sqrt;
using CGAL::orientation;
using CGAL::squared_distance;

//typedef CGAL::CK<CGAL::Cartesian<double>>	CK;
//typedef CGAL::Segment_Delaunay_graph_traits_without_intersections_2<CK>	SDGt;
typedef CGAL::Quotient<CGAL::MP_Float>  ENT;
typedef CGAL::Simple_cartesian<double>     CK;
typedef CGAL::Simple_cartesian<ENT>        EK;
typedef CGAL::Segment_Delaunay_graph_filtered_traits_2<CK,CGAL::Field_with_sqrt_tag, EK, CGAL::Field_tag>  SDGt;
typedef CGAL::Segment_Delaunay_graph_2<SDGt>	SDG_2;
typedef CGAL::Polygon_2<CK>	SDG_Polygon_2;
typedef CGAL::Segment_2<CK>	SDG_Segment_2;
typedef CGAL::Line_2<CK>		SDG_Line_2;
typedef CGAL::Vector_2<CK>		SDG_Vector_2;
typedef SDG_2::Site_2	Site_2;
typedef SDG_2::Point_2	SDG_Point_2;
typedef CGAL::Aff_transformation_2<CK>		Aff_transformation_2;
typedef SDGt::Object_2	Object;
//using CGAL::Object;

using namespace std;

namespace Geex{
// descriptor of a segment of medial axis of a simple polygon
struct MASegment 
{
	typedef enum {SEGMENT, PARABOLA} type;
	MASegment(type t, const SDG_Point_2& end0, const SDG_Point_2& end1, const Object& se0, const Object& se1);
	void setTransformation(const Aff_transformation_2& trans){transformation = trans;}
	//sample a parabolic segment to obtain a series of line segment
	list<SDG_Point_2>* sampleParabolicSegment( const unsigned int nbSampling = 100);

	/*data fields*/
	type curve_type;
	SDG_Point_2 endpoint[2];
	Object supporting_elements[2];//the original two sites defining the bisector
	Aff_transformation_2	transformation;//possible transformation
};

class SphericalBounder
{
public:
	// SphericalBounder(const SDG_Polygon_2& aPolygon);//initialized with the polygon to be approximated/bounded
	SphericalBounder(const vector<Segment_3>& aPolygon);
	~SphericalBounder();
	//CircleUnion* getCircleUnionBounder(const int nbSampling = 100);//get the circle union bounder got from sampling the medial axis
	// we change to use weighted_point to simply represent the circles now
	//const vector<Weighted_point>& getMedialBounder(const unsigned int nbSampling = 100);
	void genMarginBounder(const unsigned int nbSamping = 50);
	const vector<Weighted_point>& getMarginBounder() const { return *msu;}
	void drawPolygon();
	void drawMedialAxis();
private:
	//get the intersection point between an arbitrary parabola and a line passing through focus
	//note that only the result closer to the focus is returned
	//also, give the transformation 'afftrans' which transforms the parabola to one x-symmetric
	//internal usage
	Aff_transformation_2 intersectLineParabola(const SDG_Line_2& l, const SDG_Point_2& focus,
														const SDG_Line_2& directix, SDG_Point_2 *p0, SDG_Point_2 *p1);
	void createLinearMAList();//construct the display list of linear medial axis
	void createParabolicMAList();//construct the display list of parabolic medial axis
private:
	//SDG_Polygon_2 pgn;
	vector<Segment_3> pgn;
	SDG_2 sdg;//generalized delaunay triangulation, where sites are vertices and (open) edge of input polygon
	//CircleUnion *cu;//resulting circle union
	vector<Weighted_point> *su;
//	int nbCircles;//number of circles in this union
	list<MASegment*> *ma;//medial axis
	vector<Weighted_point> *msu;//margin spherical union
};

}//end of namespace Geex