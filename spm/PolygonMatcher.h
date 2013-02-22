#pragma  once

#include <vector>
#include <cmath>
#include <algorithm>
//#include <CGAL/Aff_transformation_2.h>
//#include <CGAL/Vector_2.h>
//#include <CGAL/Polygon_2.h>
#include <cv.h>
#include <opencv2/core/types_c.h>
#include "spm_cgal.h"
#include "CmAffine.h"
#include <CGAL/convex_hull_2.h>
#include <CGAL/bounding_box.h>

namespace Geex
{
	static const double PI = 3.14159265358979323846;
	typedef struct _ConstPiece
	{
		double l; 
		double u; // a piece for piece wise function
		double f; // at interval [l, u), the function remains f
		_ConstPiece(double funcVal, double lowerBound, double upperBound):f(funcVal), l(lowerBound), u(upperBound){}
	}ConstPiece;

	/*
	** piecewise constant function on the interval [0,1)
	** and the function satisfies the following property
	**		f(x+n) = f(x) + n*2*pi,
	** where n is an integer.
	*/
	class PiecewiseConstFunc : public std::vector<ConstPiece>
	{
	public:
		inline double leftValueAt(int i) const
		{ 
			return rightValueAt(i-1);
		}
		inline double rightValueAt(int i) const
		{
			if (i >= 0)
				return at(i%size()).f + (i/size())*2.0*PI;
			else
			{
				int reflex = i + (-i/size()+1)*size();
				return at(reflex%size()).f -(reflex/size())*2.0*PI;
			}
		}
		// compute f^2(x)
		inline void square()
		{
			for (unsigned int i = 0; i < size(); i++)
				at(i).f *= at(i).f; 
		}
		double integral() const
		{
			//accumulate the upper and lower separately, in order to minimize roundoff errors
			double upperSum = 0.0, lowerSum = 0.0;
			for (unsigned int i = 0; i < size(); i++)
			{
				upperSum += at(i).u*at(i).f;
				lowerSum += at(i).l*at(i).f;
			}
			return (upperSum - lowerSum);
		}
	};
	// compute f(x)+s
	PiecewiseConstFunc operator+(const PiecewiseConstFunc& func, double s);
	// compute f + g on their common domain
	PiecewiseConstFunc operator+(const PiecewiseConstFunc& f, const PiecewiseConstFunc& g);
	// f(x)-g(x)
	PiecewiseConstFunc operator-(const PiecewiseConstFunc& f, const PiecewiseConstFunc& g);
	struct triplet
	{
		double val;
		int left;
		int right;
		triplet(): val(0.0), left(0), right(0) {}
		triplet(double v, int l, int r) : val(v), left(l), right(r) {}
		bool operator<=(const triplet& tri) { return (val <= tri.val); }
	} ;
	inline bool cmpTriplet(const triplet& t0, const triplet& t1) { return (t0.val < t1.val); }
	
	/*
	** class for computing the similarity between two polygons
	*/
	using CGAL::sqrt;
	using CGAL::transform;
	using CGAL::orientation;
	using std::vector;
	//typedef CGAL::Polygon_2<K> Polygon_2;
	//typedef CGAL::Aff_transformation_2<K>	Transformation_2;
	//typedef CGAL::Vector_2<K>	Vector_2;
	class PolygonMatcher
	{
	public:
		PolygonMatcher(const Polygon_2& model);
		PolygonMatcher(const vector<Segment_2>& model);
		// match based on points
		PolygonMatcher(const vector<Point_2>& model);
		~PolygonMatcher();
		// insert a polygon to compare it with the model, return their similarity distance
		double insert(const Polygon_2& instance);
		/************************************************************************/
		/*							match with affine registration             */
		/************************************************************************/
		double affineMatch(const Polygon_2& instance, double *scale, double *theta, double *tx, double *ty);
		double getAffineMatchAngle() const;
		double getAffineTranslation() const;
		/************************************************************************/
		double getSimilarity(int i) const { return similarityDistances[i]; }
		double getMatchAngle(int i)	const { return matchAngles[i]; }
		bool getyAxisReflection(int i) const { return yAxisReflected[i]; }
		//size_t nbInstancePolygons()	const { return instances.size(); }
	private:
		// compute the similarity distance between the polygon with turning function "tf" and the model function
		double computeSimilarity(const PiecewiseConstFunc& tf, double *matchAngle);
		inline double radian(double x, double y);
		void constructTurningFunc(const Polygon_2& pgn, PiecewiseConstFunc *func);
		// compute the integral of the turning function on [0,1]
		//double integral(const PiecewiseConstFunc& func);
	private:
		Polygon_2 model;
		PointSet modelPtsSet;
		CvMat *modelMat;
		//vector<Polygon_2> instances;
		PiecewiseConstFunc modelTurningFunc;//piecewise turning function of model polygon
		double modelInt;//integral on [0,1] of the turning function of the model polygon
		vector<double> similarityDistances;
		vector<double> matchAngles;
		vector<bool> yAxisReflected;
		double modelArea;
};

}

