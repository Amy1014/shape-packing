#ifndef _PACKING_OBJECT_H_
#define _PACKING_OBJECT_H_

#include <GL/GL.h>
#include <GL/glut.h>
#include "spm_cgal.h"

namespace Geex
{
	struct Local_frame // represent a local frame for which a transformation is effective
	{
		Vector_3 u;
		Vector_3 v;
		Vector_3 w;
		Point_3 o;
		inline Point_2 to_uv(const Point_3& p) const
		{
			Vector_3 op(o, p);
			return Point_2(op*u, op*v);
		}
		inline Segment_2 to_uv(const Segment_3& s) const 
		{
			return Segment_2(to_uv(s.source()), to_uv(s.target()));
		}
		inline Point_3 to_xy(const Point_2& p) const
		{
			return o + ( p.x()*u + p.y()*v );
		}
	};

	struct Parameter // represent a 2D transformation parameter
	{
		double k;
		double theta;
		double tx;
		double ty;
		Parameter(void): k(1.0), theta(0.0), tx(0.0), ty(0.0) {}
		Parameter(double _k, double _theta, double _tx, double _ty) : k(_k), theta(_theta), tx(_tx), ty(_ty) {}
		Parameter& operator*=(double f)
		{
			//k = std::max(f*k, 1.0); 
			//k = 1.0 + (k - 1.0)*f;
			theta *= f;	tx *= f; ty *= f;
			return *this;
		}
		// compute the remaining transformation after applying transformation p
		Parameter& operator-=(const Parameter& p)
		{
			k /= p.k;
			theta -= p.theta;
			tx -= p.tx;
			ty -= p.ty;
			return *this;
		}
		inline Parameter operator*(double f)
		{
			//return Packer::Parameter(1.0 + (k - 1.0)*f, f*theta, f*tx, f*ty);
			return Parameter(k, f*theta, f*tx, f*ty);
		}
	};

	class Packing_object : public Polygon_3
	{
	public:
		Packing_object() : factor(1.0)/*, group_id(-1)*/ {}

		Packing_object(const Polygon_2& xoy_polygon, const Vector_3& n, const Point_3& p)
			: Polygon_3(xoy_polygon, n, p), factor(1.0) {}

		Packing_object(const Polygon_2& xoy_polygon, const Vector_3& n, const Point_3& p, double f)
			: Polygon_3(xoy_polygon, n, p), factor(f) {}

		// scale
		Packing_object& operator*=(double factor)
		{
			Polygon_3::operator*=(factor);
			this->factor *= factor;
			return *this;
		}

		inline void activate() { active = true; }

		Local_frame local_frame() const
		{
			Local_frame lf;
			lf.o = centroid();
			lf.w = norm(); // assume this vector has already been normalized
			Vector_3 u(lf.o, vertex(0));
			lf.u = u/CGAL::sqrt(u.squared_length());
			lf.v = CGAL::cross_product(lf.w, lf.u);
			return lf;
		}

	public:
		double factor; //scaling factor
		unsigned int lib_idx;
		int facet_idx; // index of triangle on which this object lies on
		std::vector<Point_2> texture_coord;
		unsigned int texture_id;
		bool active;
		bool reach_barrier;
		//int group_id;
	};

	struct Ex_polygon_2 : public Polygon_2
	{
		std::vector<Point_2> texture_coords;
		unsigned int texture_id;
		//double normalize_factor;
	};

}

#endif