#ifndef _PACKING_OBJECT_H_
#define _PACKING_OBJECT_H_

#include <GL/GL.h>
#include <GL/glut.h>
#include "spm_cgal.h"

namespace Geex
{
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

	public:
		double factor; //scaling factor
		unsigned int lib_idx;
		int facet_idx; // index of triangle on which this object lies on
		std::vector<Point_2> texture_coord;
		unsigned int texture_id;
		bool active;
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