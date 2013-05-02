#ifndef	_UTILITIES_H_
#define  _UTILITIES_H_

#include "spm_cgal.h"

namespace Geex
{
	double cur_size_map(double cur);
// check whether a point is inside or outside a polygon represented as an array of line segments
// use even-odd rule
bool inside_polygon(const Point_2& p, const std::vector<Segment_2>& bounded_edges);

}

#endif