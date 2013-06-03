#include "utilities.h"

namespace Geex
{
	double cur_size_map(double cur)
	{
		//return 0.4*std::pow(cur + 1.0, 0.35);
		//return 0.201*std::pow(cur, 0.6742) + 0.02238; // fitted with matlab
		return 0.25*std::pow(cur+1, 0.6);
		//return 1.0;
	}

	bool inside_polygon(const Point_2& p, const std::vector<Segment_2>& bounded_edges)
	{
		// construct an array from p to infinity
		Point_2 inf_p(1.0e30, 1.0e30);
		Ray_2 r(p, inf_p);
		unsigned int nb_intersetion = 0;
		for (unsigned int i = 0; i < bounded_edges.size(); i++)
		{
			if (CGAL::do_intersect(r, bounded_edges[i]))
				nb_intersetion++;
		}
		return (nb_intersetion%2 != 0);
	}
}