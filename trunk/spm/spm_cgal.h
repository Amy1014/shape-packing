#ifndef _SPM_CGAL_H_
#define _SPM_CGAL_H_

#include <ostream>
#include <vector>
#include <string>
#include <CGAL/Polygon_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAl/Regular_triangulation_filtered_traits_3.h>
#include <CGAL/Triangulation_hierarchy_vertex_base_3.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Extended_cartesian.h>
#include <CGAL/Object.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/utility.h>
#include <Geex/CVT/geometry.h>
#include "Polygon_3.h"


namespace Geex {
	
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef K::Point_2		Point_2;
	typedef K::Vector_2		Vector_2;
	typedef K::Line_2		Line_2;
	typedef K::Segment_2	Segment_2;
	typedef CGAL::Polygon_2<K>	Polygon_2;
	typedef CGAL::Aff_transformation_2<K>	Transformation_2;
	typedef CGAL::Triangle_2<K>				Triangle_2;
	typedef K::Point_3    Point_3;
	typedef K::Vector_3   Vector_3;
	typedef K::Segment_3	Segment_3;
	typedef K::Ray_3		Ray_3;
	typedef K::Direction_3	Direction_3;
	typedef CGAL::Plane_3<K>	Plane_3;
	typedef CGAL::Triangle_3<K> Triangle_3;
	typedef CGAL::Object	Object;
	typedef CGAL::Aff_transformation_3<K>	Transformation_3;
	typedef MyPolygon_3<K> Polygon_3;

	
	inline void prompt_and_exit(std::string prompt)
	{
		std::cerr<<prompt<<std::endl;
		system("pause");
		exit(0);
	}

	inline Point_3 to_cgal_pnt(const vec3& p)
	{
		return Point_3(p.x, p.y, p.z);
	}

	inline vec3 to_geex_pnt(const Point_3& p)
	{
		return vec3(p.x(), p.y(), p.z());
	}

	inline Vector_3 to_cgal_vec(const vec3& v)
	{
		return Vector_3(v.x, v.y, v.z);
	}

	inline vec3 to_geex_vec(const Vector_3& v)
	{
		return vec3(v.x(), v.y(), v.z());
	}

	template <class VectorType>
	inline void cgal_vec_normalize(VectorType& v)
	{
		v = v / CGAL::sqrt(v.squared_length());
	}

}

#endif