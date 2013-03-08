#ifndef _POLYGON_3_H_
#define  _POLYGON_3_H_

/*************************************************************************/
/**						My 3D polygon type                              **/
/**@Author: Wenchao Hu													**/
/************************************************************************/
#include <vector>
#include <algorithm>
#include <iterator>
#include <CGAL/Vector_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Dimension.h>
#include <CGAL/centroid.h>
#include <CGAL/Aff_transformation_3.h>

namespace Geex
{
template < class Kernel, class Container = std::vector<typename Kernel::Point_3> >
class MyPolygon_3
{
	
public:
	typedef typename Kernel::FT			FT;
	typedef typename Kernel::Point_3	Point_3;
	typedef typename Kernel::Segment_3	Segment_3;
	typedef typename Kernel::Vector_3	Vector_3;
	typedef CGAL::Aff_transformation_3<Kernel> Transformation_3;

	typedef typename Container::size_type size_type;
	
	template <class Kernel, class Container>
	friend MyPolygon_3<Kernel, Container> operator*(const Transformation_3& t, const MyPolygon_3<Kernel, Container>& p);

public:
	/** constructors **/
	MyPolygon_3();

	template<class PolygonTraits_2, class Container_2>
	MyPolygon_3(const CGAL::Polygon_2<PolygonTraits_2, Container_2>& xoy_polygon, const Vector_3& v, const Point_3& pos);

	template <class InputIterator> MyPolygon_3(InputIterator first, InputIterator last);

	/** destructors **/
	~MyPolygon_3();

	/** set functions **/
	void push_back(const Point_3& p);
	void clear();

	/** access functions **/
	Point_3 vertex(size_type i) const 
	{ 
		assert(i<verts.size());
		return verts[i]; 
	}

	Point_3 centroid() const { return cent; }

	Vector_3 norm() const { return normal; }

	size_type size() const { return verts.size(); }

	bool empty() const { return size() == 0; }

	Segment_3 edge(size_type i) const 
	{ 
		assert(i < verts.size());
		return Segment_3(verts[i%verts.size()], verts[(i+1)%verts.size()]);
	}

	/** computation **/
	// unsigned area
	FT area() const ;
	MyPolygon_3& operator*=(const Transformation_3& t);
	

private:
	Container verts;//vertices
	Point_3 cent;//centroid
	Vector_3 normal;
	//CGAL::Plane_3<Kernel> embed_plane;
};

template <class Kernel, class Container> MyPolygon_3<Kernel, Container>::MyPolygon_3()
{
	cent = Point_3(FT(0), FT(0), FT(0));
	normal = CGAL::NULL_VECTOR;
	//embed_plane = Plane_3(cent, normal);
}

/* 
** construct a 3d polygon from an xOy polygon, by aligning with the normal v and fixing at point pos 
*/
template <class Kernel, class Container> template <class PolygonTraits_2, class Container_2> 
MyPolygon_3<Kernel, Container>::MyPolygon_3(const CGAL::Polygon_2<PolygonTraits_2, Container_2>& xoy_polygon, const Vector_3& v, const Point_3& pos)
{
	// centroid of xOy polygon
	PolygonTraits_2::Point_2 cent_2 = CGAL::centroid(xoy_polygon.vertices_begin(), xoy_polygon.vertices_end(), CGAL::Dimension_tag<0>());
	Transformation_3 to_org(CGAL::TRANSLATION, Vector_3(-cent_2.x(), -cent_2.y(), 0.0));
	Transformation_3 to_pos(CGAL::TRANSLATION, Vector_3(CGAL::ORIGIN, pos));
	Vector_3 z(FT(0), FT(0), FT(1));
	Vector_3 n = v / CGAL::sqrt(v.squared_length());
	Vector_3 lz = CGAL::cross_product(z, n);
	Transformation_3 t(CGAL::IDENTITY);
	FT lz_len2 = lz.squared_length();
	if (lz_len2 != FT(0))
	{
		lz = lz / CGAL::sqrt(lz_len2);
		Vector_3 lx = z;
		Vector_3 ly = CGAL::cross_product(lz, lx);
		// to local frame
		Vector_3 ln(n*lx, n*ly, FT(0));
		FT len = CGAL::sqrt(ln.squared_length());
		Vector_3 ex(Vector_3(FT(1), FT(0), FT(0)));
		FT cos_theta = ln*ex/len, sin_theta = CGAL::sqrt(CGAL::cross_product(ex, ln).squared_length())/len;
		Transformation_3 rot(cos_theta, -sin_theta, FT(0), sin_theta, cos_theta, FT(0), FT(0), FT(0), FT(1));
		Transformation_3 to_global(lx.x(), ly.x(), lz.x(), lx.y(), ly.y(), lz.y(), lx.z(), ly.z(), lz.z());
		Transformation_3 to_local(lx.x(), lx.y(), lx.z(), ly.x(), ly.y(), ly.z(), lz.x(), lz.y(), lz.z());
		t = (to_pos*to_global)*rot*(to_local*to_org);
	}
	else
		t = to_pos*to_org;
	for (unsigned int i = 0; i < xoy_polygon.size(); i++)
	{
		PolygonTraits_2::Point_2 vert = xoy_polygon.vertex(i);
		verts.push_back(t(Point_3(vert.x(), vert.y(), FT(0))));
	}
	normal = n;
	cent = CGAL::centroid(verts.begin(), verts.end(), CGAL::Dimension_tag<0>());

}

template <class Kernel, class Container> MyPolygon_3<Kernel, Container>::~MyPolygon_3()
{
	verts.clear();
}

template <class Kernel, class Container> void MyPolygon_3<Kernel, Container>::push_back(const Point_3& p)
{
	verts.push_back(p);
	cent = CGAL::centroid(verts.begin(), verts.end(), CGAL::Dimension_tag<0>());
}

template <class Kernel, class Container> void MyPolygon_3<Kernel, Container>::clear()
{
	verts.clear();
	normal = CGAL::NULL_VECTOR;
	cent = Point_3(FT(0.0), FT(0.0), FT(0.0));
}

template <class Kernel, class Container> typename MyPolygon_3<Kernel, Container>::FT MyPolygon_3<Kernel, Container>::area() const 
{
	Vector_3 a(FT(0.0), FT(0.0), FT(0.0));
	for (unsigned int i = 0; i < verts.size(); i++)
	{
		Vector_3 ta = CGAL::cross_product(Vector_3(cent, verts[i]), Vector_3(cent, verts[(i+1)%verts.size()]));
		a = a + ta;
	}
	return 0.5*CGAL::sqrt(a.squared_length());
}

template <class Kernel, class Container> 
MyPolygon_3<Kernel, Container>& MyPolygon_3<Kernel, Container>::operator*=(const Transformation_3& t)
{
	std::transform(verts.begin(), verts.end(), verts.begin(), t);
	cent = CGAL::centroid(verts.begin(), verts.end(), CGAL::Dimension_tag<0>());
	normal = normal.transform(t);
	return *this;
}

/* 
** apply transformation to a polygon and return a new result 
*/
template <class Kernel, class Container> 
MyPolygon_3<Kernel, Container> operator*(const typename MyPolygon_3<Kernel, Container>::Transformation_3& t, const MyPolygon_3<Kernel, Container>& p)
{
	MyPolygon_3<Kernel, Container> res;
	std::transform(p.verts.begin(), p.verts.end(), std::back_inserter(res.verts), t);
	res.cent = CGAL::centroid(res.verts.begin(), res.verts.end(), CGAL::Dimension_tag<0>());
	res.normal = p.normal.transform(t);
	return res;
}

}
#endif
