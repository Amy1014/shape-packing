/*
 *  _____ _____ ________  _
 * /  __//  __//  __/\  \//
 * | |  _|  \  |  \   \  / 
 * | |_//|  /_ |  /_  /  \ 
 * \____\\____\\____\/__/\\
 *
 * Graphics Environment for EXperimentations.
 *  Copyright (C) 2006 INRIA - Project ALICE
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: 
 *
 *     ALICE Project - INRIA
 *     INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 *  Note that the GNU General Public License does not permit incorporating
 *  the Software into proprietary programs. 
 */

#ifndef __GEOMETRY_CCVT__ 	/*! add by dmyan*/
#define __GEOMETRY_CCVT__

#include <Geex/symbolics/expression.h>
#include <Geex/mathematics/glsl_linear.h>
//#include <ginac/flags.h>

namespace Geex {

    template <class T> inline T tetra_radius_ratio(
        const vec3g<T>& p, const vec3g<T>& q, 
        const vec3g<T>& r, const vec3g<T>& s
    ) {
		return T(3)*tetra_incirbedradius(p, q, r, s)/tetra_circumradius(p, q, r, s);
	}

    template <class T> inline T tetra_incirbedradius(
        const vec3g<T>& p, const vec3g<T>& q, 
        const vec3g<T>& r, const vec3g<T>& s
    ) {
		return T(3)*tetra_volume(p, q, r, s)/(tri_area(q, r, s)+tri_area(r, s, p)+tri_area(s, p, q)+tri_area(p, r, q)) ;
	}

    template <class T> inline T tetra_circumradius(
        const vec3g<T>& p, const vec3g<T>& q, 
        const vec3g<T>& r, const vec3g<T>& s
    ) {
		return distance(p, tetra_circumcenter(p, q, r, s) ) ;
	}

	template <class T> inline bool on_segment(const vec3g<T>& p, const vec3g<T>& a, const vec3g<T>& b) {
		if(distance(p, a) + distance(p, b) == distance(a, b))
			return true ;
		return false ; 
	}

	template <class T> inline bool inside_triangle(const vec3g<T>& p, const vec3g<T>& a, const vec3g<T>& b, const vec3g<T>& c) {
		vec3g<T> n = cross(b-a, c-a);
		//T a0 = dot(cross(a-p, b-p), n); 
		//T a1 = dot(cross(b-p, c-p), n); 
		//T a2 = dot(cross(c-p, a-p), n);
		if(dot(cross(a-p, b-p), n)<-1e-10 || dot(cross(b-p, c-p), n)<-1e-10 || dot(cross(c-p, a-p), n)<-1e-10) 
			return false;
		else 
			return true;
	}
		
	template <class T> inline bool ray_triangle_intersection(vec3g<T>& p0, vec3g<T>& dir, vec3g<T>& a, vec3g<T>& b, vec3g<T>& c, vec3g<T>& n, vec3g<T>& intp) {
		Plane<T> P(a, n);
		if(ray_plane_intersects(p0, dir, P)) {
			intp = ray_plane_intersection(p0, dir, P);
			if(inside_triangle(intp, a, b, c))
				return true;
		}
		return false;
	}

	template <class T> inline bool ray_plane_intersects(vec3g<T>& p0, vec3g<T>& dir, const Plane<T>& plane) {
		bool side = plane.side(p0) > 0.0;
		bool sign = dot(dir, plane.normal()) > 0.0;

		if(side && sign || !side && !sign)
			return false;
		return true;
	}

	template <class T> inline bool seg_plane_intersects(vec3g<T>& p0, vec3g<T>& p1, const Plane<T>& plane) {
		if(plane.side(p0)>0.0 != plane.side(p1)>0.0)
			return true;
		return false;
	}

	template <class T> inline bool seg_triangle_intersection(vec3g<T>& p0, vec3g<T>& p1, vec3g<T>& a, vec3g<T>& b, vec3g<T>& c, vec3g<T>& n, vec3g<T>& intp) {
		Plane<T> P(a, n);
		if(seg_plane_intersects(p0, p1, P)) {
			intp = seg_plane_intersection(p0, p1, P);
			if(inside_triangle(intp, a, b, c))
				return true;
		}
		return false;
	}

	template <class T> inline vec3g<T> project_to_plane(const vec3g<T>& p, const Plane<T>& plane) {
		return p-plane.side(p)*plane.normal() ;
	}

	template <class T> inline bool project_to_triangle(const vec3g<T>& p, const vec3g<T>& a, const vec3g<T>& b, const vec3g<T>& c, vec3g<T>& proj) {
		vec3g<T> pc = p - c;
		vec3g<T> ac = a - c;
		vec3g<T> bc = b - c;
		T a00 = ac.length2();
		T a01 = dot(ac, bc);
		T a11 = bc.length2();
		T b0 = dot(pc, ac);
		T b1 = dot(pc, bc);
		T mdet = a00 * a11 - a01 * a01;
		T s = (a11 * b0 - a01 * b1) / mdet;
		T t = (a00 * b1 - a01 * b0) / mdet;
		proj = s * a + t * b + (1 - s - t ) * c;
		if ( s < 0 || s > 1 || t < 0 || t > 1)
		{
			return false;
		}
		else
			return true;
		//vec3g<T>& n = cross(b-a, c-a);
		//n = normalize(n);
		//proj = project_to_plane(p, Plane<T>(a, n));
		//if(inside_triangle(proj, a, b, c))
		//	return true;
		//return false;
	}

		template <class T> inline void project_to_linesegment(const vec3g<T>& p, const vec3g<T>& a, const vec3g<T>& b, vec3g<T> & fp, T & sqdist)
	{
		//double t;

		vec3g<T> ab = a - b;
		vec3g<T> pb = p - b;
		T t = dot(pb, ab) / ab.length2();
		//t = (-a.x*b.x + b.x*b.x - a.y*b.y + b.y*b.y - a.z*b.z + b.z*b.z + a.x*p.x
		//	-b.x*p.x +a.y*p.y-b.y*p.y+a.z*p.z-b.z*p.z) 
		//	/ 
		//	(a.x*a.x+a.y*a.y+a.z*a.z-2.*a.x*b.x-2.*a.y*b.y-2.*a.z*b.z+b.x*b.x+b.y*b.y+b.z*b.z);

		if(t < 0)
		{
			fp = b;
			sqdist = pb.length2();
			return;
		}
		else if(t > 1)
		{
			fp = a;
			sqdist = (p-a).length2();
			return;
		}
		else
		{
			fp = t*a+(1-t)*b;
			sqdist = (p-fp).length2();
			return;
		}
	}

    template <class T> inline T tri_area(
        const vec3g<T>& p, const vec3g<T>& q, 
        const vec3g<T>& r
    ) {
        return length(cross(q-p, r-p)) / T(2);
    }

    template <class T> inline vec3g<T> tri_centroid(
        const vec3g<T>& p, const vec3g<T>& q, 
        const vec3g<T>& r
    ) {
        return T(1)/T(3)*(p+q+r);
    }

	template <class T> inline T tri_mass(
        const vec3g<T>& p, const vec3g<T>& q, const vec3g<T>& r,
		T a, T b, T c // a, b, c are density at p, q, r
    ) {
		return tri_area(p, q, r)/T(3)*(sqrt(fabs(a))+sqrt(fabs(b))+sqrt(fabs(c))) ;
	}

	template <class T> inline void tri_centroid(
		const vec3g<T>& p, const vec3g<T>& q, const vec3g<T>& r,
		vec3g<T>& g, T& V
		) {
			V = tri_area(p, q, r) ;
			g = T(1.)/T(3.) * ( p+q+r ) ;
	}

    template <class T> inline void tri_centroid(
        const vec3g<T>& p, const vec3g<T>& q, const vec3g<T>& r,
		T a, T b, T c, vec3g<T>& g, T& V // a, b, c are density at p, q, r
    ) {
		T abc = a+b+c ; 
		V = tri_area(p, q, r) ;
		g =  T(1) / (4*abc) * ( (a+abc)*p+(b+abc)*q+(c+abc)*r );
	}

    template <class T> inline void tetra_centroid(
        const vec3g<T>& p, const vec3g<T>& q, const vec3g<T>& r, const vec3g<T>& s,
		T a, T b, T c, T d, vec3g<T>& g, T& V // a, b, c d are density at p, q, r, s
    ) {
		T abcd = a+b+c+d ; 
		V = tetra_volume(p, q, r, s) ;
		g =  T(1) / (5*abcd) *  ( (a+abcd)*p+(b+abcd)*q+(c+abcd)*r+(d+abcd)*s );
	}

    template <class T> inline bool inside_tetra(const vec3g<T>& p0, 
        const vec3g<T>& p, const vec3g<T>& q, 
        const vec3g<T>& r, const vec3g<T>& s
    ) {
		vec3 dir[4];
		dir[0] = p - p0;
		dir[1] = q - p0;
		dir[2] = r - p0;
		dir[3] = s - p0;

		real msignvolume[4] ;
		msignvolume[0] = Geex::dot(dir[3], Geex::cross(dir[1], dir[2]));
		msignvolume[1] = -Geex::dot(dir[3], Geex::cross(dir[0], dir[2]));
		msignvolume[2] = Geex::dot(dir[3], Geex::cross(dir[0], dir[1]));
		msignvolume[3] = -Geex::dot(dir[2], Geex::cross(dir[0], dir[1]));

		if (msignvolume[0]>= 0 && 
			msignvolume[1]>= 0 && 
			msignvolume[2]>= 0 && 
			msignvolume[3]>= 0 )
		{
			return true;
		}
		return false;
	}

    template <class T> inline bool inside_tetra(const vec3g<T>& p0, 
        const vec3g<T>& p, const vec3g<T>& q, 
        const vec3g<T>& r, const vec3g<T>& s, real msignvolume[] 
    ) {
		vec3 dir[4];
		dir[0] = p - p0;
		dir[1] = q - p0;
		dir[2] = r - p0;
		dir[3] = s - p0;

		msignvolume[0] = Geex::dot(dir[3], Geex::cross(dir[1], dir[2]));
		msignvolume[1] = -Geex::dot(dir[3], Geex::cross(dir[0], dir[2]));
		msignvolume[2] = Geex::dot(dir[3], Geex::cross(dir[0], dir[1]));
		msignvolume[3] = -Geex::dot(dir[2], Geex::cross(dir[0], dir[1]));

		if (msignvolume[0]>= 0 && 
			msignvolume[1]>= 0 && 
			msignvolume[2]>= 0 && 
			msignvolume[3]>= 0 )
		{
			return true;
		}
		return false;
	}

	template<class T> inline bool inside_sphere(const vec3g<T>& p, const vec3g<T>& o, const T r) {
		if(distance2(p, o) >= r*r)
			return false ;
		return true ;
	}

	template<class T> inline bool triangle_sphere_intersects(
		const vec3g<T>& a, 
		const vec3g<T>& b, 
		const vec3g<T>& c, 
		const vec3g<T>& o,
		const T r) {
			if(inside_sphere(a, o, r) || inside_sphere(b, o, r) || inside_sphere(c, o, r)) {
				return true ;
			}
			else {
				vec3g<T> proj ;
				if(project_to_triangle(o, a, b, c, proj)) {
					if(distance2(o, proj) < r*r)
						return true ;
				}
			}
			return false; 
	}

    template<class T> inline bool inside_box(const vec3g<T>& p, vec3g<T>& bmin, vec3g<T>& bmax) {
	if(p[0]<bmin[0] || p[0]>bmax[0] || p[1]<bmin[1] || p[1]>bmax[1] || p[1]<bmin[1] || p[1]>bmax[1])
	    return false ; 
	return true ; 
    }
    
    template<class T> inline bool is_nan(vec3g<T>& p) {
	if(Numeric::is_nan(p.x) || Numeric::is_nan(p.y) || Numeric::is_nan(p.z)) {
	    std::cerr << "Nan !" << std::endl ;
	    return true;
	}
	return false ;
    }
}

#endif
