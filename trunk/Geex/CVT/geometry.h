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

#ifndef __CVT_GEOMETRY__
#define __CVT_GEOMETRY__

#include <Geex/symbolics/expression.h>
#include <Geex/mathematics/glsl_linear.h>

namespace Geex {

    template <class T> class Plane {
    public:
        Plane() { }

        Plane(T a_in, T b_in, T c_in, T d_in) :
            a(a_in), b(b_in), c(c_in), d(d_in) {
        }

        Plane(
            const vec3g<T>& P, const vec3g<T>& N
        ) : a(N.x), b(N.y), c(N.z), d(-dot(P,N)) {
        }

        T side(const vec3g<T>& p) const { 
            T result = a*p.x + b*p.y + c*p.z + d ; 
            return result ;
        }

        vec3g<T> normal() const { return vec3(a,b,c) ; }

        T a,b,c,d ;
    } ;

    template <class T> double distance2(const Plane<T>& P, const vec3g<T>& p) {
        T s = P.side(p) ;
        return s*s / (P.a * P.a + P.b * P.b + P.c * P.c) ;
    }

    template <class T> inline T det2x2(
        T a00, T a01,
        T a10, T a11
    ) {
        return a00*a11 - a10 * a01 ;
    }

    template <class T> inline T det3x3(
        T a00,  T a01,  T a02,
        T a10,  T a11,  T a12,
        T a20,  T a21,  T a22
    ) {
        T m01 = a00*a11 - a10*a01;
        T m02 = a00*a21 - a20*a01;
        T m12 = a10*a21 - a20*a11;
        return m01*a22 - m02*a12 + m12*a02;
    }

    template <class T> inline void comatrix3x3(
        T a00,  T a01,  T a02,
        T a10,  T a11,  T a12,
        T a20,  T a21,  T a22,
        T& b00,  T& b01,  T& b02,
        T& b10,  T& b11,  T& b12,
        T& b20,  T& b21,  T& b22
    ) {
        b00 = a11*a22 - a12*a21 ;
        b01 = a12*a20 - a10*a22 ;
        b02 = a10*a21 - a11*a20 ;
        b10 = a02*a21 - a01*a22 ;
        b11 = a00*a22 - a02*a20 ;
        b12 = a01*a20 - a00*a21 ;
        b20 = a01*a12 - a02*a11 ;
        b21 = a02*a10 - a00*a12 ;
        b22 = a00*a11 - a01*a10 ;
    }

    template <class T> inline void invert3x3(
        T a00,  T a01,  T a02,
        T a10,  T a11,  T a12,
        T a20,  T a21,  T a22,
        T& b00,  T& b01,  T& b02,
        T& b10,  T& b11,  T& b12,
        T& b20,  T& b21,  T& b22
    ) {
        //T d = det3x3(
        //    a00, a01, a02,
        //    a10, a11, a12,
        //    a20, a21, a22
        //)  ;

        // Compute transpose of the comatrix
        // (note that bij's are transposed in func call)
        comatrix3x3(
            a00, a01, a02,
            a10, a11, a12,
            a20, a21, a22,
            b00, b10, b20,
            b01, b11, b21,
            b02, b21, b22
        ) ;

        T d = a00*b00 + a01*b10 + a02*b20;


        b00 /= d ; b01 /= d ; b02 /= d ;
        b10 /= d ; b11 /= d ; b12 /= d ;
        b20 /= d ; b21 /= d ; b22 /= d ;
    }

    template <class T> inline T mixed_product(
        const vec3g<T>& v1, const vec3g<T>& v2, const vec3g<T>& v3
    ) {
        return det3x3(
            v1.x, v2.x, v3.x,
            v1.y, v2.y, v3.y,
            v1.z, v2.z, v3.z
        ) ;
    }

    template <class T> inline vec3g<T>
    triangle_normal_with_area(const vec3g<T>& v1, const vec3g<T>& v2, const vec3g<T>& v3) {
        return cross(v2-v1, v3-v1) ;
    }

    template <class T> inline vec3g<T> tetra_circumcenter(
        const vec3g<T>& p, const vec3g<T>& q, 
        const vec3g<T>& r, const vec3g<T>& s) {

        vec3g<T> qp = q-p;
        vec3g<T> rp = r-p;
        vec3g<T> sp = s-p;

        T tet_vol_x_12 = T(2) * mixed_product(qp,rp,sp) ;

        return p +  T(1) / tet_vol_x_12 *
			(qp.length2() * cross(rp, sp) + 
			rp.length2() * cross(sp, qp) +
			sp.length2() * cross(qp, rp) );
    }

    template <class T> inline bool tetra_circumcenter_squaredradius(
        const vec3g<T>& p, const vec3g<T>& q, 
        const vec3g<T>& r, const vec3g<T>& s,
        vec3g<T>& t_circumcenter, T& squared_radius, T *dihedral_angles = 0) {

            const vec3g<T> qp = q-p;
            const vec3g<T> rp = r-p;
            const vec3g<T> sp = s-p;

            const T tet_vol_x_12 = T(2) * mixed_product(qp,rp,sp) ;

            if (tet_vol_x_12 == 0) // it is not a safe check, one should avoid to input a degenerated tetra.
            {
                return false;
            }

            t_circumcenter = 
                qp.length2() * cross(rp, sp) + 
                rp.length2() * cross(sp, qp) +
                sp.length2() * cross(qp, rp);
            squared_radius = t_circumcenter.length2() / (tet_vol_x_12*tet_vol_x_12);
            t_circumcenter /=  tet_vol_x_12;
            t_circumcenter += p;

            if (dihedral_angles)
            {
                //compute Dihedral angle directly
                //the following computation is not optimized.
                vec3g<T> point[4];
                point[0] = p; point[1] = q;  point[2] = r; point[3] = s;
                static int pid[6][4] = {{0,1,2,3},{0,2,1,3},{0,3,1,2},{1,2,0,3},{1,3,2,0},{2,3,0,1}};
                for (int i = 0; i < 6; i++)
                {
                    vec3g<T> b1 = point[pid[i][0]]-point[pid[i][2]];
                    vec3g<T> b2 = point[pid[i][1]]-point[pid[i][0]];
                    vec3g<T> b3 = point[pid[i][3]]-point[pid[i][1]];
                    vec3g<T> b2b3 = cross(b2,b3);
                    dihedral_angles[i] = fabs(atan2( b2.length() * dot(b1, b2b3), dot(cross(b1,b2), b2b3) ));
                }
            }
            return true;
    }

    template <class T> inline bool tetra_incenter_squaredradius(
        const vec3g<T>& p, const vec3g<T>& q, 
        const vec3g<T>& r, const vec3g<T>& s,
        vec3g<T>& t_incenter, T& squared_radius) {

            const vec3g<T> qp = q-p;
            const vec3g<T> rp = r-p;
            const vec3g<T> sp = s-p;
            const vec3g<T> rq = r-q;
            const vec3g<T> sq = s-q;
			vec3g<T> N0 = cross(sq, rq);
			vec3g<T> N1 = cross(rp, sp);
			vec3g<T> N2 = cross(sp, qp);
			vec3g<T> N3 = cross(qp, rp);
			T l0 = N0.length();
			if (l0 == 0)
			{
				return false;
			}
			T l1 = N1.length();
			if (l1 == 0)
			{
				return false;
			}
			T l2 = N2.length();
			if (l2 == 0)
			{
				return false;
			}
			T l3 = N3.length();
			if (l3 == 0)
			{
				return false;
			}
			N0 /= l0; N1 /= l1; N2 /= l2; N3 /= l3;

			T b00, b01, b02 ;
			T b10, b11, b12 ;
			T b20, b21, b22 ;
			const vec3g<T> N10 = N1 - N0;
			// Note: b is transposed
			comatrix3x3(
				N10.x, N10.y, N10.z,
				N2.x - N0.x, N2.y - N0.y, N2.z - N0.z,
				N3.x - N0.x, N3.y - N0.y, N3.z - N0.z,
				b00, b10, b20,
				b01, b11, b21,
				b02, b12, b22
				) ;


			const T det3 = N10.x*b00 + N10.y*b10 + N10.z*b20;
			const T val = -dot(N3, sp);
			const vec3g<T> sol = (val/det3) * vec3g<T>(b02, b12, b22);
			t_incenter = s + sol;
			squared_radius = dot(N0, sol);
			squared_radius *= squared_radius;
            return true;
    }

	//tetrahedron in regular triangulation
	template <class T> inline vec3g<T> tetra_circumcenter(
		const vec3g<T>& p, const vec3g<T>& q, 
		const vec3g<T>& r, const vec3g<T>& s,
		const T& pw, const T& qw, const T& rw, const T& sw) {
			return three_planes_intersection(
					median_plane(p, q, pw, qw), 
					median_plane(p, r, pw, rw), 
					median_plane(p, s, pw, sw));
	}

    template <class T> inline T tetra_signed_volume(
        const vec3g<T>& p, const vec3g<T>& q, 
        const vec3g<T>& r, const vec3g<T>& s
    ) {
        vec3g<T> U = q-p ;
        vec3g<T> V = r-p ;
        vec3g<T> W = s-p ;
        
        return mixed_product(U,V,W) / T(6) ;
    }

    template <class T> inline vec3g<T> tetra_centroid(
        const vec3g<T>& p, const vec3g<T>& q, 
        const vec3g<T>& r, const vec3g<T>& s
    ) {
        return T(1)/T(4)*(p+q+r+s);
    }

    template <class T> inline bool ray_plane_intersects(
        const vec3g<T>& p0, const vec3g<T>& dir, const Plane<T>& P
        ) {
            T dside = P.side(p0);
            bool side = dside > 0.0;
            bool sign = dot(dir, P.normal()) > 0.0;
            if (dside == 0 || side == !sign )
            {
                return true;
            }
            return false;
    }

    template <class T> inline vec3g<T> ray_plane_intersection(
        const vec3g<T>& p0, const vec3g<T>& dir, const Plane<T>& P
    ) {
        //precondition: ray_plane_intersects(p0, dir, P) == true
        const vec3g<T>& N= P.normal() ;
        if (dot(dir, N) == 0)
            return p0;
        T t = gx_max(0.0, - (P.d + dot(p0,N)) / dot(dir,N));
        return p0 + t * dir ;
    }

    template <class T> inline bool seg_plane_intersects(
        const vec3g<T>& p1, const vec3g<T>& p2, const Plane<T>& P
        ) {
            return (P.side(p1) * P.side(p2)) <= T(0) ;
    }

    template <class T> inline vec3g<T> seg_plane_intersection(
        const vec3g<T>& p1, const vec3g<T>& p2, const Plane<T>& P
    ) {
        //precondition: seg_plane_intersects(p1, p2, P) == true
        T h = fabs(P.side(p1));
        T l = fabs(P.side(p2));
        T hl = h + l;
        if (hl > 0) {
            return (l / hl) * p1 + (h / hl) * p2;
        } else {
            return T(1)/T(2)*(p1+p2);
        }
    }

    template <class T> inline vec3g<T> three_planes_intersection(
        const Plane<T>& P1, const Plane<T>& P2, const Plane<T>& P3
    ) {
        T b00, b01, b02 ;
        T b10, b11, b12 ;
        T b20, b21, b22 ;
        // Note: b is transposed
        comatrix3x3(
            P1.a, P1.b, P1.c,
            P2.a, P2.b, P2.c,
            P3.a, P3.b, P3.c,
            b00, b10, b20,
            b01, b11, b21,
            b02, b12, b22
        ) ;

      
         T det3 = P1.a*b00 + P1.b*b10 + P1.c*b20;

        return vec3g<T>(
            -(P1.d * b00 + P2.d * b01 + P3.d * b02)/det3,
            -(P1.d * b10 + P2.d * b11 + P3.d * b12)/det3,
            -(P1.d * b20 + P2.d * b21 + P3.d * b22)/det3
        );

  
        //return -T(1) / det3x3(
        //    P1.a, P1.b, P1.c,
        //    P2.a, P2.b, P2.c,
        //    P3.a, P3.b, P3.c
        //) * vec3g<T>(
        //    P1.d * b00 + P2.d * b01 + P3.d * b02,
        //    P1.d * b10 + P2.d * b11 + P3.d * b12,
        //    P1.d * b20 + P2.d * b21 + P3.d * b22
        //) ;
    }


    template <class T> inline Plane<T> median_plane(const vec3g<T>& p1, const vec3g<T>& p2) {
        return Plane<T>( T(1)/T(2)*(p1 + p2), p2 - p1) ;
    }

	//bisector in power diagram
	template <class T> inline Plane<T> median_plane(const vec3g<T>& p1, const vec3g<T>& p2,
		const T& w1, const T& w2) {
			return Geex::Plane<T>(
				p1.x - p2.x, 
				p1.y - p2.y, 
				p1.z - p2.z,  
				T(1)/T(2) * ( p2.length2() - p1.length2() + w1 - w2)
				) ;
	}

	template <class T> inline vec3g<T> tri_circumcenter(
		const vec3g<T>& p, const vec3g<T>& q, const vec3g<T>& r) {
			const vec3g<T> N = cross(p-q, p-r);
			return three_planes_intersection(
				median_plane(p, q), median_plane(p, r),
				Plane<T>(p, N));

	}
    
    template <class T> inline void sample_tri(
        const vec3g<T>& p1, const vec3g<T>& p2, const vec3g<T>& p3, 
        int nb, std::vector<vec3g<T> >& pts
    ) {
        T s, t;
        for(int i=0; i<nb; ++i) {
            // Greg Turk's method 2
            s = (T)Numeric::random_float64() ;
            t = (T)Numeric::random_float64() ;
            if (s + t > 1)
            {
                s = 1 - s;
                t = 1 - t;
            }
            pts.push_back(vec3g<T>((1-s-t)*p1 + s*p2 + t*p3)) ;
        }
    }

    template <class T> inline void sample_tetra(
        const vec3g<T>& p1, const vec3g<T>& p2, const vec3g<T>& p3, const vec3g<T>& p4,
        int nb, std::vector<vec3g<T> >& pts
    ) {
        T s, t, u;
        for(int i=0; i<nb; ++i) {
            T s = (T)Numeric::random_float64() ;
            T t = (T)Numeric::random_float64() ;
            T u = (T)Numeric::random_float64() ;
            if (s+t>1.0) { // cut'n fold the cube into a prism
                s = 1.0 - s ;
                t = 1.0 - t ;
            }
            if (t+u>1.0) { // cut'n fold the prism into a tetrahedron
                T tmp = u ;
                u = 1.0 - s - t ;
                t = 1.0 - tmp ;
            } else if (s+t+u>1.0) {
                T tmp = u ;
                u = s + t + u - 1.0 ;
                s = 1 - t - tmp ;
            }
            T a=1-s-t-u ; // a,s,t,u are the barycentric coordinates of the random point.
            pts.push_back(vec3g<T>(a*p1 + s*p2 + t*p3 + u*p4)) ;
        }
    }

	template <class T> inline bool same_point(const vec3g<T>& p0, const vec3g<T>& p1) {
		if(distance(p0, p1)==0.0)
			return true ;
		return false ; 
	}

}

#endif
