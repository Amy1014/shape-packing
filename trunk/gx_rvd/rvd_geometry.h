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

#ifndef __RVD_GEOMETRY__
#define __RVD_GEOMETRY__

#include <Geex/CVT/delaunay.h>
#include <Geex/CVT/topo_poly_mesh.h>
#include <Geex/CVT/RVD.h>
//#include "custom_rvd.h"
#include <GL/gl.h>
#include <map>

namespace Geex {

// Geometry functions, copied from gx_ccvt

    template <class T> inline T tri_area(
        const vec3g<T>& p, const vec3g<T>& q, 
        const vec3g<T>& r
    ) {
        return T(1) / T(2) * length(cross(q-p, r-p)) ;
    }

    template <class T> inline vec3g<T> tri_centroid(
        const vec3g<T>& p, const vec3g<T>& q, 
        const vec3g<T>& r
    ) {
        return T(1) / T(3) * (p+q+r) ;
    }
    
    template <class T> inline T tri_mass(
        const vec3g<T>& p, const vec3g<T>& q, const vec3g<T>& r,
        T a, T b, T c // a, b, c are density at p, q, r
    ) {
        return tri_area(p, q, r)/T(3)*(sqrt(fabs(a))+sqrt(fabs(b))+sqrt(fabs(c))) ;
    }

    template <class T> inline void tri_centroid(
        const vec3g<T>& p, const vec3g<T>& q, const vec3g<T>& r,
        T a, T b, T c, vec3g<T>& Vg, T& V // a, b, c are density at p, q, r
    ) {
        T abc = a+b+c ; 
        T area = tri_area(p, q, r) ;
        V = area/T(3)*abc ;
        Vg = (area/T(12))*((a+abc)*p+(b+abc)*q+(c+abc)*r) ;
    }

    class RVDGeometry {
    public:
        RVDGeometry(const std::string& del_algo="CGAL_regular") ;
        ~RVDGeometry() ;
        
        void load_boundary(const std::string& filename, double gamma=1.0) ;
        void init_vertices(unsigned int nb = 0) ;
        void init_vertices(const std::string& filename) ;

        void get_bbox(
            real& x_min, real& y_min, real& z_min,
            real& x_max, real& y_max, real& z_max
        ) ;

        const std::string& delaunay_algo() { return del_algo_ ; }
        Delaunay* delaunay() { return delaunay_ ; }
        TopoPolyMesh* boundary() { return boundary_ ; }
        RestrictedVoronoiDiagram_poly& RVD() { 
            RVD_.set_exact(exact_) ;
            return RVD_ ; 
        }

        const RestrictedVoronoiDiagram_poly& RVD() const { 
            const_cast<RVDGeometry*>(this)->RVD_.set_exact(exact_) ;
            return RVD_ ; 
        }
        
	unsigned int size_of_rvd() { return RVD_.size() ; }
	unsigned int triangle_id() { return RVD_.triangle_id() ; }
        std::vector<vec3>& vertices() { return vertices_ ; }
        std::vector<double>& weights() { return weights_ ; }
        std::vector<unsigned int>& facets() { return facets_ ; }
        std::vector<double>& energy() { return energy_ ; }
        std::vector<vec4>& planes() { return planes_ ; }
	GLfloat& lloyd_factor() {return lloyd_factor_ ;}
        void update(bool redraw = false, const double* data=nil) ; 
        GLboolean& exact() { return exact_ ; }

        void save_primal(const std::string& filename = "out.obj") ;
        void load_points(const std::string& filename = "out.obj") ;

        // ---------- Anisotropy

        void load_anisotropy(const std::string& filename) ;

        void get_anisotropy(const vec3& p, vec3& U, vec3& V, vec3& W) const {
            unsigned int i = anisotropy_delaunay_->nearest_vertex_id(p) ;
            U = U_[i] ; V = V_[i] ; W = W_[i] ;
        }

        void get_anisotropy(const vec3& p, vec3& U, vec3& V, vec3& W, vec3& K) const {
            unsigned int i = anisotropy_delaunay_->nearest_vertex_id(p) ;
            U = U_[i] ; V = V_[i] ; W = W_[i] ; 
            K = K_[i];
        }

        void get_anisotropy(unsigned int f, vec3& U, vec3& V, vec3& W) const {
            unsigned int i = f2v_[f] ;
            U = U_[i] ; V = V_[i] ; W = W_[i] ;
        }

        const std::vector<vec3>& get_anisotropy_array() {
            return anisotropy_;
        }


        GLenum& adt_mode() {
            return aniso_dt_method_;
        }

        // ---------- Other geometric operations
        void spread_points(std::vector<int> indices) ;

        // Returns false if there where no duplicated vertices
        bool remove_duplicated_vertices() ;

        // Isolated vertices are the vertices that have an empty
        // Restricted Vornoi cell.
        // Returns false if there where no isolated vertices
        bool remove_isolated_vertices() ;

        vec3 boundary_normal(unsigned int vertex) ;

	void facet_barycentric_coords(
	  unsigned int f, 
	  const vec3& point, 
	  std::vector<double>& coords) const ;

	void facet_barycentric_coords_diff(
	  unsigned int f, 
	  const vec3& point, 
	  std::vector<double>& coords,
	  std::vector<vec3>& gradients) const ;

	vec3 facet_point_normal(unsigned int f, const vec3& point) ;
        void nearest_on_boundary(
	  std::vector<double>& new_vertices, 
	  std::vector<unsigned int>& new_faces,
	  std::multimap< double,
			 std::pair<vec3,unsigned int>,
			 std::greater<double> > *candidates = NULL
	) ;
        void project_on_boundary(
	  std::multimap< double,
			 std::pair<vec3,unsigned int>,
			 std::greater<double> > *candidates = NULL
	) ;
	void compute_vertex_normals() ;
        void split_non_manifold_primal_triangles() ;
	void setup_vertices_with_information() ;

    private:
        std::string del_algo_ ;
        std::vector<vec3> vertices_ ;
        std::vector<double> weights_ ;
        std::vector<unsigned int> facets_ ;
	std::vector<vec3> boundary_normals_ ;
        std::vector<double> energy_ ;
        std::vector<vec4> planes_ ;
        Delaunay* delaunay_ ;
        TopoPolyMesh* boundary_ ;
        RestrictedVoronoiDiagram_poly RVD_ ;
        GLboolean exact_ ;

        std::vector<vec3> anisotropy_vertices_ ;
        Delaunay* anisotropy_delaunay_ ;
        std::vector<vec3> anisotropy_ ;
        std::vector<vec3> U_ ;
        std::vector<vec3> V_ ;
        std::vector<vec3> W_ ;
        std::vector<vec3> K_ ;
        std::vector<unsigned int> f2v_ ;

        GLenum aniso_dt_method_;

	GLfloat lloyd_factor_;
    } ;


} 

#endif
