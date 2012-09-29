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

#include "rvd_geometry.h"
#include "rvd_primal.h"
#include <glut_viewer/glut_viewer.h>
#include <Geex/basics/line_stream.h>
#include <Geex/basics/types.h>

#include <Geex/combinatorics/map.h>
#include <Geex/combinatorics/map_properties.h>
#include <Geex/combinatorics/map_io.h>
#include <Geex/combinatorics/map_geometry.h>
#include <Geex/basics/line_stream.h>

#include <fstream>

namespace Geex {

    RVDGeometry::RVDGeometry(
        const std::string& del_algo
    ) : del_algo_(del_algo),
        delaunay_(Delaunay::create(del_algo)),
        boundary_(new TopoPolyMesh),
        RVD_(delaunay_, boundary_),
        exact_(GL_FALSE),
	lloyd_factor_(0.1),
        anisotropy_delaunay_(Delaunay::create("CGAL"))  // This one is just used to query nearest point, 
                                                        //to read anisotropy from surface.
    {
    }

    RVDGeometry::~RVDGeometry() {
        delete boundary_ ;
        delete delaunay_ ;
        delete anisotropy_delaunay_ ;
    }
    
    void RVDGeometry::load_boundary(const std::string& filename, double gamma) {
        boundary_->load(filename, gamma) ;
    }

    struct GetTriangles {
        GetTriangles (std::vector<const TopoPolyVertexEdge*>& tri_in) : tri(tri_in) {}
        void operator()(
            const TopoPolyVertexEdge& p1, 
            const TopoPolyVertexEdge& p2, 
            const TopoPolyVertexEdge& p3
        ) const {
            tri.push_back(&p1) ; tri.push_back(&p2) ; tri.push_back(&p3) ;
        }
    private:
        std::vector<const TopoPolyVertexEdge*>& tri ;
    } ;

    static vec3 random_point(
        const TopoPolyVertexEdge& p1,
        const TopoPolyVertexEdge& p2,
        const TopoPolyVertexEdge& p3
    ) {
        double s = Numeric::random_float64() ;
        double t = Numeric::random_float64() ;
        //use Turk's second method
        if (s + t > 1)
            {
                s = 1 - s;
                t = 1 - t;
            }
        return (1 - s - t) * p1 + s * p2 + t * p3;
    }

    static vec3 random_point(
        const std::vector<const TopoPolyVertexEdge*>& triangles,
        const std::vector<double>& tribound
    ) {
        double u = Numeric::random_float64() ;
        unsigned int i = 0 ;
        while(tribound[i+1] < u && i+1 < tribound.size() - 1) {
            i++ ;
        }
        return random_point(*triangles[3*i], *triangles[3*i+1], *triangles[3*i+2]) ;
    }

    static unsigned int random_facet(
        const std::vector<double>& tribound
    ) {
        double u = Numeric::random_float64() ;
        unsigned int i = 0 ;
        while(tribound[i+1] < u && i+1 < tribound.size() - 1) {
            i++ ;
        }
        return i ;
    }

    void RVDGeometry::setup_vertices_with_information()	{
        if(delaunay_algo() == "CGAL_regular" || delaunay_algo() == "CGAL_regular_sort") {
            unsigned int n = (unsigned int)vertices_.size() ;
            if (weights_.size() != n)
            {
              weights_.resize(n);
	      memset(&weights_[0], 0, sizeof(double)*weights_.size());
	      //for(int i=0; i<n; ++i)
	      //  weights_[n-i-1] = Numeric::random_float64();
            }
	    delaunay()->set_weight(n,&weights_[0]) ;
        }
    }

    void RVDGeometry::spread_points(std::vector<int> indices) {
      std::vector<const TopoPolyVertexEdge*> triangles ;
      boundary_->for_each_triangle(GetTriangles(triangles)) ;
      unsigned int ntri = (unsigned int)triangles.size() / 3 ;
      std::vector<double> trimass(ntri) ;
      double totalmass = 0.0 ;
      for(unsigned int i=0; i<ntri; i++) {
          trimass[i] = tri_mass(
              *triangles[3*i], *triangles[3*i+1], *triangles[3*i+2],
              triangles[3*i]->w, triangles[3*i+1]->w, triangles[3*i+2]->w
          ) ;
          totalmass += trimass[i] ;
      }
      std::vector<double> tribound(ntri + 1) ;
      tribound[0] = 0.0 ;
      for(unsigned int i=1; i<=ntri; i++) {
          tribound[i] = tribound[i - 1] + (trimass[i-1] / totalmass) ;
      }
      for(unsigned int i = 0; i<indices.size(); ++i)
	vertices_[indices[i]] = random_point(triangles, tribound) ;
    }
    
    void RVDGeometry::init_vertices(unsigned int nb) {
        if(nb == 0) { nb = boundary_->nb_facets(); }
        nb = gx_max(nb, (unsigned int)(4)) ;
        vertices_.resize(nb) ;

        planes_.resize(nb,vec4(0.,0.,0.,0.)) ;

        if(true) {
            std::cerr << "Sampling with " << nb << " points" << std::endl ;
            std::vector<const TopoPolyVertexEdge*> triangles ;
            boundary_->for_each_triangle(GetTriangles(triangles)) ;
            unsigned int ntri = (unsigned int)triangles.size() / 3 ;
            std::vector<double> trimass(ntri) ;
            double totalmass = 0.0 ;
            for(unsigned int i=0; i<ntri; i++) {
                trimass[i] = tri_mass(
                    *triangles[3*i], *triangles[3*i+1], *triangles[3*i+2],
                    triangles[3*i]->w, triangles[3*i+1]->w, triangles[3*i+2]->w
                ) ;
                totalmass += trimass[i] ;
            }
            std::vector<double> tribound(ntri + 1) ;
            tribound[0] = 0.0 ;
            for(unsigned int i=1; i<=ntri; i++) {
                tribound[i] = tribound[i - 1] + (trimass[i-1] / totalmass) ;
            }
            for(unsigned int i=0; i<vertices_.size(); i++) {
                vertices_[i] = random_point(triangles, tribound) ;
            }
        } else {
            nb = gx_min(nb, boundary_->nb_facets()) ;
            nb = gx_max(nb, (unsigned int)(4)) ;
            vertices_.resize(nb) ;
            std::vector<unsigned int> mrandom_num(boundary_->nb_facets());
            for (unsigned int i = 0; i < boundary_->nb_facets(); i++) {
                mrandom_num[i] = i;
            }
            std::random_shuffle(mrandom_num.begin(),mrandom_num.end());
            for(unsigned int f=0; f<nb; f++) {
                vertices_[f] = boundary_->facet_center(mrandom_num[f]) ;
            }
        }

        setup_vertices_with_information();

	std::cout << "sending vertices to the delaunay structure ..." << std::endl ;
        delaunay_->set_vertices(vertices_) ;
	std::cout << "done" << std::endl ;
        energy_.resize(vertices_.size()) ;
    }

    void RVDGeometry::init_vertices(const std::string& filename) {
        std::ifstream input(filename.c_str()) ;
        if(!input) {
            std::cerr << "could not open file." << std::endl ;
            return ;
        }
        vertices_.clear() ;
        LineInputStream in(input) ;
        while(!in.eof()) {
            in.get_line() ;
            std::string keyword ;
            in >> keyword ;
            if(keyword == "v") {
                vec3 p ;
                in >> p ;
                vertices_.push_back(p) ;
            }
        }  
        if(delaunay_algo() == "CGAL_aniso" || delaunay_algo() == "CGAL_aniso_sort") {
            unsigned int n = vertices_.size() ;
            anisotropy_.resize(3*n) ;
            for(unsigned int i=0; i<n; i++) {
                get_anisotropy(
                    vertices_[i],
                    anisotropy_[3*i],
                    anisotropy_[3*i+1],
                    anisotropy_[3*i+2]
                ) ;                
            }
            delaunay_->set_anisotropy(n, &(anisotropy_[0].x), aniso_dt_method_) ;
        }
        delaunay_->set_vertices(vertices_) ;
        energy_.resize(vertices_.size()) ;
    }

    void RVDGeometry::get_bbox(
        real& x_min, real& y_min, real& z_min,
        real& x_max, real& y_max, real& z_max
    ) {
        boundary_->get_bbox(x_min, y_min, z_min, x_max, y_max, z_max) ;
    }


    void RVDGeometry::update(bool redraw, const double* data) {
        if(data != 0 && data != &(vertices_[0].x)) {
            for(unsigned int i=0; i<vertices_.size(); i++) {
                vertices_[i].x = data[4*i] ;
                vertices_[i].y = data[4*i+1] ;
                vertices_[i].z = data[4*i+2] ;
		weights_[i] = data[4*i+3] ;
            }
        }
        if(delaunay_algo() == "CGAL_aniso" || delaunay_algo() == "CGAL_aniso_sort") {
            unsigned int n = vertices_.size() ;
            anisotropy_.resize(3*n) ;
            for(unsigned int i=0; i<n; i++) {
                get_anisotropy(
                    vertices_[i],
                    anisotropy_[3*i],
                    anisotropy_[3*i+1],
                    anisotropy_[3*i+2]
                ) ;                
            }
            delaunay_->set_anisotropy(n, &(anisotropy_[0].x), aniso_dt_method_) ;
        }
        delaunay_->set_vertices(vertices_) ;
	delaunay_->set_weight(vertices_.size(),&weights_[0]) ;
        energy_.resize(vertices_.size()) ;
        if(redraw) {
            glut_viewer_redraw() ;
        }
    }


    class SavePrimalTriangle {
    public:
        SavePrimalTriangle(
            std::ofstream& out
        ) : out_(&out) { 
        }
        void operator()(unsigned int i, unsigned int j, unsigned int k) const {
            (*out_) << "f " << i+1 << " " << j+1 << " " << k+1 << std::endl ;
        }

    private:
        std::ofstream* out_ ;
    } ;

    void RVDGeometry::save_primal(const std::string& filename) {
        std::cerr << "Writing primal mesh to " << filename << std::endl ;
        std::ofstream out(filename.c_str()) ;
        for(unsigned int i=0; i<vertices_.size(); i++) {
            out << "v " << vertices_[i] << std::endl ;
        }
        RVD_.for_each_primal_triangle(SavePrimalTriangle(out)) ;
        std::cerr << "Done." << std::endl ;
    }


    void RVDGeometry::load_points(const std::string& filename) {
        vertices_.clear() ;
        std::ifstream input(filename.c_str()) ;
        if(!input) {
            std::cerr << "could not open file." << std::endl ;
            return ;
        }
        LineInputStream in(input) ;
        while(!in.eof()) {
            in.get_line() ;
            std::string keyword ;
            
            in >> keyword ;
            
            if(keyword == "v") {
                vec3 p ;
                in >> p ;
                vertices_.push_back(p) ;
            }
        }        
	weights_.resize(vertices_.size(),0) ;
        update(true) ;
    }

    void RVDGeometry::load_anisotropy(const std::string& filename) {

        std::ifstream in(filename.c_str()) ;
        if(!in) {
            std::cerr << "Could not open anisotropy file: " << filename << std::endl ;
            exit(-1) ;
        }
        Map* map = new Map ;
        load_map(in, map, false, false) ;

        MapVertexProperty<vec3> U ;
        U.bind_if_defined(map, "K1") ;

        MapVertexProperty<vec3> V ;
        V.bind_if_defined(map, "K2") ;

        MapVertexProperty<vec3> W ;
        W.bind_if_defined(map, "N") ;

        bool bound = true ;
        if(!U.is_bound()) {
            std::cerr << "Missing U property, using geometry" << std::endl ;
            bound = false ;
        }
        if(!V.is_bound()) {
            std::cerr << "Missing V property, using geometry" << std::endl ;
            bound = false ;
        }

        anisotropy_vertices_.clear() ;
        U_.clear() ;
        V_.clear() ;
        W_.clear() ;

        FOR_EACH_VERTEX(Map, map, it) {
            anisotropy_vertices_.push_back(it->point()) ;
            vec3 cur_U,cur_V,cur_W ; 
            if(bound) {
                cur_U = normalize(U[it]) ;
                cur_V = normalize(V[it]) ;
                if(false && W.is_bound()) {
                    cur_W = normalize(W[it]) ;
                } else {
                    cur_W = cross(cur_U, cur_V) ;
                }
            } else {
                vec3 X = Geom::vertex_normal(it) ;
                if(Numeric::has_nan(X)) {
                    std::cerr << "Zero normal" << std::endl ;
                    X = vec3(1.0, 0.0, 0.0) ;
                }
                vec3 Y ;
                vec3 Y1(0.0, 1.0, 0.0) ;
                vec3 Y2(1.0, 0.0, 0.0) ;
                if(::fabs(dot(Y1,X)) > ::fabs(dot(Y2,X))) {
                    Y = Y2 ;
                } else {
                    Y = Y1 ; 
                }
                vec3 Z = cross(X,Y) ; 
                Z = normalize(Z) ;
                Y = cross(Z,X) ;
                if(Numeric::has_nan(X) || Numeric::has_nan(Y) || Numeric::has_nan(Z)) {
                    std::cerr << "Nan aniso" << std::endl ;
                    X = vec3(1.0, 0.0, 0.0) ;
                    Y = vec3(0.0, 1.0, 0.0) ;
                    Z = vec3(0.0, 0.0, 1.0) ;
                }
                cur_U = X ; cur_V = Y ; cur_W = Z ;
            }
            U_.push_back(cur_U) ;
            V_.push_back(cur_V) ;
            W_.push_back(cur_W) ;
        }

        U.unbind() ;
        V.unbind() ;
        W.unbind() ;

        delete map ;

        anisotropy_delaunay_->set_vertices(anisotropy_vertices_) ;
        for(unsigned int i=0; i<boundary_->nb_facets(); i++) {
            f2v_.push_back(anisotropy_delaunay_->nearest_vertex_id(boundary_->facet_center(i))) ;
        }
    }


    //===============================================================================================

    inline vec3 projection_on_line(const vec3& p1, const vec3& p2, const vec3& q) {
        vec3 v = p2 - p1 ;
        if(length2(v) < 1e-10) {
            return p1 ;
        }
        v = normalize(v) ;
        return p1 + dot((q-p1),v)*v ;
    }

    inline bool in_segment(const vec3& p1, const vec3& p2, const vec3& h) {
        return (  dot(p2 - h, p1 - h) <= 0 ) ;
    }

    inline vec3 projection_on_plane(const vec3& p1, const vec3& p2, const vec3& p3, const vec3& q) {
        vec3 n = cross(p2-p1, p3-p1) ;
        if(length2(n) < 1e-10) {
            double l1 = length2(p2-p1) ;
            double l2 = length2(p3-p2) ;
            double l3 = length2(p1-p3) ;
            if(l1 > l2 && l1 > l3) {
                return projection_on_line(p1, p2, q) ;
            } else if(l2 > l1 && l2 < l3) {
                return projection_on_line(p2, p3, q) ;
            } else {
                return projection_on_line(p3, p1, q) ;
            }
        }
        n = normalize(n) ;
        return q - dot((q - p1), n) * n ;
    }

    inline bool in_triangle(const vec3& p1, const vec3& p2, const vec3& p3, const vec3& h) {
        vec3 n = cross(p2-p1, p3-p1) ;
        
        double s1 = mixed_product(p1-p2, p1-h, n) ;
        double s2 = mixed_product(p2-p3, p3-h, n) ;
        double s3 = mixed_product(p3-p1, p3-h, n) ;

        return (
            (s1 >= 0 && s2 >= 0 && s3 >= 0) ||
            (s1 <= 0 && s2 <= 0 && s3 <= 0)
        ) ;
    }

    inline vec3 nearest_in_triangle(const vec3& p1, const vec3& p2, const vec3& p3, const vec3& q) {
        vec3 result = projection_on_plane(p1, p2, p3, q) ;
        if(in_triangle(p1, p2, p3, result)) { return result ; }
        
        result = p1 ;
        if(length2(p2 - q) < length2(result - q)) { result = p2 ; }
        if(length2(p3 - q) < length2(result - q)) { result = p3 ; }

        vec3 h = projection_on_line(p1, p2, q) ;
        if(in_segment(p1, p2, h) && length2(h-q) < length2(result - q)) 
	  result = h ; 
        h = projection_on_line(p2, p3, q) ;
        if(in_segment(p2, p3, h) && length2(h-q) < length2(result - q))
          result = h ;
        h = projection_on_line(p3, p1, q) ;
        if(in_segment(p3, p1, h) && length2(h-q) < length2(result - q))
          result = h ;
        
        return result ;
    }

    
    class ProjectOnBoundary {
    public:
        ProjectOnBoundary(
	    const RVDGeometry *geometry,
            const std::vector<vec3>& vertices, 
            std::vector<double>& new_vertices, 
            std::vector<bool>& visited,
	    std::vector<unsigned int>& new_facets
        ) : geometry_(geometry), vertices_(vertices), new_vertices_(new_vertices), 
	visited_(visited), new_facets_(new_facets) {
        }
        void operator() (
            unsigned int v,
            const vec3& v1,
            const vec3& v2,
            const vec3& v3
        ) const {
            vec3 P = nearest_in_triangle(v1, v2, v3, vertices_[v]) ;
	    vec3 current_vertex(new_vertices_[3*v], 
				new_vertices_[3*v+1], 
				new_vertices_[3*v+2]) ;
            if(
                !visited_[v] || (
                    length2(P - vertices_[v]) < 
                    length2(current_vertex - vertices_[v])
                )
            ) {
              visited_[v] = true ;
	      for(int i=0; i<3; ++i)
		new_vertices_[3*v+i] = P[i] ;
	      new_facets_[v] = geometry_->RVD().current_facet();
            }
        }
    private:
	const RVDGeometry *geometry_;
        const std::vector<vec3>& vertices_ ;
        mutable std::vector<double>& new_vertices_ ;
        mutable std::vector<unsigned int>& new_facets_ ;
        mutable std::vector<bool>& visited_ ;
    } ;

    void RVDGeometry::nearest_on_boundary(
      std::vector<double>& new_vertices, 
      std::vector<unsigned int>& new_facets,
      std::multimap< double,
		     std::pair<vec3,unsigned int>,
		     std::greater<double> > *candidates
    ) {
      if(boundary_normals_.size() == 0)
	compute_vertex_normals() ;
      std::vector<bool> visited ;
      new_vertices.resize(3*vertices().size()) ;
      new_facets.resize(vertices().size()) ;
      visited.resize(vertices().size()) ;
      std::fill(visited.begin(), visited.end(), false) ;
      RVD().for_each_triangle( 
        ProjectOnBoundary( this,
          vertices(), new_vertices, 
          visited, new_facets)
      ) ;
      for(int i = 0; i<vertices().size(); ++i) {
	if(!visited[i]) {
	  if(candidates && candidates->size() > 0) {
	    vec3 candidate = candidates->begin()->second.first ;
	    for(int j=0; j<3; ++j)
	      new_vertices[3*i+j] = candidate[j] ;
	    new_facets[i] = candidates->begin()->second.second ;
	    candidates->erase(candidates->begin()) ;
	  } else {
	    std::vector<const TopoPolyVertexEdge*> triangles ;
            boundary_->for_each_triangle(GetTriangles(triangles)) ;
            unsigned int ntri = (unsigned int)triangles.size() / 3 ;
            std::vector<double> trimass(ntri) ;
            double totalmass = 0.0 ;
            for(unsigned int j=0; j<ntri; j++) {
                trimass[j] = tri_mass(
                    *triangles[3*j], *triangles[3*j+1], *triangles[3*j+2],
                    triangles[3*j]->w, triangles[3*j+1]->w, triangles[3*j+2]->w
                ) ;
                totalmass += trimass[j] ;
            }
            std::vector<double> tribound(ntri + 1) ;
            tribound[0] = 0.0 ;
            for(unsigned int j=1; j<=ntri; j++) {
                tribound[j] = tribound[j - 1] + (trimass[j-1] / totalmass) ;
	    }
	    unsigned int facet = random_facet(tribound) ;
	    vec3 rp = random_point(*triangles[3*facet],
	          		 *triangles[3*facet+1],
	          		 *triangles[3*facet+2]) ;
	    for(int j=0; j<3; ++j)
	      new_vertices[3*i+j] = rp[j] ;
	    new_facets[i] = facet ;
	  }
	  //std::cout << "new position for " << i << " [ " ;
	  //for(int j=0; j<3; ++j)
	  //  std::cout << new_vertices[3*i+j] << " " ;
	  //std::cout << "]" << std::endl ;
	}
      }
    }

    void RVDGeometry::project_on_boundary(
      std::multimap< double,
		     std::pair<vec3,unsigned int>,
		     std::greater<double> > *candidates
    ) {
        std::vector<double> new_vertices ;
	nearest_on_boundary(new_vertices, facets_, candidates) ;
        update(false, &new_vertices[0]) ;
    }

    vec3 RVDGeometry::boundary_normal(unsigned int vertex) {
      if(boundary_normals_.size() == 0)
	compute_vertex_normals() ;

      return boundary_normals_[boundary_->vertex_index(vertex)] ;
    }

    void RVDGeometry::compute_vertex_normals() {
      //compute the amount of non duplicated vertices
      unsigned int nb_unique_vertices = 0 ;

      for(unsigned int v = 0 ; v<boundary_->nb_vertices(); ++v)
	nb_unique_vertices = gx_max( nb_unique_vertices, 
				     boundary_->vertex_index(v)) ;
      //resize the normals vector
      boundary_normals_.resize(nb_unique_vertices+1, vec3(0.,0.,0.)) ;

      //compute per vertex normal
      for(unsigned int f = 0 ; f<boundary_->nb_facets(); ++f) {
	for(
	  unsigned int v = boundary_->facet_begin(f) ; 
	  v < boundary_->facet_end(f) ;
	  ++v
	) {
	  boundary_normals_[boundary_->vertex_index(v)] 
	    += boundary_->facet_normal(f) ;
	}
      }
    }

    bool RVDGeometry::remove_duplicated_vertices() {
        unsigned int count = 0 ;
        const DelaunaySkeleton* skl = delaunay_->skeleton() ;
        std::vector<vec3> new_vertices ;
        std::vector<unsigned int> new_facets ;
        bool has_duplicated_vertices = false ;
        for(unsigned int i=0; i<vertices_.size(); i++) {
            if(skl->nb_neighbors(i) == 0) {
                count++ ;
                has_duplicated_vertices = true ;
            } else {
                new_vertices.push_back(vertices_[i]) ;
		if(facets_.size() > 0 )
		  new_facets.push_back(facets_[i]) ; 
            }
        }

        if(has_duplicated_vertices) {
            vertices_ = new_vertices ;
            facets_ = new_facets ;
            update() ;
            std::cerr << "Removed " << count << " duplicated vertices" << std::endl ;
        }

        return has_duplicated_vertices ;
    }

    class VisitVertices {
    public:
        VisitVertices(std::vector<bool>& visited) : visited_(visited) {
        }
        void operator()(unsigned int c, TopoPolyMesh* M) const {
            visited_[c] = true ;
        }        
    private:
        std::vector<bool>& visited_ ;
    } ;

    bool RVDGeometry::remove_isolated_vertices() {
        std::vector<bool> visited(vertices_.size()) ;
        std::fill(visited.begin(), visited.end(), false) ;
        RVD().for_each_facet(VisitVertices(visited)) ;
        std::vector<vec3> new_vertices ;
        std::vector<unsigned int> new_facets ;
        for(unsigned int i=0; i<vertices_.size(); i++) {
            if(visited[i]) { 
	      new_vertices.push_back(vertices_[i]) ; 
	      new_facets.push_back(facets_[i]) ; 
	    }
        }
        if(new_vertices.size() != vertices_.size()) {
            std::cerr << "Removed " << vertices_.size() - new_vertices.size() << " isolated vertices"
                      << std::endl ;
            vertices_ = new_vertices ;
            facets_ = new_facets ;
            update() ;
        } 
        return true;
    }

    
    class SplitPrimalSingularTriangles {
    public:
        SplitPrimalSingularTriangles(
            std::vector<vec3>& seeds_in, std::map<trindex, int>& triangles_in
        ) : seeds_(seeds_in), triangles_(triangles_in) {
        }
        void operator()(unsigned int i, unsigned int j, unsigned int k) const {
            trindex t(i,j,k) ;
            if(triangles_[t] > 1) {
                const vec3& v1 = seeds_[i] ;
                const vec3& v2 = seeds_[j] ;
                const vec3& v3 = seeds_[k] ;
                vec3 N = cross(v2-v1, v3-v1) ;
                vec3 G = (1.0 / 3.0) * (v1 + v2 + v3) ;
                seeds_.push_back(G + 0.01 * N) ;
            }
        }
    private:
        std::vector<vec3>& seeds_ ;
        std::map<trindex, int>& triangles_ ;
    } ;


    void RVDGeometry::split_non_manifold_primal_triangles() {
        std::map<trindex, int> triangles ;        
        RVD().for_each_primal_triangle(CountPrimalSingularTriangles(triangles)) ;
        RVD().for_each_primal_triangle(SplitPrimalSingularTriangles(vertices(), triangles)) ;
        update(false) ;
    }

    void RVDGeometry::facet_barycentric_coords(
      unsigned int f, 
      const vec3& point, 
      std::vector<double>& coords
    ) const {
      //clear the result vector
      coords.clear() ;

      //normalizing coefficient in the end
      double total = 0.0 ;

      //size of the polygon
      unsigned int size = boundary_->facet_size(f) ;

      //iteration over the vertices of the facet
      for ( 
	unsigned int j = boundary_->facet_begin(f) ; 
	j < boundary_->facet_end(f) ; 
	++j
      )  {
	//take the surrounding vertices
        unsigned int i = boundary_->prev_around_facet(f,j) ;
        unsigned int k = boundary_->next_around_facet(f,j) ;

        vec3 pI = boundary_->vertex(i);
        vec3 pJ = boundary_->vertex(j);
        vec3 pK = boundary_->vertex(k);
          
        vec3 vectJI = pI - pJ;
        vec3 vectJP = point - pJ;
        vec3 vectJK = pK - pJ;

	if(size > 3) {
	  double dotIJP = dot(vectJP,vectJI) ;
          vec3 cIJP = cross(vectJP,vectJI) ;
	  double l2IJP = length2(cIJP) ;
	  double crossIJP = sqrt(l2IJP) ;
          double dotPJK = dot(vectJK,vectJP) ;
          vec3 cPJK = cross(vectJK,vectJP) ;
	  double l2PJK = length2(cPJK) ;
	  double crossPJK = sqrt(l2PJK) ;

          // cotangent of angles
          double cotIJP = dotIJP/crossIJP;
          double cotPJK = dotPJK/crossPJK;

          double factor = (cotIJP+cotPJK)/length2(vectJP);
	  //if(Numeric::is_nan(factor)) {
	  //  std::cout << "NAN !! cot IJP: " << cotIJP << ", cot PJK: " << cotPJK << std::endl ;
	  //  std::cout << "       JI: " << vectJI ;
	  //  std::cout << " JK: " << vectJK ;
	  //  std::cout << " JP: " << vectJP << std::endl ;
	  //  std::cout << "       dotIJP: " << dotIJP ;
	  //  std::cout << " crossIJP: " << crossIJP << std::endl ;
	  //  std::cout << "       dotPJK: " << dotPJK ;
	  //  std::cout << " crossPJK: " << crossPJK << std::endl ;
	  //}
          coords.push_back(factor);
          total += factor;
	} else {
	  double factor = length(cross(pK-point,pI-point)) ;
	  coords.push_back(factor);
	  total += factor;
	}
      }

      bool is_valid = true ;
      //normalize the result
      for(unsigned int i=0; i < size; ++i) {
	coords[i] /= total ;
	is_valid = is_valid && coords[i]>=0 ;
      }

      if(!is_valid) {
	/*std::cout << "invalid barycentric coords " ;
	for(unsigned int i=0; i < boundary_->facet_size(f); ++i)
	  std::cout << coords[i]*total << " " ;
	std::cout << std::endl ;
	std::cout << "point : " << point << std::endl ;
	std::cout << "polygon : " ;
	for ( 
      	  unsigned int j = boundary_->facet_begin(f) ; 
      	  j < boundary_->facet_end(f) ; 
      	  ++j
      	)  {
	  std::cout << "[" << boundary_->vertex(j) << "] " ;
	}
	std::cout << std::endl ;*/
      }
    } 

    void RVDGeometry::facet_barycentric_coords_diff(
      unsigned int f, 
      const vec3& point, 
      std::vector<double>& coords,
      std::vector<vec3>& gradients
    ) const {
      //clear the result vector
      coords.clear() ;

      //      /                                        
      //     /                                         
      //    /                                          
      // I /                                           
      //   \                         |                 
      //    \          .P            |                 
      //     \        /              |                 
      //      \      /               |                 
      //       \    /                |                 
      //        \  /                 |                 
      //         \/__________________|                 
      //         J                   K                 

      //normalizing coefficient in the end
      double total = 0.0 ;
      vec3 total_grad(0.,0.,0.) ;

      //size of the polygon
      unsigned int size = boundary_->facet_size(f) ;
      coords.clear() ;
      gradients.clear() ;
    
      vec3 normal = boundary_->facet_normal(f) ;
      normal = normalize(normal) ;
      vec3 coord_gradients[] = {vec3(1.,0.,0.),
				vec3(0.,1.,0.),
				vec3(0.,0.,1.)} ;

      //iteration over the vertices of the facet
      for ( 
	unsigned int j = boundary_->facet_begin(f) ; 
	j < boundary_->facet_end(f) ; 
	++j
      )  {
	//take the surrounding vertices
        unsigned int i = boundary_->prev_around_facet(f,j) ;
        unsigned int k = boundary_->next_around_facet(f,j) ;

        vec3 pI = boundary_->vertex(i);
        vec3 pJ = boundary_->vertex(j);
        vec3 pK = boundary_->vertex(k);
          
        vec3 vectJI = pI - pJ;
        vec3 vectJP = point - pJ;
        vec3 vectJK = pK - pJ;

	if(size>3) {
	  double dotIJP = dot(vectJP,vectJI) ;
          double crossIJP = length(cross(vectJP,vectJI)) ;
          double dotPJK = dot(vectJK,vectJP) ;
          double crossPJK = length(cross(vectJK,vectJP)) ; 

          // cotangent of angles
          double cotIJP = dotIJP/crossIJP;
          double cotPJK = dotPJK/crossPJK;

          double factor = (cotIJP+cotPJK)/length2(vectJP);
          coords.push_back(factor);
          total += factor;

	  //derivatives
	  vec3 tmpdot = 2*(cotIJP+cotPJK)*vectJP ;
	  vec3 tmpcross = (cotIJP*cotIJP-cotPJK*cotPJK)*normal ;
	  double len4 = length2(vectJP) ;
	  len4 *= len4 ;
	  vec3 grad(0.,0.,0.) ;
	  for(int m=0; m<3; ++m){
	    grad[m] = -dot(cross(coord_gradients[m],vectJP),tmpcross) ;
	    grad[m] += dot(coord_gradients[m],tmpdot) ;
	    grad[m] /= len4 ;
	  }
	  total_grad += grad ;
	  gradients.push_back(grad) ;
	} else {
	  //barycentric coordinates
	  vec3 height_vect = normalize(cross(normal,pK-pI)) ;
	  double edge_len = length(pK-pI) ;
	  double factor = edge_len*dot(height_vect,point-pK) ;
	  coords.push_back(factor);
	  total += factor;

	  //derivatives
	  vec3 grad(0.,0.,0.) ;
	  for(int m=0; m<3; ++m)
	    grad[m] = edge_len*dot(height_vect,coord_gradients[m]) ;

	  total_grad += grad ;
	  gradients.push_back(grad) ;
	}

      }

      //normalize the weights
      bool is_valid = true ;
      for(unsigned int i=0; i < size; ++i) {
	gradients[i] = (total*gradients[i] - coords[i]*total_grad) ;
	gradients[i] /= total*total ;
	coords[i] /= total ;
	is_valid = is_valid && coords[i]>=0 ;
      }

      if(!is_valid) {
	//std::cout << "invalid barycentric coords diff " ;
	//for(unsigned int i=0; i < boundary_->facet_size(f); ++i)
	//  std::cout << coords[i] << " " ;
	//std::cout << std::endl ;
	//std::cout << "point : " << point << std::endl ;
	//std::cout << "polygon : " ;
	//for ( 
      	//  unsigned int j = boundary_->facet_begin(f) ; 
      	//  j < boundary_->facet_end(f) ; 
      	//  ++j
      	//)  {
	//  std::cout << "[" << boundary_->vertex(j) << "] " ;
	//}
	//std::cout << std::endl ;
      }
    } 

    vec3 RVDGeometry::facet_point_normal(unsigned int f, const vec3& point) {
      // barycentric coordinates
      std::vector<double> coords ;
      facet_barycentric_coords(f, point, coords) ;

      // interpolate normal
      vec3 normal(0.,0.,0.) ;
      int i=0 ;
      for(
	unsigned int v = boundary_->facet_begin(f) ;
	v < boundary_->facet_end(f) ;
	++v, ++i
      ) {
	double coord = coords[i] ;
	normal += coord*boundary_normal(v) ;
      }
      return normal ;
    }

} 


