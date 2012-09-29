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

#ifndef __DELAUNAY__
#define __DELAUNAY__

//#define NOMINMAX 

//#include <queue>

#include "delaunay_cgal.h"
#include "meshes.h"
#include "rvd_ccvt.h"
#include <Geex/CVT/topo_poly_mesh.h>
#include <Geex/graphics/opengl.h>


namespace Geex {

    class DelaunayGraphics ;
    class DelaunayCVT ;
	class RestrictedVoronoiDiagram ;

	class MyPoint
	{
	public:
		MyPoint()
		:id(-1)
		{

		}

	public:
		Point pt;
		int id;
	};

    //______________________________________________________________________________________
    class MyDelaunay : public DelaunayBase{
        typedef DelaunayBase baseclass ;

    public:
        MyDelaunay() ;
        virtual ~MyDelaunay() ;

        GLboolean& is_clip_volume() { return clip_volume_ ; }
        GLboolean& is_clip_feature() { return clip_feature_ ; }
        GLenum&    clip_mode() { return rvd_->rvd_mode() ; }
		RestrictedVoronoiDiagram* rvd() { return rvd_ ; }
		void partition_mesh() ;

        void save(const std::string& filename) ;
        void load(const std::string& filename) ;
        void load_boundary(const std::string& filename, double perturb=0.0, bool disable_normalize_=false) ;
        void set_non_convex_mode(bool x) { non_convex_mode_ = x ; }
        void save_clipped_mesh(const std::string& filename) ;
        void save_remesh(const std::string& filename) ;
		
        void save_normalized(const std::string& filename) ;
        void load_remesh(const std::string& filename) ;

        bool load_tetmesh(const std::string& filename) ;
        void save_tetmesh(const std::string& filename) ;

        void export_eps(const std::string& file_name) ;

        void get_bbox(
            real& x_min, real& y_min, real& z_min,
            real& x_max, real& y_max, real& z_max
            ) ;

        GLfloat& voronoi_fitering_angle() {return voronoi_fitering_angle_; }

        //bool in_boundary(const vec3& p) { 
        //    return non_convex_mode_ ? boundary_plucker_.contains(p) : boundary_planes_.contains(p) ; 
        //}
		TriMesh& boundary() { return boundary_ ; }

        // ------------------------------------ MyDelaunay 

        int nb_vertices() const { return (int)all_vertices_.size() ; }
		int nb_corners()  const { return (int)boundary_.corners_.size() ; }
        void clear() ;
        void begin_insert() ;
        Vertex_handle insert(const vec3& p) ;
        void end_insert(bool redraw = true) ;

        void insert_corners() ;
        void insert_vertices(int nb, int vertex_type=1) ;
        void insert_vertices(std::string& ptsfilename) ;
		void set_vertices(int N, const double* x, bool tag = false) ;
		void set_vertices(std::vector<vec3>& new_points, bool tag = false) ;
		//new insertion code
		void prepare_sort_points(int N, const double *x, std::vector<MyPoint>& mypointvec);
		void prepare_sort_points(const std::vector<vec3>& new_points, std::vector<MyPoint>& mypointvec);
		void insert_sort_points(const std::vector<MyPoint>& mypointvec, bool tag);

        void insert_boundary_vertices(int nb) ;
        void insert_random_vertex() ;
        void insert_random_vertices(int nb) ;
        void insert_inner_vertices(int nb) ;
        void insert_random_vertices_in_box(int nb) ;
		/*insert_facet_vertices is used to test the convergence rate of newton lloyd vs lloyd*/
		void insert_facet_vertices(int nb) ;
        void sample_boundary(int nb) ;

        inline std::vector<Vertex_handle>& vertices() { return all_vertices_ ; }
		inline std::vector<std::vector<Edge> >& stars() { return stars_ ; } 
		inline std::vector<std::vector<Geex::Plane<real> > >& cells() { return cells_ ; } 

        // ------------------------------------ MyDelaunay combinatorics ----------------

        //        vec3 dual(Cell_handle c) { return c->dual ; }

        void get_incident_cells(Vertex_handle v, std::set<Cell_handle>& c) ;

        /** returns false if this contains infinite or cropped cells */
        enum CellStatus { CS_INTERIOR, CS_CROSS, CS_INFINITE } ;
        CellStatus get_incident_edges(Vertex_handle v, std::vector<Edge>& E) ;
        const std::vector<Edge>& get_incident_edges(Vertex_handle v) {
            if(stars_dirty_) { const_cast<MyDelaunay*>(this)->update_stars() ; }
            gx_debug_assert(v->idx() >= 0 && v->idx() < stars_.size()) ;
            return stars_[v->idx()] ;
        }
        int dual_facet_degree(Edge e) ;

        Cell_handle edge_key(Edge e) ;

        /** make sure that the same edge is always represented by the same tetrahedron */
        void normalize_oriented_edge(Edge& e) ;

        // ----------------------------------- Geometry --------------------------------

        inline Geex::Plane<real> get_dual_plane(Edge e) {
            gx_assert(!is_infinite(e.first->vertex(e.third))) ;
            vec3 p1 = to_geex(e.first->vertex(e.second)->point()) ;
            vec3 p2 = to_geex(e.first->vertex(e.third)->point()) ;
            vec3 N = p2 - p1 ; N = - normalize(N) ;
            return Geex::Plane<real>(0.5*(p1 + p2), N) ;
        }

        TriMesh* dual_convex_clip(std::vector<Edge>& edges, bool close = true) ;   

        //--------------------------------------------------------------------------------

        void save_delaunay(std::string& filename) ;
        void load_delaunay(std::string& filename) ;
        bool load_density(std::string& filename, real gamma=1.0) ; 
        bool load_feature(std::string& filename) ; 

        inline vec3 project_to_mesh(const vec3& pt, vec3& norm) { return boundary_.project_to_mesh(pt, norm) ; } 

        void voronoi_clip_timing(int max_seed_num=5000000) ;

		void compute_rvd() ;

		void perturb_boundary_vertices() ; 
        void extract_mesh() ;
		void extract_dual() ;
        void voronoi_filtering();
        void alpha_shape();

		void compute_lfs() ;

        //void init_vcells() ; // initialize Voronoi cells
        //void build_kdtree() ;

        //// linear voronoi clip by propagation
        //void voronoi_clip_linear() ;
        //void clip_surface_linear() ;
        //void clip_tri_by_cell(int fidx, int cell_idx, std::vector<Geex::Facet>& from) ;

        //// voronoi clip by kd-tree query and using incident cells
        //void clip_surface_mlogn() ;

        //// voronoi clip by kd-tree query and triangle splitting
        //void voronoi_clip_splitting() ;
        //void clip_surface_splitting() ;
        //void subdivide_tri(const Geex::Facet& F, int fidx) ;
        //void subdivide_mesh(VoroCell& cell, Geex::TriMesh& mesh, int fidx) ;

        //// voronoi clip by kd-tree query and using incident cells
        //void check_clip() ;

        //void voronoi_clip() ;
        //void clip_surface() ;
        ////void subdivide_tri(int fidx) ;
        //void convex_clip_tri(int fidx, int vidx) ;
        //void clip_tri_by_plane(Geex::Facet& F, const Geex::Plane<real>& P, int oppvidx, std::vector<Geex::Facet>& inlist, std::vector<Geex::Facet>& outlist) ;

        ////void voronoi_grow_clip() ;
        ////void compute_cell_voronoi(std::pair<int, int>& clippair, std::queue<std::pair<int, int> >& clip_queue) ;
        ////void clip_face(int cellidx, int fidx) ;

        //// meshing and other functions
        //void add_clip_point(vec3& pt, int fidx, Vertex_handle& vh, Edge& e0, Edge& e1) ;

        //void extract_tetmesh() ;
        //double weight(vec3& p) ; // user defined weighting function
        //void compute_sizing_field(real graduation=0.5) ;
        ////        void set_vertex_tex() ;

        ////------------------------ out-of-core Voronoi computation ---------------------------------------------
        //void convex_clip_tri(Geex::Facet& F, int fidx, int vidx, std::set<int> face_cells, std::ofstream& ofile) ;
        //void stream_clip(std::string fname, std::string ptsfname) ;

        ////------------------------ feature edge clipping algorithm ----------------------
        //void clip_edge() ;
        //void clip_edge_by_cell(int eidx, int vidx,  std::set<int> edge_cells) ;
        //void clip_edge_by_plane(Geex::LineSeg& E, const Geex::Plane<real>& P, int oppvidx, std::vector<Geex::LineSeg>& inlist) ;

        void refine_feature() ; // insert new seeds if one seed occupied more than 1 feature lineseg

        ////------------------------ corner clipping algorithm ----------------------
        //void clip_corner() ; // used for identify corners after load data
    protected:
        void update_stars() ;

    protected:
        std::vector<Vertex_handle> all_vertices_ ;
//		std::vector<Vertex_handle> corner_vertices_ ;
        std::vector<std::vector<Edge> > stars_ ;
		std::vector<std::vector<Geex::Plane<real> > > cells_ ;
        bool stars_dirty_ ;

		std::vector<MyPoint> mypointvec;

		static int faces_to_vertices_[4][3] ;
        static int vertices_to_edge_[4][4] ;
        static int edge_to_vertices_[6][2] ;
        static int facet_to_edges_[4][3] ;
        static int edge_to_faces_[6][2] ;
        //static int configs[8][3] ;
        //static int configv[16][4] ;

        //TriMesh remesh_ ;
        //TriMesh clipped_mesh_ ;
        TriMesh boundary_ ;
        TriMesh boundary_save_ ;
        Convex boundary_planes_ ;
//        PluckerMesh boundary_plucker_ ;

		TriMesh ping_ ;
        TriMesh pong_ ;

        TopoPolyMesh boundary_topo_ ;
        Array1d<TopoPolyMesh> boundary_topo_parts_ ;

        bool non_convex_mode_ ;
        GLboolean clip_volume_ ;
        GLboolean clip_feature_ ;
        GLenum    clip_mode_ ;
        GLfloat voronoi_fitering_angle_;

        ANNKdTree_3           *kdtree_ ;
//        std::vector<VoroCell>  vcells_ ;

//        PolyMesh boundary_poly_ ;
//        std::vector< Plane<real> > boundary_poly_plane_ ;
//        PolyMesh ping_poly_ ;
//        PolyMesh pong_poly_ ;
//
//        // TopoTriMesh boundary_topo_ ;
//#ifdef POLY_RVD
//        TopoPolyMesh boundary_topo_ ;
//#else
//        TopoTriMesh boundary_topo_ ;
//#endif
//
//        vec3 normalize_g_ ;
//        double normalize_r_ ;

		RestrictedVoronoiDiagram*  rvd_ ;
		std::vector<Geex::TriMesh> parts_ ; // submeshes for parallel rvd computation

        friend class DelaunayGraphics ;
        friend class DelaunayCVT ;
		friend class RestrictedVoronoiDiagram ;
    } ;

    //______________________________________________________________________________________

}

#endif
