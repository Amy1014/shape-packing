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
 *  You should have received a copy of the GNU General Public License<
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

#include "delaunay.h"
#include "rvd_ccvt.h"
#include "geometry_ccvt.h"

#include <Geex/basics/file_system.h>
#include <Geex/basics/processor.h>
#include <Geex/numerics/matrix_partition.h>
#include <Geex/CVT/ann_kdtree.h>
#include <Geex/CVT/CVT_spatial_sort.h>
#include <Geex/CVT/MedialAxis.h>
#include <glut_viewer/glut_viewer.h>

#include <fstream>

//#include "sizing_field.h" 

namespace Geex {

    int MyDelaunay::faces_to_vertices_[4][3] = {
	{1, 2, 3},
	{2, 0, 3}, 
	{3, 0, 1},
	{0, 2, 1}
    } ;
    
    int MyDelaunay::vertices_to_edge_[4][4] = {
        {-1,  0,  2,  5},
        { 0, -1,  1,  3},
        { 2,  1, -1,  4},
        { 5,  3,  4, -1}
    } ;
    
    int MyDelaunay::edge_to_vertices_[6][2] = {
        {0,1},
        {0,2},
        {0,3},
        {1,2},
        {1,3},
        {2,3}
    } ;
    
    int MyDelaunay::edge_to_faces_[6][2] = {
        {2,3},
        {1,3},
        {1,2},
        {0,3},
        {0,2},
        {0,1}
	} ;
    
    int MyDelaunay::facet_to_edges_[4][3] = {
	{5, 4, 3},//{1, 4, 3},
	{2, 1, 5}, //{2, 5, 4},
	{0, 4, 2},//{0, 3, 5},
	{2, 1, 0}//{2, 1, 0}
    } ;
    
    MyDelaunay::MyDelaunay() {
        non_convex_mode_ = false ;
		clip_feature_ = false ;
		clip_volume_ = false ;
		kdtree_ = nil ; 
		stars_dirty_ = true ;
        voronoi_fitering_angle_ = 90;
		rvd_ = new RestrictedVoronoiDiagram() ;
    }

    MyDelaunay::~MyDelaunay() {
		if(kdtree_ != nil)      delete kdtree_;	
		delete rvd_ ;
    }

    void MyDelaunay::save(const std::string& filename) {
        std::ofstream out(filename.c_str()) ;
        for(Vertex_iterator it = vertices_begin() ; it != vertices_end() ; it++) {
            if(is_infinite(it)) { continue ; }
            out << to_geex(it->point()) ;
        }
    }
    
    void MyDelaunay::load(const std::string& filename) {
        std::cerr << "loading " << filename << std::endl ;
        clear() ;
        std::ifstream in(filename.c_str()) ;
        if(!in) {
            std::cerr << "could not open file" << std::endl ;
            return ;
        }
		std::vector<vec3> pts;
		vec3 p ;
		while(in) {
			in >> p ;
			if(in) {  // we need to do this additional check else we got the last point twice !
				pts.push_back(p); 
			}
		} ;

		set_vertices(pts);

        //vec3 p ;
        //while(in) {
        //    in >> p ;
        //    if(in) {  // we need to do this additional check else we got the last point twice !
        //        insert(p) ; 
        //    }
        //} ;
    }

    void MyDelaunay::load_boundary(const std::string& filename, double perturb, bool no_normalize) {
		vec3 g ;
		double R ;
		int len = filename.length() ;
		bool is_obj = filename[len-1]=='j' ; 

		if(is_obj) {
			boundary_.load(filename, &boundary_planes_) ;

			if(perturb>0.0) {
				std::string perturb_fname = filename ;
				int len = filename.length() ;
				perturb_fname[len-3] = 't' ;
				perturb_fname[len-2] = 'r' ;
				perturb_fname[len-1] = 'i' ;
				boundary_.perturb_facet_vertices(perturb) ; 
				boundary_.save_tri(perturb_fname) ;
			}
		} else {
			boundary_.load_tri(filename) ;
		}

        boundary_.normalize(g, R, no_normalize) ;
        boundary_.measure_quality() ; // temp use
		boundary_.build_kdtree() ;

		boundary_save_.clear() ;
        for(unsigned int i=0; i<boundary_.size(); i++) {
            boundary_save_.push_back(boundary_[i]) ;
        }

        //boundary_topo_.clear() ;
        //boundary_topo_.load(filename) ;
        //boundary_topo_.normalize(g, R, no_normalize) ;
        //int nb_threads_ = Geex::Processor::number_of_cores() ;
//		
//        if(nb_threads_ > 1) {
//            std::cerr << "Partitioning into " << nb_threads_ << " parts" << std::endl ;
//            boundary_topo_.partition(nb_threads_, boundary_topo_parts_) ;
//            for(int i=0; i<nb_threads_; i++) {
//                std::cerr << "part " << i << " : nb_facets = " 
//                          << boundary_topo_parts_[i].nb_facets() << std::endl ;
////                boundary_topo_parts_[i].save("part" + String::to_string(i) + ".obj") ;
//			}
//		}
    }

	void MyDelaunay::perturb_boundary_vertices() { 
		boundary_.perturb_facet_vertices() ; 
		boundary_.build_kdtree() ;

		boundary_save_.clear() ;
        for(unsigned int i=0; i<boundary_.size(); i++) {
            boundary_save_.push_back(boundary_[i]) ;
        }

		unsigned nv = nb_vertices() ;
		this->clear() ;
		insert_vertices(nv) ;
		compute_rvd() ;
	}

	void MyDelaunay::save_normalized(const std::string& filename) {
		boundary_.save_obj(filename) ;
	}

    void MyDelaunay::get_bbox(
        real& x_min, real& y_min, real& z_min,
        real& x_max, real& y_max, real& z_max
	) {
        x_min = y_min = z_min =  1e30 ;
        x_max = y_max = z_max = -1e30 ;
        for(unsigned int i=0; i<boundary_.size(); i++) {
            for(unsigned int j=0; j<3; j++) {
                const vec3& p = boundary_[i].vertex[j] ;
                x_min = gx_min(x_min, p.x) ;
                y_min = gx_min(y_min, p.y) ;
                z_min = gx_min(z_min, p.z) ;
                x_max = gx_max(x_max, p.x) ;
                y_max = gx_max(y_max, p.y) ;
                z_max = gx_max(z_max, p.z) ;
            }
        }
    }
    
    // ------------------------------------ MyDelaunay 

    void MyDelaunay::clear() {
        baseclass::clear() ;
        all_vertices_.clear() ;
		if(kdtree_ != nil) {
		    delete kdtree_;
		    kdtree_ = nil; 
		}
    }

    MyDelaunay::Vertex_handle MyDelaunay::insert(const vec3& p) {
		std::cerr << "insert function should not be used since now!!! \n";
		return 0;
  //      if(Numeric::is_nan(p.x) || Numeric::is_nan(p.y) || Numeric::is_nan(p.z)) {
  //          std::cerr << "Nan ! delaunay insert vertex" << std::endl ;
  //          return 0 ;
  //      }
  //      Vertex_handle result = baseclass::insert(to_cgal(p)) ;
		//result->index = (int)all_vertices_.size() ; // dmyan
  //      all_vertices_.push_back(result) ;
  //      return result ;

    }

    void MyDelaunay::begin_insert() { 
	}

    void MyDelaunay::end_insert(bool redraw) {
		//moved insert_corners function to set_vertices
		//insert_corners() ; // insert constrained corner vertices
        stars_dirty_ = true ; 
		update_stars() ;
        compute_rvd() ;

        if(redraw) {
            glut_viewer_redraw() ;            
        }
    }

    inline vec3 random_v() {
        return vec3( 
            Numeric::random_float64(),
            Numeric::random_float64(),
            Numeric::random_float64()
	    ) ;
    }

	void MyDelaunay::insert_facet_vertices(int nb) {
		int fidx[4] ;
		fidx[0] = ::rand() % boundary_.size() ; 
		for(int i=0; i<3; ++i) fidx[i+1] = boundary_[fidx[0]].face_index[i] ;
		std::vector<vec3> pts ;
		for(int i=0; i<nb; ++i) {
			const Geex::Facet& F = boundary_[fidx[0]] ; //boundary_[fidx[rand()%4]] ;
			double r0 = Numeric::random_float64() ;
			double r1 = (1-r0)*Numeric::random_float64() ;
			pts.push_back(r0*F.vertex[0]+r1*F.vertex[1]+(1-r0-r1)*F.vertex[2]);
		}		

		set_vertices(pts) ;
	}

    void MyDelaunay::insert_random_vertex() {
		std::cerr << "insert_random_vertex function should not be used since now!!! \n";
		return;

        //Geex::vec3 p = random_v() ;
        //int nb_tries = 0 ;
        //while(!in_boundary(p)) {
        //    nb_tries++ ;
        //    if(!non_convex_mode_ && nb_tries > 1000) {
        //        std::cerr << "Could not insert point, probably missing +non_convex flag" << std::endl ;
        //        exit(-1) ;
        //    }
        //    p = random_v() ;
        //}
        //insert(p) ;
    }

    void MyDelaunay::insert_random_vertices(int nb) {
////        begin_insert() ;
//        
//		std::vector<vec3> pts(nb);
//		
//		for (int i = 0; i < nb; i++)
//		{
//			pts[i] = random_v() ;
//			int nb_tries = 0 ;
//			while(!in_boundary(pts[i])) {
//				nb_tries++ ;
//				if(!non_convex_mode_ && nb_tries > 1000) {
//					std::cerr << "Could not insert point, probably missing +non_convex flag" << std::endl ;
//					exit(-1) ;
//				}
//				pts[i] = random_v() ;
//			}
//		}
//
//		set_vertices(pts);
//		//for(int i=0; i<nb; i++) {
//  //          insert_random_vertex() ;
//  //          std::cerr << (i+1) << '/' << nb << std::endl ;
//  //      }
//
//
////        end_insert(false) ;
    }

	void MyDelaunay::insert_random_vertices_in_box(int nb) {
		vec3 minb, maxb ;		
		boundary_.get_bounding_box(minb, maxb) ;

		std::vector<vec3> pts(nb);
		
		for(int i=0; i<nb; ++i) {
			pts[i] = random_v() ;
			pts[i][0] = minb[0]+pts[i][0]*(maxb[0]-minb[0]) ;
			pts[i][1] = minb[1]+pts[i][1]*(maxb[1]-minb[1]) ;
			pts[i][2] = minb[2]+pts[i][2]*(maxb[2]-minb[2]) ;
		}

		set_vertices(pts);
	}

    void MyDelaunay::insert_boundary_vertices(int nb) {
		if(nb == 0) {
			std::vector<vec3> pts(boundary_.nb_vertices());

			for(int i=0; i<boundary_.nb_vertices(); ++i) {
				pts[i] = boundary_.vertex(i).pos_;
				//insert(boundary_.vertex(i).pos_) ;	
			}
			set_vertices(pts);
			for(unsigned int i=0; i<boundary_.corners_.size(); ++i) {
				all_vertices_[boundary_.corners_[i]]->type() = V_CORNER ;
			}

		} else {
			sample_boundary(nb-nb_corners()) ;
		}
    }

	void MyDelaunay::sample_boundary(int nb) {
		double totalw = 0. ; 
		double res = 0 ;
        int    check = 0  ;
        
        std::vector<double> face_mass(boundary_.size());
		for(unsigned int i=0; i<boundary_.size(); ++i) {
			const Geex::Facet& F = boundary_[i] ;
			face_mass[i] = Geex::tri_mass(F.vertex[0], F.vertex[1], F.vertex[2], F.vertex_weight[0], F.vertex_weight[1], F.vertex_weight[2]) ;
            totalw += face_mass[i];
		}

		std::vector<vec3> pts;

		for(unsigned int i=0; i<boundary_.size(); ++i) {
			const Geex::Facet& F = boundary_[i] ;
			double nf = nb*face_mass[i]/totalw ;
			int n = (int)(nf+res) ;
			if(n >= 1) {
				res = nf+res - n ;
				check += n ;
				//for(int j=0; j<n; ++j) {
				//	double r0 = Numeric::random_float64() ;
				//	double r1 = (1-r0)*Numeric::random_float64() ;
				//	pts.push_back(r0*F.vertex[0]+r1*F.vertex[1]+(1-r0-r1)*F.vertex[2]);
				//	//insert(r0*F.vertex[0]+r1*F.vertex[1]+(1-r0-r1)*F.vertex[2]) ;			
				//}
				sample_tri(F.vertex[0], F.vertex[1], F.vertex[2], n, pts) ;
			}
			else
				res += nf ;
		}

		set_vertices(pts);

//		std::cout << "required sample: " << nb << " number of samples: " << check << std::endl ;
	}

	void MyDelaunay::insert_inner_vertices(int nb) {
		//if(sizingfield_!=nil) {
		//	std::vector<vec3> pts ;
		//	sizingfield_->sample_grid(pts, nb) ;
		//	sizingfield_->sample_grid_uniform(pts, nb) ;
		//	for(int i=0; i<pts.size(); ++i)
		//		insert(pts[i]) ;
		//}
		////if(!volume_.is_valid()) {
		////	std::cerr << "volume tetmesh is invalid..." << std::endl ;
		////	return ;
		////}
		////for(int i=0; i<nb; ++i){ 
		////	int              tidx = abs(Numeric::random_int32())%volume_.ntets_ ; // choose a random tet
		////	const Geex::Tet& tet = volume_.tets_[tidx] ;
		////	double r0 = Numeric::random_float64() ;
		////	double r1 = (1-r0)*Numeric::random_float64() ;
		////	double r2 = (1-r0-r1)*Numeric::random_float64() ;
		////	insert(r0*tet.verts_[0]+r1*tet.verts_[1]+r2*tet.verts_[2]+(1-r0-r1-r2)*tet.verts_[3]) ;
		////}
	}

	void MyDelaunay::insert_corners() {
		std::cerr <<"insert_corners function should not be used since now !!!\n";
//		std::vector<int>& corners = boundary_.corners_ ; 
////		corner_vertices_.resize(cons.size()) ;
//		for(unsigned int i=0; i<corners.size(); ++i) {
//			insert(boundary_.vertex(corners[i]).point()) ;
//			all_vertices_[all_vertices_.size()-1]->type() = V_CORNER ;
//		}
	}

	void MyDelaunay::insert_vertices(int nb, int vertex_type) {
		//begin_insert() ;
		switch(vertex_type) {
			case 0:
				insert_random_vertices_in_box(nb) ;
				break ;
			case 1:
				insert_boundary_vertices(nb) ;
				break ;
			case 2:
				insert_inner_vertices(nb) ;
				break ;
			case 3:
				insert_facet_vertices(nb) ;
				break ;
		}
//		insert_random_vertices(nb) ;
		//end_insert(false) ;
	}

	void MyDelaunay::insert_vertices(std::string& ptsfilename) {
	    std::cout << "insert vertices from file: " << ptsfilename << std::endl ;
	    std::ifstream pts(ptsfilename.c_str()) ;
	    int num ;
	    //double x, y, z ;
	    //begin_insert() ;
	    pts >> num ;
	    std::cout << "number of points: " << num << std::endl ;
		std::vector<vec3> mpoints(num);
	    for(int i=0; i<num; ++i) {
			pts >> mpoints[i];
		//pts >> x >> y >> z ;
		//insert(vec3(x, y, z)) ;
	    }
		set_vertices(mpoints);
	    //end_insert(false) ;
	    pts.close() ;
	    std::cout << "end insert vertices from file " << std::endl ;
	}

	void MyDelaunay::set_vertices(std::vector<vec3>& new_points, bool tag) {
		prepare_sort_points(new_points, mypointvec);
		insert_sort_points(mypointvec, tag);
	}
	void MyDelaunay::set_vertices(int N, const double* x, bool tag) {
		prepare_sort_points(N, x, mypointvec);
		insert_sort_points(mypointvec, tag);
	}

	// new insertion code
	void MyDelaunay::prepare_sort_points(int N, const double *x, std::vector<MyPoint>& mypointvec)
	{
		int n = N/3;

		std::vector<int>& corners = boundary_.corners_ ; 
		int n_corner = (int)boundary_.corners_.size();
		bool first = false;
		if ((int)mypointvec.size() != n+n_corner)
		{
			mypointvec.clear();
			mypointvec.resize(n+n_corner);
			for (int i = 0; i < n+n_corner; i++)
			{
				mypointvec[i].id = i;
			}
			for (int i = n; i < n+n_corner; i++)
			{
				mypointvec[i].pt = to_cgal(boundary_.vertex(corners[i - n]).point());
			}
			first = true;
		}
		int count = 0;
		for(int i = 0; i < n+n_corner; i++)
		{
			if ( mypointvec[i].id < n)
			{
				count = 3*mypointvec[i].id;
				mypointvec[i].pt = Point(x[count], x[count+1], x[count+2]);
			}
		}
		K k;
		if (first)
		{
			std::random_shuffle(mypointvec.begin(), mypointvec.end());
		}
		CGAL::CVT_spatial_sort_3(mypointvec.begin(), mypointvec.end(), k);
	}
	void MyDelaunay::prepare_sort_points(const std::vector<vec3>& new_points, std::vector<MyPoint>& mypointvec)
	{
		int N = (int)new_points.size();

		std::vector<int>& corners = boundary_.corners_ ; 
		int n_corner = (int)boundary_.corners_.size();

		bool first = false;
		if ((int)mypointvec.size() != N+n_corner)
		{
			mypointvec.resize(N+n_corner);
			for (int i = 0; i < N+n_corner; i++)
			{
				mypointvec[i].id = i;
			}
			for (int i = N; i < N+n_corner; i++)
			{
				mypointvec[i].pt = to_cgal(boundary_.vertex(corners[i - N]).point());
			}
			first = true;
		}

		for(int i = 0; i < N+n_corner; i++)
		{
			if (mypointvec[i].id < N)
			{
				mypointvec[i].pt = to_cgal(new_points[mypointvec[i].id]);
			}
		}

		K k;
		if (first)
		{
			std::random_shuffle(mypointvec.begin(), mypointvec.end());
		}
		CGAL::CVT_spatial_sort_3(mypointvec.begin(), mypointvec.end(), k);

	}
	void MyDelaunay::insert_sort_points(const std::vector<MyPoint>& mypointvec, bool tag)
	{
		baseclass::clear() ;
		
		if(kdtree_ != nil) {
			delete kdtree_;
			kdtree_ = nil; 
		}

		begin_insert() ;

		int N = (int)mypointvec.size();

		if ((int)all_vertices_.size() != N)
		{
			all_vertices_.resize(N);
		}
		
		int n_corners = (int) boundary_.corners_.size();


		Cell_handle ch;
		Vertex_handle vh;
		//std::cout << "======================\n";
		for (std::vector<MyPoint>::const_iterator p = mypointvec.begin(), end = mypointvec.end();
			p != end; ++p)
		{
			//std::cout << p->pt <<"\n";
			//std::cout.flush();
			vh = baseclass::insert(p->pt, ch);
			ch = vh->cell();
			vh->idx() = p->id;
			all_vertices_[p->id] = vh;	
			if (p->id >= N - n_corners)
			{
				all_vertices_[p->id]->type() = V_CORNER;
			}
		}
		//std::cout << "=====================\n";
		//std::cout.flush();
		end_insert(tag) ;
		
	}



    // ------------------------------------ MyDelaunay combinatorics

    void MyDelaunay::get_incident_cells(Vertex_handle v, std::set<Cell_handle>& star) {
        std::stack<Cell_handle> S ;
        S.push(v->cell()) ;
        star.insert(v->cell()) ;
        while(!S.empty()) {
            Cell_handle c = S.top() ;
            S.pop() ;
            for(unsigned int i=0; i<4; i++) {
                Cell_handle d = c->neighbor(i) ;
                if(d->has_vertex(v) && (star.find(d) == star.end())) {
                    S.push(d) ;
                    star.insert(d) ;
                }
            }
        }
    }


    MyDelaunay::CellStatus MyDelaunay::get_incident_edges(
        Vertex_handle v, std::vector<Edge>& E
	) {
        CellStatus result = CS_INTERIOR ;

        std::stack<Edge> S ;
        std::set<Vertex_handle> visited ;
        Cell_handle c = v->cell() ;
        int i = c->index(v) ;
        for(int j=0; j<4; j++) {
            if(j != i) {
                S.push(Edge(c,i,j)) ;
                E.push_back(S.top()) ;
                visited.insert(c->vertex(j)) ;
            }
        }
        while(!S.empty()) {
            Edge e = S.top() ;
            S.pop() ;
            Cell_circulator ccir = incident_cells(e) ;
            Cell_circulator cdone = ccir ;
            do {
                if(ccir->infinite) { 
                    result = CS_INFINITE ; 
                } else if(result == CS_INTERIOR) {
                    if(ccir->dual_outside) {
                        result = CS_CROSS ;
                    }
                }
                int i = ccir->index(v) ;
                for(int j=0; j<4; j++) {
                    if((j != i) && (visited.find(ccir->vertex(j)) == visited.end())) {
                        S.push(Edge(ccir,i,j)) ;
                        E.push_back(S.top()) ;
                        visited.insert(ccir->vertex(j)) ;
                    }
                }
                ccir++ ;
            } while(ccir != cdone) ;
        }
        return result ;
    }
    
    int MyDelaunay::dual_facet_degree(Edge e) {
        int result = 0 ;
        Cell_handle prev = nil ;
        Cell_circulator ccir = incident_cells(e) ;
        Cell_circulator cdone = ccir ;
        do {
            Cell_handle p = ccir ;
            prev = ccir ;
            ccir++ ;
            result++ ;
        } while(ccir != cdone) ;
        return result ;
    }

    MyDelaunay::Cell_handle MyDelaunay::edge_key(Edge e) {
        Cell_handle result = nil ;
        Cell_circulator ccir = incident_cells(e) ;
        Cell_circulator cdone = ccir ;
        do {
            Cell_handle cur = ccir ;
            if(&*cur > &*result) { result = cur ; }
            ccir++ ;
        } while(ccir != cdone) ;
        return result ;
    }

    void MyDelaunay::normalize_oriented_edge(Edge& e) {
        Vertex_handle v1 = e.first->vertex(e.second) ;
        Vertex_handle v2 = e.first->vertex(e.third) ;
        e.first = edge_key(e) ;
        e.second = e.first->index(v1) ;
        e.third = e.first->index(v2) ;
    }

    TriMesh* MyDelaunay::dual_convex_clip(std::vector<Edge>& edges, bool close) {

        TriMesh* from = &boundary_ ;
        TriMesh* to = &pong_ ;

        // No need to have special cases for infinite vertices: they
        // do not yield any clipping plane !
        for(unsigned int i=0; i<edges.size(); i++) {
            Edge e = edges[i] ;
//          normalize_oriented_edge(e) ;  // Not needed ?
            if(is_infinite(e.first->vertex(e.third))) { continue ; }
            Geex::Plane<real> P = get_dual_plane(e) ;
            from->convex_clip(*to, P, close ? TriMesh::CLIP_CONVEX : TriMesh::CLIP_BNDRY ) ;
            if(from == &boundary_) {
                from = &pong_ ;
                to = &ping_ ;
            } else {
                gx_swap(from, to) ;
            }
        }
        return from ;
    }

	// -----------------------------------------------------------------------------------------
	// ------------------------- Voronoi Diagram Computation and Remeshing ---------------------
	// -----------------------------------------------------------------------------------------

	void MyDelaunay::compute_lfs() {
		std::vector<vec3> pts ;
		for(unsigned i=0; i<boundary_.size(); ++i) {
			pts.push_back(boundary_[i].center()) ;
		}
		Geex::MedialAxis lfs(pts) ;
		for(unsigned i=0; i<boundary_.size(); ++i) {
			for(int j=0; j<3; ++j) {
				boundary_[i].vertex_weight[j] = 1.0/lfs.squared_lfs(boundary_[i].vertex[j]) ;
			}
		}

		//// output lfs to text file
		//std::ofstream out("out.lfs") ;
		//out << boundary_.nb_vertices() << std::endl ;
		//for(unsigned i=0; i<boundary_.nb_vertices(); ++i) {
		//	out << 1.0/lfs.squared_lfs(boundary_.vertex(i).point()) << std::endl ;
		//}
		//out.close() ; 
	}

	void MyDelaunay::compute_rvd() {
		CGAL::Timer timer ;
		double ts;
		timer.start() ;
		ts = timer.time() ;
		rvd_->compute(&boundary_, this) ;
		timer.stop() ;
//		std::cout << "RVD computation time:  " << timer.time()-ts << std::endl ;
	}

	void MyDelaunay::refine_feature() {
		rvd_->refine_feature() ;
	}

	void MyDelaunay::extract_mesh() {
		rvd_->extract_mesh() ;
	}
	void MyDelaunay::extract_dual() {
		rvd_->extract_dual() ;
	}

    void MyDelaunay::voronoi_filtering()
    {
        rvd_->voronoi_filtering((double)voronoi_fitering_angle_) ;
    }

    void MyDelaunay::alpha_shape()
    {
        rvd_->alpha_shape();
    }

    void MyDelaunay::save_remesh(const std::string& filename) {
		rvd_->dual_mesh().save_obj(filename) ;
    }
    
	void MyDelaunay::load_remesh(const std::string& filename) {
		std::cout << "load remesh: " << filename << std::endl ;
		rvd_->dual_mesh().load(filename) ;
		//vec3 g ;
		//double R ;
		//rvd_->dual_mesh().normalize(g, R) ;
    }

	void MyDelaunay::export_eps(const std::string& file_name) {
//		    std::string ext = FileSystem::extension(file_name) ;
//            GLint format = 0 ;
//            if(ext == "ps") {
//                format = GL2PS_PS ;
//            } else if(ext == "eps") {
//                format = GL2PS_EPS ;
//            } else if(ext == "pdf") {
//                format = GL2PS_PDF ;
//            } else {
////                ogf_assert_not_reached ;
//            }
//
//            GLint options =  
//                GL2PS_DRAW_BACKGROUND      |  
//      //          GL2PS_SIMPLE_LINE_OFFSET   | 
//				GL2PS_SILENT               | 
//				GL2PS_NO_BLENDING          |
//                GL2PS_OCCLUSION_CULL       |
//                GL2PS_USE_CURRENT_VIEWPORT | 
//                GL2PS_COMPRESS ;
//
//       //     make_current() ;
//            FILE* snapshot_output_ = fopen(file_name.c_str(),"wb") ;
//            gl2psBeginPage(
//                "snapshot", "GEEX",
//                nil,
//                format, 
//                GL2PS_SIMPLE_SORT,
//                options,
//                GL_RGBA,
//                0,
//                0,
//                16,16,4,
//                16*1024*1024,
//                snapshot_output_, file_name.c_str()
//            ) ;
//
//            //gl2psPointSize(4.0) ;
//            //gl2psLineWidth(1.0) ;
//            gl2psEnable(GL2PS_POLYGON_OFFSET_FILL) ; 
//
//			glut_viewer_redraw() ; // draw scene
//			gl2psEndPage() ;
//			fclose(snapshot_output_) ;
	}

    void MyDelaunay::save_delaunay(std::string& filename) {
	std::ofstream ofile(filename.c_str()) ;
	ofile << number_of_vertices() << std::endl ;
	for(Finite_vertices_iterator vit=finite_vertices_begin(); vit!=finite_vertices_end(); ++vit) {
	    ofile << vit->point() <<std::endl;			
	}
	ofile.close() ;
    }

    void MyDelaunay::load_delaunay(std::string& filename) {
	std::ifstream ifile(filename.c_str()) ;
	int size ;
	Point p ;
	CGAL::Timer timer ;

	clear() ;
	//begin_insert() ;
	timer.start() ; 
	ifile >> size ;
	std::vector<vec3> pts(size);
	for(int i=0; i<size; ++i) {
	    ifile >> pts[i] ;
	    //insert(to_geex(p)) ;
	}
	ifile.close() ;
	set_vertices(pts, true);
	end_insert(true) ;

	timer.stop() ;
	std::cout << "triangulation and clip time: " << timer.time() << std::endl ;
	//timer.start() ;
	//voronoi_clip() ;
	//clip_corner() ; // identify fixed corner vertices
	//timer.stop() ;
	//std::cout << "voronoi clip  time: " << timer.time() << std::endl ;
    }

    bool MyDelaunay::load_density(std::string& filename, real gamma) {
	if(boundary_.load_density(filename, gamma)) {
//	    set_vertex_tex() ;
	    return true ;
	}
	return false; 
    }
	
	//void MyDelaunay::set_vertex_tex() {
	//	ramp_.begin_range(); 
	//	for(int i=0; i<boundary_.nb_vertices(); ++i) {
	//		ramp_.add_value_to_range(boundary_.vertex(i).weight_) ;
	//		boundary_.vertex(i).tex_coord_ = boundary_.vertex(i).weight_ ;
	//	}
	//	ramp_.end_range() ;
	//}

	bool MyDelaunay::load_feature(std::string& filename) {
		return boundary_.load_feature(filename) ;
	}

	void MyDelaunay::voronoi_clip_timing(int max_seed_num) {
	    CGAL::Timer timer ;
	    double ts, te, dttime, vdtime ;
	    int n = 10 ;
	    double p = 1. ;
	    
	    while(n < max_seed_num) {
		n = (int)pow(10., p) ;
		for(int i=1; i<10; ++i) {
		    clear() ;
		    timer.start() ;
		    ts = timer.time() ;
		    insert_boundary_vertices(n) ;
		    timer.stop() ;
		    te = timer.time() ;
		    dttime = te - ts ;
		    timer.start() ;
		    ts = timer.time() ;
//		    voronoi_clip() ;
			compute_rvd() ;
		    timer.stop() ;
		    glut_viewer_redraw() ;            
		    te = timer.time() ;
		    vdtime = te - ts ;
		    std::cout << n << "\t" << dttime <<"\t" << vdtime << std::endl ;
		    n += (int)pow(10., p) ;
		}
		p += 1. ;
	    }
	} 

    void MyDelaunay::update_stars() {
        if(stars_.size() != all_vertices_.size()) {
            stars_.resize(all_vertices_.size()) ;
			cells_.resize(all_vertices_.size()) ;
        }
        for(unsigned int i=0; i<all_vertices_.size(); i++) {
			std::vector<Edge> edges ;
            get_incident_edges(all_vertices_[i], edges) ;
            stars_[i].clear() ;
			cells_[i].clear() ;

			int nb_edges =(int)edges.size() ;
			for(int j=0; j<nb_edges; ++j) {
				if(!is_infinite(edges[j].first->vertex(edges[j].third))){
					stars_[i].push_back(edges[j]) ;
					cells_[i].push_back(get_dual_plane(edges[j])) ;
				}
			}
        }
        stars_dirty_ = false ;
    }

	void MyDelaunay::partition_mesh() {
		int                   nb_cores = Geex::Processor::number_of_cores() ;
		int                   nf = (int)boundary_.size() ;
		SparseMatrix          adjmat(nf, nf) ; 
		MatrixPartition_METIS partition ; 
		Array1d<int>          result ;
		result.allocate(nf);

		std::cout << "number of processors: " << nb_cores << std::endl ;
		for(int i=0; i<nf; ++i) {
			for(int j=0; j<3; ++j) 
				if(boundary_[i].face_index[j]>=0)
					adjmat.add(i, boundary_[i].face_index[j], 1) ;
		}

		partition.partition_matrix(&adjmat, nb_cores, result) ;

		parts_.resize(nb_cores) ;
		std::vector<int> new_indices(nf) ;

		for(int i=0; i<nf; ++i) {
//			std::cout << result[i] << std::endl ;
			new_indices[i] = (int)parts_[result[i]].size() ;
			parts_[result[i]].push_back(boundary_[i]) ;
		}

		// reset neighbor faces of each part
		for(int i=0; i<nb_cores; ++i) {
			for(int j=0; j<(int)parts_[i].size(); ++j) {
				for(int k=0; k<3; ++k) {
					int nid = parts_[i][j].face_index[k] ;
                    if (nid >= 0)
                    {
					if(result[nid] != i) // in different parts
						parts_[i][j].face_index[k] = -1 ; 
					else
						parts_[i][j].face_index[k] = new_indices[nid] ;
                    }
				}
			}
		}
	}
}
