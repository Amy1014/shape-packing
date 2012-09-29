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

#include "rvd_ccvt.h"
#include "delaunay.h"

namespace Geex {
//---------------------- Edge clipping algorithm -----------------------
	void RestrictedVoronoiDiagram::clip_feature() {
//		for(unsigned int i=0; i<rvcells_.size(); ++i) 
//			rvcells_[i].clip_facet_list_.clear() ; // used for tets

		for(int i=0; i<(int)mesh_->features_.size(); ++i) {
			Geex::LineSeg& seg = mesh_->features_[i] ;
			std::set<int> edge_cells ;
			vec3 g = 0.5*(seg.v0_+seg.v1_) ;
			
			if(Numeric::is_nan(g.x) || Numeric::is_nan(g.y) || Numeric::is_nan(g.z)) {
				std::cerr << "Nan ! edge clip" << std::endl ;
				return ;
			}
			
			int vidx = kdtree_->queryNearestNeighbor(g[0], g[1], g[2]) ;

			if(rvcells_[vidx].insert_clipped_seg(i) )	
			    //if(edge_cells.insert(vidx).second) 
				clip_edge_by_cell(i, vidx, edge_cells) ;
		}
	}

	void RestrictedVoronoiDiagram::clip_edge_by_cell(int eidx, int vidx, std::set<int> edge_cells) {
		std::vector<Geex::LineSeg> ping, pong, outlist;
		std::vector<Geex::LineSeg> *from, *to ;
		Geex::LineSeg& S = mesh_->features_[eidx] ;
		std::vector<TDS::Edge>& edges = dt_->stars_[vidx] ;
		std::vector<Geex::Plane<real> >& planes = dt_->cells_[vidx] ;
		
		ping.push_back(S) ;
		from = &ping ;
		to = &pong ;

		for(unsigned int j=0; j<planes.size(); ++j) {
			const Geex::Plane<real>& P = planes[j] ;
			const TDS::Edge& e = edges[j] ;
			int oppvidx = e.first->vertex(e.second)==rvcells_[vidx].vh()?e.first->vertex(e.third)->idx():e.first->vertex(e.second)->idx() ;

			to->clear() ;
			for(unsigned int k=0; k<from->size(); ++k) {
				clip_edge_by_plane((*from)[k], P, oppvidx, *to) ;
			}			
			gx_swap(from, to) ;

			if(from->size() == 0) 
				break ; 
		}

		//
		if(from->size() > 0) {
			for(unsigned int i=0; i<from->size(); ++i)
				rvcells_[vidx].insert_lineseg((*from)[i]) ;

			// recursive clip this tet by neighbor tet
			for(unsigned int i=0; i<from->size(); ++i) {
				for(int j=0; j<2; ++j) {
					int nidx = (*from)[i].nclip_[j] ;
					if( nidx >= 0) {
						if(rvcells_[nidx].insert_clipped_seg(eidx)) {
						    //if(edge_cells.insert(nidx).second) {
							clip_edge_by_cell(eidx, nidx, edge_cells) ;
						}
					}
				}
			} // end clip by neighbor
		}
	}

	void RestrictedVoronoiDiagram::clip_edge_by_plane(Geex::LineSeg& E, const Geex::Plane<real>& P, int oppvidx, std::vector<Geex::LineSeg>& inlist) {
		int config = 
		    ((P.side(E.v0_) > 1e-10) ? 1 : 0) |
		    ((P.side(E.v1_) > 1e-10) ? 2 : 0) ;
	            
		switch(config) {
		case 0: 
		    break ; // Completely outside
		case 3: 
			inlist.push_back(E) ;
			break ;
		case 1: {
		    vec3 P01 = seg_plane_intersection(E.v0_, E.v1_, P) ;
			inlist.push_back(Geex::LineSeg(E.v0_, P01, E.nclip_[0], oppvidx, E.parentidx())) ;
			} break ; 
		case 2: {
		    vec3 P01 = seg_plane_intersection(E.v0_, E.v1_, P) ;
			inlist.push_back(Geex::LineSeg(P01, E.v1_, oppvidx, E.nclip_[1], E.parentidx())) ;
			} break ;
		}
	}

		
	//------------------------ clip corners --------------------------------
	void RestrictedVoronoiDiagram::clip_corner() {
		for(unsigned int i=0; i<mesh_->corners_.size(); ++i) {
			vec3& g = mesh_->vertex(mesh_->corners_[i]).point() ;
			int vidx = kdtree_->queryNearestNeighbor(g[0], g[1], g[2]) ;
			dt_->all_vertices_[vidx]->type() = V_CORNER ;
			dt_->all_vertices_[vidx]->point() = to_cgal(g) ;
		}
	}

	////----- insert new seeds if one seed occupied more than 1 feature lineseg ------
	void RestrictedVoronoiDiagram::refine_feature() {
	//        std::vector<vec3> new_points ;
	//    
	//		for(unsigned int i=0; i<rvcells_.size(); ++i) {
	//			if(!rvcells_[i].is_crease()) {
	//				double V ;
	//				vec3 g ;
	//				rvcells_[i].get_centroid(g, V) ;
	//				new_points.push_back(g) ; //to_geex(rvcells_[i].vh_->point())) ;
	//			}
	//			else {
	//				std::vector<vec3> gs ;
	//				rvcells_[i].get_crease_centroids(gs) ;
	//				new_points.insert(new_points.end(), gs.begin(), gs.end()) ;
	//				if(gs.size() > 0 ) std::cout << "split feature..." << std::endl ;
	//			}				
	//		}
	//    
	//	    clear() ;
 //           begin_insert() ;
	//		//insert_constrains() ;
 //           for(unsigned int j=0; j<new_points.size(); j++) {
 //               insert(new_points[j]) ;
 //           }
 //           end_insert(false) ;
	//    
	//		std::cout <<"finish refining feature. " << std::endl ;
	//		// DT should be updated before voronoi clip
	//		voronoi_clip() ; // compute Voronoi diagram
	}


}
