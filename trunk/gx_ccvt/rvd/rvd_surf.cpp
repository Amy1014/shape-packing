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

#include "delaunay.h" 
#include "rvc_ccvt.h"
#include "rvd_ccvt.h"

#include <Geex/CVT/RVD.h>
#include <Geex/CVT/delaunay_skel_CGAL.h>
#include <Geex/CVT/RVD_poly_kdtree.h>
#include <Geex/CVT/VoronoiFiltering.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_hierarchy_3.h>
#include <CGAL/Alpha_shape_3.h>



#include <omp.h>

namespace Geex {

	int RestrictedVoronoiDiagram::configs[8][3] = {
		{0, 1, 2},
		{0, 1, 2},
		{1, 2, 0},
		{0, 1, 2},
		{2, 0, 1},
		{2, 0, 1},
		{1, 2, 0},
		{0, 1, 2},
	} ;

	void RestrictedVoronoiDiagram::compute(Geex::TriMesh* mesh, MyDelaunay* dt) {
		mesh_ = mesh ;
		dt_   = dt ;
		init_vcells() ;

		if(rvd_mode_==TRI_KD_INCIDENT || rvd_mode_==TRI_KD_SPLIT || rvd_mode_==POLY_KD_INCIDENT) { 
			build_kdtree() ;
		} 

		switch(rvd_mode_) {
			case TRI_LINEAR:
				clip_surface_linear() ;
				//std::cout << "compute RVD by triangle clipping propagation." << std::endl ;
				break ;
			case TRI_KD_INCIDENT:
				clip_surface_mlogn() ;
				//std::cout << "compute RVD by triangle clipping kd-query" << std::endl ;
				break ;
			case TRI_KD_SPLIT:
				clip_surface_splitting() ;
				//std::cout << "compute RVD by triangle splitting" << std::endl ;
				break ;
			case POLY_LINEAR: {
				DelaunaySkeleton skel ;
				copy_CGAL_Delaunay_to_skel(dt_, &skel, dt->vertices()) ;
				clip_poly_linear(&dt->boundary_topo_, dt_, &skel) ;
				//std::cout << "compute RVD by polygon clipping propagation." << std::endl ;
			}
				break ;
			case POLY_KD_INCIDENT: {
				DelaunaySkeleton skel ;
				copy_CGAL_Delaunay_to_skel(dt_, &skel, dt->vertices()) ;
				clip_poly_kdtree(&dt->boundary_topo_, dt_, &skel, kdtree_) ;
								   }
				break ;
			case PARA_TRI:
				compute(dt_->parts_, dt_) ;
				break;
			case PARA_POLY:
				compute(dt_->boundary_topo_parts_, dt_) ;
				break;
			default:
				break ;
		}

		if(dt_->is_clip_feature())
			collect_crease_segments() ;
	}

	void RestrictedVoronoiDiagram::compute(Array1d<Geex::TopoPolyMesh>& parts, MyDelaunay* dt) {
        std::vector<RestrictedVoronoiDiagram> threads ;
		int nb_threads = parts.size() ;
		DelaunaySkeleton skel ;
		copy_CGAL_Delaunay_to_skel(dt, &skel) ;

		for(int i=0; i<nb_threads; ++i)
			threads.push_back(RestrictedVoronoiDiagram(dt, nil, POLY_LINEAR)) ;

		int i ;
#pragma omp parallel for shared(nb_threads) private(i)
		for(i=0; i<nb_threads; ++i) 
		//	threads[i].compute(&parts[i], dt) ;
			threads[i].clip_poly_linear(&parts[i], dt, &skel) ;

		// merge results of threads
		for(int i=0; i<nb_threads; ++i) {
			int nv = (int)rvcells_.size() ;
			std::vector<RestrictedVoronoiCell>& subcells = threads[i].rvcells() ;
			for(int j=0; j<nv; ++j) {
				rvcells_[j].append(subcells[j]) ;
			}
		}

		std::cout << "finish parallel RVD computation by linear polygon clipping." << std::endl ;
	}

	void RestrictedVoronoiDiagram::compute(std::vector<Geex::TriMesh>& parts, MyDelaunay* dt) {
                std::vector<RestrictedVoronoiDiagram> threads ;
		int nb_threads = (int)parts.size() ;

		for(int i=0; i<nb_threads; ++i)
			threads.push_back(RestrictedVoronoiDiagram(dt, &parts[i], TRI_LINEAR)) ;

		int i ;
//		omp_set_num_threads(nb_cores) ;
#pragma omp parallel for shared(nb_threads) private(i)
		for(i=0; i<nb_threads; ++i) 
			threads[i].compute(&parts[i], dt) ;

		// merge results of threads
		for(int i=0; i<nb_threads; ++i) {
			int nv = (int)rvcells_.size() ;
			std::vector<RestrictedVoronoiCell>& subcells = threads[i].rvcells() ;
			for(int j=0; j<nv; ++j) {
				rvcells_[j].append(subcells[j]) ;
			}
		}

		std::cout << "finish parallel RVD computation by linear triangle clipping." << std::endl ;
	}

	//------------------------------------------------------------------------------------
	//----------------------- collect feature edges from RVD -----------------------------
	//------------------------------------------------------------------------------------

	void RestrictedVoronoiDiagram::collect_crease_segments() {
		int i, nv = dt_->nb_vertices() ;
#pragma omp parallel for shared(nv) private(i)
		for(i=0; i<nv; ++i)
			rvcells_[i].collect_crease_segments() ;
	}

	void RestrictedVoronoiDiagram::voronoi_clip_linear() {
		init_vcells() ;   // init voronoi cells 
		clip_surface_linear() ;
	}

	//------------------------------------------------------------------------------------
	//----------------------- kdtree based clipping algorithm ----------------------------
	//------------------------------------------------------------------------------------

	void RestrictedVoronoiDiagram::clip_surface_mlogn() {
		typedef std::pair<int, int> CellFacePair ;
//		std::queue<CellFacePair>    clip_queue ;
		std::vector<int>            nearest_seeds ;
		int emptyinter =  0 ;
		int nf = (int)(*mesh_).size(), idx ;
		int i ;

		nearest_seeds.resize(nf) ;
		for(i=0; i<nf; ++i) {
			//		int idx = 0 ;
			//		int cidx = nearest_vertex(to_cgal(boundary_[idx].center()))->idx() ;
			const vec3 g = (*mesh_)[i].center() ;
			nearest_seeds[i] = kdtree_->queryNearestNeighbor(g[0], g[1], g[2]) ;
		}

//#pragma omp parallel for shared(nf) private(idx)
		for(idx=0; idx<nf; ++idx) {
			// init queue
			std::queue<CellFacePair>    clip_queue ;
			clip_queue.push(CellFacePair(nearest_seeds[idx], idx)) ;
			rvcells_[nearest_seeds[idx]].insert_clipped_facet(idx) ;
			
			//const vec3 g = (*mesh_)[idx].center() ;
			//int nearest = kdtree_->queryNearestNeighbor(g[0], g[1], g[2]) ;
			//clip_queue.push(CellFacePair(nearest, idx)) ;
			//rvcells_[nearest].insert_clipped_facet(idx) ;
			
			while(!clip_queue.empty()) {
				CellFacePair cfp = clip_queue.front() ;
				clip_queue.pop() ;
				//			if(vcells_[cfp.first].insert_clipped_facet(cfp.second)) {
				std::vector<Geex::Facet> L ; 
				clip_tri_by_cell(cfp.second, cfp.first, L) ;
				//if(L.size() <= 0) emptyinter++ ;
				for(unsigned int i=0; i<L.size(); ++i) {
					rvcells_[cfp.first].insert_facet(L[i]) ;
				}

				for(int i=0; i<(int)L.size(); ++i) {
					for(int j=0; j<3; ++j) {
						int nidx = L[i].nclip_[j] ;
						if( nidx >= 0) { // clip by incident cell
							if(rvcells_[nidx].insert_clipped_facet(cfp.second)) 
								clip_queue.push(CellFacePair(nidx, cfp.second));
						} 
						//else if(nidx==-1) { // clip cell by incident face
						//	for(int j=0; j<3; ++j) {
						//		int nei = L[i].face_index[j] ;
						//		if(nei<0) continue ;
						//		if(vcells_[cfp.first].insert_clipped_facet(nei)) 
						//			clip_queue.push(CellFacePair(cfp.first, nei));
						//	}
						//}
					}					
				}
			}
		}
	}


	//------------------------------------------------------------------------------------
    //----------------------- linear clipping algorithm ----------------------------------
	//------------------------------------------------------------------------------------

	void RestrictedVoronoiDiagram::clip_surface_linear() {
		typedef std::pair<int, int> CellFacePair ;
		std::queue<CellFacePair>    clip_queue ;
		int idx = 0 ;
		vec3 g = (*mesh_)[idx].center() ;
		Vertex_handle vh = dt_->nearest_vertex(to_cgal(g)) ;
		int cidx = vh->idx() ;
		clip_queue.push(CellFacePair(cidx, idx)) ;
		rvcells_[cidx].insert_clipped_facet(idx) ;
		
		while(!clip_queue.empty()) {
			CellFacePair cfp = clip_queue.front() ;
			clip_queue.pop() ;
			std::vector<Geex::Facet> L ; 
			clip_tri_by_cell(cfp.second, cfp.first, L) ;
			for(unsigned int i=0; i<L.size(); ++i) {
			    rvcells_[cfp.first].insert_facet(L[i]) ;
			}

			for(int i=0; i<(int)L.size(); ++i) {
				for(int j=0; j<3; ++j) {
					int nidx = L[i].nclip_[j] ;
					if( nidx >= 0) { // clip by incident cell
					 	if(rvcells_[nidx].insert_clipped_facet(cfp.second)) 
					 		clip_queue.push(CellFacePair(nidx, cfp.second));
					} 
					else if(nidx==-1) { // clip cell by incident face
						for(int j=0; j<3; ++j) {
							int nei = L[i].face_index[j] ;
							if(nei<0) continue ;
							if(rvcells_[cfp.first].insert_clipped_facet(nei)) 
								clip_queue.push(CellFacePair(cfp.first, nei));
						}
					}
				}					
			}
		}
	}

	void RestrictedVoronoiDiagram::clip_tri_by_cell(int fidx, int vidx, std::vector<Geex::Facet>& L) {
		std::vector<Geex::Facet> ping, pong, outlist;
		std::vector<Geex::Facet> *from, *to ;
		Geex::Facet& F = (*mesh_)[fidx] ;
		std::vector<TDS::Edge>& edges = dt_->stars()[vidx] ;
		std::vector<Geex::Plane<real> >& planes = dt_->cells()[vidx] ;

		ping.push_back(F) ;
		from = &ping ;
		to = &pong ;

		for(unsigned int j=0; j<planes.size(); ++j) {
			const Geex::Plane<real>& P = planes[j] ;
			const TDS::Edge& e = edges[j] ;
			int oppvidx = e.first->vertex(e.second)==rvcells_[vidx].vh()?e.first->vertex(e.third)->idx():e.first->vertex(e.second)->idx() ;

			to->clear() ;
			for(unsigned int k=0; k<from->size(); ++k) {
				clip_tri_by_plane((*from)[k], P, oppvidx, *to, outlist) ;
			}			
			gx_swap(from, to) ;

			if(from->size() == 0) 
				break ; 
		}

		L = *from ;
	}

	void RestrictedVoronoiDiagram::clip_tri_by_plane(Geex::Facet& F, const Geex::Plane<real>& P, int oppvidx, std::vector<Geex::Facet>& inlist, std::vector<Geex::Facet>& outlist) {
		int config = 
		    ((P.side(F.vertex[0]) > 1e-10) ? 1 : 0) |
		    ((P.side(F.vertex[1]) > 1e-10) ? 2 : 0) |
		    ((P.side(F.vertex[2]) > 1e-10) ? 4 : 0) ;
	            
		switch(config) {
		case 0: 
		    break ; // Completely outside
		case 7: 
			inlist.push_back(F) ;
			break ;
		case 1: 
		case 2:
		case 4: {
		    vec3 P01 = seg_plane_intersection(F.vertex[configs[config][0]], F.vertex[configs[config][1]], P) ;
		    vec3 P20 = seg_plane_intersection(F.vertex[configs[config][2]], F.vertex[configs[config][0]], P) ;
		    inlist.push_back(Geex::Facet(F.vertex[configs[config][0]], P01, P20, 2, F.edge_flag[configs[config][1]], F.edge_flag[configs[config][2]])) ;
		    inlist[inlist.size()-1].set_clip_neighbor(oppvidx, F.nclip_[configs[config][1]], F.nclip_[configs[config][2]]) ;
		    inlist[inlist.size()-1].set_face_neighbor(-2, F.face_index[configs[config][1]], F.face_index[configs[config][2]]) ;
			real w01 = F.weight(configs[config][0], configs[config][1], P01) ;
			real w20 = F.weight(configs[config][2], configs[config][0], P20) ;
			inlist[inlist.size()-1].set_vertex_weight(F.vertex_weight[configs[config][0]], w01, w20) ;
		} break ;
		case 3:
		case 5:
		case 6: {
		    vec3 P12 = seg_plane_intersection(F.vertex[configs[config][1]], F.vertex[configs[config][2]], P) ;
		    vec3 P20 = seg_plane_intersection(F.vertex[configs[config][2]], F.vertex[configs[config][0]], P) ;
			real w12 = F.weight(configs[config][1], configs[config][2], P12) ;
			real w20 = F.weight(configs[config][2], configs[config][0], P20) ;
		    inlist.push_back(Geex::Facet(F.vertex[configs[config][0]], P12, P20, 2, F.edge_flag[configs[config][1]], 0)) ;
			inlist[inlist.size()-1].set_vertex_weight(F.vertex_weight[configs[config][0]], w12, w20) ;
		    inlist[inlist.size()-1].set_clip_neighbor(oppvidx, F.nclip_[configs[config][1]], -2) ;
		    inlist[inlist.size()-1].set_face_neighbor(-1, F.face_index[configs[config][1]], -2) ;
		    inlist.push_back(Geex::Facet(F.vertex[configs[config][0]], F.vertex[configs[config][1]], P12,  F.edge_flag[configs[config][0]], 0, F.edge_flag[configs[config][2]])) ;
			inlist[inlist.size()-1].set_vertex_weight(F.vertex_weight[configs[config][0]], F.vertex_weight[configs[config][1]], w12) ;
		    inlist[inlist.size()-1].set_clip_neighbor(F.nclip_[configs[config][0]], -2, F.nclip_[configs[config][2]]) ;
		    inlist[inlist.size()-1].set_face_neighbor(F.face_index[configs[config][0]], -2, F.face_index[configs[config][2]]) ;
		} break ;
		}
	}

	// -------------------------------------------------------------------------------------
	// -------------- voronoi clip by kd-tree query and triangle splitting -----------------
	// -------------------------------------------------------------------------------------

	void RestrictedVoronoiDiagram::clip_surface_splitting() {

		for(unsigned int i=0; i<mesh_->size(); ++i) {
			subdivide_tri((*mesh_)[i], i) ;
		}		

		//for(unsigned int i=0; i<vcells_.size(); ++i) {
		//	if(vcells_[i].mesh_.size() > 0) 
		//		vcells_[i].set_edge_flag() ;
		//}
	}

	void RestrictedVoronoiDiagram::subdivide_tri(const Geex::Facet& F, int fidx) {
			vec3 c = 1.0 / 3.0 * (F.vertex[0]+F.vertex[1]+F.vertex[2]) ;
			int vid = kdtree_->queryNearestNeighbor(c[0], c[1], c[2]) ;
			TriMesh submesh ;
			submesh.push_back(F) ;
			subdivide_mesh(vid, submesh, fidx) ;				
	}

	void RestrictedVoronoiDiagram::subdivide_mesh(int vid, Geex::TriMesh& mesh, int fidx) {
		RestrictedVoronoiCell& rvc = rvcells_[vid] ;
		std::vector<Geex::Plane<real> >& planes = dt_->cells()[vid] ;
		std::vector<TDS::Edge>&          edges = dt_->stars()[vid] ;
		Geex::TriMesh pong;
		Geex::TriMesh* from = &mesh ;
        Geex::TriMesh* to = &pong ;
		double toler = 1e-8 ;
		bool clipped = false ;

		for(unsigned int i=0; i<planes.size(); ++i) {
			const Geex::Plane<real>& P = planes[i] ;
			int   oppidx = edges[i].first->vertex(edges[i].third)->idx() ;
			Geex::TriMesh m1 ; // the remain trimmed mesh;
			std::vector<Geex::Facet> outlist ;
			to->clear() ; 

			for(unsigned int j=0; j<from->size(); ++j) {
				Geex::Facet& F = (*from)[j] ;
				int config = 
					((P.side(F.vertex[0]) > -1e-8) ? 1 : 0) |
					((P.side(F.vertex[1]) > -1e-8) ? 2 : 0) |
					((P.side(F.vertex[2]) > -1e-8) ? 4 : 0) ;
	            
				switch(config) {

				case 0: {
					// Completely outside
					outlist.push_back(F) ;
				} break ;

				case 7: {
					// Completely inside
					to->push_back(F) ;
				} break ;

				case 1: 
				case 2:
				case 4: {
					vec3 P01 = seg_plane_intersection(F.vertex[configs[config][0]], F.vertex[configs[config][1]], P) ;
					vec3 P20 = seg_plane_intersection(F.vertex[configs[config][2]], F.vertex[configs[config][0]], P) ;
					to->push_back(Geex::Facet(F.vertex[configs[config][0]], P01, P20, 2, F.edge_flag[configs[config][1]], F.edge_flag[configs[config][2]])) ;
					(*to)[to->size()-1].set_clip_neighbor(oppidx, F.nclip_[configs[config][1]], F.nclip_[configs[config][2]]) ;
					m1.push_back(Geex::Facet(P01, F.vertex[configs[config][1]], F.vertex[configs[config][2]], F.edge_flag[configs[config][0]], 0, F.edge_flag[configs[config][2]])) ; 
					m1[m1.size()-1].set_clip_neighbor(F.nclip_[configs[config][0]], -1, F.nclip_[configs[config][2]]) ;
					m1.push_back(Geex::Facet(P01, F.vertex[configs[config][2]], P20, F.edge_flag[configs[config][1]], 1, false)) ;   
					m1[m1.size()-1].set_clip_neighbor(F.nclip_[configs[config][1]], oppidx, -1) ;
				} break ;

				case 3:
				case 5:
				case 6:{
					vec3 P12 = seg_plane_intersection(F.vertex[configs[config][1]], F.vertex[configs[config][2]], P) ;
					vec3 P20 = seg_plane_intersection(F.vertex[configs[config][2]], F.vertex[configs[config][0]], P) ;
					to->push_back(Geex::Facet(F.vertex[configs[config][0]], P12, P20, 2, F.edge_flag[configs[config][1]], 0)) ;
					(*to)[to->size()-1].set_clip_neighbor(oppidx, F.nclip_[configs[config][1]], -1) ;
					to->push_back(Geex::Facet(F.vertex[configs[config][0]], F.vertex[configs[config][1]], P12,  F.edge_flag[configs[config][0]], 0, F.edge_flag[configs[config][2]])) ;
					(*to)[to->size()-1].set_clip_neighbor(F.nclip_[configs[config][0]], -1, F.nclip_[configs[config][2]]) ;
					m1.push_back(Geex::Facet(P20, P12, F.vertex[configs[config][2]], F.edge_flag[configs[config][0]], F.edge_flag[configs[config][1]], 1)) ;   
					m1[m1.size()-1].set_clip_neighbor(F.nclip_[configs[config][0]], F.nclip_[configs[config][1]], oppidx) ;
				} break ;
				}
			}

			gx_swap(from, to) ;

			for(unsigned int j=0; j<outlist.size(); ++j) {
				subdivide_tri(outlist[j], fidx) ;
			}

			if(m1.size() > 0) {
				const MyDelaunay::Edge& e = edges[i] ;
				assert(e.first->vertex(e.second) == rvc.vh()) ;
				int nei = (e.first->vertex(e.third))->idx() ;
				subdivide_mesh(nei, m1, fidx) ;
			}
		}

//		assert(from->size()!=0) ;
		for(unsigned int i=0; i<from->size(); ++i) {
			rvc.insert_facet((*from)[i]) ;
		}
	}

	// -------------------------------------------------------------------------------------
	// ------------------------------ rvd poly computation ---------------------------------
	// -------------------------------------------------------------------------------------
	class Accumulation_poly {
	protected:
		RestrictedVoronoiDiagram *rvd_ ;
	public:
		Accumulation_poly(RestrictedVoronoiDiagram* rvd) : rvd_(rvd) {} 
        void operator()(unsigned int v, TopoPolyMesh* M) const {
			for(int i=0; i<(int)M->nb_facets(); ++i) {
                unsigned int facet_base = M->facet_begin(i) ;
                unsigned int facet_n = M->facet_size(i) ;
                for(unsigned int i=1; i<facet_n-1; i++) {
					rvd_->rvcells()[v].insert_facet(Geex::Facet(M->vertex(0), M->vertex(i), M->vertex(i+1))) ;
                }
			}
		}
	} ;

	void RestrictedVoronoiDiagram::clip_poly_linear(TopoPolyMesh* poly, MyDelaunay* dt, DelaunaySkeleton* skel) {
		Accumulation_poly action(this) ;
		RestrictedVoronoiDiagram_poly RVD_;
		RVD_.set_exact(false);
		RVD_.compute(poly, dt, skel, action) ;
//		RVD_.compute(&(dt->boundary_topo_parts_[1]), dt, skel, action) ;
	}

	void RestrictedVoronoiDiagram::clip_poly_kdtree(TopoPolyMesh* poly, MyDelaunay* dt, DelaunaySkeleton* skel, ANNKdTree_3* kdtree) {
		Accumulation_poly action(this) ;
		RestrictedVoronoiDiagram_poly_kdtree RVD_;
		RVD_.set_exact(false);
		RVD_.compute(poly, dt, skel, kdtree_, action) ;
	}


	// -------------------------------------------------------------------------------------
	// ----------------------------  Extract mesh from RVD ---------------------------------
	// -------------------------------------------------------------------------------------

	void RestrictedVoronoiDiagram::extract_mesh() {
		std::vector<std::set<int> > adj_cells ;
		Geex::TriMesh               newmesh ;
		Geex::TriMesh               clipped_mesh ;
		std::vector<int>            VoroCell_idx ;
		int                         nv = (int)rvcells_.size() ;


		std::cout << "extract_mesh() " << std::endl ;

		remesh_.clear() ;
		clipped_mesh.clear() ;

		if(dt_->is_clip_feature()) {
			this->collect_crease_segments() ;
		//	clip_corner() ;
		//	clip_feature() ;
		}

		for(int i=0; i<nv; ++i) {
			for(int j=0; j<(int)rvcells_[i].mesh().size(); ++j) {
				newmesh.push_back(rvcells_[i].mesh()[j]) ; 
				VoroCell_idx.push_back(i) ;
			}
		}

		std::cout << "begin merging same vertices..." << std::endl ;

		newmesh.merge_same_vertices(64) ;
		newmesh.merge_same_vertices(63) ;

		std::cout << "finish merging same vertices..." << std::endl ;

		for(int i=0; i<newmesh.nb_vertices(); ++i) {
			std::set<int> v_adj_cells ; 	
			MeshVertex& v = newmesh.vertex(i) ;
			for(unsigned int j=0; j<v.faces_.size(); ++j)
				v_adj_cells.insert(VoroCell_idx[v.faces_[j]]) ; //newmesh[v.faces_[j]].cell_idx) ;
			adj_cells.push_back(v_adj_cells) ;
		}

		clipped_mesh = newmesh ; 

		//		remesh_ = newmesh ;
		//	remesh_.clear() ;
		for(unsigned int i=0; i<adj_cells.size(); ++i) {
			if(adj_cells[i].size() ==3) {
				vec3 pts[3] ;
				int  j = 0 ;
				for(std::set<int>::iterator it=adj_cells[i].begin(); it!=adj_cells[i].end(); ++it) {
					vec3 g ;
					double V ;
					rvcells_[*it].get_centroid(g, V) ;
					pts[j++] = g ; 
					//to_geex(vcells_[*it].vh_->point()) ;
				}
				// compute average normal of current vertex
				vec3 aven(0, 0, 0) ;
				std::vector<int>& vfacets = clipped_mesh.vertex(i).faces_ ;
				for(unsigned int j=0; j<vfacets.size(); ++j) {
					Geex::Facet& F = clipped_mesh[vfacets[j]] ;
					aven = aven + tri_area(F.vertex[0], F.vertex[1], F.vertex[2]) * F.normal() ;
				}

				//		Geex::Facet& F = vcells_[*adj_cells[i].begin()].mesh_[0] ;
				//		if(Geex::dot(Geex::cross(pts[1]-pts[0], pts[2]-pts[0]), F.normal()) < 0 )
				vec3 center = 1.0/3.0*(pts[0]+pts[2]+pts[1]) ;
				vec3 normal ;
				dt_->project_to_mesh(center, normal) ;

				if(Geex::dot(Geex::cross(pts[1]-pts[0], pts[2]-pts[0]), aven) <= 0 )
					remesh_.push_back(Geex::Facet(pts[0], pts[1], pts[2], 1, 1, 1)) ;
				else 
					remesh_.push_back(Geex::Facet(pts[0], pts[2], pts[1], 1, 1, 1)) ;		
			}
			else if (adj_cells[i].size()>3) {
				std::cout << "degenerate cases ... " << adj_cells[i].size() << std::endl ;
			}
		}
		remesh_.merge_same_vertices(64) ;

		std::cout << "finish extracting mesh..." << std::endl ;

		remesh_.measure_quality() ;
		clipped_mesh.clear() ;
	}


	class PrimalEdge {
	public:
		PrimalEdge(int v0=-1, int v1=-1) {
			v0_ = v0 ;
			v1_ = v1 ;
			ref_ = 0 ;
		}

		bool operator == (const PrimalEdge& rhs) const {
			return v0_==rhs.v0_ && v1_==rhs.v1_ ;
		}

		bool operator < (const PrimalEdge& rhs) const {
			if(v0_ < rhs.v0_) 
				return true ;
			else if(v0_ > rhs.v0_) 
				return false ;
			else if(v1_ < rhs.v1_) 
				return true ;
			else 
				return false ;
		}

		void increase() { ref_++ ; }
		void decrease() { ref_-- ; }

	public:
		int v0_ ;
		int v1_ ;
		int ref_ ;
	} ;

	void RestrictedVoronoiDiagram::extract_dual() {

        std::cout << "extract_dual() " << std::endl ;

		typedef PrimalEdge Edge ;
        std::set< Edge > graph ;
        std::map<int, int> vmap ;
        int vidx = 0 ;

        remesh_.clear_all() ;

        std::vector<vec3> vertices_normal(rvcells_.size()); //collect vertices normal for correcting orientation
        // build connectivity graph of mesh
        for(unsigned int i=0; i<rvcells_.size(); ++i) {
            if(!rvcells_[i].is_boundary())
            {
                vmap[i] = - 1; 
                continue ;
            }

            Geex::TriMesh& m = rvcells_[i].mesh_ ;
            for(unsigned int j=0; j<m.size(); ++j) {
                Geex::Facet& F = m[j] ;
                for(int k=0; k<3; ++k) {
                    if(F.nclip_[k]>=0) 
                        if(i < F.nclip_[k])
                            graph.insert(Edge(i, F.nclip_[k])) ;
                        else
                            graph.insert(Edge(F.nclip_[k], i)) ;
                }
            }
            dt_->project_to_mesh(rvcells_[i].point(), vertices_normal[vidx]) ;
            remesh_.add_vertex(rvcells_[i].point()) ;
            vmap.insert(std::pair<int,int>(i, vidx++)) ;
            
        }

        int v0,v1,v2;
        Edge edges[3] ;
        
        for (MyDelaunay::Finite_facets_iterator fit=dt_->finite_facets_begin(); fit!=dt_->finite_facets_end(); fit++)
        {
            MyDelaunay::Cell_handle ch= fit->first;
            v0 = ch->vertex(MyDelaunay::faces_to_vertices_[fit->second][0])->idx();
            v1 = ch->vertex(MyDelaunay::faces_to_vertices_[fit->second][1])->idx();
            v2 = ch->vertex(MyDelaunay::faces_to_vertices_[fit->second][2])->idx();

            if (vmap[v0] == -1 || vmap[v1] == -1 || vmap[v2] == -1)
            {
                continue;
            }
            

            edges[0] = v0<v1 ? Edge(v0, v1) : Edge(v1, v0) ;
            if ( graph.find(edges[0])==graph.end() )
                continue;
            edges[1] = v1<v2 ? Edge(v1, v2) : Edge(v2, v1) ;
            if ( graph.find(edges[1])==graph.end() )
                continue;
            edges[2] = v2<v0 ? Edge(v2, v0) : Edge(v0, v2) ;
            if ( graph.find(edges[2])==graph.end() )
                continue;
          
            if ( Geex::dot(
                Geex::cross(rvcells_[v1].point()-rvcells_[v0].point(), rvcells_[v2].point()-rvcells_[v0].point()), 
                vertices_normal[vmap[v0]]+vertices_normal[vmap[v1]]+vertices_normal[vmap[v2]]) <= 0
                )
            {
                remesh_.push_back(Geex::Facet(rvcells_[v0].point(), rvcells_[v1].point(), rvcells_[v2].point(), 1, 1, 1, vmap[v0], vmap[v1], vmap[v2])) ;
            }
            else
            {
                remesh_.push_back(Geex::Facet(rvcells_[v1].point(), rvcells_[v0].point(), rvcells_[v2].point(), 1, 1, 1, vmap[v1], vmap[v0], vmap[v2])) ;
            }

        }
        
 		// fill holes
        std::cout <<"done \n";
		remesh_.measure_quality() ;
    }

    void RestrictedVoronoiDiagram::voronoi_filtering(double input_theta_degree)
    {
        std::cout << "voronoi_filtering() " << std::endl ;

        remesh_.clear_all();

        std::vector<vec3> pts;
        pts.reserve(rvcells_.size());
        for(unsigned int i=0; i<rvcells_.size(); ++i) {
            if (rvcells_[i].is_boundary())
            {
                vec3 tv = rvcells_[i].point();
                pts.push_back(tv);
                remesh_.add_vertex(tv);
            }
        }

        Voronoi_Filtering vfilter(pts, input_theta_degree);
        std::vector<vec3> vertices_normal(pts.size());
        for (unsigned int i=0; i<pts.size(); ++i)
        {
            dt_->project_to_mesh(pts[i], vertices_normal[i]) ;
        }

        const std::list<VF_Face*>& vfacelist = vfilter.get_face_indice();
        int v0, v1, v2;
        for (std::list<VF_Face*>::const_iterator fiter = vfacelist.begin();
            fiter != vfacelist.end(); fiter++)
        {
            VF_Face* vf = *fiter;
            if (!vf->valid)
            {
                continue;
            }
			
            v0 = vf->indices[0];
            v1 = vf->indices[1];
            v2 = vf->indices[2];

            if ( Geex::dot(
                Geex::cross(pts[v1]-pts[v0], pts[v2]-pts[v0]), 
                vertices_normal[v0]+vertices_normal[v1]+vertices_normal[v2]) <= 0
                )
            {
                remesh_.push_back(Geex::Facet(pts[v0], pts[v1], pts[v2], 1, 1, 1, v0, v1, v2)) ;
            }
            else
            {
                remesh_.push_back(Geex::Facet(pts[v1], pts[v0], pts[v2], 1, 1, 1, v1, v0, v2)) ;
            }
            

            //if (vf->is_outside)
            //    remesh_.push_back(Geex::Facet(pts[v0], pts[v1], pts[v2], 1, 1, 1, v0, v1, v2)) ;
            //else
            //    remesh_.push_back(Geex::Facet(pts[v1], pts[v0], pts[v2], 1, 1, 1, v1, v0, v2)) ;
        }
        
        if (!remesh_.empty())
            remesh_.measure_quality() ;
        else
            std::cerr <<"no face has been added !" << std::endl;
    }

    void RestrictedVoronoiDiagram::alpha_shape()
    {
        typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

        typedef CGAL::Alpha_shape_vertex_base_3<K>               Vb;
        typedef CGAL::Triangulation_hierarchy_vertex_base_3<Vb>  Vbh;
        typedef CGAL::Alpha_shape_cell_base_3<K>                 Fb;
        typedef CGAL::Triangulation_data_structure_3<Vbh,Fb>     Tds;
        typedef CGAL::Delaunay_triangulation_3<K,Tds>            Delaunay;
        typedef CGAL::Triangulation_hierarchy_3<Delaunay>        Delaunay_hierarchy;
        typedef CGAL::Alpha_shape_3<Delaunay_hierarchy>          Alpha_shape_3;

        typedef K::Point_3                                  Point;
        typedef Alpha_shape_3::Alpha_iterator               Alpha_iterator;
        typedef Alpha_shape_3::NT                           NT;

        std::cout << "alpha_shape() " << std::endl ;

        std::vector<Point> pts;
        pts.reserve(rvcells_.size());
        for(unsigned int i=0; i<dt_->number_of_vertices(); ++i) {
            if (rvcells_[i].is_boundary())
               pts.push_back(rvcells_[i].vh()->point());
        }
        Delaunay_hierarchy dt;
        dt.insert(pts.begin(), pts.end());

        // compute alpha shape
        Alpha_shape_3 as(dt);
        std::cerr << "Alpha shape computed in REGULARIZED mode by default."  << std::endl;

        // find optimal alpha values
        Alpha_iterator opt = as.find_optimal_alpha(1);
        std::cerr << "Optimal alpha value to get one connected component is "
            <<  *opt    << std::endl;

        as.set_alpha(*opt);

        remesh_.clear_all();
        
        std::vector<vec3> vertices_normal(pts.size());
        std::map<Alpha_shape_3::Vertex_handle, int> vertex_id;
        int count = 0;
        for (Alpha_shape_3::Finite_vertices_iterator viter = as.finite_vertices_begin();
            viter != as.finite_vertices_end(); viter++)
        {
            vec3 p = to_geex( viter->point() );
            dt_->project_to_mesh(p, vertices_normal[count]) ;
            remesh_.add_vertex(p);
            vertex_id[viter] = count++;
        }

        for( Alpha_shape_3::Finite_facets_iterator fit = as.finite_facets_begin(); 
              fit != as.finite_facets_end() ; ++fit)
          {
                if (as.classify(*fit, *opt) == Alpha_shape_3::REGULAR  )
                {
                    const Alpha_shape_3::Cell_handle& ch = fit->first;
                    const int index = fit->second;

                    Alpha_shape_3::Vertex_handle a = ch->vertex((index+1)&3);
                    Alpha_shape_3::Vertex_handle b = ch->vertex((index+2)&3);
                    Alpha_shape_3::Vertex_handle c = ch->vertex((index+3)&3);

                    vec3 p0 = to_geex(a->point());
                    vec3 p1 = to_geex(b->point());
                    vec3 p2 = to_geex(c->point());

                    if ( Geex::dot(
                        Geex::cross(p1-p0, p2-p0), 
                        vertices_normal[vertex_id[a]]+vertices_normal[vertex_id[b]]+vertices_normal[vertex_id[c]]) <= 0
                        )
                    {
                        remesh_.push_back(Geex::Facet(p0, p1, p2, 1, 1, 1, vertex_id[a], vertex_id[b], vertex_id[c])) ;
                    }
                    else
                    {
                        remesh_.push_back(Geex::Facet(p1, p0, p2, 1, 1, 1, vertex_id[b], vertex_id[a], vertex_id[c])) ;
                    }
                }
          }

        std::cerr << " Done !" << std::endl; 
        if (!remesh_.empty())
            remesh_.measure_quality() ;
        else
            std::cerr <<"no face has been added !" << std::endl;



    }
} // end of namespace Geex

