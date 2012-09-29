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
#include <Geex/CVT/delaunay_skel_CGAL.h>
#include <Geex/CVT/RVD.h>
#include <omp.h>

namespace Geex {

	RestrictedVoronoiDiagram::RestrictedVoronoiDiagram(MyDelaunay* dt, Geex::TriMesh* mesh, RVDMode mode) 
	: dt_(dt)
	, mesh_(mesh)
	, rvd_mode_(mode) {
		kdtree_ = nil ;
		init_vcells() ;
	}

	RestrictedVoronoiDiagram::~RestrictedVoronoiDiagram() {
		if(kdtree_ != nil) delete kdtree_ ;
	}

	void RestrictedVoronoiDiagram::init_vcells() {
		int nv = (int)dt_->vertices().size() ; 
		rvcells_.clear() ;
		rvcells_.resize(nv) ; 

		for(int i=0; i<nv; ++i) {
			rvcells_[i].vh() = dt_->vertices()[i] ; // = RestrictedVoronoiCell(dt_->vertices()[i]) ;
		}
	}

    void RestrictedVoronoiDiagram::build_kdtree() 
    {
		int i, nv = dt_->nb_vertices() ;
		double * data = new double[nv*3];

#pragma omp parallel for shared(nv) private(i)
		for(i=0; i<nv; ++i) {
			data[i*3+0] = dt_->all_vertices_[i]->point().x();
			data[i*3+1] = dt_->all_vertices_[i]->point().y();
			data[i*3+2] = dt_->all_vertices_[i]->point().z();
		}

		if(kdtree_ != nil) delete kdtree_;
		kdtree_ = new ANNKdTree_3(data, nv, 3);
		delete [] data;
    }	

	// -------------------------------------------------------------------------------------
	// --------------------------------- Interface for CVT ---------------------------------
	// -------------------------------------------------------------------------------------

	void RestrictedVoronoiDiagram::get_centroids(std::vector<vec3>& gs, std::vector<double>& Vs, bool weight) {
		int i, nv = (int)rvcells_.size() - dt_->nb_corners();
		gs.resize(nv) ;
		Vs.resize(nv) ;

		if(weight) {
#pragma omp parallel for shared(nv) private(i)
			for(i=0; i<nv; ++i) {
				rvcells_[i].get_centroidw(gs[i], Vs[i]) ;
			}
		} else {
#pragma omp parallel for shared(nv) private(i)
			for(i=0; i<nv; ++i) {
				rvcells_[i].get_centroid(gs[i], Vs[i]) ;
			}
		}
	}

	void RestrictedVoronoiDiagram::get_constrained_centroids(std::vector<vec3>& gs, std::vector<double>& Vs, bool weight) {
		int i, nv = (int)rvcells_.size() - dt_->nb_corners() ; // corner will not be changed in optimization
		vec3 n ; 
		gs.resize(nv) ;
		Vs.resize(nv) ;

//#pragma omp parallel for shared(nv) private(i, n)
		if(weight) {
			for(i=0; i<nv; ++i) {
				rvcells_[i].get_centroidw(gs[i], Vs[i]) ;
				if(rvcells_[i].is_corner()) {}
				else if(rvcells_[i].is_crease()) {
					rvcells_[i].project_to_crease(gs[i], n) ;
				} 
				else {
					gs[i] = project_to_mesh(gs[i], n) ;
				}
			}
		} else {
			for(i=0; i<nv; ++i) {
				rvcells_[i].get_centroid(gs[i], Vs[i]) ;
				if(rvcells_[i].is_corner()) {}
				else if(rvcells_[i].is_crease()) {
					rvcells_[i].project_to_crease(gs[i], n) ;
				} 
				else {
					gs[i] = project_to_mesh(gs[i], n) ;
				}
			}
		}
	}

//	void RestrictedVoronoiDiagram::get_gfs(double* f, std::vector<vec3>& grad_fs, bool constrained) {
//		int i, nv = rvcells_.size() ;
//		grad_fs.resize(nv) ;
//
//		*f = 0.0 ;
//
//		if(constrained) {
//		} 
//		else {
////#pragma omp parallel for shared(nv) private(i)
//			for(i=0; i<nv; ++i) {
//				double cur_f ;
//				rvcells_[i].get_fg(cur_f, grad_fs[i]) ;
//				*f+=cur_f ;
//			}
//		}
//	}

	void RestrictedVoronoiDiagram::get_constrained_fgs(int N, double* x_in, double* f, double* g, bool weight) {
		int nv = N/3 ;
		std::vector<vec3> normals ;
		normals.resize(nv) ;
		
		// TODO 
		// The new Optimizer framework in Geex does not allow to modify x
		//double* x = const_cast<double * >(x_in);
		double* x = x_in;

		*f = 0.0 ;

		// re-compute rvd to identify features, project feature seeds onto feature
		if(dt_->is_clip_feature()) { 

			// identify features
//			dt_->set_vertices(N, x) ;

			// project vertices
			for(int i=0; i<nv; ++i) {
				vec3 p(x[i*3], x[i*3+1], x[i*3+2]) ;

				//if(rvcells_[i].is_corner()) {
				//	// no corner will be optimized, remove later. --dmyan 2009-02-21
				//	//project_to_mesh(p, normals[i]) ;
				//	p = rvcells_[i].point() ;
				//}
				//else if(rvcells_[i].is_crease()) {
				//	p = rvcells_[i].project_to_crease(p, normals[i]) ;
				//	// gradient is already adjusted in get_fg
				//	//cur_grad = dot(cur_grad, surfacenormal[i])*surfacenormal[i] ;
				//}
				//else { // smooth or out layer.
					p= project_to_mesh(p, normals[i]) ;
				//}

				// update x
				x[i*3] = p[0] ; 
				x[i*3+1] = p[1] ; 
				x[i*3+2] = p[2] ;
			}

			// re-compute RVD based on new constrained
			dt_->set_vertices(N, x) ;

			// compute gradient
			double cur_f ;
			Geex::vec3 cur_grad ;
			for(int i=0; i<nv; ++i) {
				if(weight) 
					rvcells_[i].get_fgw(cur_f, cur_grad) ;
				else
					rvcells_[i].get_fg(cur_f, cur_grad) ;
				*f += cur_f ;

				// gradient is already computed in get_fg
				if(rvcells_[i].is_corner()) {
					//cur_grad = vec3(0, 0, 0) ;
				}
				else if(rvcells_[i].is_crease()) {
					// gradient is already adjusted in get_fg
					//cur_grad = dot(cur_grad, surfacenormal[i])*surfacenormal[i] ;
				}
				else { // smooth or out layer.
					cur_grad = cur_grad - (dot(cur_grad,normals[i])) * normals[i];
				}

				// update x and g
				g[i*3] = cur_grad[0] ;
				g[i*3+1] = cur_grad[1] ;
				g[i*3+2] = cur_grad[2] ;
			}
		} 
		else { // no feature constrains smooth case
			int i=0 ;
			while(i+2<N) {
				vec3 p(x[i], x[i+1], x[i+2]) ;
				p = project_to_mesh(p, normals[i/3]) ;
				x[i] = p[0] ;
				x[i+1] = p[1] ;
				x[i+2] = p[2] ;
				i+=3 ;
			}

			dt_->set_vertices(N, x) ;

			i=0 ;
			double cur_f ;
			vec3   cur_grad ;
			while(i+2<N) {
				//vec3 p(x[i], x[i+1], x[i+2]) ;
				if(weight) 
					rvcells_[i/3].get_fgw(cur_f, cur_grad) ;
				else
					rvcells_[i/3].get_fg(cur_f, cur_grad) ;
				*f+=cur_f ;
				cur_grad = cur_grad - (dot(cur_grad,normals[i/3])) * normals[i/3];
				g[i] = cur_grad[0] ;
				g[i+1] = cur_grad[1] ;
				g[i+2] = cur_grad[2] ;
				i+=3 ;
			}
		}

		// accumulate energy of corners
		double cur_f ;
		vec3   cur_grad ;
		for(int i=nv; i<(int)rvcells_.size(); ++i) {
			//rvcells_[i/3].get_fg(cur_f, cur_grad) ;
			if(weight)
				rvcells_[i].get_fgw(cur_f, cur_grad) ;
			else
				rvcells_[i].get_fg(cur_f, cur_grad) ;
			*f+=cur_f ;
		}
	}

	void RestrictedVoronoiDiagram::get_fgs(double* f, double* grad_fs, bool weight) {
		int i, nv = (int)rvcells_.size() - dt_->nb_corners() ;

		*f = 0.0 ;
		vec3 grad ;
		double cur_f ;
		std::vector<double> fs ;
		fs.resize(nv) ;
//#pragma omp parallel for shared(nv) private(i, grad, cur_f)
		if(weight) {
			for(unsigned i=0; i<nv; ++i) {
				rvcells_[i].get_fgw(cur_f, grad) ;
				grad_fs[i*3+0] = grad[0] ;
				grad_fs[i*3+1] = grad[1] ;
				grad_fs[i*3+2] = grad[2] ;
				*f += cur_f ;
			}

			// accumulate energy of corners
			for(unsigned i=nv; i<rvcells_.size(); ++i) {
				rvcells_[i].get_fgw(cur_f, grad) ;
				*f += cur_f ;
			}
		} else {
			for(i=0; i<nv; ++i) {
				rvcells_[i].get_fg(cur_f, grad) ;
				grad_fs[i*3+0] = grad[0] ;
				grad_fs[i*3+1] = grad[1] ;
				grad_fs[i*3+2] = grad[2] ;
				*f += cur_f ;
			}

			// accumulate energy of corners
			for(unsigned i=nv; i<rvcells_.size(); ++i) {
				rvcells_[i].get_fg(cur_f, grad) ;
				*f += cur_f ;
			}
		}
	}
} // end of namespace Geex

