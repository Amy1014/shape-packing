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

#include <Geex/CVT/geometry.h>
#include <glut_viewer/glut_viewer.h>
#include "delaunay_cvt.h"
#include "lloyd_energy.h"
#include "rvd_ccvt.h"
//#include <Geex/basics/stopwatch.h>

bool first = true;
int  curriter = 0 ;
int  call_nb_iter = 0;
int  currfundgraditer = 0 ;

////////////////////////////////////////////////////////////////////////// 
void funcgrad_ccvt(int N, double* x, double& f, double* g)
{
	Geex::DelaunayCVT* cvt = Geex::DelaunayCVT::instance() ;

	if(cvt->constrained()) {
		Geex::RestrictedVoronoiDiagram* rvd = cvt->delaunay()->rvd() ;
		rvd->get_constrained_fgs(N, x, &f, g, cvt->weighted()) ;
	}
	else {
		cvt->delaunay()->set_vertices(N, x, true) ;
		Geex::RestrictedVoronoiDiagram* rvd = cvt->delaunay()->rvd() ;
		rvd->get_fgs(&f, g, cvt->weighted()) ;
	}

	call_nb_iter++;
}

//////////////////////////////////////////////////////////////////////////  
void newiteration_ccvt(int n, const double* x, double f, const double* g, double gnorm)
{
	std::cout.precision(16);
	std::cout << curriter <<": " << call_nb_iter <<" " << std::scientific << f <<" " << gnorm  << std::endl;
	curriter ++;
}



//////////////////////////////////////////////////////////////////////////
namespace Geex {

	DelaunayCVT* DelaunayCVT::instance_ = nil ;

	DelaunayCVT::DelaunayCVT(MyDelaunay* delaunay) : delaunay_(delaunay) {
		gx_assert(instance_ == nil) ;
		instance_ = this ;
		weighted_  = GL_FALSE ;
		constrained_ = GL_FALSE ;
        optimizer_mode_ = HLBFGS;
        opt_m_ = 7;
	}

	DelaunayCVT::~DelaunayCVT() { 
		instance_ = nil ; 
	}

	// ---------------------------------------------------------------------
	// ----------------------- Lloyd CVT iteration -------------------------
	// ---------------------------------------------------------------------

	void DelaunayCVT::lloyd(int nb_iter) {
		for(int k=0; k<nb_iter; k++) {
			std::vector<vec3> new_points ;
			real f = 0, gnorm2 = 0 ;
			RestrictedVoronoiDiagram* rvd = delaunay_->rvd() ;
			std::vector<double> Vs ; 
			if(constrained_) 
				rvd->get_constrained_centroids(new_points, Vs, weighted()) ;
			else {
				int nn = delaunay_->nb_vertices()-delaunay_->nb_corners();
				double * g = new double[nn*3] ;
				rvd->get_fgs(&f, g, weighted()) ;
				rvd->get_centroids(new_points, Vs, weighted()) ;
				for(int i=0; i<nn*3; ++i)
					gnorm2 += g[i]*g[i] ;
				delete [] g ;
			}
			delaunay_->set_vertices(new_points, true) ;
			std::cout.precision(16) ;
			std::cout << k <<": " <<std::scientific << f <<" " << std::scientific << ::sqrt(gnorm2) << std::endl ;
		}
	}

	//----------------------------------------------------------------------------------------------------

	void DelaunayCVT::newton_lloyd(int nb_iter) {
		int nv = delaunay_->nb_vertices() - delaunay_->nb_corners() ;
		int n = nv * 3 ;
		int m = (int)opt_m_ ;
        if (optimizer_mode_ == LBFGSB && m == 0)
        {
            m = 1;
        }

		double *x = new double[n];

		for(unsigned int i=0; i< nv; i++) {
			MyDelaunay::Vertex_handle it = delaunay_->all_vertices_[i] ;
			x[3*i  ] = it->point().x();
			x[3*i+1] = it->point().y();
			x[3*i+2] = it->point().z();
		}

		curriter = 1;
		call_nb_iter = 0;

		double epsg = 0, epsf=0, epsx=0;

        Optimizer* opt = nil;
        std::cerr << "Starting Newton (warming up...)" << std::endl ;
        switch(optimizer_mode_) {
        case LBFGSB:
            opt = new LBFGSBOptimizer();
            break;
        case HLBFGS:
            opt = new HLBFGSOptimizer() ;
            break ;
        case HM1QN3:
            opt = new HLBFGSOptimizer() ;
            static_cast<HLBFGSOptimizer*>(opt)->set_m1qn3(true);
            break ;
        case HCG:
            opt = new HLBFGSOptimizer() ;
            static_cast<HLBFGSOptimizer*>(opt)->set_cg(true);
            break ;
        default:
            gx_assert_not_reached ;
        }

        if (optimizer_mode_ == LBFGSB) {
            int *rhs = new int[n];
            memset(rhs, 0, sizeof(int)*n);
            double *lu = new double[n];
            memset(lu, 0, sizeof(double)*n);
            static_cast<LBFGSBOptimizer*>(opt)->set_nbd(n, rhs);
            static_cast<LBFGSBOptimizer*>(opt)->set_l(n, lu);
            static_cast<LBFGSBOptimizer*>(opt)->set_u(n, lu);
            delete[] rhs;
            delete[] lu;
        }

		opt->set_epsg(epsg) ;
		opt->set_epsf(epsf) ;
		opt->set_epsx(epsx) ;
		
		opt->set_M(m) ;
		opt->set_N(n) ;
		opt->set_max_iter(nb_iter) ;

		opt->set_newiteration_callback(newiteration_ccvt) ;
		opt->set_funcgrad_callback(funcgrad_ccvt) ;
		
		opt->optimize(x) ;
		
		//set_vertices(x) ;
		delete opt;
		delete [] x;
	}


//	void DelaunayCVT::set_vertices(const double* x) {
//		int N = ((int)delaunay_->all_vertices_.size() - (int)delaunay_->nb_corners()) * 3 ;
//		delaunay_->set_vertices(N, x, true);
////		delaunay_->end_insert(true) ;
//	}
	//----------------------------------------------------------------------------------------------------
}
