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

#ifndef __DELAUNAY_CVT__
#define __DELAUNAY_CVT__

#include "delaunay.h"
#include <Geex/graphics/opengl.h> // use GLboolean

#include <Geex/numerics/optimizer.h>
#include <Geex/numerics/lbfgs_interface.h>
#include <ostream>

namespace Geex {

    enum OptimizerMode { LBFGSB, HLBFGS, HM1QN3, HCG, HLBFGS_HESS } ;
#define OptimizerModeNames  "LBFGSB" | "HLBFGS" | "HM1QN3" | "HCG" | "HESSIAN"

    class DelaunayCVT {
    public:
        DelaunayCVT(MyDelaunay* delaunay) ;
        ~DelaunayCVT() ;
	
	// CVT
        void lloyd(int nb_iter = 1) ;
        void newton_lloyd(int nb_iter = 1) ;

	// ODT
        void lloyd_odt(int nb_iter = 1) ;
        void lloyd_weighted_odt(int nb_iter = 1) ;
            
	// Properties
        static DelaunayCVT* instance() { return instance_ ; }
        std::vector<MyDelaunay::Vertex_handle>& all_vertices() { return delaunay_->all_vertices_ ; }
		
        void set_vertices(const double* x) ;    
		//new insertion code
		void prepare_sort_points(int N, const double *x, std::vector<MyPoint>& mypointvec);
		void inserte_sort_points(const std::vector<MyPoint>& mypointvec, bool tag);
		
        inline GLboolean& weighted()    { return weighted_ ; }
        inline GLboolean& constrained() { return constrained_ ; }
		
        MyDelaunay* delaunay(){return delaunay_;}
	
        GLenum& optimizer_mode() { return optimizer_mode_ ; }
        float &set_opt_m() {return opt_m_ ; }
    private:
        MyDelaunay* delaunay_ ;
        static DelaunayCVT* instance_ ;
        GLboolean weighted_ ;
        GLboolean constrained_ ;

        GLenum optimizer_mode_ ;
        float opt_m_;
    } ;
    

}

#endif
