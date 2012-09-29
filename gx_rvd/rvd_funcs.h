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

#ifndef __CVT_FUNC__
#define __CVT_FUNC__

#include "generated/P2P.h"
#include "generated/P0.h"
#include "generated/P1.h"
#include "generated/P2.h"
#include "generated/P3.h"
#include "generated/P2tet.h"

namespace Geex {

    class CVTFuncs {
    protected:
        Function* new_P_func(int point_type) {
            Function* result = nil ;
            switch(point_type) {
            case  0: result = new P0() ;  break ;
            case  1: result = new P1() ;  break ;
            case  2: result = new P2() ;  break ;
            case  3: result = new P3() ;  break ;
            case -1: result = new P2P() ; break ;
            default: gx_assert_not_reached ; break ;
            }
            return result ;
        }
    } ;

    //---------------------------------------------------------------------------

    class CVTDualTetFuncs : public CVTFuncs {
    public:
        Function* get_template(unsigned int i, unsigned int j, unsigned int k) {
            gx_assert(i<4 && j<4 && k<4) ;
            return func_table_[i][j][k] ;
        }
    protected:
        Function_var func_table_[4][4][4] ;
    } ;

    //---------------------------------------------------------------------------

    template <class TETFUNC> class GenericCVTDualTetFuncs : public CVTDualTetFuncs {
    public:
        GenericCVTDualTetFuncs(int nb_additional_parameters = 0) { 
            for(int i=0; i<4; i++) {
                for(int j=0; j<4; j++) {
                    for(int k=0; k<4; k++) {
                        Function_var F = new TETFUNC() ;
                        std::vector<Function_var> G ;
                        G.push_back(new P2P()) ;
                        G.push_back(new_P_func(i)) ;
                        G.push_back(new_P_func(j)) ;
                        G.push_back(new_P_func(k)) ;
                        if(nb_additional_parameters != 0) {
                            G.push_back(new Identity(nb_additional_parameters)) ;
                        }
                        func_table_[i][j][k] = new CompositeFunction(F, G) ;
                    }
                }
            }
        }
    } ;

    //---------------------------------------------------------------------------

    class CVTDualEdgeFuncs : public CVTFuncs {
    public:
        Function* get_template(unsigned int i, unsigned int j) {
            gx_assert(i<4 && j<4) ;
            return func_table_[i][j] ;
        }
    protected:
        Function_var func_table_[4][4] ;
    } ;

    //---------------------------------------------------------------------------

    template <class EDGEFUNC> class GenericCVTDualEdgeFuncs : public CVTDualEdgeFuncs {
    public:
        GenericCVTDualEdgeFuncs(int nb_additional_parameters = 0) { 
            for(int i=0; i<4; i++) {
                for(int j=0; j<4; j++) {
                    Function_var F = new EDGEFUNC() ;
                    std::vector<Function_var> G ;
                    G.push_back(new P2P()) ;
                    G.push_back(new_P_func(i)) ;
                    G.push_back(new_P_func(j)) ;
                    if(nb_additional_parameters != 0) {
                        G.push_back(new Identity(nb_additional_parameters)) ;
                    }
                    func_table_[i][j] = new CompositeFunction(F, G) ;
                }
            }
        }
    } ;

    //---------------------------------------------------------------------------

}

#endif
