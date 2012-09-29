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

#ifndef __RVD_PRIMAL__
#define __RVD_PRIMAL__

#include <map>
#include <algorithm>

namespace Geex {

    struct trindex {
        trindex(
            unsigned int i,
            unsigned int j,
            unsigned int k
        ) {
            indices[0] = i ; indices[1] = j ; indices[2] = k ;
            std::sort(indices, indices+3) ;
        }
        bool operator<(const trindex& rhs) const {
            for(unsigned int i=0; i<3; i++) {
                if(indices[i] < rhs.indices[i]) { return true ;  }
                if(indices[i] > rhs.indices[i]) { return false ; }
            }
            return false ;
        }
        unsigned int indices[3] ;
    } ;

    struct bindex {
        bindex(
            unsigned int i, 
            unsigned int j
        ) {
            if(i < j) {
                indices[0] = i ;
                indices[1] = j ;
            } else {
                indices[0] = j ;
                indices[1] = i ;
            }
        }
        bool operator<(const bindex& rhs) const {
            if(indices[0] < rhs.indices[0]) { return true ;  }
            if(indices[0] > rhs.indices[0]) { return false ; }
            if(indices[1] < rhs.indices[1]) { return true ;  }
            return false ;
        }
        unsigned int indices[2] ;
    } ;

    //=============================================================

    class CountPrimalSingularTriangles {
    public:
        CountPrimalSingularTriangles(std::map< trindex, int >& triangles) : triangles_(triangles) {
        }
        void operator()(unsigned int i, unsigned int j, unsigned int k) const {
            triangles_[trindex(i,j,k)]++ ;
        }
        std::map< trindex, int >& triangles_ ;        
    } ;

    class CountPrimalSingularEdges {
    public:
        CountPrimalSingularEdges(std::map< bindex, unsigned int >& edges) : edges_(edges) {
        }
        void operator()(unsigned int i, unsigned int j, unsigned int k) const {
            edges_[bindex(i,j)]++ ;
            edges_[bindex(j,k)]++ ;
            edges_[bindex(k,i)]++ ;
        }
        std::map< bindex, unsigned int >& edges_ ;        
    } ;

}

#endif
