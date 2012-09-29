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


#include "cvt.h"

namespace Geex {
    
    CVT::CVT() : DelaunayGraphics(this), DelaunayCVT(this) {
//        set_culling(GL_FRONT_AND_BACK) ;
        set_culling(GL_BACK) ;
    }

    CVT::~CVT() {
        
    }

    void CVT::do_draw() {
        DelaunayGraphics::draw() ;
        // glDisable(GL_LIGHTING) ;
        // glLineWidth(3) ;
        // glBegin(GL_LINES) ;
        // glColor3f(1, 0, 0) ;
        // glVertex3f(0, 0, 0) ;
        // glVertex3f(1, 0, 0) ;
        // glColor3f(0, 1, 0) ;
        // glVertex3f(0, 0, 0) ;
        // glVertex3f(0, 1, 0) ;
        // glColor3f(0, 0, 1) ;
        // glVertex3f(0, 0, 0) ;
        // glVertex3f(0, 0, 1) ;
        // glEnd() ;
    }

    void CVT::get_bbox(
            real& x_min, real& y_min, real& z_min,
            real& x_max, real& y_max, real& z_max
    ) {
        MyDelaunay::get_bbox(x_min, y_min, z_min, x_max, y_max, z_max) ;
    }


}
