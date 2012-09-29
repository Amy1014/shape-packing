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

#ifndef __RVD_GRAPHICS__
#define __RVD_GRAPHICS__


#include <Geex/graphics/opengl.h>

namespace Geex {

    class RVDGeometry ;

    enum SurfaceMode { None, Plain, Cells, Density, Energy, Linear, Connected} ;
#define SurfaceModeNames "None" | "Plain" | "Cells" | "Density" | "Energy" | "Linear" | "Connected"

    enum SliceAxis { AXIS_X, AXIS_Y, AXIS_Z } ;
#define SliceAxisNames "X" | "Y" | "Z" 

    class RVDGraphics {
    public:
        RVDGraphics(RVDGeometry* geometry) ;
        ~RVDGraphics() ;
        void draw(bool geometry_is_dirty = false) ;
        GLboolean& use_display_lists() { return use_display_lists_ ; }
        GLenum& surface_mode()         { return surface_mode_ ; }
        GLboolean& show_domain_mesh()  { return show_domain_mesh_ ; }
        GLfloat& vertices_size()       { return vertices_size_ ; }
        GLboolean& show_rvd_mesh()     { return show_rvd_mesh_ ; }
        GLboolean& show_rvd_submesh()  { return show_rvd_submesh_ ; }
        GLboolean& show_primal()       { return show_primal_ ; }
        GLboolean& show_primal_mesh()  { return show_primal_mesh_ ; }
        GLboolean& show_non_manifold() { return show_non_manifold_ ; }
        GLfloat& shrink()              { m_geometry_dirty_ = true; return shrink_ ; }
        GLfloat& id_threshold()         { return id_threshold_ ; }
        GLboolean& shiny()             { return shiny_ ; }
		GLboolean& show_vertices()     { return show_vertices_ ; }

        RVDGeometry* geometry()        { return geometry_ ; }


        static void gl_random_color() ;
        static void gl_randomize_colors(int index = 0) ;
        static void gl_random_color(int index) ;

    protected:
        void draw_domain() ; 
        void draw_domain_mesh() ;
        void draw_vertices() ;
        void draw_rvd() ;
        void draw_rvd_mesh() ;
        void draw_primal() ;
        void draw_primal_mesh() ;
        void create_colormap_texture() ;

    protected:

        inline void begin_list(GLuint& L) {
            if(use_display_lists_) {
                delete_list(L) ;
                L = glGenLists(1);
                glNewList(L, GL_COMPILE_AND_EXECUTE);
            }
        }

        inline void end_list() {
            if(use_display_lists_) {
                glEndList() ;
            }
        }

        inline void delete_list(GLuint& L) {
            if (glIsList(L)) { 
                glDeleteLists(L, 1); 
            }
            L = 0 ;
        }

        inline bool playback_list(GLuint& L) {
            if(!use_display_lists_) {
                return false ;
            }
            if(!glIsList(L)) {
                return false ;
            }
            glCallList(L) ;
            return true ;
        }

        void delete_lists() ;

    private:
        RVDGeometry* geometry_ ;
        GLenum    surface_mode_ ;
        GLboolean show_domain_mesh_ ;
        GLfloat vertices_size_ ;
        GLboolean show_rvd_mesh_ ;
        GLboolean show_rvd_submesh_ ;
        GLboolean show_primal_ ;
        GLboolean show_primal_mesh_ ;
        GLfloat shrink_ ;
	GLfloat id_threshold_ ;
        GLboolean show_non_manifold_ ;
        GLuint colormap_texture_ ;
        GLboolean shiny_ ;
		GLboolean show_vertices_;

        GLboolean use_display_lists_ ;
        GLuint u_primal_;
        GLuint u_primal_mesh_;
        GLuint u_vertices_;
        GLuint u_singular_;
        GLuint u_singular_mesh_;
        GLuint u_rvd_none_, u_rvd_cell_, u_rvd_connected_, u_rvd_energy_, u_rvd_linear_, u_rvd_shrink_, u_rvd_mesh_, u_rvd_shrink_mesh_;
        GLboolean m_geometry_dirty_; //only used for singular triangles
    } ;

}

#endif
