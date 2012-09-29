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


#ifndef __DELAUNAY_GRAPHICS__
#define __DELAUNAY_GRAPHICS__

#include "delaunay.h"
#include <Geex/graphics/opengl.h>

namespace Geex {
    
    class DelaunayGraphics {
    public:
        DelaunayGraphics(MyDelaunay* delaunay) ;
        void draw() ;

        GLboolean& show_domain_mesh() { return show_domain_mesh_ ; }
        GLboolean& show_domain() { return show_domain_ ; }
        GLfloat& vertices_size() { return vertices_size_ ; }
        GLboolean& show_mesh() { return show_mesh_ ; }
        GLboolean& colorize() { return colorize_ ; }
        GLboolean& show_primal() { return show_primal_ ; }
        GLboolean& show_inner_cells() { return show_inner_cells_ ; }
        GLboolean& show_clipped_cells() { return show_clipped_cells_ ; }
        GLboolean& clipped_boundary_only() { return clipped_boundary_only_ ; }
		GLboolean& show_surface_clipping() { return show_surface_clipping_ ; }
        GLboolean& show_volume_clipping() { return show_volume_clipping_ ; }
		GLboolean& show_remesh() { return show_remesh_ ; } 
		GLboolean& show_volume() { return show_volume_ ; }
		GLboolean& show_retet() { return show_retet_ ; }
        float& slice_x() { return slice_x_ ; }
        float& slice_y() { return slice_y_ ; }
        float& slice_z() { return slice_z_ ; }
        float& selection() { return selection_ ; }
		GLboolean& show_rvd() { return show_rvd_ ; }
		GLboolean& set_cull_back() {return show_cull_back_; }

    protected:
        void draw_domain_mesh() ;
        void draw_domain() ;
        void draw_vertices() ;
        void draw_selection() ;

        void draw_cells() ;

        void draw_dual_cell(MyDelaunay::Vertex_handle v, bool mesh = false) ;
        void draw_primal_cell(MyDelaunay::Cell_handle c, bool mesh = false) ;
		void draw_primal_cell(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3, bool mesh = false) ;

        //----------------------------------------------------------------------------------------------
        void draw_clip_points() ;
		void draw_seeds() ;
        void draw_surface_clipping() ;        
		void draw_remesh() ;
		//void draw_volume(TetMesh& tetmesh) ;
		//void draw_tetmesh_mesh(std::vector<Geex::Tet>& tets) ;
		//void draw_volume_clipping() ;
		void draw_feature() ;

		//void draw_grid() ;
		//void draw_colormap() ;
        //----------------------------------------------------------------------------------------------

        inline void draw_dual_face(MyDelaunay::Edge e, bool mesh=false) {
			if(!mesh) {
				// The normal corresponds to the direction of the primal edge
				MyDelaunay::Vertex_handle v1 = e.first->vertex(e.second) ;
				MyDelaunay::Vertex_handle v2 = e.first->vertex(e.third) ;
				CGAL_Vector V = v2->point() - v1->point() ;
				glNormal3f(
					float(V.cartesian(0)), float(V.cartesian(1)), float(V.cartesian(2))
				) ;

				glBegin(GL_POLYGON) ;
				MyDelaunay::Cell_circulator ccir = delaunay_->incident_cells(e) ;
				MyDelaunay::Cell_circulator cdone = ccir ;
				do {
					glVertex(ccir->dual) ;
					ccir++ ;
				} while(ccir != cdone) ;
				glEnd() ;
			} else {
				glBegin(GL_LINE_LOOP) ;
				MyDelaunay::Cell_circulator ccir = delaunay_->incident_cells(e) ;
				MyDelaunay::Cell_circulator cdone = ccir ;
				do {
					glVertex(ccir->dual) ;
					ccir++ ;
				} while(ccir != cdone) ;
				glEnd() ;
			}
        }

        void draw_trimesh(TriMesh& m) {
            glBegin(GL_TRIANGLES) ;
            for(unsigned int i=0; i<m.size(); i++) {
                glNormal(-m[i].normal()) ;
                glVertex(m[i].vertex[0]) ;
                glVertex(m[i].vertex[1]) ;
                glVertex(m[i].vertex[2]) ;
            }
            glEnd() ;
        }

        bool check_edge_flag(const Facet& F, int e) {
            switch(F.edge_flag[e]) {
            case 1:
                glColor4f(0.2f, 0.2f, 0.2f, 0.5f) ;
		//		glColor4f(0., 0.9, 0.0, 0.6) ;
                return true ;
            case 2:
                glColor3f(0.5f, 0.5f, 0.5f) ;
                return true ;
			case 3:
				glColor3f(0.1f, 0.1f, 0.9f) ;
				return true ;
            }
            return false ;
        }

        void draw_trimesh_mesh(TriMesh& m) {
            glBegin(GL_LINES) ;
            for(unsigned int i=0; i<m.size(); i++) {
                const Facet& F = m[i] ;
                if(check_edge_flag(F,2)) { glVertex(F.vertex[0]) ; glVertex(F.vertex[1]) ; }
                if(check_edge_flag(F,0)) { glVertex(F.vertex[1]) ; glVertex(F.vertex[2]) ; }
                if(check_edge_flag(F,1)) { glVertex(F.vertex[2]) ; glVertex(F.vertex[0]) ; }
            }
            glEnd() ;
        }

	void draw_vfacet(const Facet& F) {
            glNormal(-F.normal()) ;
            glVertex(F.vertex[0]) ;
            glVertex(F.vertex[1]) ;
            glVertex(F.vertex[2]) ;	
	}
	void draw_facet_mesh(const Facet& F) {
            if(check_voronoi_edge_flag(F,2)) { glVertex(F.vertex[0]) ; glVertex(F.vertex[1]) ; }
            if(check_voronoi_edge_flag(F,0)) { glVertex(F.vertex[1]) ; glVertex(F.vertex[2]) ; }
            if(check_voronoi_edge_flag(F,1)) { glVertex(F.vertex[2]) ; glVertex(F.vertex[0]) ; }
	}
	
	bool check_voronoi_edge_flag(const Facet& F, int e) {
            switch(F.edge_flag[e]) {
            case 0:
				return false ;             // aux edge
                //glColor3f(0.5f, 0.5f, 0.5f) ;
                //return true ;
            case 1:
                glColor3f(0.5f, 0.5f, 0.5f) ; // orignial edge
                return false ;
            case 2:
                glColor3f(0.1f, 0.1f, 0.1f) ; // clipped edge
                return true ;
            case 3:
                glColor3f(0.1f, 0.1f, 0.9f) ; // feature edge
                return true ;
            }
            return false ;
        }

        void draw_voronoi_trimesh_mesh(TriMesh& m) {
            glBegin(GL_LINES) ;
            for(unsigned int i=0; i<m.size(); i++) {
                const Facet& F = m[i] ;
                if(check_voronoi_edge_flag(F,2)) { glVertex(F.vertex[0]) ; glVertex(F.vertex[1]) ; }
                if(check_voronoi_edge_flag(F,0)) { glVertex(F.vertex[1]) ; glVertex(F.vertex[2]) ; }
                if(check_voronoi_edge_flag(F,1)) { glVertex(F.vertex[2]) ; glVertex(F.vertex[0]) ; }
            }
            glEnd() ;
        }
	
	bool check_aux_edge_flag(const Facet& F, int e) {
            switch(F.edge_flag[e]) {
	    case 0: // aux edge
		glColor3f(0.2f, 0.2f, 0.2f) ;
		return true ;
            }
            return false ;
        }
        void draw_voronoi_split_mesh(TriMesh& m) {
            glBegin(GL_LINES) ;
            for(unsigned int i=0; i<m.size(); i++) {
                const Facet& F = m[i] ;
                if(check_aux_edge_flag(F,2)) { glVertex(F.vertex[0]) ; glVertex(F.vertex[1]) ; }
                if(check_aux_edge_flag(F,0)) { glVertex(F.vertex[1]) ; glVertex(F.vertex[2]) ; }
                if(check_aux_edge_flag(F,1)) { glVertex(F.vertex[2]) ; glVertex(F.vertex[0]) ; }
            }
            glEnd() ;
        }
        
		void draw_rvd() const ;
		void draw_rvd_poly() const ;

    private:
        MyDelaunay* delaunay_ ;
        GLboolean show_domain_mesh_ ;
        GLboolean show_domain_ ;
        GLfloat vertices_size_ ;
        GLboolean show_mesh_ ;
        GLboolean colorize_ ;
        GLboolean show_primal_ ;
        GLboolean show_inner_cells_ ;
        GLboolean show_clipped_cells_ ;
        GLboolean clipped_boundary_only_ ;
        float slice_x_ ;
        float slice_y_ ;
        float slice_z_ ;
        float selection_ ;
        bool non_convex_ ;
		GLboolean show_rvd_ ;

        //---------------------------------------------------------------------------------------
        GLboolean show_surface_clipping_ ;
		GLboolean show_volume_clipping_ ;
		GLboolean show_remesh_ ;
		GLboolean show_volume_ ;
		GLboolean show_retet_ ;
		GLboolean show_cull_back_;
    } ;

} 

#endif
