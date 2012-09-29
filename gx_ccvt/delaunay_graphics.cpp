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

#include "delaunay_graphics.h"
#include <glut_viewer/glut_viewer.h>
#include <Geex/basics/stopwatch.h>
#include <Geex/graphics/geexob.h>
#include "rvd_ccvt.h"
#include "rvc_ccvt.h"

namespace Geex {

//	static const GLfloat surf_color[3] = {1.0, 0.835, 0.72} ;
	static const GLfloat surf_color[3] = {.9f, 0.9f, 0.6f} ;
    static const double c1 = 0.35 ;
    static const double c2 = 0.5 ;
    static const double c3 = 1.0 ;

    static double color_table[12][3] = 
    {
        {c3, c2, c2},
        {c2, c3, c2},
        {c2, c2, c3},
        {c2, c3, c3},
        {c3, c2, c3},
        {c3, c3, c2},

        {c1, c2, c2},
        {c2, c1, c2},
        {c2, c2, c1},
        {c2, c1, c1},
        {c1, c2, c1},
        {c1, c1, c2}

    } ;

    static int random_color_index_ = 0 ;


    static void gl_random_color() {
		glColor3dv(color_table[random_color_index_]);
        random_color_index_ = (random_color_index_ + 1) % 12 ;
    }
    
    static void gl_randomize_colors(int index = 0) {
        random_color_index_ = (index % 12) ;
    }

    
    DelaunayGraphics::DelaunayGraphics(MyDelaunay* delaunay) : delaunay_(delaunay) {
        show_domain_ = GL_FALSE ;
        show_domain_mesh_ = GL_FALSE ;
        vertices_size_ = 0.0 ;
        show_mesh_ = GL_TRUE ;
        colorize_ = GL_TRUE ;
        show_primal_ = GL_FALSE ;
        show_inner_cells_ = GL_TRUE ;
        show_clipped_cells_ = GL_FALSE ;
        clipped_boundary_only_ = GL_TRUE ;
		show_surface_clipping_ = GL_TRUE ;
		show_volume_clipping_  = GL_TRUE ;
		show_remesh_ = GL_FALSE ;
		show_volume_ = GL_FALSE ;
		show_retet_  = GL_FALSE ;
        slice_x_ = delaunay_->boundary_.xmin() ; //0.0 ;
        slice_y_ = delaunay_->boundary_.ymin() ; //0.0 ;
		slice_z_ = delaunay_->boundary_.zmax() ; //1.0 ;
        selection_ = 0.0 ;
		show_rvd_ = GL_FALSE ;
		show_cull_back_ = GL_TRUE ;
    }
    
    void DelaunayGraphics::draw() {
        non_convex_ = delaunay_->non_convex_mode_ ;

        if (Geexob::program_objects_supported()) {
            glUseProgramObjectARB(0) ;
        }

        if(show_domain_) {
            glColor3f(1.0, 1.0, 1.0) ;
            glCullFace(GL_FRONT) ;
            draw_domain() ;
			if (show_cull_back_)
				glCullFace(GL_BACK) ;
			else 
				glCullFace(GL_FRONT) ;
        }

        if(show_domain_mesh_) {
//	    glEnable(GL_BLEND) ;
//	    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA) ;
            draw_domain_mesh() ;
//	    glDisable(GL_BLEND) ;
        }
        if(vertices_size_ != 0.0) {
	    draw_seeds() ;
//            draw_vertices() ;
        }
        if(selection_ != 0.0) {
            draw_selection() ;
        }
        if(show_inner_cells_ || show_clipped_cells_) {
            draw_cells() ;
        }

        if(show_domain_) {
			if (show_cull_back_)
				glCullFace(GL_BACK) ;
			else 
				glCullFace(GL_FRONT) ;
//            glColor4f(1.0, 1.0, 1.0, 0.2) ;
 //           glEnable(GL_BLEND) ;
            draw_domain() ;
//            glDisable(GL_BLEND) ;
        }
		if(show_rvd_) {
            glDisable(GL_CULL_FACE) ;
//            draw_rvd() ;
            glEnable(GL_CULL_FACE) ;		
		}
	if(show_surface_clipping_ && !show_primal_) {
	    draw_surface_clipping() ;
	}
	//if(show_volume_clipping_ && show_clipped_cells_) 
	//	draw_volume_clipping() ;
	
	if(show_remesh_) {
        glEnable(GL_CULL_FACE) ;
	    draw_remesh() ;
	}

	if(show_volume_)
		draw_feature() ; 
////		draw_volume(*delaunay_->sizingfield_->ctlmesh_) ;
//	if(show_retet_)
//		draw_volume(delaunay_->retet_);
    }

    void DelaunayGraphics::draw_domain_mesh() {
        glDisable(GL_LIGHTING) ;
		glDepthRange(0, 0.999) ;
        glLineWidth(1.0) ;
        if(glut_viewer_is_enabled(GLUT_VIEWER_BACKGROUND) /* || colorize_*/) {
            glColor3f(1.f,1.f,0.f) ;
        } else {
            glColor3f(0.8f,0.8f,0.8f) ;      
        }
//        draw_trimesh_mesh(delaunay_->boundary_save_) ;
        draw_trimesh_mesh(delaunay_->boundary_) ;
		glDepthRange(0, 1) ;
    }
    
    void DelaunayGraphics::draw_domain() {
        glEnable(GL_LIGHTING) ;
		if (show_cull_back_)
			glCullFace(GL_BACK) ;
		else 
			glCullFace(GL_FRONT) ;
//	    glColor3f(.9, .9, .6) ;
		glColor3fv(surf_color) ;
        draw_trimesh(delaunay_->boundary_save_) ;
		//glDisable(GL_LIGHTING) ;
		//draw_colormap() ;
    }

    void DelaunayGraphics::draw_remesh() {
		Geex::TriMesh& m = delaunay_->rvd()->dual_mesh() ;
	glEnable(GL_LIGHTING) ;
    //glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0f);
	if (show_cull_back_)
		glCullFace(GL_BACK) ;
	else 
		glCullFace(GL_FRONT) ;
	//glColor3f(.9, .9, .6) ;
	glColor3fv(surf_color) ;
    glBegin(GL_TRIANGLES) ;
    for(unsigned int i=0; i<m.size(); i++) {
		const vec3 p = m[i].center() ;
//		if((p[0] <= slice_x_ || p[1] <= slice_y_ || p[2] <= slice_z_)) {
			glNormal(-m[i].normal()) ;
			glVertex(m[i].vertex[0]) ;
			glVertex(m[i].vertex[1]) ;
			glVertex(m[i].vertex[2]) ;
//		}
    }
    glEnd() ;
	//draw_trimesh(mesh) ;
	
	if(show_mesh_) {
	    glDisable(GL_LIGHTING) ;
		glDepthRange(0, 0.999) ;
	    glLineWidth(2) ;
//	    draw_trimesh_mesh(m) ;
        glBegin(GL_LINES) ;
        for(unsigned int i=0; i<m.size(); i++) {
            const Facet& F = m[i] ;
			const vec3 p = m[i].center() ;
//			if((p[0] <= slice_x_ || p[1] <= slice_y_ || p[2] <= slice_z_)) {
				if(check_edge_flag(F,2)) { glVertex(F.vertex[0]) ; glVertex(F.vertex[1]) ; }
				if(check_edge_flag(F,0)) { glVertex(F.vertex[1]) ; glVertex(F.vertex[2]) ; }
				if(check_edge_flag(F,1)) { glVertex(F.vertex[2]) ; glVertex(F.vertex[0]) ; }
//			}
		}
        glEnd() ;
//	    draw_trimesh_mesh(delaunay_->remesh_) ;
		glDepthRange(0, 1) ;
	}
	glPolygonMode(GL_FRONT_AND_BACK,GL_FILL) ;
    }

    void DelaunayGraphics::draw_vertices() {
        glDisable(GL_LIGHTING) ;
        glPointSize(int(vertices_size_ * 20)) ;
        glColor3f(1.0f,1.0f,0.3f) ;
        int w,h ;
        glut_viewer_get_screen_size(&w, &h) ;
        notify_screen_resize(w,h) ;
        begin_spheres() ;
        for(unsigned int i=0; i<delaunay_->all_vertices_.size(); i++) {
            MyDelaunay::Vertex_handle v = delaunay_->all_vertices_[i] ;
            glVertex3d(v->point().x(), v->point().y(), v->point().z()) ;
        }            
        end_spheres() ;
    }

    void DelaunayGraphics::draw_selection() {
        int iv = int(selection_ * (delaunay_->all_vertices_.size() - 1)) ;
        MyDelaunay::Vertex_handle v = delaunay_->all_vertices_[iv] ;
        glDisable(GL_LIGHTING) ;
        glPointSize(gx_max(int(vertices_size_ * 40), 4)*1.0f) ;
        glColor3f(1.0f,0.0f,0.3f) ;
        int w,h ;
        glut_viewer_get_screen_size(&w, &h) ;
        notify_screen_resize(w,h) ;
        begin_spheres() ;
        glVertex3d(v->point().x(), v->point().y(), v->point().z()) ;
        end_spheres() ;
        glLineWidth(2) ;
        glBegin(GL_LINES) ;
        glVertex3d(v->point().x(), v->point().y(), v->point().z()) ;
        glVertex3d(1e3, 1e3, 1e3) ;
        glEnd() ;

        //PluckerSegment E(to_geex(v->point()), to_geex(v->point()) + vec3(0, 0, 1e3)) ;
        //glDisable(GL_CULL_FACE) ;
        //glBegin(GL_TRIANGLES) ;
        //int isect_count = 0 ;
        //for(unsigned int i=0; i<delaunay_->boundary_plucker_.size(); i++) {
        //    const PluckerFacet& PF = delaunay_->boundary_plucker_[i] ;
        //    const Facet& F = delaunay_->boundary_[i]; 
        //    if(PF.intersects(E)) {
        //        isect_count++ ;
        //        glVertex(F.vertex[0]) ;
        //        glVertex(F.vertex[1]) ;
        //        glVertex(F.vertex[2]) ;
        //    }
        //}
        //glEnd() ;
        //glEnable(GL_CULL_FACE) ;
        //std::cerr << "nb isect = " << isect_count << std::endl ;
    }

    void DelaunayGraphics::draw_cells() {
        glEnable(GL_COLOR_MATERIAL) ;
        gl_randomize_colors() ;
        glEnable(GL_LIGHTING) ;
        glColor3f(1,1,1) ;

        if(show_primal_) {
            for(MyDelaunay::Cell_iterator it = delaunay_->cells_begin(); it != delaunay_->cells_end(); it++) {
                if(delaunay_->is_infinite(it)/* || it->aux*/) { continue ; }
                if(colorize_) { gl_random_color() ; }
                draw_primal_cell(it) ;
            }
            if(show_mesh_) {
                glLineWidth(2) ;
                glDisable(GL_LIGHTING) ;
                glColor3f(0,0,0) ;
                glPolygonMode(GL_FRONT_AND_BACK,GL_LINE) ;
                for(MyDelaunay::Cell_iterator it = delaunay_->cells_begin(); it != delaunay_->cells_end(); it++) {
                    if(delaunay_->is_infinite(it)/* || it->aux*/) { continue ; }
                    draw_primal_cell(it, true) ;
                }
            }
            glPolygonMode(GL_FRONT_AND_BACK,GL_FILL) ;
        } else {
            //for(unsigned int i=0; i<delaunay_->all_vertices_.size(); i++) {
            //    if(colorize_) { gl_random_color() ; }
            //    draw_dual_cell(delaunay_->all_vertices_[i]) ;
            //}
            //if(show_mesh_) {
            //    glLineWidth(2) ;
            //    glDisable(GL_LIGHTING) ;
            //    glColor3f(0,0,0) ;
            //    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE) ;
            //    for(unsigned int i=0; i<delaunay_->all_vertices_.size(); i++) {
            //        draw_dual_cell(delaunay_->all_vertices_[i], true) ;
            //    }
            //}
            //glPolygonMode(GL_FRONT_AND_BACK,GL_FILL) ;
        }
    }

    void DelaunayGraphics::draw_dual_cell(MyDelaunay::Vertex_handle v, bool mesh) {
    //    const Point& p = v->point() ;        
    //    if( /* (non_convex_ && slice_z_ > 0.99) || */ (p.x() <= slice_x_ || p.y() <= slice_y_ || p.z() <= slice_z_)) {
    //        std::vector<MyDelaunay::Edge> edges ;
    //        MyDelaunay::CellStatus stat = delaunay_->get_incident_edges(v, edges) ;
    //        if(stat == MyDelaunay::CS_INTERIOR && show_inner_cells_) {
    //            //for(unsigned int i=0; i<edges.size(); i++) {
    //            //    draw_dual_face(edges[i]) ;
    //            //}
    //        } else if(show_clipped_cells_ && stat != MyDelaunay::CS_INTERIOR) {
    //            //TriMesh* m = delaunay_->dual_convex_clip(edges, !clipped_boundary_only_) ;
    //            //if(mesh) {
    //            //    draw_trimesh_mesh(*m) ;
    //            //} else {
    //            //    draw_trimesh(*m) ;
    //            //}
				//// draw_fast_clipping() ;
    //        }
    //    }
    }


/*
                3       3
                    0
                1       2

                    3
*/

    static inline void draw_facet(const vec3& p1, const vec3& p2, const vec3& p3, bool mesh) {
        if(!mesh) {
            glNormal(cross(p2-p1, p3-p1)) ;
        }
        glVertex(p1) ; glVertex(p2) ; glVertex(p3) ;
    }

    void DelaunayGraphics::draw_primal_cell(MyDelaunay::Cell_handle c, bool mesh) {
        vec3 p0 = to_geex(c->vertex(0)->point()) ;
        vec3 p1 = to_geex(c->vertex(1)->point()) ;
        vec3 p2 = to_geex(c->vertex(2)->point()) ;
        vec3 p3 = to_geex(c->vertex(3)->point()) ;

        vec3 g = 0.25 * (p0 + p1 + p2 + p3) ;
        if( /* (non_convex_ && slice_z_ > 0.99) || */ (g.x <= slice_x_ || g.y <= slice_y_ || g.z <= slice_z_)) {

            glBegin(GL_TRIANGLES) ;
            
            draw_facet(p2, p1, p0, mesh) ;
            draw_facet(p0, p1, p3, mesh) ;
            draw_facet(p3, p2, p0, mesh) ;
            draw_facet(p2, p3, p1, mesh) ;
            
            glEnd() ;
        }
    }

	void DelaunayGraphics::draw_primal_cell(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3, bool mesh)  {
//        vec3 g = 0.25 * (p0 + p1 + p2 + p3) ;

//		if( /* (non_convex_ && slice_z_ > 0.99) || */ (g.x <= slice_x_ || g.y <= slice_y_ || g.z <= slice_z_)) {

            glBegin(GL_TRIANGLES) ;
            
            draw_facet(p2, p1, p0, mesh) ;
            draw_facet(p0, p1, p3, mesh) ;
            draw_facet(p3, p2, p0, mesh) ;
            draw_facet(p2, p3, p1, mesh) ;
            
            glEnd() ;
//        }
	}

	void DelaunayGraphics::draw_seeds() {
		glDisable(GL_LIGHTING) ;
		glColor3f(.8f, 0.f, 0.f) ;
		glPointSize(float(vertices_size_ * 20)) ;
		glDepthRange(0, 0.998) ;
		glEnable(GL_POINT_SMOOTH) ;
		glBegin(GL_POINTS) ;
		for(int i=0; i<delaunay_->all_vertices_.size(); ++i) {
			glVertex(to_geex(delaunay_->all_vertices_[i]->point())) ;
		}
		glEnd() ;
		glDisable(GL_POINT_SMOOTH) ;
		glPointSize(1.0) ;
		glDepthRange(0, 1) ;
	}

    void DelaunayGraphics::draw_surface_clipping() {
        glEnable(GL_COLOR_MATERIAL) ;
        gl_randomize_colors() ;
        glEnable(GL_LIGHTING) ;
		glEnable(GL_CULL_FACE) ;
		if (show_cull_back_)
			glCullFace(GL_BACK) ;
		else 
			glCullFace(GL_FRONT) ;
		
		int nv = (int)delaunay_->vertices().size() ;
		std::vector<RestrictedVoronoiCell>& rvcells = delaunay_->rvd()->rvcells() ;
	
		// draw Voronoi facets
		glBegin(GL_TRIANGLES) ;
		for(int i=0; i<nv; ++i) {
			if(colorize_) { gl_random_color() ; }
			else glColor3fv(surf_color) ; 

			const Point& p = rvcells[i].vh()->point() ; 
//			if((p[0] <= slice_x_ || p[1] <= slice_y_ || p[2] <= slice_z_)) {
				if(rvcells[i].is_boundary()) {
					TriMesh& m = rvcells[i].mesh() ;
					for(int j=0; j<m.size(); ++j) 
						draw_vfacet(m[j]) ;
				}
//			}
        }
		glEnd() ;

		glDisable(GL_CULL_FACE) ;

		// draw Voronoi edges
	    if(show_mesh_) {
			glLineWidth(2) ;
			glDisable(GL_LIGHTING) ;
			glDepthRange(0, 1.0) ;

			glBegin(GL_LINES) ;
			for(int i=0; i<nv; ++i) {
				const Point& p = rvcells[i].vh()->point() ; 
//				if((p[0] <= slice_x_ || p[1] <= slice_y_ || p[2] <= slice_z_)) {
					if(rvcells[i].is_boundary()) {
						TriMesh& m = rvcells[i].mesh() ;
						for(int j=0; j<m.size(); ++j) 
							draw_facet_mesh(m[j]) ;
					}
//				}
			}
			glEnd() ;
			glLineWidth(1) ;
			glEnable(GL_LIGHTING) ;
			glDepthRange(0, 1) ;
			glPolygonMode(GL_FRONT_AND_BACK,GL_FILL) ;
		}
    }

	void DelaunayGraphics::draw_feature() {
		std::vector<LineSeg>& edges = delaunay_->boundary_.features_ ;
		glLineWidth(4) ;
		glDisable(GL_LIGHTING) ;
		glColor3f(0., 0., 0.) ;
		glDepthRange(0, 0.998) ;

		glBegin(GL_LINES) ;
		for(unsigned int i=0; i<edges.size(); ++i) {
			glVertex3dv(&edges[i].v0_[0]);
			glVertex3dv(&edges[i].v1_[0]);
			//glVertex3f(edges[i].v0_[0], edges[i].v0_[1], edges[i].v0_[2]) ;
			//glVertex3f(edges[i].v1_[0], edges[i].v1_[1], edges[i].v1_[2]) ;
		}
		glEnd() ;
		glDepthRange(0, 1) ;
		glEnable(GL_LIGHTING) ;
	    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL) ;
	}

	void DelaunayGraphics::draw_rvd() const {
	}

	void DelaunayGraphics::draw_rvd_poly() const {
	}

}
