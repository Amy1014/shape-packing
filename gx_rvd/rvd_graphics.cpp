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

#include "rvd_graphics.h"
#include "rvd_geometry.h"
#include "rvd_primal.h"
#include "glut_viewer/glut_viewer.h"
#include <Geex/graphics/opengl.h>
#include <Geex/graphics/geexob.h>

namespace Geex {

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


  static float white[4] = {
      0.8f, 0.8f, 0.3f, 1.0f
  } ;

  void RVDGraphics::gl_random_color() {
      glColor3f(
          color_table[random_color_index_][0], 
          color_table[random_color_index_][1], 
          color_table[random_color_index_][2]
      ) ;
      random_color_index_ = (random_color_index_ + 1) % 12 ;
  }
  
  void RVDGraphics::gl_randomize_colors(int index) {
      random_color_index_ = (index % 12) ;
  }

  void RVDGraphics::gl_random_color(int index) {
      if(index >= 0) {
          gl_randomize_colors(index) ;
          gl_random_color() ;
      } else {
          glColor4fv(white) ;
      }
  }

  RVDGraphics::RVDGraphics(RVDGeometry* geometry): geometry_(geometry) {
      surface_mode_ = Plain ;
      show_domain_mesh_ = GL_FALSE ;
      vertices_size_ = 0.2f ;
      show_rvd_mesh_ = GL_FALSE ;
      show_rvd_submesh_ = GL_FALSE ;
      show_primal_ = GL_FALSE ;
      show_primal_mesh_ = GL_FALSE ;
      shrink_ = 0.0f ;
      id_threshold_ = 1.0f ;
      show_non_manifold_ = GL_FALSE ;
      colormap_texture_ = 0 ;
      shiny_ = GL_TRUE ;
	  show_vertices_ = GL_TRUE;
      u_primal_ = u_primal_mesh_ = u_vertices_ = u_singular_ = u_singular_mesh_ = 0;
      u_rvd_none_ = u_rvd_cell_ = u_rvd_energy_ = u_rvd_linear_ = u_rvd_connected_ = u_rvd_shrink_ = u_rvd_mesh_ = u_rvd_shrink_mesh_ = 0; 
      use_display_lists_ = GL_FALSE ;
  }

  void RVDGraphics::delete_lists() {
      delete_list(u_primal_) ;
      delete_list(u_primal_mesh_);
      delete_list(u_vertices_);
      delete_list(u_singular_) ;
      delete_list(u_singular_mesh_) ;
      delete_list(u_rvd_none_);
      delete_list(u_rvd_cell_) ;
      delete_list(u_rvd_connected_) ;
      delete_list(u_rvd_energy_) ;
      delete_list(u_rvd_shrink_) ;
      delete_list(u_rvd_mesh_) ;
      delete_list(u_rvd_shrink_mesh_) ;
  }

  RVDGraphics::~RVDGraphics() {
      delete_lists() ;
  }
  
  inline void draw_triangle(const vec3& v1, const vec3& v2, const vec3& v3) {
      glVertex(v1) ;
      glVertex(v2) ;
      glVertex(v3) ;
  }

  inline void draw_triangle_with_normal(
      const vec3& v1, const vec3& v2, const vec3& v3
  ) {
      glNormal(cross(v2-v1, v3-v1)) ;
      glVertex(v1) ;
      glVertex(v2) ;
      glVertex(v3) ;
  }

  class draw_triangle_with_normal_and_colors {
  public:
      draw_triangle_with_normal_and_colors(
          double wmin, double wmax
      ) : bias_(wmin), scale_(1.0/(wmax - wmin)) {
      }
      void operator() (
          const TopoPolyVertexEdge& v1, 
          const TopoPolyVertexEdge& v2, 
          const TopoPolyVertexEdge& v3
      ) const {
          glNormal(cross(v2-v1, v3-v1)) ;
          send_vertex(v1) ;
          send_vertex(v2) ;
          send_vertex(v3) ;
      }
  protected:
      void send_vertex(const TopoPolyVertexEdge& v) const {
          double w = 1.0 - scale_ * (1.0 / v.w - bias_) ;
          glTexCoord1d(w) ;
          glVertex(v) ;
      }
  private:
      double bias_ ;
      double scale_ ;
  } ;


  static float spec[4] = {
      0.7f, 0.7f, 0.7f, 1.0f
  } ;

  static float shininess = 20.0f ;

  void RVDGraphics::draw(bool geometry_is_dirty) {
      m_geometry_dirty_ = geometry_is_dirty;
      show_primal_mesh_ = show_primal_ ;
      double x_min, y_min, z_min ;
      double x_max, y_max, z_max ;
      geometry()->get_bbox(x_min, y_min, z_min, x_max, y_max, z_max) ;
      if(!use_display_lists_ || m_geometry_dirty_) {
          delete_lists() ;
      }

      if(Geexob::program_objects_supported()) {
          glUseProgramObjectARB(0) ;
      }

      if(shiny_) {
          white[0] = 0.8f ;
          white[1] = 0.8f ;
          white[2] = 0.3f ;
          white[3] = 1.0f ;
          spec[0] = 0.7f ;
          spec[1] = 0.7f ;
          spec[2] = 0.7f ;
          spec[3] = 1.0f ;
          shininess = 20.0f ;
      } else {
          white[0] = 1.0f ;
          white[1] = 1.0f ;
          white[2] = 1.0f ;
          white[3] = 1.0f ;
          spec[0] = 0.0f ;
          spec[1] = 0.0f ;
          spec[2] = 0.0f ;
          spec[3] = 1.0f ;
          shininess = 0.0f ;
      }

      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,   spec) ;
      glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, &shininess) ;

      if(colormap_texture_ == 0) {
          create_colormap_texture() ;
      }

      glEnable(GL_COLOR_MATERIAL) ;

      bool show_dom = (
          surface_mode_ == Plain || surface_mode_ == Density
      ) && !show_primal_ && (shrink_ == 0.0) ;

      if(show_dom) {
          glColor4fv(white) ;
          glCullFace(GL_FRONT) ;
          draw_domain() ;
          glCullFace(GL_BACK) ;
      }

      if(show_domain_mesh_ && shrink_ == 0.0) {
          draw_domain_mesh() ;
      }

      if(
          !show_primal_ &&
          (surface_mode_ == Cells || surface_mode_ == Energy 
	    || surface_mode_ == Linear || surface_mode_ == Connected 
	    ||(shrink_ != 0.0)
          )
      ) {
          draw_rvd() ;
      }

      if(show_rvd_mesh_ || (show_domain_mesh_ && shrink_ != 0.0)) {
          draw_rvd_mesh() ;
      }
      
      if(show_primal_) {
          draw_primal() ;
      }

      if(show_primal_mesh_) {
          draw_primal_mesh() ;
      }

      if(vertices_size_ > 0.0) {
          draw_vertices() ;
      }

      if(show_dom) {
          glCullFace(GL_BACK) ;
          glColor4f(1.0f, 1.0f, 1.0f, 0.2f) ;
          glEnable(GL_BLEND) ;
          draw_domain() ;
          glDisable(GL_BLEND) ;
      }
  }

  void RVDGraphics::draw_domain_mesh() {
      glDisable(GL_LIGHTING) ;
      glLineWidth(1) ;
      glColor3f(0.0, 0.5, 1.0) ;
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) ;
      glBegin(GL_TRIANGLES) ;
      geometry_->boundary()->for_each_triangle(draw_triangle) ;
      glEnd() ;
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;
  }

  void RVDGraphics::draw_domain() {
      if(surface_mode_ == Density) {
          glDisable(GL_LIGHTING) ;
          glEnable(GL_TEXTURE_1D) ;
          glBindTexture(GL_TEXTURE_1D, colormap_texture_) ;
      } else {
          glEnable(GL_LIGHTING) ;
      }
      glBegin(GL_TRIANGLES) ;
      if(surface_mode_ == Density) {
          double wmin = 1.0 / geometry_->boundary()->vertex(0).w ;
          double wmax = wmin ;
          for(unsigned int i=0; i<geometry_->boundary()->nb_vertices(); i++) {
              wmin = gx_min(wmin, 1.0 / geometry_->boundary()->vertex(i).w) ;
              wmax = gx_max(wmax, 1.0 / geometry_->boundary()->vertex(i).w) ;
          }
          glColor3f(1.0, 1.0, 1.0) ;
          geometry_->boundary()->for_each_triangle(
              draw_triangle_with_normal_and_colors(wmin, wmax)
          ) ;
      } else {
          glColor4fv(white) ;
          geometry_->boundary()->for_each_triangle(draw_triangle_with_normal) ;
      }
      glEnd() ;
      if(surface_mode_ == Density) {
          glDisable(GL_TEXTURE_1D) ;
      }
  }

  void RVDGraphics::draw_vertices() {
	  if (!show_vertices()) return;

      glDisable(GL_LIGHTING) ;
      glDepthRange(0.0, 0.998) ;
      glEnable(GL_POINT_SMOOTH) ;
      glPointSize(vertices_size_ * 20.0f) ;
      glLineWidth(vertices_size_ * 10.0f) ;
      glColor3f(0.0, 0.0, 1.0) ;	  
      if(!playback_list(u_vertices_)) {
          begin_list(u_vertices_) ;
          glBegin(GL_POINTS) ;
          for(unsigned int i=0; i<geometry_->vertices().size(); i++) {
              glVertex(geometry_->vertices()[i]);
          }
          glEnd() ;
          //glBegin(GL_LINES) ;
	  //if(geometry_->facets().size() != 0) {
	  //  for(unsigned int i=0; i<geometry_->delaunay()->nb_vertices(); i++) {
          //      const vec3& vertex = geometry_->vertices()[i] ;
	  //      unsigned int facet = geometry_->facets()[i] ;
	  //      vec3 normal = normalize(geometry_->facet_point_normal(facet, vertex)) ;
	  //      glVertex(vertex) ;
	  //      glVertex(vertex + normal) ;
	  //}
          //}
          //glEnd() ;
          end_list() ;
      }
      glDepthRange(0.0, 1.0) ;
  }


  inline void c_draw_triangle(
      unsigned int c, const vec3& v1, const vec3& v2, const vec3& v3
  ) {
      glVertex(v1) ;
      glVertex(v2) ;
      glVertex(v3) ;
  } 

  inline void c_draw_triangle_with_normal(
      unsigned int c, const vec3& v1, const vec3& v2, const vec3& v3
  ) {
      glNormal(cross(v2-v1, v3-v1)) ;
      glVertex(v1) ;
      glVertex(v2) ;
      glVertex(v3) ;
  } 

  inline void c_draw_triangle_with_normal_and_color(
      unsigned int c, const vec3& v1, const vec3& v2, const vec3& v3
  ) {
      RVDGraphics::gl_randomize_colors(c) ;
      RVDGraphics::gl_random_color() ;
      glNormal(cross(v2 - v1, v3 - v1)) ;
      glVertex(v1) ;
      glVertex(v2) ;
      glVertex(v3) ;
  }

  class DrawThresholdTriangle{
  public:
    DrawThresholdTriangle(
      RVDGraphics& graphics
    ): threshold_(graphics.id_threshold()),
       geometry_(graphics.geometry()){
    }

    inline void operator()(
          unsigned int c, const vec3& v1, const vec3& v2, const vec3& v3
      ) const {
      if(
          (geometry_->size_of_rvd() == 0) ||
          (threshold_*(geometry_->size_of_rvd()) > geometry_->triangle_id())
	)
	c_draw_triangle_with_normal_and_color(c,v1,v2,v3) ;
    }

    double threshold_ ;
    RVDGeometry* geometry_ ;
  } ;

  class DrawShrinkTriangle {
  public:
      DrawShrinkTriangle(
          RVDGraphics& graphics
      ) : graphics_(graphics),
          shrink_(graphics.shrink()), vertices_(graphics.geometry()->vertices()), 
          colorize_(graphics.surface_mode() == Cells 
		    || graphics.surface_mode() == Energy) {
      }
      inline void operator()(
          unsigned int c, const vec3& v1, const vec3& v2, const vec3& v3
      ) const {
        if(colorize_) {
            RVDGraphics::gl_randomize_colors(c) ;
            RVDGraphics::gl_random_color() ;
        }
        glNormal(cross(v2 - v1, v3 - v1)) ;
        const vec3& p = vertices_[c] ;
        glVertex(shrink_ * p + (1.0 - shrink_) * v1) ;
        glVertex(shrink_ * p + (1.0 - shrink_) * v2) ;
        glVertex(shrink_ * p + (1.0 - shrink_) * v3) ;
      }
      RVDGraphics& graphics_ ;
      double shrink_ ;
      const std::vector<vec3>& vertices_ ;
      bool colorize_ ;
  } ;
  
  class DrawTriangleWithEnergy {
  public:
      DrawTriangleWithEnergy(
          const std::vector<double>& E_in
      ) : E(E_in), c1_(0.0, 0.0, 0.0), c2_(1.0, 1.0, 1.0) {
          double emin = 0. ;
          double emax = -1e30 ;
          for(unsigned int i=0; i<E.size(); i++) {
              emin = gx_min(emin, E[i]) ;
              emax = gx_max(emax, E[i]) ;
          }
          bias_ = emin ;
          if(emax - emin > 1e-30) {
              scale_ = 1.0 / (emax - emin) ;
          } else {
              scale_ = 1.0 ;
          }
      }

      bool is_boundary_edge(
        const SymbolicVertex s1,
        const SymbolicVertex s2
      ) const {
        SymbolicVertex I ;
        sets_intersect(s1,s2,I) ;
        return (I.size()>= 2 && I[0] < 0 && I[1] > 0) ;
      }

      void draw_edge(const TopoPolyVertexEdge& v1, const TopoPolyVertexEdge& v2) const {
	if(is_boundary_edge(v1.sym,v2.sym)){
	  glEnable(GL_LINE_SMOOTH) ;
	  glLineWidth(1) ;
	  glColor3f(0.,0.,0.) ;
	  glBegin(GL_LINES) ;
	  glVertex(v1) ;
	  glVertex(v2) ;
	  glEnd() ;
	}
      }

      void operator()(
          unsigned int c, const TopoPolyVertexEdge& v1, const TopoPolyVertexEdge& v2, const TopoPolyVertexEdge& v3
      ) const {
          double e = (E[c] - bias_)*scale_ ;
          //double e = gx_min(E[c],1.) ;
	  //glBegin(GL_TRIANGLES) ;
          glColor(mix(c1_, c2_, e)) ;
	  //glNormal(cross(v2 - v1, v3 - v1)) ;
          glVertex(v1) ;
          glVertex(v2) ;
          glVertex(v3) ;
	  //glEnd() ;

	  //draw_edge(v1,v2) ;
	  //draw_edge(v2,v3) ;
	  //draw_edge(v3,v1) ;
	    
      }
  private:
      const std::vector<double>& E ;
      double bias_ ;
      double scale_ ;
      const vec3 c1_ ; 
      const vec3 c2_ ;
  } ;

  class DoNothing {
    public:
      void operator()(
          unsigned int c, const vec3& v1, const vec3& v2, const vec3& v3
      ) const { }
  } ;
  
  class DrawTriangleWithLinear {
  public:
      DrawTriangleWithLinear(
          const std::vector<vec4>& planes_in
      ) : planes(planes_in), c1_(0.0, 0.0, 0.0), c2_(1.0, 1.0, 1.0) {
      }
      void operator()(
          unsigned int c, const vec3& v1, const vec3& v2, const vec3& v3
      ) const {
	  vec4 plane = planes[c] ;
          double e = dot(plane,vec4(v1[0],v1[1],v1[2],1.)) ;
          glColor(mix(c1_, c2_, e)) ;
          glVertex(v1) ;
          e = dot(plane,vec4(v2[0],v2[1],v2[2],1.)) ;
          glColor(mix(c1_, c2_, e)) ;
          glVertex(v2) ;
          e = dot(plane,vec4(v3[0],v3[1],v3[2],1.)) ;
          glColor(mix(c1_, c2_, e)) ;
          glVertex(v3) ;
      }
  private:
      const std::vector<vec4>& planes ;
      const vec3 c1_ ; 
      const vec3 c2_ ;
  } ;

  void RVDGraphics::draw_rvd() {

      if(surface_mode_ == None) { return ; }

      glColor4fv(white) ;

      if(shrink_ != 0.0) {
          glEnable(GL_LIGHTING) ;
          if(!playback_list(u_rvd_shrink_)) {
              begin_list(u_rvd_shrink_) ;
              glBegin(GL_TRIANGLES) ;
              geometry_->RVD().for_each_triangle(DrawShrinkTriangle(*this)) ;
              glEnd() ;
              end_list() ;
          }
      } else if(surface_mode_ == Energy) {
          glDisable(GL_LIGHTING) ;
          if(!playback_list(u_rvd_energy_)) {
              begin_list(u_rvd_energy_) ;
              glBegin(GL_TRIANGLES) ;
	      geometry_->RVD().set_symbolic(false) ;
              //geometry_->RVD().for_each_triangle(DrawTriangleWithEnergy(geometry()->weights())) ;
              geometry_->RVD().for_each_triangle(DrawTriangleWithEnergy(geometry()->energy())) ;
              glEnd() ;
              end_list() ;
          }
      } else if(surface_mode_ == Linear) {
          glDisable(GL_LIGHTING) ;
          if(!playback_list(u_rvd_linear_)) {
              begin_list(u_rvd_linear_) ;
              glBegin(GL_TRIANGLES) ;
              geometry_->RVD().set_symbolic(false) ;
              geometry_->RVD().for_each_triangle(DrawTriangleWithLinear(geometry()->planes())) ;
              glEnd() ;
              end_list() ;
          }
      } else if(surface_mode_ == Cells) {
          glEnable(GL_LIGHTING) ;
          if(!playback_list(u_rvd_cell_)) {
              begin_list(u_rvd_cell_) ;
              glBegin(GL_TRIANGLES) ;
              geometry_->RVD().for_each_triangle(DrawThresholdTriangle(*this)) ;
              glEnd() ;
              end_list() ;
          }
      //} else if(surface_mode_ == Connected) {
      //    glDisable(GL_LIGHTING) ;
      //    if(!playback_list(u_rvd_connected_)) {
      //        begin_list(u_rvd_connected_) ;
      //        glBegin(GL_TRIANGLES) ;
      //        geometry_->RVD().for_each_triangle(DoNothing()) ;
      //        geometry_->RVD().for_each_triangle_connected(DrawTriangleWithEnergy(geometry()->energy()),&(geometry()->facets())) ;
      //        glEnd() ;
      //        end_list() ;
      //    }
      } else {
          glEnable(GL_LIGHTING) ;
          if(!playback_list(u_rvd_none_)) {
              begin_list(u_rvd_none_) ;
              glBegin(GL_TRIANGLES) ;
              geometry_->RVD().for_each_triangle(c_draw_triangle_with_normal) ;
              glEnd() ;
              end_list() ;
          }
      }
  }

  void c_draw_edge(unsigned int c, const vec3& v1, const vec3& v2) {
      glVertex(v1) ;
      glVertex(v2) ;
  }

  class c_draw_shrink_edge {
  public:
      c_draw_shrink_edge(
          double shrink, const std::vector<vec3>& seeds
      ) : shrink_(shrink), seeds_(seeds) {
      }
      void operator()(unsigned int c, const vec3& v1, const vec3& v2) const {
          const vec3& p = seeds_[c] ;
          glVertex(shrink_ * p + (1.0 - shrink_) * v1) ;
          glVertex(shrink_ * p + (1.0 - shrink_) * v2) ;
      }
  private:
      double shrink_ ;
      const std::vector<vec3>& seeds_ ;
  } ;


  void draw_original_mesh(unsigned int v, TopoPolyMesh* M) {
      for(unsigned int f=0; f<M->nb_facets(); f++) {
          unsigned int facet_base = M->facet_begin(f) ;
          unsigned int facet_n = M->facet_size(f) ;
          for(unsigned int i=0; i<facet_n; i++) {
              if(M->vertex(i).flags == TopoPolyVertexEdge::ORIGINAL) {
                  glVertex(M->vertex(facet_base + i)) ;
                  glVertex(M->vertex(facet_base + ((i+1)%facet_n))) ;
              }
          }
      }
  } 

  class draw_shrink_original_mesh {
  public:
      draw_shrink_original_mesh(
          double shrink, const std::vector<vec3>& seeds
      ) : shrink_(shrink), seeds_(seeds) { 
      }
      void operator()(unsigned int v, TopoPolyMesh* M) const {
          const vec3& p = seeds_[v] ;
          for(unsigned int f=0; f<M->nb_facets(); f++) {
              unsigned int facet_base = M->facet_begin(f) ;
              unsigned int facet_n = M->facet_size(f) ;
              for(unsigned int i=0; i<facet_n; i++) {
                  if(M->vertex(i).flags == TopoPolyVertexEdge::ORIGINAL) {
                      glVertex(
                          shrink_ * p + (1.0 - shrink_) * M->vertex(facet_base + i)
                      ) ;
                      glVertex(
                          shrink_ * p + (1.0 - shrink_) * M->vertex(facet_base + ((i+1)%facet_n))
                      ) ;
                  }
              }
          }
      }
  private:
      double shrink_ ;
      const std::vector<vec3>& seeds_ ;
  } ;


  void RVDGraphics::draw_rvd_mesh() {
      glDisable(GL_LIGHTING) ;
      glLineWidth(1) ;

      if(shrink_ == 0.0) {
          if(!playback_list(u_rvd_mesh_)) {
              begin_list(u_rvd_mesh_) ;
              glBegin(GL_LINES) ;
              glColor3f(0.,0.,0.) ;
              geometry_->RVD().for_each_halfedge(c_draw_edge) ;
              glEnd() ;
              end_list() ;
          }
          if(show_rvd_submesh_) {
              glColor3f(0.0, 0.5, 1.0) ;
              glBegin(GL_LINES) ;
              geometry_->RVD().for_each_facet(draw_original_mesh) ;
              glEnd() ;
          }
      } else {
          if(!playback_list(u_rvd_shrink_mesh_)) {
              begin_list(u_rvd_shrink_mesh_) ;
              glBegin(GL_LINES) ;
              glColor3f(0.2f, 0.2f, 0.2f) ;
              geometry_->RVD().for_each_halfedge(c_draw_shrink_edge(shrink_, geometry()->vertices())) ;
              glEnd() ;
              end_list() ;
          }
          if(show_rvd_submesh_ || show_domain_mesh_) {
              glColor3f(0.0f, 0.5f, 1.0f) ;
              glBegin(GL_LINES) ;
              geometry_->RVD().for_each_facet(draw_shrink_original_mesh(shrink_, geometry()->vertices())) ;
              glEnd() ;
          }
      }
  }

  class DrawPrimal {
  public:
      DrawPrimal(const std::vector<vec3>& seed_in) : seed_(seed_in) { }
      void operator()(unsigned int i, unsigned int j, unsigned int k) const {
          const vec3& v1 = seed_[i] ;
          const vec3& v2 = seed_[j] ;
          const vec3& v3 = seed_[k] ;
          glNormal(cross(v2-v1, v3-v1)) ;
          glVertex(v1) ;
          glVertex(v2) ;
          glVertex(v3) ;
      }
  private:
      const std::vector<vec3>& seed_ ;
  } ;


  class DrawPrimalSingularTriangles {
  public:
      DrawPrimalSingularTriangles(
          const std::vector<vec3>& seed_in, std::map< trindex, int >& triangles_in
      ) : seed_(seed_in), triangles_(triangles_in) { }
      void operator()(unsigned int i, unsigned int j, unsigned int k) const {
          trindex t(i,j,k) ;
          if(triangles_[t] > 1) {
              glColor3f(1.0, 0.0, 0.0) ;
          } else {
              glColor4fv(white) ;
          }
          const vec3& v1 = seed_[i] ;
          const vec3& v2 = seed_[j] ;
          const vec3& v3 = seed_[k] ;
          glNormal(cross(v2-v1, v3-v1)) ;
          glVertex(v1) ;
          glVertex(v2) ;
          glVertex(v3) ;
      }
  private:
      const std::vector<vec3>& seed_ ;
      std::map< trindex, int >& triangles_ ;        
  } ;


  void RVDGraphics::draw_primal() {
      glColor4fv(white) ;
      glEnable(GL_LIGHTING) ;
      if(show_non_manifold_) {
          if(!playback_list(u_singular_)) {
              begin_list(u_singular_) ;
              std::map<trindex, int> triangles ;        
              geometry_->RVD().for_each_primal_triangle(CountPrimalSingularTriangles(triangles)) ;
              glBegin(GL_TRIANGLES) ;
              geometry_->RVD().for_each_primal_triangle(
                  DrawPrimalSingularTriangles(geometry_->vertices(), triangles)
              ) ;
              glEnd() ;
              end_list() ;
          }
      } else {
          if(!playback_list(u_primal_)) {
              begin_list(u_primal_) ;
              glBegin(GL_TRIANGLES) ;
              geometry_->RVD().for_each_primal_triangle(DrawPrimal(geometry_->vertices())) ;
              glEnd() ;
              end_list() ;
          }
      }
  }


  class DrawPrimalMesh {
  public:
      DrawPrimalMesh(RVDGeometry* geometry) : vertices_(geometry->vertices()) {
      }

      void operator()(unsigned int c, TopoPolyMesh* M) const {
          for(unsigned int i=0; i<M->nb_vertices(); i++) {
              draw(M->vertex(i)) ;
          }
      }        
  protected:
      void draw(const TopoPolyVertexEdge& v) const {
          // if v has two bisectors, then v.sym[0] is a boundary facet (<0)
          // and the two bisectors are v.sym[1] and v.sym[2]. 
          // Note: this draws the same edge several times, we'd probably
          // better use draw_primal() in wireframe mode ...
          // -> This is what we do now, this class is no longer used, kept
          // here for reference. 
          if(v.sym.nb_bisectors() >= 2) {
              unsigned int i1 = v.sym.bisector(0) ;
              unsigned int i2 = v.sym.bisector(1) ;
              glVertex(vertices_[i1]) ;
              glVertex(vertices_[i2]) ;
          }
      }
  private:
      const std::vector<vec3>& vertices_ ;
  } ;


  class DrawPrimalSingularEdges {
  public:
      DrawPrimalSingularEdges(
          const std::vector<vec3>& seeds,
          std::map<bindex, unsigned int>& edges
      ) : seeds_(seeds), edges_(edges) {
      }

      void operator()(unsigned int c, TopoPolyMesh* M) const {
          for(unsigned int i=0; i<M->nb_vertices(); i++) {
              if(M->vertex(i).sym.nb_bisectors() > 0) {
                  std::vector<unsigned int> v ;
                  v.push_back(c) ;
                  for(unsigned int j = 0; j<M->vertex(i).sym.nb_bisectors(); j++) {
                      v.push_back(M->vertex(i).sym.bisector(j)) ;
                  }
                  for(unsigned int j1=0; j1<v.size(); j1++) {
                      for(unsigned int j2=j1+1; j2<v.size(); j2++) {
                          unsigned int v1 = v[j1] ;
                          unsigned int v2 = v[j2] ;
                          bindex E(v1, v2) ;
                          std::map<bindex, unsigned int>::const_iterator it = edges_.find(E) ;
                          if(it == edges_.end() || it->second > 2) {
                              glColor3f(1.0, 0.0, 0.0) ;
                          } else {
                              glColor3f(0.0, 0.0, 0.0) ;
                          }
                          glVertex(seeds_[v1]) ;
                          glVertex(seeds_[v2]) ;
                      }
                  }
              }
          }
      }        
  private:
      const std::vector<vec3>& seeds_ ;
      std::map<bindex, unsigned int>& edges_ ;
  } ;



  void RVDGraphics::draw_primal_mesh() {
      glDepthRange(0.0, 0.999) ;
      glDisable(GL_LIGHTING) ;
      glLineWidth(2.0) ;
      glColor3f(0.0, 0.0, 0.0) ;

      if(show_non_manifold_) {
          if(!playback_list(u_singular_mesh_)) {
              begin_list(u_singular_mesh_) ;
              geometry_->RVD().set_symbolic(true) ;
              std::map<bindex, unsigned int> edges ;
              geometry_->RVD().for_each_primal_triangle(CountPrimalSingularEdges(edges)) ;
              glBegin(GL_LINES) ;
              geometry_->RVD().for_each_facet(
                  DrawPrimalSingularEdges(geometry_->vertices(), edges)
              ) ;
              glEnd() ;
              geometry_->RVD().set_symbolic(false) ;
              end_list() ;
          }
      } else {
          glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) ;
          if(!playback_list(u_primal_mesh_)) {
              begin_list(u_primal_mesh_) ;
              glBegin(GL_TRIANGLES) ;
              geometry_->RVD().for_each_primal_triangle(DrawPrimal(geometry_->vertices())) ;
              glEnd() ;
              end_list() ;
          }
          glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;
      }
      glDepthRange(0.0, 1.0) ;
  }

  static void color_ramp(
      GLfloat* cmap, unsigned int i1, vec3 c1, unsigned int i2, vec3 c2
  ) {
      for(unsigned int i=i1; i<=i2; i++) {
          double s = double(i-i1) / double(i2-i1) ;
          vec3 c = s*c2 + (1.0 - s)*c1 ;
          cmap[3*i]   = GLfloat(c.x) ;
          cmap[3*i+1] = GLfloat(c.y) ;
          cmap[3*i+2] = GLfloat(c.z) ;
      }
  }

  void RVDGraphics::create_colormap_texture() {
      GLfloat RGB[256*3] ;


      color_ramp(RGB, 0, vec3(0.0, 0.0, 1.0), 127, vec3(0.7, 0.0, 0.0)) ;
      color_ramp(RGB, 127, vec3(0.7, 0.0, 0.0), 255, vec3(1.0, 1.0, 0.0)) ;

      glGenTextures(1, &colormap_texture_) ;
      glEnable(GL_TEXTURE_1D) ;
      glBindTexture(GL_TEXTURE_1D, colormap_texture_) ;
      gluBuild1DMipmaps(GL_TEXTURE_1D, GL_RGB, 256, GL_RGB, GL_FLOAT, RGB) ; 
      glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);            
      glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR) ;
      glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR) ;
      glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP) ;
      glDisable(GL_TEXTURE_1D) ;
  }


}
