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

#include <Geex/basics/file_system.h>
#include <Geex/basics/processor.h>
#include "rvd_app.h"
#include <AntTweakBar.h>

namespace Geex {

  RVDApp::RVDApp(int argc, char** argv) : GeexApp(argc, argv) { 
    hdr_ = false ;

    aniso_filename_ = get_file_arg("eobj") ;
    boundary_filename_ = get_file_arg("obj") ;
    if(boundary_filename_.length() == 0) {
        boundary_filename_ = aniso_filename_ ;
    }
    if(boundary_filename_.length() == 0) {
        boundary_filename_ = "C.obj" ;
    }
    if(boundary_filename_.length() > 0) {
        if(!Geex::FileSystem::is_file(boundary_filename_)) {
            boundary_filename_ = Geex::FileSystem::get_project_root() + 
                "/gx_rvd" + boundary_filename_ ;
        }
    }
    if(aniso_filename_.length() > 0) {
        if(!Geex::FileSystem::is_file(aniso_filename_)) {
            aniso_filename_ = Geex::FileSystem::get_project_root() + 
                "/gx_rvd/" + aniso_filename_ ;
        }
    }

    nb_points_ = 0 ;
    get_arg("nb_pts", nb_points_) ;
    nb_iter_ = 30 ;
    get_arg("nb_iter", nb_iter_) ;            
    use_precond_ = GL_FALSE ;
    nb_threads_ = Geex::Processor::number_of_cores() ;
    std::cerr << "Number of cores = " << nb_threads_ << std::endl ;
    get_arg("nb_threads", nb_threads_) ;
    delaunay_algo_ = "CGAL_regular" ;
    get_arg("delaunay", delaunay_algo_) ;
    gamma_ = 1.0 ;
    get_arg("gamma", gamma_) ;
    exact_ = GL_FALSE ;
    get_arg("exact", exact_) ;
    if (exact_) {
        std::cerr << "Using exact mode"<< std::endl ;
    } else {
        std::cerr << "Using inexact mode"<< std::endl ;
    }
    output_filename_ = "out.obj" ;
    get_arg("out", output_filename_) ;
    edit_ = GL_FALSE ;
  }

  void RVDApp::save() {
    rvd()->save_primal(output_filename_) ;
  }

  void RVDApp::load() {
    rvd()->load_points(output_filename_) ;
  }

  void RVDApp::init_scene() {
    scene_ = new RVDob(delaunay_algo_) ;
    // rvd()->set_nb_threads(nb_threads_) ;
    if(boundary_filename_.length() > 0) {
        rvd()->load_boundary(boundary_filename_, gamma_) ;
    }
    if(aniso_filename_.length() > 0) {
        rvd()->load_anisotropy(aniso_filename_) ;                
    } else {
        rvd()->load_anisotropy(boundary_filename_) ;                
    }
    rvd()->exact() = exact_ ;
    if(nb_points_ >= 0) {
        rvd()->init_vertices(nb_points_) ;
    } else {
        rvd()->init_vertices(boundary_filename_) ;
    }
  }

  void RVDApp::init_gui() {
    GeexApp::init_gui() ;

    // New-style GUI =====================================================================================
    TwBar* graphics_bar = TwNewBar("Graphics") ;
    TwAddVarRW(graphics_bar, "DLists", TW_TYPE_BOOL8, &rvd()->use_display_lists(), "") ;
    TwAddVarRW(graphics_bar, "Shiny", TW_TYPE_BOOL8, &rvd()->shiny(), "") ;
    TwAddSeparator(graphics_bar, "Dual graphics", "") ;
    TwEnumVal surface_mode_def[] = {
        {None, "None"}, {Plain, "Plain"}, {Cells, "Cells"}, {Density, "Density"}, {Energy, "Energy"}, {Linear, "Linear"}, {Connected, "Connected"}
    } ;
    TwType tw_surface_mode = TwDefineEnum("SurfaceType", surface_mode_def, 7) ;
    TwAddVarRW(graphics_bar, "Surface", tw_surface_mode, &rvd()->surface_mode(), "") ;
    TwAddVarRW(graphics_bar, "RVD mesh", TW_TYPE_BOOL8, &rvd()->show_rvd_mesh(), "") ;            
    TwAddVarRW(graphics_bar, "Orig. mesh", TW_TYPE_BOOL8, &rvd()->show_domain_mesh(), "") ;            
    TwAddVarRW(graphics_bar, "Vertices", TW_TYPE_FLOAT, &rvd()->vertices_size(), "min=0 max=1 step=0.01") ;
    TwAddVarRW(graphics_bar, "Shrink", TW_TYPE_FLOAT, &rvd()->shrink(), "min=0 max=1 step=0.01") ;
    TwAddVarRW(graphics_bar, "Id Threshold", TW_TYPE_FLOAT, &rvd()->id_threshold(), "min=0 max=1 step=0.0001") ;
    TwAddVarRW(graphics_bar, "Lloyd factor", TW_TYPE_FLOAT, &rvd()->lloyd_factor(), "min=0 max=1 step=0.01") ;
    TwAddSeparator(graphics_bar, "Primal graphics", "") ;
    TwAddVarRW(graphics_bar, "Primal", TW_TYPE_BOOL8, &rvd()->show_primal(), "") ;            
	TwAddVarRW(graphics_bar, "Show vertices", TW_TYPE_BOOL8, &rvd()->show_vertices(), "") ;            
    TwAddVarRW(graphics_bar, "N-manifold", TW_TYPE_BOOL8, &rvd()->show_non_manifold(), "") ;            
    TwAddSeparator(graphics_bar, "Editing graphics", "") ;
    //TwAddVarRW(graphics_bar, "Sel. only", TW_TYPE_BOOL8, &rvd()->selected_only(), "") ;            
    TwAddVarRW(graphics_bar, "Edit", TW_TYPE_BOOL8, &edit_, "") ;            

    //TwBar* numerics_bar = TwNewBar("Numerics") ;
    //TwDefine("Numerics position='16 400'") ;

    //TwAddVarRW(numerics_bar, "Volume CVT", TW_TYPE_BOOL8, &rvd()->three_D(), "") ;            

    //TwAddVarRW(numerics_bar, "Bd. weight", TW_TYPE_FLOAT, &rvd()->boundary_weight(), "min=0 max=1 step=0.01") ;

    //TwEnumVal adt_mode_def[] = {
    //    {ANISO_DT_SPHERE, "sphere"}, {ANISO_DT_ELLIPSOID, "ellipsoid"}, {ANISO_DT_CUBE, "cube"}
    //} ;
    //TwType tw_adt_mode = TwDefineEnum("AdtMode", adt_mode_def, 3) ;
    //if(delaunay_algo_ == "CGAL_aniso" || delaunay_algo_ == "CGAL_aniso_sort") {
    //    TwAddVarRW(numerics_bar, "ADT mode", tw_adt_mode, &rvd()->adt_mode(), "") ;
    //}
    //if(aniso_) {
    //    TwAddVarRW(numerics_bar, "Lp order", TW_TYPE_FLOAT, &rvd()->Lp(), "min=0, max=6, step=1") ;
    //} ;

    //TwEnumVal center_mode_def[] = {
    //    {CENTROID, "centroid"}, {QUASI_CIRCUMCENTER, "quasi-circumcenter3D"}, 
    //    {QUASI_INCENTER, "quasi-incenter3D"}, {QUASI_MIDCENTER, "quasi-midcenter3D"}
    //};
    //TwType tw_center_mode = TwDefineEnum("CenterMode", center_mode_def, 4) ;
    //TwAddVarRW(numerics_bar, "CenterMode", tw_center_mode, &rvd()->center_mode(), "") ;

    //TwEnumVal optimizer_def[] = {
    //    {LBFGSB,"LBFGSB"}, {HLBFGS,"HLBFGS"}, {HM1QN3,"HM1QN3"}, {HCG,"HCG"}, {HLBFGS_HESS,"HESSIAN"}
    //} ;
    //TwType tw_optimizer = TwDefineEnum("Optimizer", optimizer_def, 5) ;
    //TwAddVarRW(numerics_bar, "Optimizer", tw_optimizer, &rvd()->optimizer_mode(), "") ;
    //TwAddVarRW(numerics_bar, "LBFGS:M", TW_TYPE_FLOAT, &rvd()->set_opt_m(), "min=0 max=20 step=1") ;
    //TwAddVarRW(numerics_bar, "LBFGS:T", TW_TYPE_FLOAT, &rvd()->set_opt_t(), "min=0 max=20 step=1") ;

    //TwAddVarRW(numerics_bar, "iteration", TW_TYPE_INT32, &rvd()->nb_iter(), "min=1 max=500 step=10") ;

    //======================================================================================================

    // Old-stype GUI
    viewer_properties_->hide() ;
    viewer_properties_->add_separator("Graphics") ;
    viewer_properties_->add_toggle("DLists", rvd()->use_display_lists()) ;
    viewer_properties_->add_enum(
        "Surface", rvd()->surface_mode(), GlutViewerGUI::LabelList() | SurfaceModeNames
    ) ;
    viewer_properties_->add_toggle("RVD  mesh", rvd()->show_rvd_mesh()) ;
    viewer_properties_->add_toggle("Orig mesh", rvd()->show_domain_mesh()) ;
    viewer_properties_->add_slider("Vertices", rvd()->vertices_size()) ;
    viewer_properties_->add_slider("Shrink",  rvd()->shrink()) ;
    viewer_properties_->add_slider("Id threshold",  rvd()->id_threshold()) ;
    viewer_properties_->add_separator("Primal") ;
    viewer_properties_->add_toggle("Primal", rvd()->show_primal()) ;
    viewer_properties_->add_toggle("Non manifold", rvd()->show_non_manifold()) ;
    viewer_properties_->add_toggle("Exact", rvd()->exact()) ;
    viewer_properties_->add_toggle("Edit", edit_) ;
    toggle_skybox_CB() ;
    
    glut_viewer_add_toggle(
        'B', glut_viewer_is_enabled_ptr(GLUT_VIEWER_BACKGROUND), "switch Color/BW"
    ) ;
    glut_viewer_add_toggle(
        'T', &viewer_properties_->visible(), "viewer properties"
    ) ;

  }

  void RVDApp::reset() {
    rvd()->init_vertices(nb_points_) ;
  }

  GLboolean RVDApp::mouse(float x, float y, int button, enum GlutViewerEvent event) {
    static int timestamp = 0 ;
    static int last_timestamp = 0 ;
    timestamp++ ;
    static int mode = 0 ;
    static bool down = false ;

    if(GeexApp::mouse(x, y, button, event)) { return GL_TRUE ; }
    if(edit_) {
        
        bool change = false ;
        if(event == GLUT_VIEWER_DOWN) { mode = button ; change = true ; down = true ; }
        if(event == GLUT_VIEWER_UP) { down = false ; }
        if(!down) { return GL_FALSE ; }

        vec3 p ;
        GLboolean hit_background ;
        glut_viewer_get_picked_point(&(p.x), &hit_background) ;
        if(hit_background) {
            return GL_FALSE ;
        }

        // Add vertex
        if(mode == 0 && (change || timestamp - last_timestamp > 5)) {
            last_timestamp = timestamp ;
            rvd()->vertices().push_back(p) ;
            rvd()->update(true) ;
        }

        // Move vertex
        if(mode == 1 && (change || timestamp - last_timestamp > 5)) {
            last_timestamp = timestamp ;
            unsigned int v = rvd()->delaunay()->nearest_vertex_id(p) ;
            rvd()->vertices()[v] = p ;
            rvd()->update(true) ;
        }

        // Remove vertex
        if(mode == 2 && (change || timestamp - last_timestamp > 5)) {
            last_timestamp = timestamp ;
            if(rvd()->vertices().size() > 4) {
                unsigned int v = rvd()->delaunay()->nearest_vertex_id(p) ;
                rvd()->vertices().erase(rvd()->vertices().begin() + v) ;
                rvd()->update(true) ;
            }
        }

        return GL_TRUE ;
    }
    return GL_FALSE ;
  }
}
