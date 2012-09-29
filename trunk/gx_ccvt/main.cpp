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

#include <Geex/graphics/geexapp.h>
#include <Geex/basics/file_system.h>
#include <Geex/basics/processor.h>
#include "cvt.h"
#include <ctime>
#include <cstdlib>
#include <AntTweakBar.h>
#include "meshes.h"

namespace Geex {

    class CVTApp : public GeexApp {
    public:
        CVTApp(int argc, char** argv) : GeexApp(argc, argv) { 
            hdr_ = false ;
            boundary_filename_ = get_file_arg("obj") ;
            if(boundary_filename_.length()==0)
                boundary_filename_ = get_file_arg("tri") ;
            if(boundary_filename_.length() > 0) {
                if(!Geex::FileSystem::is_file(boundary_filename_)) {
                    boundary_filename_ = Geex::FileSystem::get_project_root() + "/gx_ccvt/" + boundary_filename_ ;
                }
            }
            points_filename_ = get_file_arg("pts") ;
            if(points_filename_.length() > 0) {
                if(!Geex::FileSystem::is_file(points_filename_)) {
                    points_filename_ = Geex::FileSystem::get_project_root() + "/gx_ccvt/" + points_filename_ ;
                }
            }
            nb_points_ = 100 ;
            get_arg("nb_pts", nb_points_) ;
            nb_iter_ = 10 ;
            get_arg("nb_iter", nb_iter_) ;            
            non_convex_ = GL_FALSE ;
            get_arg("non_convex", non_convex_) ;
            vertex_type_ = 1 ;
            get_arg("vertex_type", vertex_type_) ;
            gamma_ = 1. ;
            get_arg("gamma", gamma_) ;
			perturb_ = 0.0f ;
			get_arg("perturb", perturb_) ;
            disable_normalize_ = GL_FALSE;
            get_arg("non_normalize", disable_normalize_);
        }

        ~CVTApp() {
            std::cout << "desconstruct cvt app. " << std::endl ;
        }

        CVT* cvt() { return static_cast<CVT*>(scene()) ; }

        virtual void init_scene() {
            scene_ = new CVT ;
            std::cerr << "Non convex = " << (non_convex_ ? "true" : "false") << "(use +non_convex to set)" << std::endl ;

            cvt()->set_non_convex_mode(non_convex_) ;
            if(boundary_filename_.length() > 0) {
                cvt()->load_boundary(boundary_filename_, perturb_, disable_normalize_) ;
            }

            cvt()->clipped_boundary_only() = non_convex_ ;
            if(LoadFeature()) {
                std::cout << "feature is loaded..." << std::endl ;
            } 
            if(LoadDensity()) {
                std::cout << "density is loaded..." << std::endl ;
            }
            int nb_cores = Geex::Processor::number_of_cores() ;
			if(nb_cores > 1)
				cvt()->partition_mesh() ; // for parallel computation


            CGAL::Timer timer ;
            std::cout << "begin inserting vertices: " << std::endl ;
            timer.start() ;
            if(points_filename_.size() > 0) 
                cvt()->insert_vertices(points_filename_) ;
            else
                cvt()->insert_vertices(nb_points_, vertex_type_) ;
            timer.stop() ;
            std::cout << timer.time() <<"s" << std::endl ;
        }

        void Lloyd() {
            CGAL::Timer timer ;
            timer.start() ;
            cvt()->lloyd(nb_iter_) ;
            timer.stop() ;
            if(cvt()->weighted())
                std::cout << nb_iter_ << " weighted lloyd iterations: " << timer.time() << 's' << std::endl ;
            else
                std::cout << nb_iter_ << " lloyd iterations: " << timer.time() << 's' << std::endl ;
        }

        void NewtonLloyd() {
            CGAL::Timer timer ;

            if(cvt()->is_clip_feature()) { // project feature seeds
				cvt()->compute_rvd() ;
                cvt()->lloyd(1) ;
            }

            timer.start() ;
            cvt()->newton_lloyd(nb_iter_) ;
            timer.stop() ;
            std::cout << nb_iter_ << " newton lloyd iterations: " << timer.time() << 's' << std::endl ;
        }

        void reset() {
            cvt()->clear() ;
            if(points_filename_.size() > 0) 
                cvt()->insert_vertices(points_filename_) ;
            else
                cvt()->insert_vertices(nb_points_, vertex_type_) ;
        }

        void RefineFeature() {
            cvt()->refine_feature() ;
        }

        void ExtractMesh() {
            cvt()->extract_mesh() ;
        }

        void ExtractDual() {
            cvt()->extract_dual() ;
        }

        void SaveDelaunay() {
            std::string data_filename = boundary_filename_ ;
            int size = data_filename.size() ;
            data_filename[size-1] = 's';
            data_filename[size-2] = 't';
            data_filename[size-3] = 'p';
            cvt()->save_delaunay(data_filename) ;
        }

        void LoadDelaunay() {
            std::string data_filename = boundary_filename_ ;
            int size = data_filename.size() ;
            data_filename[size-1] = 's';
            data_filename[size-2] = 't';
            data_filename[size-3] = 'p';
            cvt()->load_delaunay(data_filename) ;
        }

        bool LoadFeature() {
            std::string data_filename = boundary_filename_ ;
            int size = data_filename.size() ;
            data_filename[size-1] = 'a';
            data_filename[size-2] = 'e';
            data_filename[size-3] = 'f';
            return cvt()->load_feature(data_filename) ;            
        }

        bool LoadDensity() {
            std::string data_filename = boundary_filename_ ;
            int size = data_filename.size() ;
            data_filename[size-1] = 'r';
            data_filename[size-2] = 'u';
            data_filename[size-3] = 'c';
            if(cvt()->load_density(data_filename, gamma_) ) {
                cvt()->weighted() = GL_TRUE ;
                //        cvt()->voronoi_clip() ;
                return true ;
            }
            return false ;
        }

        void SaveRemesh() {
            std::string clip_filename = boundary_filename_ ;            
            std::string remesh_filename = boundary_filename_ ;
            std::string normalized_filename = boundary_filename_ ;
            //        clip_filename.insert(clip_filename.size()-4, "_clip") ;
            //    cvt()->save_clipped_mesh(clip_filename) ;        
            remesh_filename.insert(remesh_filename.size()-4, "_remesh") ;
            cvt()->save_remesh(remesh_filename) ;
            // save a normalized mesh for computing L^2 dist
            normalized_filename.insert(normalized_filename.size()-4, "_normalized") ;
            cvt()->save_normalized(normalized_filename) ;        
        }

        void LoadRemesh() {
            std::string remesh_filename = boundary_filename_ ;
            remesh_filename.insert(remesh_filename.size()-4, "_remesh") ;
            cvt()->load_remesh(remesh_filename) ;
        }

        void Timing() {
            cvt()->voronoi_clip_timing(nb_points_) ;
        }

        void ExportEPS() {
            std::string filename = boundary_filename_ ;
            int len = filename.length() ;
            filename[len-1] = 'f' ;
            filename[len-2] = 'd' ;
            filename[len-3] = 'p' ;
            std::cout << "export eps: " << filename << std::endl ;
            return cvt()->export_eps(filename) ;
        }

		void PerturbVertices() {
			cvt()->perturb_boundary_vertices() ;
		}

        void VoronoiFiltering()
        {
            cvt()->voronoi_filtering();
        }

        void AlphaShape()
        {
            cvt()->alpha_shape();
        }
		void ComputeLFS() {
			cvt()->compute_lfs() ;
		}

		void MergeSameVertices() {
			std::string filename = boundary_filename_ ;
			filename.insert(filename.size()-4, "_merged") ;
			std::cout << "merging same vertices of boundary mesh " << filename << std::endl ;
			cvt()->boundary().merge_same_vertices(64, 1e-5) ;
			cvt()->boundary().merge_same_vertices(63, 1e-5) ;
			cvt()->boundary().save_obj(filename) ;
		}

        virtual void init_gui() {
            GeexApp::init_gui() ;

            // New-style GUI =====================================================================================
            TwBar* graphics_bar = TwNewBar("Graphics") ;
            TwAddVarRW(graphics_bar, "Cull Back", TW_TYPE_BOOL8, &cvt()->set_cull_back(), "") ;
            TwAddVarRW(graphics_bar, "Domain mesh", TW_TYPE_BOOL8, &cvt()->show_domain_mesh(), "") ;
			TwAddVarRW(graphics_bar, "Domain", TW_TYPE_BOOL8, &cvt()->show_domain(), "") ;
			TwAddVarRW(graphics_bar, "Features", TW_TYPE_BOOL8, &cvt()->show_volume(), "") ;
            TwAddVarRW(graphics_bar, "Vertices", TW_TYPE_FLOAT, &cvt()->vertices_size(), "min=0 max=1 step=0.01") ;

            TwAddSeparator(graphics_bar, "Dual graphics", "") ;
			TwAddVarRW(graphics_bar, "Mesh", TW_TYPE_BOOL8, &cvt()->show_mesh(), "") ;
			TwAddVarRW(graphics_bar, "Colorize", TW_TYPE_BOOL8, &cvt()->colorize(), "") ;
			TwAddVarRW(graphics_bar, "Primal", TW_TYPE_BOOL8, &cvt()->show_primal(), "") ;
			TwAddVarRW(graphics_bar, "Remesh", TW_TYPE_BOOL8, &cvt()->show_remesh(), "") ;
            TwAddVarRW(graphics_bar, "Slice", TW_TYPE_FLOAT, &cvt()->slice_z(), "min=0 max=1 step=0.01") ;
			TwEnumVal rvd_mode_def[] = {
                {TRI_KD_INCIDENT,"TriKD"}, {POLY_KD_INCIDENT,"PolyKD"}, 
				{TRI_LINEAR,"TriLinear"}, {POLY_LINEAR,"PolyLinear"}, 
				{TRI_KD_SPLIT,"TriSplit"}, {PARA_TRI, "TriParallel"}, {PARA_POLY, "PolyParallel"}
            } ;
			TwType    tw_rvd_mode = TwDefineEnum("RVDMode", rvd_mode_def, 7) ;
			TwAddVarRW(graphics_bar, "RVD Mode", tw_rvd_mode, &cvt()->clip_mode(), "") ;
			TwAddVarRW(graphics_bar, "Surface RVD", TW_TYPE_BOOL8, &cvt()->show_surface_clipping(), "") ;
			TwAddVarRW(graphics_bar, "Density", TW_TYPE_BOOL8, &cvt()->weighted(), "") ;
			TwAddVarRW(graphics_bar, "Surface Constrain", TW_TYPE_BOOL8, &cvt()->constrained(), "") ;
			TwAddVarRW(graphics_bar, "Feature Constrain", TW_TYPE_BOOL8, &cvt()->constrained(), "") ;
			TwAddVarRW(graphics_bar, "VoronoiFiltering", TW_TYPE_FLOAT, &cvt()->voronoi_fitering_angle(), "min=0 max=90 step=0.1") ;
            
			TwAddSeparator(graphics_bar, "Numerics", "") ;
            TwEnumVal optimizer_def[] = {
                {LBFGSB,"LBFGSB"}, {HLBFGS,"HLBFGS"}, {HM1QN3,"HM1QN3"}, {HCG,"HCG"}, {HLBFGS_HESS,"HESSIAN"}
            } ;
            TwType tw_optimizer = TwDefineEnum("Optimizer", optimizer_def, 5) ;
            TwAddVarRW(graphics_bar, "Optimizer", tw_optimizer, &cvt()->optimizer_mode(), "") ;
            TwAddVarRW(graphics_bar, "LBFGS:M", TW_TYPE_FLOAT, &cvt()->set_opt_m(), "min=0 max=20 step=1") ;
//            TwAddVarRW(graphics_bar, "LBFGS:T", TW_TYPE_FLOAT, &cvt()->set_opt_t(), "min=0 max=20 step=1") ;

            //TwEnumVal surface_mode_def[] = {
            //    {None, "None"}, {Plain, "Plain"}, {Cells, "Cells"}, {Density, "Density"}, {Energy, "Energy"}
            //} ;
//            TwType tw_surface_mode = TwDefineEnum("SurfaceType", surface_mode_def, 5) ;
//            TwAddVarRW(graphics_bar, "Surface", tw_surface_mode, &cvt()->surface_mode(), "") ;
//            TwAddVarRW(graphics_bar, "RVD mesh", TW_TYPE_BOOL8, &cvt()->show_rvd_mesh(), "") ;            
//            TwAddVarRW(graphics_bar, "Orig. mesh", TW_TYPE_BOOL8, &cvt()->show_domain_mesh(), "") ;            
//            TwAddVarRW(graphics_bar, "Vertices", TW_TYPE_FLOAT, &cvt()->vertices_size(), "min=0 max=1 step=0.01") ;
//            TwAddVarRW(graphics_bar, "Shrink", TW_TYPE_FLOAT, &cvt()->shrink(), "min=0 max=1 step=0.01") ;
//            TwAddSeparator(graphics_bar, "Primal graphics", "") ;
//            TwAddVarRW(graphics_bar, "Primal", TW_TYPE_BOOL8, &cvt()->show_primal(), "") ;            
//            TwAddVarRW(graphics_bar, "N-manifold", TW_TYPE_BOOL8, &cvt()->show_non_manifold(), "") ;            
//            TwAddSeparator(graphics_bar, "Cells graphics", "") ;
//            TwEnumVal axis_def[] = { {AXIS_X, "X"}, {AXIS_Y, "Y"}, {AXIS_Z, "Z"} } ;
//            TwType tw_axis = TwDefineEnum("AxisType", axis_def, 3) ;
//            TwAddVarRW(graphics_bar, "Axis", tw_axis, &cvt()->slice_axis(), "") ;
//            TwAddVarRW(graphics_bar, "Slice", TW_TYPE_FLOAT, &cvt()->slice(), "min=0 max=1 step=0.01") ;
//            TwAddSeparator(graphics_bar, "Editing graphics", "") ;
//            TwAddVarRW(graphics_bar, "Sel. only", TW_TYPE_BOOL8, &cvt()->selected_only(), "") ;            
//            TwAddVarRW(graphics_bar, "Edit", TW_TYPE_BOOL8, &edit_, "") ;          

            //======================================================================================================
			
			// Old-stype GUI
            viewer_properties_->add_separator("Graphics") ;
			viewer_properties_->add_toggle("Cull back", cvt()->set_cull_back()) ;
            viewer_properties_->add_toggle("Domain mesh", cvt()->show_domain_mesh()) ;
            viewer_properties_->add_toggle("Domain", cvt()->show_domain()) ;
            viewer_properties_->add_toggle("Features", cvt()->show_volume()) ;
            viewer_properties_->add_slider("Vertices", cvt()->vertices_size()) ;
            viewer_properties_->add_toggle("Mesh", cvt()->show_mesh()) ;
            viewer_properties_->add_toggle("Colorize", cvt()->colorize()) ;
            viewer_properties_->add_toggle("Primal", cvt()->show_primal()) ;
            viewer_properties_->add_toggle("Remesh", cvt()->show_remesh()) ;
            viewer_properties_->add_slider("Slice", cvt()->slice_z(), -1.2, 1.2) ;
            viewer_properties_->add_enum("ClipMode", cvt()->clip_mode(), GlutViewerGUI::LabelList() | ClipModeNames) ;
            viewer_properties_->add_toggle("Surface RVD", cvt()->show_surface_clipping()) ;
            viewer_properties_->add_toggle("Density", cvt()->weighted()) ;
            viewer_properties_->add_toggle("Surface Constrain", cvt()->constrained()) ;
            viewer_properties_->add_toggle("Feature Constrain", cvt()->is_clip_feature()) ;
            viewer_properties_->add_slider("VoronoiFiltering", cvt()->voronoi_fitering_angle(), 0.0, 90.0) ;
            viewer_properties_->add_enum(
                "Optimizer", cvt()->optimizer_mode(), GlutViewerGUI::LabelList() | OptimizerModeNames
                ) ;
            viewer_properties_->add_slider("LBFGS:M", cvt()->set_opt_m(), 0.0, 45.0) ;
            toggle_skybox_CB() ;

            glut_viewer_add_toggle('B', glut_viewer_is_enabled_ptr(GLUT_VIEWER_BACKGROUND), "switch Color/BW") ;
            glut_viewer_add_toggle('T', &viewer_properties_->visible(), "viewer properties") ;
        }

    private:
        std::string boundary_filename_ ;
        std::string points_filename_ ; // input points for Delaunay Triangulation
        int nb_points_ ;
        int nb_iter_ ;
        int vertex_type_ ;
        GLboolean non_convex_ ;
        GLfloat gamma_ ; // grading factor of density function
		GLfloat perturb_ ;
        GLboolean disable_normalize_;
    } ;
}

Geex::CVTApp* cvt_app() { return static_cast<Geex::CVTApp*>(Geex::GeexApp::instance()) ; }

void Lloyd() {
    cvt_app()->Lloyd() ;
}

void NewtonLloyd() {
    cvt_app()->NewtonLloyd() ;
}

void reset() {
    cvt_app()->reset() ;
}

void RefineFeature() {
    cvt_app()->RefineFeature() ;
}

void ExtractMesh() {
    cvt_app()->ExtractMesh() ;
}

void ExtractDual() {
    cvt_app()->ExtractDual() ;
}

void LoadRemesh() {
    cvt_app()->LoadRemesh() ;
}

void SaveDelaunay() {
    cvt_app()->SaveDelaunay() ;  
}

void LoadDelaunay() {
    cvt_app()->LoadDelaunay() ;  
}

void LoadDensity() {
    cvt_app()->LoadDensity() ;
}

void SaveRemesh() {
    cvt_app()->SaveRemesh() ;
}

void ExportEPS() {
    cvt_app()->ExportEPS() ;
}
void Timing() {
    cvt_app()->Timing() ;
}
void PerturbVertices() {
	cvt_app()->PerturbVertices() ;
}

void VoronoiFiltering()
{
    cvt_app()->VoronoiFiltering();
}

void AlphaShape()
{
    cvt_app()->AlphaShape();
}
void ComputeLFS() {
	cvt_app()->ComputeLFS() ;
}
void MergeSameVertices() {
	cvt_app()->MergeSameVertices() ;
}

#if (defined(WIN32) && defined(_DEBUG))
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#endif

#include <iostream>
#include <fstream>

int main(int argc, char** argv) {

#if (defined(WIN32) && defined(_DEBUG))
    _CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
#endif

    srand((unsigned int)time(0));

    Geex::initialize() ;
    glut_viewer_add_key_func('k', Lloyd, "Lloyd iterations") ;
    glut_viewer_add_key_func('m', NewtonLloyd, "Newton-Lloyd iterations") ;
    glut_viewer_add_key_func('Z', reset, "reset") ;
    glut_viewer_add_key_func('s', SaveDelaunay, "save delaunay data") ;
    glut_viewer_add_key_func('o', LoadDelaunay, "load delaunay data") ;
    glut_viewer_add_key_func('d', LoadDensity, "load density function") ;
    glut_viewer_add_key_func('r', ExtractMesh, "extract remesh") ;
	glut_viewer_add_key_func('e', ExtractDual, "extract dual mesh") ;
    glut_viewer_add_key_func('R', LoadRemesh, "load remesh") ;
    glut_viewer_add_key_func('c', SaveRemesh, "save clipped mesh") ;
    glut_viewer_add_key_func('f', RefineFeature, "refine features") ;
    glut_viewer_add_key_func('v', Timing, "voronoi clipping timing") ;
    glut_viewer_add_key_func('1', ExportEPS, "export eps") ;
	glut_viewer_add_key_func('p', PerturbVertices, "perturb vertices") ;
    glut_viewer_add_key_func('V', VoronoiFiltering, "Voronoi filtering") ;
    glut_viewer_add_key_func('A', AlphaShape, "Alpha Shape") ;
	glut_viewer_add_key_func('C', ComputeLFS, "Compute LFS") ;
	glut_viewer_add_key_func('x', MergeSameVertices, "Merge same boundary vertices") ;
    glut_viewer_add_toggle('2', glut_viewer_is_enabled_ptr(GLUT_VIEWER_3D), "toggle 2D display") ;

    Geex::CVTApp app(argc, argv) ;
    app.main_loop() ;
    Geex::terminate() ;
}
