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
#include <glut_viewer/tweak_bar.h>
#include <cstdlib>
#include "SPM.h"
#include "meshes.h"

namespace Geex {

	void TW_CALL tw_reset(void *clientData);
	void TW_CALL tw_adjust(void *clientData);
	
    class SPMApp : public GeexApp {
    public:
        SPMApp(int argc, char** argv) : GeexApp(argc, argv) { 
            hdr_ = false ;
            boundary_filename_ = get_file_arg("obj") ;
            if(boundary_filename_.length()==0)
                boundary_filename_ = get_file_arg("tri") ;
            if(boundary_filename_.length() > 0) {
                if(!Geex::FileSystem::is_file(boundary_filename_)) {
                    boundary_filename_ = Geex::FileSystem::get_project_root() + "/gx_ccvt/" + boundary_filename_ ;
                }
            }
            pgn_filename_ = get_file_arg("poly") ;
            if(pgn_filename_.length() > 0) {
                if(!Geex::FileSystem::is_file(pgn_filename_)) {
                    pgn_filename_ = Geex::FileSystem::get_project_root() + "/gx_ccvt/" + pgn_filename_ ;
                }
            }
            nb_iter_ = 10 ;
            get_arg("nb_iter", nb_iter_) ;            
			scalingFactor = 1.0;
        }

        SPM* spm() { return static_cast<SPM*>(scene()) ; }

        void init_scene() {
            scene_ = new SPM ;
			if (boundary_filename_.length()>0)
				spm()->loadMeshDomain(boundary_filename_);
			else
			{
				std::cerr << " \nFailure to load mesh domain\n";
				exit(1);
			}
			if (pgn_filename_.length()>0)
				spm()->loadPolygons(pgn_filename_);
			else
			{
				std::cerr<<"\nFailure to load polygons\n";
				exit(1);;
			}
			//glut_viewer_redraw();
        }

        //void Lloyd() {
        //    spm()->lloyd(nb_iter_) ;
        //}

		void adjust()
		{
			spm()->adjustPolygons(scalingFactor);
			glut_viewer_redraw();
		}
        void reset() {
            spm()->clear() ;
            spm()->loadPolygons(pgn_filename_);
			spm()->loadMeshDomain(boundary_filename_);
			glut_viewer_redraw();
            //else
            //    spm()->insert_vertices(nb_points_, vertex_type_) ;
        }

        virtual void init_gui() {
            GeexApp::init_gui() ;
			
            TwBar* graphics_bar = TwNewBar("Graphics");
			TwDefine("Graphics position='16 20' size='200 100' alpha=200"); 
			TwAddVarRW(graphics_bar, "Show Domain Mesh", TW_TYPE_BOOL8, &spm()->show_domain_mesh(), "");
			TwAddVarRW(graphics_bar, "Show Mesh Skeleton", TW_TYPE_BOOL8, &spm()->show_mesh_skeleton(), "");
			TwAddVarRW(graphics_bar, "Show Bounding Spheres", TW_TYPE_BOOL8, &spm()->show_bounding_spheres(), "");
			TwAddVarRW(graphics_bar, "Show Triangulation", TW_TYPE_BOOL8, &spm()->show_regular_triangulation(), "");
			TwBar* function_bar = TwNewBar("Functions");
			TwDefine("Functions position='16 220' size='200 100' alpha=200");
			TwAddVarRW(function_bar, "Scale Factor", TW_TYPE_DOUBLE, &scalingFactor, "min=0.0 step=0.01");
			TwAddButton(function_bar, "Scale", tw_adjust, NULL, "key=a");
			TwAddButton(function_bar, "Reset", tw_reset, NULL, "key=Z");
			toggle_skybox_CB() ;
            glut_viewer_add_toggle('B', glut_viewer_is_enabled_ptr(GLUT_VIEWER_BACKGROUND), "switch Color/BW") ;
            glut_viewer_add_toggle('T', &viewer_properties_->visible(), "viewer properties") ;
        }

    private:
        std::string boundary_filename_ ;
        std::string pgn_filename_ ; // input polygons
		double scalingFactor;
        int nb_iter_ ;
    } ;


Geex::SPMApp* spm_app() { return static_cast<Geex::SPMApp*>(Geex::GeexApp::instance()) ; }

void TW_CALL tw_reset(void *clientData) { spm_app()->reset(); }
void TW_CALL tw_adjust(void *clientData) { spm_app()->adjust(); }

}

int main(int argc, char** argv) 
{
    Geex::initialize();
    Geex::SPMApp app(argc, argv) ;
    app.main_loop() ;
    Geex::terminate() ;
}
