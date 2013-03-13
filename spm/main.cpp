
#include <cstdlib>
#include <ctime>
#include <string>
#include <Geex/graphics/geexapp.h>
#include <Geex/basics/file_system.h>
#include <glut_viewer/tweak_bar.h>
#include <CGAL/Timer.h>
#include "spm_cgal.h"
#include "SPM.h"
#include "project_io.h"

namespace Geex {
	//void TW_CALL tw_reset(void *clientData);
	//void TW_CALL tw_adjust(void *clientData);
	//void TW_CALL tw_lloyd(void *clientData);
	//void TW_CALL tw_pack(void *clientData);
	//void TW_CALL tw_rpack(void*);
	//void TW_CALL tw_replace(void *clientData);
	//void TW_CALL tw_detect_holes(void*);
	//void TW_CALL tw_save(void*);
	//void TW_CALL tw_fill(void*);
	//void TW_CALL tw_merge(void*);
	//void TW_CALL tw_feature_set_callback(const void*, void *);
	//void TW_CALL tw_feature_get_callback(void*, void*);
	//void TW_CALL tw_curvature_get_callback(void*, void*);
	//void TW_CALL tw_curvature_set_callback(const void*, void*);
	//void TW_CALL tw_correct_normals(void*);
	//void TW_CALL tw_smooth_normals(void*);
	//void TW_CALL tw_flip_normals(void*);
	//void TW_CALL tw_dump(void*);
	//void TW_CALL tw_restore(void*);
	//void TW_CALL tw_cluster_detect_holes(void*);
	//void TW_CALL tw_affine_fill(void*);
	
    class SPMApp : public GeexApp 
	{
    public:
        SPMApp(int argc, char** argv) : GeexApp(argc, argv) 
		{ 
            hdr_ = false ;
			if (argc < 2)
				prompt_and_exit("Error: No project configuration file input!");
			prj_config_file = argv[1];
       }

        SPM* spm() { return static_cast<SPM*>(scene()) ; }

        void init_scene() 
		{
            scene_ = new SPM;
			spm()->load_project(prj_config_file);
        }

        void init_gui() 
		{
            GeexApp::init_gui() ;
			
            TwBar* graphics_bar = TwNewBar("Graphics");
			TwDefine("Graphics position='16 10' size='200 250' alpha=200"); 
			TwAddVarRW(graphics_bar, "Domain Mesh", TW_TYPE_BOOL8, &spm()->show_mesh(), "");
			TwAddVarRW(graphics_bar, "Triangulation", TW_TYPE_BOOL8, &spm()->show_triangulation(), "");
			TwAddVarRW(graphics_bar, "Voronoi Cell", TW_TYPE_BOOL8, &spm()->show_voronoi_cell(), "");
			TwAddVarRW(graphics_bar, "Highlight", TW_TYPE_INT32, &spm()->highlighted_group_id(), "");
			TwBar* function_bar = TwNewBar("Functions");
			TwDefine("Functions position='16 250' size='200 250' alpha=200");
			//TwAddVarRW(function_bar, "Iter Number.", TW_TYPE_INT32, &spm()->setPackIterLim(), "min=1");
			//TwAddVarRW(function_bar, "Min Scale", TW_TYPE_DOUBLE, &spm()->setMinScalor(), ""/*"min=0.01 max=0.99"*/);
			//TwAddVarRW(function_bar, "Max Scale", TW_TYPE_DOUBLE, &spm()->setMaxScalor(), ""/*"min=1.0 max=1.99"*/);
			//TwAddVarRW(function_bar, "Area Coverage", TW_TYPE_DOUBLE, &spm()->setAreaCoverage(), "min=0.01 max=1.0");
			//TwAddButton(function_bar, "Lloyd", tw_lloyd, NULL, "key=k");
			//TwAddVarRW(function_bar, "Scale Factor", TW_TYPE_DOUBLE, &scalingFactor, "min=0.0 step=0.01");
			//TwAddButton(function_bar, "Scale", tw_adjust, NULL, "key=a");
			toggle_skybox_CB() ;
            glut_viewer_add_toggle('b', glut_viewer_is_enabled_ptr(GLUT_VIEWER_BACKGROUND), "switch Color/BW") ;
            glut_viewer_add_toggle('T', &viewer_properties_->visible(), "viewer properties") ;
        }

    private:
		std::string prj_config_file;
    } ;


Geex::SPMApp* spm_app() { return static_cast<Geex::SPMApp*>(Geex::GeexApp::instance()) ; }

//void TW_CALL tw_reset(void *clientData) { spm_app()->reset(); }
//void TW_CALL tw_adjust(void *clientData) { spm_app()->adjust(); }
//void TW_CALL tw_lloyd(void *clientData) {spm_app()->Lloyd();}
//void TW_CALL tw_pack(void *clientData) { spm_app()->pack(); }
//void TW_CALL tw_replace(void *clientData) { spm_app()->replace(); }
//void TW_CALL tw_detect_holes(void *clientData) { spm_app()->detect_holes(); }
//void TW_CALL tw_save( void *clientData ) { spm_app()->save(); }
//void TW_CALL tw_fill( void *clientData ) { spm_app()->fill(); }
//void TW_CALL tw_merge(void *clientData) { spm_app()->merge(); /*spm_app()->fitPlanes();*/ }
//void TW_CALL tw_rpack(void *clientData) { spm_app()->rpack(); }
//void TW_CALL tw_feature_set_callback(const void* fc, void *clientData) 
//{ 
//	spm_app()->feature_set_callback(*(double*)fc); 
//}
//void TW_CALL tw_feature_get_callback(void *value, void *clientData) 
//{ 
//	spm_app()->feature_get_callback((double*)value); 
//}
//void TW_CALL tw_curvature_get_callback(void *value, void *clientData)
//{
//	spm_app()->curvature_get_callback((double*)value);
//}
//void TW_CALL tw_curvature_set_callback(const void *value, void *clientData)
//{
//	spm_app()->curvature_set_callback(*(double*)value);
//}
////void TW_CALL tw_correct_normals(void* clientData) { spm_app()->correct_normals(); }
////void TW_CALL tw_smooth_normals(void* clientData) { spm_app()->smooth_normals(); }
//void TW_CALL tw_flip_normals(void* clientData) { spm_app()->flip_normals(); }
//void TW_CALL tw_dump(void *clientData) { spm_app()->dump(); }
//void TW_CALL tw_restore(void *clientData) { spm_app()->restore(); }
//void TW_CALL tw_cluster_detect_holes(void* clientData) { spm_app()->cluster_detect_holes(); }
//void TW_CALL tw_affine_fill(void* clientData) { spm_app()->affine_fill(); }
}

int main(int argc, char** argv) 
{
	//Geex::ProjectIO pio(argv[1]);
	//pio.debug_print();
    Geex::initialize();
    Geex::SPMApp app(argc, argv) ;
	app.main_loop() ;
    Geex::terminate() ;
	//system("pause");
	return 0;
}
