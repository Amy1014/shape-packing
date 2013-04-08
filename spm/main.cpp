
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
	void TW_CALL tw_lloyd(void *clientData);
	void TW_CALL tw_pack(void *clientData);
	void update();
	void TW_CALL tw_idt_update(void *clientData);
	void TW_CALL tw_detect_holes(void*);
	void TW_CALL tw_hole_size_get_callback(void*, void*);
	void TW_CALL tw_hole_size_set_callback(const void*, void *);
	void TW_CALL tw_front_len_get_callback(void*, void*);
	void TW_CALL tw_front_len_set_callback(const void*, void*);
	void TW_CALL tw_replace(void*);
	void TW_CALL tw_enlarge(void*);
	void TW_CALL tw_save_triangulation(void*);
	//void TW_CALL tw_save(void*);
	void TW_CALL tw_fill(void*);
	//void TW_CALL tw_merge(void*);
	//void TW_CALL tw_feature_set_callback(const void*, void *);
	//void TW_CALL tw_feature_get_callback(void*, void*);
	//void TW_CALL tw_curvature_get_callback(void*, void*);
	//void TW_CALL tw_curvature_set_callback(const void*, void*);
	//void TW_CALL tw_dump(void*);
	//void TW_CALL tw_restore(void*);
	//void TW_CALL tw_cluster_detect_holes(void*);
	
    class SPMApp : public GeexApp 
	{
    public:
        SPMApp(int argc, char** argv) : GeexApp(argc, argv) 
		{ 
            hdr_ = false ;
			if (argc < 2)
				prompt_and_exit("Error: No project configuration file input!");
			prj_config_file = argv[1];

			enlarge_id = -1;
			enlarge_factor = 1.0;
			enlarge_theta = 0.0;
			enlarge_tx = 0.0;
			enlarge_ty = 0.0;
       }

        SPM* spm() { return static_cast<SPM*>(scene()) ; }

        void init_scene() 
		{
            scene_ = new SPM;
			spm()->load_project(prj_config_file);
        }

		/** optimization **/
		void lloyd()
		{
			//spm()->lloyd();
			//spm()->redraw_triangulation();
			//spm()->redraw_voronoi_cell();
			//glut_viewer_redraw();
			spm()->lloyd(&update, false);
			//spm()->redraw_voronoi_cell();
			//spm()->redraw_triangulation();
			glut_viewer_redraw();
		}

		void post_update()
		{
			spm()->redraw_triangulation();
			spm()->redraw_voronoi_cell();
			glut_viewer_redraw();
		}

		void pack()
		{
			//spm()->pack(&update);
			spm()->pack(&update);
		}

		void idt_update()
		{

			spm()->update_iDT();
			spm()->redraw_voronoi_cell();
			spm()->redraw_triangulation();

			glut_viewer_redraw();
		}

		void detect_holes()
		{
			spm()->detect_holes();
		}
		void fill_holes()
		{
			spm()->fill_holes();
			post_update();
		}
		void hole_size_set_callback(double holesize)
		{
			spm()->hole_size(holesize);
			spm()->detect_holes();
		}
		void hole_size_get_callback(double *holesize)
		{
			*holesize = spm()->hole_size();
		}
		void front_len_set_callback(double front_len)
		{
			spm()->front_edge_len(front_len);
			spm()->detect_holes();
		}
		void front_len_get_callback(double *front_len)
		{
			*front_len = spm()->front_edge_len();
		}

		void replace()
		{
			spm()->replace();
			post_update();
		}

		void enlarge_polygon()
		{
			if (enlarge_id >= 0 && enlarge_factor >= 1.0)
				spm()->enlarge_one_polygon(enlarge_id, enlarge_factor, enlarge_theta, enlarge_tx, enlarge_ty);
			post_update();
		}
		void save_triangulation()
		{
			spm()->rpvd.save_triangulation("enforced_enlarge.obj");
		}
        void init_gui() 
		{
            GeexApp::init_gui() ;
			
            TwBar* graphics_bar = TwNewBar("Graphics");
			TwDefine("Graphics position='16 10' size='200 250' alpha=200"); 
			TwAddVarRW(graphics_bar, "Domain Mesh", TW_TYPE_BOOL8, &spm()->show_mesh(), "group = 'Surface' ");
			TwAddVarRW(graphics_bar, "Curvature", TW_TYPE_BOOL8, &spm()->show_curvatures(), "group = 'Surface' ");
			
			
			TwEnumVal draw_polygon_mode[] = {
				{spm()->OUTLINE_DRAW, "Outline"}, {spm()->FILL_DRAW, "Fill"}, {spm()->TEXTURE_DRAW, "Texture"}
			};
			TwType tw_draw_polygon_mode = TwDefineEnum("DrawPolygonMode", draw_polygon_mode, 3);
			TwAddVarRW(graphics_bar, "Draw Switch", TW_TYPE_BOOL8, &spm()->show_polygons(), "group = 'Polygon' ");
			TwAddVarRW(graphics_bar, "Draw Mode", tw_draw_polygon_mode, &spm()->get_polygon_draw_type(), "group = 'Polygon' ");
		
			TwAddVarRW(graphics_bar, "Triangulation", TW_TYPE_BOOL8, &spm()->show_triangulation(), "group = 'Geometry' ");
			TwAddVarRW(graphics_bar, "Voronoi Cell", TW_TYPE_BOOL8, &spm()->show_voronoi_cell(), "group = 'Geometry' ");
			TwAddVarRW(graphics_bar, "Vertices", TW_TYPE_BOOL8, &spm()->show_vertices(), "group = 'Geometry' ");

			TwAddVarRW(graphics_bar, "Highlight", TW_TYPE_INT32, &spm()->highlighted_group_id(), "group = 'Debug' ");
			TwAddVarRW(graphics_bar, "Local frame", TW_TYPE_BOOL8, &spm()->show_local_frame(), "group = 'Debug' ");

			TwAddVarRW(graphics_bar, "Hole Tri", TW_TYPE_BOOL8, &spm()->show_hole_triangles(), "group = 'Geometry' ");
			TwAddVarRW(graphics_bar, "Holes", TW_TYPE_BOOL8, &spm()->show_holes(), "group = 'Geometry' ");
			TwBar* function_bar = TwNewBar("Functions");
			TwDefine("Functions position='16 250' size='200 250' alpha=200");
			TwAddButton(function_bar, "Lloyd", tw_lloyd, NULL, "key=l");
			TwAddButton(function_bar, "Pack", tw_pack, NULL, "key=p");
			TwAddButton(function_bar, "iDT", tw_idt_update, NULL, "key=i");
			TwAddButton(function_bar, "Detect Holes", tw_detect_holes, NULL, "key=d");
			TwAddButton(function_bar, "Fill holes", tw_fill, NULL, "key=f");
			TwAddVarCB(function_bar, "Hole Size", TW_TYPE_DOUBLE, tw_hole_size_set_callback,
						tw_hole_size_get_callback, NULL, "min=0.0 step=0.001");
			TwAddVarCB(function_bar, "Front Edge", TW_TYPE_DOUBLE, tw_front_len_set_callback,
						tw_front_len_get_callback, NULL, "min=0.0 step=0.001");
			TwAddButton(function_bar, "replace", tw_replace, NULL, "key=r");
			TwAddVarRW(function_bar, "epsilon", TW_TYPE_DOUBLE, &spm()->get_epsilon(), "min=0.0 max=0.999999999");
			//TwAddVarRW(function_bar, "enlarge id", TW_TYPE_INT32, &enlarge_id, "");
			//TwAddVarRW(function_bar, "enlarge factor", TW_TYPE_DOUBLE, &enlarge_factor, "");
			//TwAddVarRW(function_bar, "rot angle", TW_TYPE_DOUBLE, &enlarge_theta, "");
			//TwAddVarRW(function_bar, "trans_x", TW_TYPE_DOUBLE, &enlarge_tx, "");
			//TwAddVarRW(function_bar, "trans_y", TW_TYPE_DOUBLE, &enlarge_ty, "");
			//TwAddButton(function_bar, "enlarge", tw_enlarge, NULL, "key=e");
			TwAddButton(function_bar, "save tri", tw_save_triangulation, NULL, "key=t");

			toggle_skybox_CB() ;
            glut_viewer_add_toggle('b', glut_viewer_is_enabled_ptr(GLUT_VIEWER_BACKGROUND), "switch Color/BW") ;
            glut_viewer_add_toggle('T', &viewer_properties_->visible(), "viewer properties") ;
        }

    private:
		std::string prj_config_file;

		//debug
		int enlarge_id;
		double enlarge_factor;
		double enlarge_theta;
		double enlarge_tx;
		double enlarge_ty;
    } ;


Geex::SPMApp* spm_app() { return static_cast<Geex::SPMApp*>(Geex::GeexApp::instance()) ; }

//void TW_CALL tw_reset(void *clientData) { spm_app()->reset(); }
//void TW_CALL tw_adjust(void *clientData) { spm_app()->adjust(); }
void TW_CALL tw_lloyd(void *clientData) {spm_app()->lloyd();}
void update() { spm_app()->post_update(); }
void TW_CALL tw_pack(void *clientData) { spm_app()->pack(); }
void TW_CALL tw_idt_update(void *clientData) { spm_app()->idt_update(); }
void TW_CALL tw_replace(void *clientData) { spm_app()->replace(); }
void TW_CALL tw_detect_holes(void *clientData) { spm_app()->detect_holes(); }
void TW_CALL tw_hole_size_set_callback(const void* size, void *clientData)
{
	spm_app()->hole_size_set_callback(*(double*)size);
}
void TW_CALL tw_hole_size_get_callback(void *size, void *clientData)
{
	spm_app()->hole_size_get_callback((double*)size);
}
void TW_CALL tw_front_len_set_callback(const void* size, void* clientData)
{
	spm_app()->front_len_set_callback(*(double*)size);
}
void TW_CALL tw_front_len_get_callback(void *size, void *clientData)
{
	spm_app()->front_len_get_callback((double*)size);
}
void TW_CALL tw_enlarge(void* clientData)
{
	spm_app()->enlarge_polygon();
}
void TW_CALL tw_save_triangulation(void* clientData)
{
	spm_app()->save_triangulation();
}
//void TW_CALL tw_save( void *clientData ) { spm_app()->save(); }
void TW_CALL tw_fill( void *clientData ) { spm_app()->fill_holes(); }
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
