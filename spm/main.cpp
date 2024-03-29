
#include <cstdlib>
#include <ctime>
#include <string>
#include <Geex/graphics/geexapp.h>
#include <Geex/basics/file_system.h>
#include <glut_viewer/tweak_bar.h>
#include <Geex/third_party/misc/gl2ps.h>
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
	//void TW_CALL tw_enlarge(void*);
	void TW_CALL tw_save_triangulation(void*);
	//void TW_CALL tw_save_cur_area(void*);
	void TW_CALL tw_fill(void*);
	//void TW_CALL tw_remove(void*);
	//void TW_CALL tw_ex_replace(void*);
	void TW_CALL tw_con_replace(void *);
	void TW_CALL tw_pack_next(void*);
	void TW_CALL tw_report(void*);
	//void TW_CALL tw_discretize_tiles(void*);
	void TW_CALL tw_save_tiles(void *);
	void TW_CALL tw_save_mat(void *);
	void TW_CALL tw_swap(void*);
	void TW_CALL tw_selective_swap(void*);
	void TW_CALL tw_adjust(void *);
	void TW_CALL tw_activate(void*);
	void TW_CALL tw_auto_scale(void*);
	void TW_CALL tw_print_pdf(void *clientData);
	void TW_CALL tw_save_scale(void*);
	void TW_CALL tw_save_CAT(void*);

    class SPMApp : public GeexApp 
	{
    public:
        SPMApp(int argc, char** argv) : GeexApp(argc, argv) 
		{ 
            hdr_ = false ;
			if (argc < 2)
				prompt_and_exit("Error: No project configuration file input!");
			prj_config_file = argv[1];

			enlarge_factor = 0.98;

			lloyd_iter_times = 5;
			pack_iter_times = 5;
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
			for (int i = 0; i < lloyd_iter_times; i++)
			{
				spm()->interface_lloyd(&update);
				glut_viewer_redraw();
			}

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
			for (int i = 0; i < pack_iter_times; i++)
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
		
		void con_replace()
		{
			spm()->con_replace();
			post_update();
		}

		void swap()
		{
			spm()->swap();
			post_update();
		}
		void selective_swap()
		{
			spm()->selective_swap();
			post_update();
		}
		void adjust()
		{
			spm()->adjust(enlarge_factor);
			spm()->redraw_triangulation();
			spm()->redraw_voronoi_cell();
			glut_viewer_redraw();
		}

		void auto_adjust()
		{
			spm()->auto_adjust();
			post_update();
		}
		void activate_all()
		{
			spm()->activate_all();
		}
		void save_triangulation()
		{
			spm()->get_rpvd().save_triangulation("enforced_enlarge.obj");
		}
		void pack_next()
		{
			spm()->start_new_packing();
		}
		void report()
		{
			spm()->report();
		}

		void save_tiles()
		{
			spm()->save_tiles();
		}

		void save_scale()
		{
			spm()->save_scale_factors();
		}

		void save_mat()
		{
			spm()->save_materials();
		}
		
		void save_CAT()
		{
			spm()->save_CAT();
		}

        void init_gui() 
		{
            GeexApp::init_gui() ;
			
            TwBar* graphics_bar = TwNewBar("Graphics");
			TwDefine("Graphics position='16 10' size='200 300' alpha=200"); 
			TwAddVarRW(graphics_bar, "Domain Mesh", TW_TYPE_BOOL8, &spm()->show_mesh(), "group = 'Surface' ");
			TwAddVarRW(graphics_bar, "Curvature", TW_TYPE_BOOL8, &spm()->show_curvatures(), "group = 'Surface' ");
			TwAddVarRW(graphics_bar, "Feature", TW_TYPE_BOOL8, &spm()->show_feature_lines(), "group = 'Surface' ");
			TwAddVarRW(graphics_bar, "Vector Field", TW_TYPE_BOOL8, &spm()->show_vector_field(), "group = 'Surface'");
			//TwAddVarRW(graphics_bar, "Red", TW_TYPE_FLOAT, &spm()->tunable_red, "min = 0.0 max = 1.0 step = 0.01 group = 'Color' ");
			//TwAddVarRW(graphics_bar, "Green", TW_TYPE_FLOAT, &spm()->tunable_green, "min = 0.0 max = 1.0 step = 0.01 group = 'Color' ");
			//TwAddVarRW(graphics_bar, "Blue", TW_TYPE_FLOAT, &spm()->tunable_blue, "min = 0.0 max = 1.0 step = 0.01 group = 'Color' ");
			
			TwEnumVal draw_polygon_mode[] = {
				{spm()->OUTLINE_DRAW, "Outline"}, {spm()->FILL_DRAW, "Fill"}, {spm()->TEXTURE_DRAW, "Texture"}, { spm()->DISCRETE_SCALE_DRAWE, "Discrete"}
			};
			TwType tw_draw_polygon_mode = TwDefineEnum("DrawPolygonMode", draw_polygon_mode, 4);
			TwAddVarRW(graphics_bar, "Draw Switch", TW_TYPE_BOOL8, &spm()->show_polygons(), "group = 'Polygon' ");
			TwAddVarRW(graphics_bar, "Draw Mode", tw_draw_polygon_mode, &spm()->get_polygon_draw_type(), "group = 'Polygon' ");
			//TwAddVarRW(graphics_bar, "In/Active", TW_TYPE_BOOL8, &spm()->show_inactive(), "group = 'Polygon' ");
		
			TwAddVarRW(graphics_bar, "Triangulation", TW_TYPE_BOOL8, &spm()->show_triangulation(), "group = 'Geometry' ");
			TwAddVarRW(graphics_bar, "Voronoi Cell", TW_TYPE_BOOL8, &spm()->show_voronoi_cell(), "group = 'Geometry' ");
			TwAddVarRW(graphics_bar, "Smoothed Voronoi", TW_TYPE_BOOL8, &spm()->show_smoothed_voronoi_cell(), "group = 'Geometry' ");
			TwAddVarRW(graphics_bar, "Midpoint Cell", TW_TYPE_BOOL8, &spm()->show_midpoint_cell(), "group = 'Geometry' ");
			//TwAddVarRW(graphics_bar, "Vertices", TW_TYPE_BOOL8, &spm()->show_vertices(), "group = 'Geometry' ");

			//TwAddVarRW(graphics_bar, "Highlight", TW_TYPE_INT32, &spm()->highlighted_group_id(), "group = 'Debug' ");
			//TwAddVarRW(graphics_bar, "Local frame", TW_TYPE_BOOL8, &spm()->show_local_frame(), "group = 'Debug' ");
			TwAddVarRW(graphics_bar, "Vicinity", TW_TYPE_BOOL8, &spm()->show_vicinity(), "group = 'Debug' ");
			TwAddVarRW(graphics_bar, "Swapped", TW_TYPE_BOOL8, &spm()->show_swapped_polygons(), "group = 'Debug'");
			//TwAddVarRW(graphics_bar, "Hole Tri", TW_TYPE_BOOL8, &spm()->show_hole_triangles(), "group = 'Geometry' ");
			TwAddVarRW(graphics_bar, "Holes", TW_TYPE_BOOL8, &spm()->show_holes(), "group = 'Geometry' ");
			
			if (spm()->contain_multi_packing())
			{
				TwAddVarRW(graphics_bar, "Seg-mesh", TW_TYPE_BOOL8, &spm()->show_multi_submeshes(), "group = 'Multi-mesh' ");
				TwAddVarRW(graphics_bar, "Seg-tiles", TW_TYPE_BOOL8, &spm()->show_multi_tiles(), "group = 'Multi-mesh' ");
			}

			TwBar* function_bar = TwNewBar("Functions");
			TwDefine("Functions position='16 320' size='200 600' alpha=200");
			TwAddVarRW(function_bar, "Lloyd Iter", TW_TYPE_INT32, &lloyd_iter_times, "min=1 group = 'Optimization' ");
			TwAddButton(function_bar, "Lloyd", tw_lloyd, NULL, "key=l group = 'Optimization' ");
			TwAddVarRW(function_bar, "Pack Iter", TW_TYPE_INT32, &pack_iter_times, "min=1 group = 'Optimization'");
			TwAddButton(function_bar, "Pack", tw_pack, NULL, "key=p group = 'Optimization' ");
			TwAddVarRW(function_bar, "Synchronize", TW_TYPE_BOOL8, &spm()->sync_optimization(), "group = 'Optimization' ");
			TwAddVarRW(function_bar, "Use Voronoi Cell", TW_TYPE_BOOL8, &spm()->use_voronoi_cell(), "group = 'Optimization' ");
			TwAddVarRW(function_bar, "Align Constraint", TW_TYPE_BOOL8, &spm()->toggle_align_constraint(), "group = 'Optimization' ");
			TwAddVarRW(function_bar, "Allow Rotation", TW_TYPE_BOOL8, &spm()->toggle_allow_rotation(), "group = 'Optimization'");
			
			//TwAddButton(function_bar, "vPack", tw_vpack, NULL, "key=v group = 'Optimization'");
			//TwAddButton(function_bar, "iDT", tw_idt_update, NULL, "key=i group = 'Geometry' ");
			TwAddButton(function_bar, "Detect Holes", tw_detect_holes, NULL, "key=d group = 'Hole' ");
			TwAddButton(function_bar, "Fill holes", tw_fill, NULL, "key=f group = 'Hole' ");
			TwAddButton(function_bar, "Activate", tw_activate, NULL, "group = 'Hole'");
			TwAddVarCB(function_bar, "Hole Size", TW_TYPE_DOUBLE, tw_hole_size_set_callback,
						tw_hole_size_get_callback, NULL, "min=0.0 step=0.001 group = 'Hole' ");
			TwAddVarCB(function_bar, "Front Edge", TW_TYPE_DOUBLE, tw_front_len_set_callback,
						tw_front_len_get_callback, NULL, "min=0.0 step=0.001 group = 'Hole' ");
			TwAddVarRW(function_bar, "Only Scale", TW_TYPE_BOOL8, &spm()->toggle_only_scale(), "group = 'Hole'");
			TwAddButton(function_bar, "Auto Scale", tw_auto_scale, NULL, "group = 'Hole'");
			TwAddVarRW(function_bar, "weight", TW_TYPE_DOUBLE, &spm()->get_match_weight(), "min=0.0 group = 'Replace' ");
			TwAddVarRW(function_bar, "factor", TW_TYPE_DOUBLE, &spm()->replace_shrink_factor(), "min=0.05 group = 'Replace' ");
			//TwAddButton(function_bar, "replace", tw_replace, NULL, "key=r group = 'Replace' ");
			TwAddButton(function_bar, "swap", tw_swap, NULL, "key=w group = 'Replace' ");
			TwAddButton(function_bar, "sel swap", tw_selective_swap, NULL, "group = 'Replace' ");
			TwAddButton(function_bar, "ext-replace", tw_con_replace, NULL, "key = x group = 'Replace' ");

			TwAddVarRW(function_bar, "k", TW_TYPE_DOUBLE, &enlarge_factor, "min=0.01 group = 'Adjustment'");
			TwAddButton(function_bar, "adjust", tw_adjust, NULL, "key = d group = 'Adjustment'");

			TwAddVarRW(function_bar, "epsilon", TW_TYPE_DOUBLE, &spm()->get_epsilon(), "min=0.0 max=0.999999999 group = 'Optimization' ");
			if (spm()->contain_multi_packing())
				TwAddButton(function_bar, "pack next", tw_pack_next, NULL, "key=n group = 'multimesh' ");
			TwAddButton(function_bar, "Report", tw_report, NULL, "group = 'Discretize' ");
			
			TwAddButton(function_bar, "save tiles", tw_save_tiles, NULL, "key=S group='File'");
			TwAddButton(function_bar, "save mat", tw_save_mat, NULL, "key=M group='File'");
			TwAddButton(function_bar, "save eps", tw_print_pdf, NULL, "group='File'");
			TwAddButton(function_bar, "save scale", tw_save_scale, NULL, "group='File'");
			TwAddButton(function_bar, "save CAT", tw_save_CAT, NULL, "group='File'");
			toggle_skybox_CB() ;
            glut_viewer_add_toggle('b', glut_viewer_is_enabled_ptr(GLUT_VIEWER_BACKGROUND), "switch Color/BW") ;
            glut_viewer_add_toggle('T', &viewer_properties_->visible(), "viewer properties") ;
        }

    private:
		std::string prj_config_file;
		int lloyd_iter_times;
		int pack_iter_times;
		//debug
 		double enlarge_factor;
    } ;


Geex::SPMApp* spm_app() { return static_cast<Geex::SPMApp*>(Geex::GeexApp::instance()) ; }

void TW_CALL tw_adjust(void *clientData) { spm_app()->adjust(); }
void TW_CALL tw_auto_scale(void* cd) { spm_app()->auto_adjust(); }
void TW_CALL tw_lloyd(void *clientData) {spm_app()->lloyd();}
void update() { spm_app()->post_update(); }
void TW_CALL tw_pack(void *clientData) { spm_app()->pack(); }

void TW_CALL tw_idt_update(void *clientData) { spm_app()->idt_update(); }
void TW_CALL tw_replace(void *clientData) { spm_app()->replace(); }
void TW_CALL tw_swap(void* cd) { spm_app()->swap(); }
void TW_CALL tw_selective_swap(void* cd) { spm_app()->selective_swap(); }
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
void TW_CALL tw_activate(void* clientData)
{
	spm_app()->activate_all();
}
void TW_CALL tw_save_triangulation(void* clientData)
{
	spm_app()->save_triangulation();
}
//void TW_CALL tw_save_cur_area( void *clientData ) { spm_app()->save_cur_area(); }
void TW_CALL tw_fill( void *clientData ) { spm_app()->fill_holes(); }
void TW_CALL tw_con_replace(void * clientData) { spm_app()->con_replace(); }
void TW_CALL tw_pack_next(void *clientData) { spm_app()->pack_next(); }
void TW_CALL tw_report(void* clientData) { spm_app()->report(); }
//void TW_CALL tw_discretize_tiles(void *clientData) { spm_app()->discretize_tiles(); }
void TW_CALL tw_save_tiles(void* clientData) { spm_app()->save_tiles(); }
void TW_CALL tw_save_mat(void *clientData) { spm_app()->save_mat(); }
void TW_CALL tw_save_scale(void *clientData) { spm_app()->save_scale(); }
void TW_CALL tw_save_CAT(void* clientData) { spm_app()->save_CAT(); }
void TW_CALL tw_print_pdf(void *clientData)
{
	std::cout << "Print to file MyFile.pdf ...";
	FILE *fp = fopen("MyFile.pdf", "wb");
	GLint buffsize = 0, state = GL2PS_OVERFLOW;
	GLint viewport[4];

	glGetIntegerv(GL_VIEWPORT, viewport);

	while( state == GL2PS_OVERFLOW ){ 
		buffsize += 1024*1024;//GL2PS_BSP_SORT 
		gl2psBeginPage ( "MyTitle", "MySoftware", viewport,
			GL2PS_PDF, GL2PS_BSP_SORT, GL2PS_SILENT |
			GL2PS_SIMPLE_LINE_OFFSET | GL2PS_NO_BLENDING |
			GL2PS_OCCLUSION_CULL | GL2PS_BEST_ROOT | 
			GL2PS_NO_PS3_SHADING | GL2PS_COMPRESS,
			GL_RGBA, 0, NULL, 2, 2, 2, buffsize,
			fp, "MyFile" );
		glut_viewer_redraw(); 
		state = gl2psEndPage();
	}

	fclose(fp);
	std::cout <<" End!" <<std::endl;
}

}

int main(int argc, char** argv) 
{
	try
	{
		Geex::initialize();
		Geex::SPMApp app(argc, argv) ;
		app.main_loop() ;
		Geex::terminate() ;
	}
	catch (std::exception& e)
	{
		std::cout<<e.what()<<std::endl;

	}
	//system("pause");
	return 0;
}
