
#include <Geex/graphics/geexapp.h>
#include <Geex/basics/file_system.h>
#include <glut_viewer/tweak_bar.h>
#include <CGAL/Timer.h>
#include <cstdlib>
#include <ctime>
#include "SPM.h"
#include "meshes.h"

void prompt_and_exit(const char* prompt)
{
	std::cerr<<prompt<<std::endl;
	system("pause");
	exit(0);
}

namespace Geex {

	void TW_CALL tw_reset(void *clientData);
	void TW_CALL tw_adjust(void *clientData);
	void TW_CALL tw_lloyd(void *clientData);
	void TW_CALL tw_pack(void *clientData);
	void TW_CALL tw_rpack(void*);
	void TW_CALL tw_replace(void *clientData);
	void TW_CALL tw_detect_holes(void*);
	void TW_CALL tw_save(void*);
	void TW_CALL tw_fill(void*);
	void TW_CALL tw_merge(void*);
	void TW_CALL tw_feature_set_callback(const void*, void *);
	void TW_CALL tw_feature_get_callback(void*, void*);
	void TW_CALL tw_curvature_get_callback(void*, void*);
	void TW_CALL tw_curvature_set_callback(const void*, void*);
	void TW_CALL tw_correct_normals(void*);
	void TW_CALL tw_smooth_normals(void*);
	void TW_CALL tw_flip_normals(void*);
	void TW_CALL tw_dump(void*);
	void TW_CALL tw_restore(void*);
	void TW_CALL tw_cluster_detect_holes(void*);
	void TW_CALL tw_affine_fill(void*);
	
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
			texture = false;
			//restore_ = false;
			get_arg("texture", texture);			
			get_arg("restore", restore_file);
			get_arg("dump", dump_dir);
			get_arg("btex", bg_tex_filename);
			get_arg("bpgn", bg_pgn_filename);
			//if (!bg_tex_filename.empty()&&!bg_pgn_filename.empty())
			//{
				mc_pts_filename = get_file_arg("pts");
			//	if (mc_pts_filename.empty())
			//		std::cout<<"Caution: no initial points input. backgroud feature disabled!\n";
			//}
			//if (!restore_file.empty())
			//{
			//	std::cout<<"Use restore file: "<<restore_file<<std::endl;
			//}
			if (texture)
			{
				std::vector<std::string> filenames;
				std::string dir_name;
				get_arg("d", dir_name);
				std::string mdir_name;//for multi-size polygons
				get_arg("md", mdir_name);
				if (!dir_name.empty())
					FileSystem::get_files(dir_name, filenames);
				else if (!mdir_name.empty())
					FileSystem::get_files(mdir_name, filenames);
				else 
					prompt_and_exit("No existing directory!");
				if (filenames.size() == 0)
					prompt_and_exit("No directory input!");
				if (!dir_name.empty())
					for (unsigned int i = 0; i < filenames.size(); i++)
					{
						if (FileSystem::extension(filenames[i]) == "tply")
							pgn_filenames.push_back(filenames[i]);
						else if (FileSystem::extension(filenames[i]) == "jpg")
							tex_filenames.push_back(filenames[i]);
					}
				else if (!mdir_name.empty())
					for (unsigned int i = 0; i < filenames.size(); i++)
					{
						if (FileSystem::extension(filenames[i]) == "tply" && filenames[i].find("small")!=string::npos)
							pgn_filenames.push_back(filenames[i]);
						else if (FileSystem::extension(filenames[i]) == "tply" && filenames[i].find("small")==string::npos)
							bg_pgn_fileset.push_back(filenames[i]);
						else if (FileSystem::extension(filenames[i]) == "jpg" && filenames[i].find("small")!=string::npos)
							tex_filenames.push_back(filenames[i]);
						else if (FileSystem::extension(filenames[i]) == "jpg" && filenames[i].find("small")==string::npos)
							bg_tex_fileset.push_back(filenames[i]);
					}
				if ( pgn_filenames.size() == 0)
					prompt_and_exit("No .tply file input for textured polygons!");
				if (tex_filenames.size() == 0)
					prompt_and_exit("No .jpg file input for texture!");
			}
			else if(pgn_filename_.length() > 0)
			{
				if(!Geex::FileSystem::is_file(pgn_filename_)) 
				{
					prompt_and_exit("File for polygon input invalid.");
				}
			}
			else
				prompt_and_exit("Absent or invalid file input!");

			cur_filename_ = get_file_arg("cur");
            nb_iter_ = 1;
            get_arg("nb_iter", nb_iter_) ;
			totalIterNb = 0;
			scalingFactor = 1.0;
       }

        SPM* spm() { return static_cast<SPM*>(scene()) ; }

        void init_scene() {
            scene_ = new SPM;
			if (boundary_filename_.length()>0)
				spm()->loadMeshDomain(boundary_filename_);
			else
				prompt_and_exit("Failure to load mesh domain");
			if (cur_filename_.length()>0)
			{
				spm()->getMeshDomain().load_density(cur_filename_);
				spm()->setCurvature();
			}
			if (texture)
			{
				if (pgn_filenames.size() > 0)
				{
					if (!bg_pgn_filename.empty() && !bg_tex_filename.empty() && !mc_pts_filename.empty())
					{
						spm()->setBackground();
						//spm()->loadPolygonsWithBackgroud(pgn_filenames, restore_file, bg_pgn_filename, mc_pts_filename);
						spm()->read_backgroud_texture(bg_tex_filename);					
					}
					else if (bg_pgn_fileset.size()>0 && bg_tex_fileset.size()>0 && !mc_pts_filename.empty())
					{
						spm()->setBackground();
						//spm()->loadPolygonsWithBackgroud(pgn_filenames, restore_file, bg_pgn_fileset, mc_pts_filename);
						spm()->read_backgroud_texture(bg_tex_fileset);
					}
					else if (!restore_file.empty()){}
						//spm()->loadPolygonsWithTexture(pgn_filenames, restore_file);
					else {}
						//spm()->loadPolygonsWithTexture(pgn_filenames);
					spm()->read_texture(tex_filenames);
				}
				else
					prompt_and_exit("no polygon files input!");
			}
			else if (pgn_filename_.length()>0)
			{
				if (pgn_filename_.substr(pgn_filename_.length()-4, pgn_filename_.length()) == "poly")
				{
					if (restore_file.empty())
						spm()->loadPolygons(pgn_filename_);
					else
					{
						std::cout<<"Use restore file: "<<restore_file<<std::endl;
						spm()->loadPolygons(pgn_filename_, restore_file);
					}
				}
				else
					prompt_and_exit("bad polygon file format!");
			}
			else
				prompt_and_exit("Failure to load polygons");
			if (!dump_dir.empty())
				spm()->immResDir = dump_dir;
			//////////////////////////////////////////////////////////////////////////
			//pack();
			//dump();
        }

        void Lloyd() {
			int i;
			CGAL::Timer timer;
			timer.start();
			double finalScaling, maxf, minf;
			spm()->getFactorStatistic(&maxf, &minf, &finalScaling);
			for (i = 0; i < nb_iter_ && finalScaling <= 1.0; i++)
			{
				std::cerr<<":::: iteration "<<totalIterNb<<" ::::"<<std::endl;
				Doc::LloydStatus s = spm()->lloyd();
				if (s != Doc::LLOYD_SUCCESS)
					break;
				spm()->getFactorStatistic(&maxf, &minf, &finalScaling);
				glut_viewer_redraw();
				totalIterNb++;
			}
			timer.stop();
			std::cout<<"Total number of polygons: "<<spm()->getPolygons().size()<<std::endl;
			std::cout<<i<<" Lloyd iteration(s) finished, time elapsed : "<<timer.time()<<" seconds"<<std::endl;
			std::cout<<"================== Enlargement relative to input: "<<finalScaling<<" ==================\n";
			glut_viewer_redraw();
        }
		void pack()
		{
			TwDefine(" Graphics/Feature visible=false ");//the feature criterion cannot be changed any more
			CGAL::Timer timer;
			timer.start();
			spm()->pack();
			timer.stop();
			std::cout<<"Time elapsed : "<<timer.time()<<" seconds"<<std::endl;
		}
		void rpack()
		{
			TwDefine(" Graphics/Feature visible=false ");//the feature criterion cannot be changed any more
			CGAL::Timer timer;
			timer.start();
			spm()->rpack();
			timer.stop();
			std::cout<<"Time elapsed : "<<timer.time()<<" seconds"<<std::endl;
		}
		void fill()
		{
			spm()->fill(0.3);

		}
		void affine_fill()
		{
			//spm()->show_holes() = false;
			nbFilled = spm()->affineFill();
			glut_viewer_redraw();
			std::cout<<nbFilled<<" holes are filled.\n";
		}
		void replace()
		{
			std::cout<<" replace polygons... "<<std::endl;
			//spm()->replace();
			//spm()->affineReplace();
			//spm()->bgreplace();
			spm()->noDupReplace();
			std::cout<< "replace finished. "<<std::endl;
		}

		void adjust()
		{
			spm()->adjustPolygons(scalingFactor);
			glut_viewer_redraw();
		}
		void detect_holes()
		{
			spm()->detectHoles();
		}
		void cluster_detect_holes()
		{
			spm()->clusterDetectHoles();
		}
		void merge()
		{
			//spm()->mergeHoles();
			// debug
			spm()->postOptimization(nbFilled);

		}
        void reset() {
            spm()->clear() ;
            spm()->loadPolygons(pgn_filename_);
			spm()->loadMeshDomain(boundary_filename_);
			glut_viewer_redraw();
            //else
            //    spm()->insert_vertices(nb_points_, vertex_type_) ;
        }
		void save()
		{
			spm()->saveCurrentPolygons("D:/research/packing/spm/data/polygons/hexagon_config.poly");
		}

		void feature_set_callback(double setValue)
		{
			featureThreshold = setValue;
			spm()->setFeatureAngle(featureThreshold);
			//std::cout<<"The feature criterion is changed to "<<featureThreshold<<std::endl;
		}
		void curvature_set_callback(double setValue)
		{
			spm()->setCurvaturePercent(setValue);
		}
		void feature_get_callback(double* getValue)
		{
			featureThreshold = spm()->getFeatureAngle();
			*getValue = featureThreshold;
		}
		void curvature_get_callback(double *getValue)
		{
			*getValue = spm()->getCurvaturePercet();
		}
		//void correct_normals()
		//{
		//	spm()->correctNormals();
		//	glut_viewer_redraw();
		//}
		//void smooth_normals()
		//{
		//	spm()->smoothNormals();
		//	glut_viewer_redraw();
		//}
		void flip_normals()
		{
			spm()->flipMeshNormals();
			glut_viewer_redraw();
		}
		void dump()
		{
			if (!dump_dir.empty())
			{
				string dump_file;
				if (dump_dir[dump_dir.length()-1] == '/' || dump_dir[dump_dir.length()-1] == '\\' )
					dump_file = dump_dir+"res.dat";
				else
					dump_file = dump_dir+"/res.dat";
				spm()->dump(dump_file);
				std::cout<<"Dump file: "<<dump_dir<<std::endl;
			}
		}
		void restore()
		{
			if (!restore_file.empty())
				spm()->restore(restore_file);
		}
		//debug
		void fitPlanes()
		{
			spm()->fitPlanes();
		}
        virtual void init_gui() {
            GeexApp::init_gui() ;
			
            TwBar* graphics_bar = TwNewBar("Graphics");
			TwDefine("Graphics position='16 10' size='200 250' alpha=200"); 
			TwAddVarRW(graphics_bar, "Domain Mesh", TW_TYPE_BOOL8, &spm()->show_domain_mesh(), "");
			//TwAddVarRW(graphics_bar, "Mesh Skeleton", TW_TYPE_BOOL8, &spm()->show_mesh_skeleton(), "");
			TwAddVarRW(graphics_bar, "Triangulation", TW_TYPE_BOOL8, &spm()->show_regular_triangulation(), "");
			TwAddVarRW(graphics_bar, "Fill Polygon", TW_TYPE_BOOL8, &spm()->fill_polygon(), "");
			TwAddVarRW(graphics_bar, "Bounding Spheres", TW_TYPE_BOOL8, &spm()->show_bounding_spheres(), "");
			//TwAddVarRW(graphics_bar, "Show Triangulation", TW_TYPE_BOOL8, &spm()->show_regular_triangulation(), "");
			TwAddVarRW(graphics_bar, "Bisectors", TW_TYPE_BOOL8, &spm()->show_bisectors(), "");
			TwAddVarRW(graphics_bar, "Highlight", TW_TYPE_INT32, &spm()->setHighlightedGroup(), "");
			TwAddVarRW(graphics_bar, "Holes", TW_TYPE_BOOL8, &spm()->show_holes(), "");
			TwAddVarRW(graphics_bar, "Boundary", TW_TYPE_BOOL8, &spm()->show_domain_boundary(), "");
			TwAddVarCB(graphics_bar, "Feature Threshold", TW_TYPE_DOUBLE, tw_feature_set_callback, tw_feature_get_callback, NULL, "min=0.0 max=180.0");
			if (spm()->hasCurvature())
			{
				TwAddVarCB(graphics_bar, "Curvature Threshold", TW_TYPE_DOUBLE, tw_curvature_set_callback, tw_curvature_get_callback, NULL, "min=0.0 max=1.0");
				TwAddVarRW(graphics_bar, "Show Curvature", TW_TYPE_BOOL8, &spm()->show_high_curvature(), "");
			}
			TwAddVarRW(graphics_bar, "Show Feature", TW_TYPE_BOOL8, &spm()->show_features(), "");
			TwAddButton(graphics_bar, "Flip Normals", tw_flip_normals, NULL, "key=f");
			if (texture)
				TwAddVarRW(graphics_bar, "Texture", TW_TYPE_BOOL8, &spm()->apply_texture(), "");
			TwBar* function_bar = TwNewBar("Functions");
			TwDefine("Functions position='16 250' size='200 250' alpha=200");
			TwAddVarRW(function_bar, "Iter Number.", TW_TYPE_INT32, &spm()->setPackIterLim(), "min=1");
			TwAddVarRW(function_bar, "Min Scale", TW_TYPE_DOUBLE, &spm()->setMinScalor(), ""/*"min=0.01 max=0.99"*/);
			TwAddVarRW(function_bar, "Max Scale", TW_TYPE_DOUBLE, &spm()->setMaxScalor(), ""/*"min=1.0 max=1.99"*/);
			TwAddVarRW(function_bar, "Area Coverage", TW_TYPE_DOUBLE, &spm()->setAreaCoverage(), "min=0.01 max=1.0");
			TwAddButton(function_bar, "Lloyd", tw_lloyd, NULL, "key=k");
			TwAddVarRW(function_bar, "Scale Factor", TW_TYPE_DOUBLE, &scalingFactor, "min=0.0 step=0.01");
			TwAddButton(function_bar, "Scale", tw_adjust, NULL, "key=a");
			//TwAddButton(function_bar, "Reset", tw_reset, NULL, "key=Z");
			TwAddButton(function_bar, "Replace", tw_replace, NULL, "key=r");
			TwAddButton(function_bar, "Pack", tw_pack, NULL, "key=p");
			//TwAddButton(function_bar, "R Pack", tw_rpack, NULL, "key=P");
			//TwAddButton(function_bar, "Detect", tw_detect_holes, NULL, "key=d");
			//TwAddButton(function_bar, "Fill", tw_fill, NULL, "key=f");
			TwAddButton(function_bar, "Post Optimization", tw_merge, NULL, "key=m");
			TwAddButton(function_bar, "Save", tw_save, NULL, "key=s");
			//TwAddButton(function_bar, "Correct Normals", tw_correct_normals, NULL, "key=c");
			//TwAddButton(function_bar, "Smooth Normals", tw_smooth_normals, NULL, "key=n");
			TwAddButton(function_bar, "Cluster Detect", tw_cluster_detect_holes, NULL, "key=c");
			TwAddButton(function_bar, "Affine Fill", tw_affine_fill, NULL, "key=a");
			if (!dump_dir.empty())
				TwAddButton(function_bar, "Dump", tw_dump, NULL, "key=d");
			//if (restore_)
			//	TwAddButton(function_bar, "Restore", tw_restore, NULL, "");
			toggle_skybox_CB() ;
            glut_viewer_add_toggle('B', glut_viewer_is_enabled_ptr(GLUT_VIEWER_BACKGROUND), "switch Color/BW") ;
            glut_viewer_add_toggle('T', &viewer_properties_->visible(), "viewer properties") ;
        }

    private:
        std::string boundary_filename_ ;
        std::string pgn_filename_ ; // input polygons
		std::vector<std::string> pgn_filenames;//for texture, multiple files input
		std::vector<std::string> tex_filenames;//for texture, multiple texture images
		std::string bg_pgn_filename;//back ground polygon file name
		std::string bg_tex_filename;//back ground texture file name
		std::vector<string> bg_pgn_fileset;
		std::vector<string> bg_tex_fileset;
		std::string mc_pts_filename;
		std::string ptr_filename_;//point distribution file name
		std::string cur_filename_;
		std::string dump_dir;//directory for intermediate result dumping
		std::string restore_file;//directory for restore
        int nb_iter_ ;
		int totalIterNb;
		double featureThreshold;
		GLboolean texture;
		//debug
		int nbFilled;
		double scalingFactor;
    } ;


Geex::SPMApp* spm_app() { return static_cast<Geex::SPMApp*>(Geex::GeexApp::instance()) ; }

void TW_CALL tw_reset(void *clientData) { spm_app()->reset(); }
void TW_CALL tw_adjust(void *clientData) { spm_app()->adjust(); }
void TW_CALL tw_lloyd(void *clientData) {spm_app()->Lloyd();}
void TW_CALL tw_pack(void *clientData) { spm_app()->pack(); }
void TW_CALL tw_replace(void *clientData) { spm_app()->replace(); }
void TW_CALL tw_detect_holes(void *clientData) { spm_app()->detect_holes(); }
void TW_CALL tw_save( void *clientData ) { spm_app()->save(); }
void TW_CALL tw_fill( void *clientData ) { spm_app()->fill(); }
void TW_CALL tw_merge(void *clientData) { spm_app()->merge(); /*spm_app()->fitPlanes();*/ }
void TW_CALL tw_rpack(void *clientData) { spm_app()->rpack(); }
void TW_CALL tw_feature_set_callback(const void* fc, void *clientData) 
{ 
	spm_app()->feature_set_callback(*(double*)fc); 
}
void TW_CALL tw_feature_get_callback(void *value, void *clientData) 
{ 
	spm_app()->feature_get_callback((double*)value); 
}
void TW_CALL tw_curvature_get_callback(void *value, void *clientData)
{
	spm_app()->curvature_get_callback((double*)value);
}
void TW_CALL tw_curvature_set_callback(const void *value, void *clientData)
{
	spm_app()->curvature_set_callback(*(double*)value);
}
//void TW_CALL tw_correct_normals(void* clientData) { spm_app()->correct_normals(); }
//void TW_CALL tw_smooth_normals(void* clientData) { spm_app()->smooth_normals(); }
void TW_CALL tw_flip_normals(void* clientData) { spm_app()->flip_normals(); }
void TW_CALL tw_dump(void *clientData) { spm_app()->dump(); }
void TW_CALL tw_restore(void *clientData) { spm_app()->restore(); }
void TW_CALL tw_cluster_detect_holes(void* clientData) { spm_app()->cluster_detect_holes(); }
void TW_CALL tw_affine_fill(void* clientData) { spm_app()->affine_fill(); }
}

int main(int argc, char** argv) 
{
    Geex::initialize();
    Geex::SPMApp app(argc, argv) ;


	//app.init_scene();
	app.main_loop() ;
    Geex::terminate() ;
}
