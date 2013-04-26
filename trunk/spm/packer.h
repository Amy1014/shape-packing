#ifndef _PACKER_H_
#define _PACKER_H_

#include <vector>
#include <queue>
#include <string>
#include <cfloat>
#include <cmath>
#include <set>
#include <numeric>
#include <sstream>


#include <CGAL/Timer.h>

#include <glut_viewer/glut_viewer.h>

#ifdef _CILK_
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>
#include <cilk/reducer_opand.h>
#endif

#include "project_io.h"
#include "packing_object.h"
#include "meshes.h"
#include "polygon_diagram.h"
#include "Containment.h"
#include "knitro.h"
#include "polygon_matcher.h"
#include "utilities.h"

namespace Geex
{
	class Packer
	{

		typedef enum {SUCCESS, FAILED} Optimization_res;
		typedef enum {MORE_ENLARGEMENT, NO_MORE_ENLARGEMENT, BARRIER_PARTIAL_ACHIEVED, BARRIER_ALL_ACHIEVED} Lloyd_res;

		typedef RestrictedPolygonVoronoiDiagram::Vertex_handle Vertex_handle;
		typedef RestrictedPolygonVoronoiDiagram::Facet_iterator Facet_iterator;
		typedef RestrictedPolygonVoronoiDiagram::Facet_handle Facet_handle;
		typedef RestrictedPolygonVoronoiDiagram::Halfedge_handle Halfedge_handle;
		typedef RestrictedPolygonVoronoiDiagram::Halfedge_iterator Halfedge_iterator;

		class Discrete_barries;

	public:

		typedef std::vector<std::pair<Vertex_handle, Vertex_handle>> Hole; // representing a hole, consisting of non-border edges
		
		Packer();

		/** initialization **/
		void load_project(const std::string& prj_config_file);

		void initialize(); // top-level initialization function. distribute polygons in some way	

		/** access functions **/
		const TriMesh& mesh_domain() const { return mesh; }
		const vector<Packing_object>& get_tiles() const { return pack_objects; }
		RestrictedPolygonVoronoiDiagram& get_rpvd() const { return rpvd; }
		double hole_size() const { return hole_face_size; }
		void hole_size(double sz) { hole_face_size = sz; }
		double front_edge_len() const { return frontier_edge_size; }
		void front_edge_len(double len) { frontier_edge_size = len; }
		std::vector<Hole>& get_holes() const { return holes; }
		double& get_epsilon() { return epsilon; }
		double& get_match_weight() { return match_weight; }
		ProjectIO& get_project_ioer() { return pio; }

		const std::vector<std::vector<Packing_object>>& get_multigroup_tiles() const { return res_pack_objects; }
		const std::vector<TriMesh>& get_multigroup_submeshes() const { return mesh_segments; }

		static Packer* instance() { return instance_; }

		void interface_lloyd(void (*post_action)() = NULL);

		/** miscellaneous **/
		void get_bbox(real& x_min, real& y_min, real& z_min, real& x_max, real& y_max, real& z_max);

		void print_area_coverage();

		void save_tiles();

		// hole detection
		void detect_holes();
		// hole filling
		void fill_holes();

		// replace
		void replace();
		void remove_polygons();// deleted for simplicity, see revision 67
		void ex_replace(); // extended replace, deleted for simplicity, see revision 67
		void con_replace(); 

		// discretize
		void discretize_tiles();
		double& max_scale_factor() { return max_scale; }
		double& min_scale_factor() { return min_scale; }
		int& discrete_levels() { return levels; }
		bool& discrete_scale() { return discrete_scaling; }
		void split_large_tiles();
		// driver
		void pack(void (*post_action)() = NULL); 
		// debug
		void update_iDT() { rpvd.iDT_update(); compute_clipped_VD();}
		void save_curvature_and_area(); // deleted for simplicity, see revision 67
		void enlarge_one_polygon(unsigned int id, double f, double theta, double tx, double ty); // deleted for simplicity, see revision 67
		CDT& get_cdt() {  return cdt; }

		void report();

	protected:
		/** multi meshes **/
		void pack_next_submesh();
		void save_sub_result();
		
	private:

		/** initialization **/
		void random_init_tiles(unsigned int nb_init_polygons);	// put polygons according to curvature (if any) and in a uniform way
		void put_polygons_on_facet(unsigned int facet_id, int n, unsigned int& pgn_lib_index);
		/** geometry **/
		void generate_RDT();

		void compute_clipped_VD(bool approx = false);

		void compute_clipped_VD(std::vector<bool> use_approx); // mixed computation of clipped voronoi diagram

		vec3 approx_normal(unsigned int facet_idx);
		
		//Local_frame compute_local_frame(const Packing_object& tile);

		/** optimization **/
		
		/** optimization **/
		void lloyd(void (*post_action)() = NULL, bool enlarge = false);

		bool discrete_lloyd(void (*post_action)() = NULL, bool enlarge = false);

		void split(std::set<unsigned int>& indices);

		// one Lloyd iteration
		Lloyd_res one_lloyd(bool enlarge, std::vector<Parameter>& solutions, std::vector<Local_frame>& lfs);

		bool discrete_one_lloyd(bool enlarge, std::vector<Parameter>& solutions, std::vector<Local_frame>& lfs, double barrier_factor, std::vector<bool>& approx_vd);

		static int callback(const int evalRequestCode, const int n, const int m, const int nnzJ, const int nnzH,
							const double * const x,	const double * const lambda, double * const obj, double * const c,
							double * const objGrad,	double * const jac,	double * const hessian,	double * const hessVector, void * userParams);
		// optimize by Knitro
		int KTR_optimize(double* io_k, double* io_theta, double* io_t1, double* io_t2, unsigned int idx);

		// get the transformation parameter for one polygon in a local frame
		Optimization_res optimize_one_polygon(unsigned int id, Local_frame& lf, Parameter& solution);

		// transform polygons inside local frames
		void transform_one_polygon(unsigned int id, Local_frame& lf, Parameter& param);

		// constraint transformation to avoid triangle flip
		void constraint_transformation(vector<Parameter>& parameters, vector<Local_frame>& lfs, bool constrain_scale);

		// mapping from curvature to transformation
		void curv_constrained_transform(Parameter& para, int fid, unsigned int pgn_id);

		/** replace and hole filling**/
		// remove one polygon and leave a hole
		void remove_one_polygon(unsigned int id, Hole& hole, std::set<Facet_handle>& removed_facets); // deleted for simplicity, see revision 67
		// fill one hole
		void fill_one_hole(Hole& hl, Packing_object& filler); 

		bool replace_one_polygon(unsigned int id, Hole& region, std::set<Facet_handle>& removed_facets); // deleted for simplicity, see revision 67

	private:

		static Packer *instance_;
		
		ProjectIO pio;
		std::vector<Packing_object> pack_objects;
		unsigned int samp_nb; // polygon sample number
		std::vector<Ex_polygon_2> pgn_lib;
		TriMesh mesh;
		double mesh_area;

		/** for segmented mesh**/
		//std::vector<std::vector<Packing_object>> res_pack_objects; // results on different sub-meshes
		std::vector<TriMesh> mesh_segments; // multiple sub-meshes from segmentation
		unsigned int sub_pack_id;
		std::vector<std::vector<Ex_polygon_2>> pgn_lib_sets; // we have multiple polygon libraries

		/** holes **/
		std::vector<Hole> holes;
		double frontier_edge_size;
		double hole_face_size;

		/** discretize **/
		double max_scale;
		double min_scale;
		int levels;
		bool discrete_scaling;

		/** optimization **/
#ifdef _CILK_
		std::vector<Containment> containments;
#else
		Containment containment;
#endif
		static double PI;
		double rot_lower_bd;
		double rot_upper_bd;
		double epsilon;

		/** control variable **/
		bool stop_update_DT;
		double match_weight;

		/** helper classes **/
	private:

		class Discrete_barriers
		{
		public:
			Discrete_barriers() {}
			Discrete_barriers(double min_factor, double max_factor, int levels) 
			{
				set(min_factor, max_factor, levels);
			}
			Discrete_barriers(const Discrete_barriers& barriers0, const Discrete_barriers& barriers1)
			{
				set(barriers0, barriers1);
			}
			void set(double min_factor, double max_factor, int levels)
			{
				grades.resize(levels+1);
				for (int i = 0; i < levels+1; i++)
					grades[i] = ( (levels-i)*min_factor + i*max_factor ) / levels;
				current_barrier = grades.begin();
			}
			void set(const Discrete_barriers& barriers0, const Discrete_barriers& barriers1)
			{
				std::merge(barriers0.grades.begin(), barriers0.grades.end(), barriers1.grades.begin(), barriers1.grades.end(), 
					std::back_insert_iterator<std::vector<double>>(this->grades));
				std::vector<double>::iterator it = std::unique(this->grades.begin(), this->grades.end());
				this->grades.resize(std::distance(this->grades.begin(), it));
				current_barrier = this->grades.begin();
			}
			bool beyond_range() { return current_barrier == grades.end(); }
			bool get_current_barrier(double& res)
			{
				if (beyond_range())
					return false;
				res = *current_barrier;
				return true;
			}
			void go_to_next_barrier()
			{
				if (!beyond_range())
					++current_barrier;
			}
			bool get_prev_barrier(double& res)
			{
				if (current_barrier == grades.begin())
					return false;
				res = *(current_barrier-1);
				return true;
			}
			void set_current_barrier(double lower_val)
			{
				current_barrier = std::lower_bound(grades.begin(), grades.end(), lower_val);
			}
			inline double get_min() const { return grades.front(); }
			inline double get_max() const { return grades.back(); }
			inline unsigned int nb_levels() const { return grades.size(); }
		private:
			std::vector<double> grades;
			std::vector<double>::iterator current_barrier;
		} phy_disc_barr, opt_disc_barr, disc_barr;

		// apply a 2D transformation to a 3D polygon inside the plane where the polygon is embedded
		struct Apply_transformation
		{
			typedef RestrictedPolygonVoronoiDiagram::Vertex_handle Vertex_handle;
			const Parameter& parameter;
			const Local_frame& local_frame;

			Transformation_2 t;

			Apply_transformation(const Parameter& _parameter, const Local_frame& _local_frame)
				: parameter(_parameter), local_frame(_local_frame)
			{
				double ksin = parameter.k*std::sin(parameter.theta), kcos = parameter.k*std::cos(parameter.theta);
				t = Transformation_2(kcos, -ksin, parameter.tx, ksin, kcos, parameter.ty);
			}

			inline Point_3 operator()(const Point_3& p)
			{
				Point_2 tp = t(local_frame.to_uv(p));
				return local_frame.to_xy(tp);
			}

			inline RestrictedPolygonVoronoiDiagram::Vertex_handle operator()(Vertex_handle v)
			{
				Point_3 p = operator()(v->point());
				v->point() = p;
				return v;
			}
		};



	public: // for debug
		RestrictedPolygonVoronoiDiagram rpvd;
		CDT cdt;
		//std::vector<Local_frame> local_frames;
		//Packing_object backup;
		std::vector<std::vector<Packing_object>> res_pack_objects;
	};


}

#endif