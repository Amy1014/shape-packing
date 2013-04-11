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
#endif

#include "project_io.h"
#include "packing_object.h"
#include "meshes.h"
#include "polygon_diagram.h"
#include "Containment.h"
#include "knitro.h"
#include "polygon_matcher.h"

namespace Geex
{
	class Packer
	{

		typedef enum {SUCCESS, FAILED} Optimization_res;
		typedef enum {MORE_ENLARGEMENT, NO_MORE_ENLARGEMENT} Lloyd_res;

		typedef RestrictedPolygonVoronoiDiagram::Vertex_handle Vertex_handle;
		typedef RestrictedPolygonVoronoiDiagram::Facet_iterator Facet_iterator;
		typedef RestrictedPolygonVoronoiDiagram::Facet_handle Facet_handle;
		typedef RestrictedPolygonVoronoiDiagram::Halfedge_handle Halfedge_handle;
		typedef RestrictedPolygonVoronoiDiagram::Halfedge_iterator Halfedge_iterator;

		struct Local_frame;
		struct Parameter;
		
	public:

		typedef std::vector<Segment_3> Hole; // representing a hole, consisting of non-border edges
		
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

		static Packer* instance() { return instance_; }

		/** optimization **/
		void lloyd(void (*post_action)() = NULL, bool enlarge = false);

		/** miscellaneous **/
		void get_bbox(real& x_min, real& y_min, real& z_min, real& x_max, real& y_max, real& z_max);

		// hole detection
		void detect_holes();
		// hole filling
		void fill_holes();

		// replace
		void replace();
		void remove_polygons();
		void replace_one_polygon(unsigned int id, Hole& region); // public for debug

		// driver
		void pack(void (*post_action)() = NULL); 
		// debug
		void rpack(void (*post_action)() = NULL);
		void update_iDT() { rpvd.iDT_update(); compute_clipped_VD();}
		void save_curvature_and_area(); 
		void enlarge_one_polygon(unsigned int id, double f, double theta, double tx, double ty);

	private:

		/** initialization **/
		void random_init_tiles(unsigned int nb_init_polygons);	// put polygons according to curvature (if any) and in a uniform way
		void put_polygons_on_facet(unsigned int facet_id, int n, unsigned int& pgn_lib_index);
		/** geometry **/
		void generate_RDT();

		void compute_clipped_VD();

		vec3 approx_normal(unsigned int facet_idx);
		
		Local_frame compute_local_frame(const Packing_object& tile);

		/** optimization **/
		// one Lloyd iteration
		Lloyd_res one_lloyd(bool enlarge, std::vector<Parameter>& solutions, std::vector<Local_frame>& lfs);

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
		void constraint_transformation(vector<Parameter>& parameters, vector<Local_frame>& lfs, vector<double>& min_factors);

		// mapping from curvature to transformation
		void curv_constrained_transform(Parameter& para, int fid, unsigned int pgn_id);

		/** replace and hole filling**/
		// remove one polygon and leave a hole
		void remove_one_polygon(unsigned int id, Hole& hole);
		// fill one hole
		void fill_one_hole(Hole& hl, Packing_object& filler);

	private:

		static Packer *instance_;
		
		ProjectIO pio;
		std::vector<Packing_object> pack_objects;
		std::vector<Ex_polygon_2> pgn_lib;
		TriMesh mesh;
		double mesh_area;

		/** holes **/
		std::vector<Hole> holes;
		double frontier_edge_size;
		double hole_face_size;
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
		struct Local_frame // represent a local frame for which a transformation is effective
		{
			Vector_3 u;
			Vector_3 v;
			Vector_3 w;
			Point_3 o;
			inline Point_2 to_uv(const Point_3& p) const
			{
				Vector_3 op(o, p);
				return Point_2(op*u, op*v);
			}
			inline Segment_2 to_uv(const Segment_3& s) const 
			{
				return Segment_2(to_uv(s.source()), to_uv(s.target()));
			}
			inline Point_3 to_xy(const Point_2& p) const
			{
				return o + ( p.x()*u + p.y()*v );
			}
		};

		struct Parameter // represent a 2D transformation parameter
		{
			double k;
			double theta;
			double tx;
			double ty;
			Parameter(void): k(1.0), theta(0.0), tx(0.0), ty(0.0) {}
			Parameter(double _k, double _theta, double _tx, double _ty) : k(_k), theta(_theta), tx(_tx), ty(_ty) {}
			Parameter& operator*=(double f)
			{
				//k = std::max(f*k, 1.0); 
				//k = 1.0 + (k - 1.0)*f;
				theta *= f;	tx *= f; ty *= f;
				return *this;
			}
			// compute the remaining transformation after applying transformation p
			Parameter& operator-=(const Parameter& p)
			{
				k /= p.k;
				theta -= p.theta;
				tx -= p.tx;
				ty -= p.ty;
				return *this;
			}
			inline Parameter operator*(double f)
			{
				//return Packer::Parameter(1.0 + (k - 1.0)*f, f*theta, f*tx, f*ty);
				return Packer::Parameter(k, f*theta, f*tx, f*ty);
			}
			inline bool is_identity() const 
			{
				return /*k == 1.0 && */theta == 0.0 && tx == 0.0 && ty == 0.0;
			}
		};

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

		//std::vector<Local_frame> local_frames;


	};


}

#endif