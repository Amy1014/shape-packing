#ifndef _PACKER_H_
#define _PACKER_H_

#include <vector>
#include <string>
#include <cfloat>
#include <cmath>
#include <set>
#include <numeric>
#include <sstream>

#include <CGAL/Timer.h>

#ifdef _CILK_
#include <cilk/cilk.h>
#endif

#include "project_io.h"
#include "packing_object.h"
#include "meshes.h"
#include "polygon_diagram.h"
#include "Containment.h"
#include "knitro.h"

namespace Geex
{
	class Packer
	{

		typedef enum {SUCCESS, FAILED} Optimization_res;
		typedef enum {MORE_ENLARGEMENT, NO_MORE_ENLARGEMENT} Lloyd_res;

		typedef RestrictedPolygonVoronoiDiagram::Vertex_handle Vertex_handle;
		
	public:
		
		Packer();

		/** initialization **/
		void load_project(const std::string& prj_config_file);

		void initialize(); // top-level initialization function. distribute polygons in some way	

		/** access functions **/
		const TriMesh& mesh_domain() const { return mesh; }
		const vector<Packing_object>& get_tiles() const { return pack_objects; }
		RestrictedPolygonVoronoiDiagram& get_rpvd() const { return rpvd; }

		static Packer* instance() { return instance_; }

		/** optimization **/
		// one Lloyd iteration
		Lloyd_res lloyd(bool enlarge = true);

		/** miscellaneous **/
		void get_bbox(real& x_min, real& y_min, real& z_min, real& x_max, real& y_max, real& z_max);

		// driver
		void pack(void (*post_action)() = NULL); 
		// debug
		void rpack(void (*post_action)() = NULL);
		void update_iDT() { for (int i = 0; i < 1; i++) rpvd.iDT_update(); compute_clipped_VD();}

	private:

		/** initialization **/
		void random_init_tiles(unsigned int nb_init_polygons);	// put polygons according to curvature (if any) and in a uniform way

		/** geometry **/
		void generate_RDT();

		void compute_clipped_VD();

		vec3 approx_normal(unsigned int facet_idx);

		/** optimization **/
		static int callback(const int evalRequestCode, const int n, const int m, const int nnzJ, const int nnzH,
							const double * const x,	const double * const lambda, double * const obj, double * const c,
							double * const objGrad,	double * const jac,	double * const hessian,	double * const hessVector, void * userParams);
		// optimize by Knitro
		int KTR_optimize(double* io_k, double* io_theta, double* io_t1, double* io_t2, unsigned int idx);

		// get the transformation parameter for one polygon in a local frame
		struct Local_frame;
		struct Parameter;
		Optimization_res optimize_one_polygon(unsigned int id, Local_frame& lf, Parameter& solution);

		// constraint transformation to avoid triangle flip
		void constraint_transformation(vector<Parameter>& parameters, vector<Local_frame>& lfs);


	private:

		static Packer *instance_;
		
		ProjectIO pio;
		std::vector<Packing_object> pack_objects;
		std::vector<Polygon_2> pgn_lib;
		TriMesh mesh;

		/** optimization **/
#ifdef _CILK_
		std::vector<Containment> containments;
#else
		Containment containment;
#endif

		static double PI;
		double rot_lower_bd;
		double rot_upper_bd;

		/** helper classes **/

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
			void print()
			{
				std::cout<<"u = ("<<u.x()<<','<<u.y()<<','<<u.z()<<"), ";
				std::cout<<"v = ("<<v.x()<<','<<v.y()<<','<<v.z()<<"), ";
				std::cout<<"w = ("<<w.x()<<','<<w.y()<<','<<w.z()<<"), ";
				std::cout<<"o = ("<<o.x()<<','<<o.y()<<','<<o.z()<<"). \n";
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
				k = std::max(f*k, 1.0); 
				theta *= f;	tx *= f; ty *= f;
				return *this;
			}
			inline Parameter operator*(double f)
			{
				return Packer::Parameter(std::max(f*k, 1.0), f*theta, f*tx, f*ty);
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
				//v->Embedded_point_2::p = operator()(v->point_3());
				Point_3 p = operator()(v->point());
				v->point() = p;
				return v;
			}
		};

	public: // for debug
		RestrictedPolygonVoronoiDiagram rpvd;


	};


}

#endif