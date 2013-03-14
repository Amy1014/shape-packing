#ifndef _PACKER_H_
#define _PACKER_H_

#include <vector>
#include <string>
#include <cfloat>
#include <cmath>
#include <set>
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

		typedef enum {LLOYD_SUCCESS, LLOYD_FAILED} Lloyd_res;

	public:
		
		Packer();

		/** initialization **/
		void load_project(const std::string& prj_config_file);

		void initialize(); // top-level initialization function. distribute polygons in some way	

		/** access functions **/
		const TriMesh& mesh_domain() const { return mesh; }
		const vector<Packing_object>& get_tiles() const { return pack_objects; }

		static Packer* instance() { return instance_; }

		/** miscellaneous **/
		void get_bbox(real& x_min, real& y_min, real& z_min, real& x_max, real& y_max, real& z_max);

		/** debug **/
		void draw_RDT() { rpvd.draw_DT(); }
		void draw_clipped_VD() {rpvd.draw_clipped_VD();}

	private:

		/** initialization **/
		void random_init_tiles(unsigned int nb_init_polygons);	// put polygons according to curvature (if any) and in a uniform way

		/** geometry **/
		void generate_RDT();

		vec3 approx_normal(unsigned int facet_idx);

		/** optimization **/
		static int callback(const int evalRequestCode, const int n, const int m, const int nnzJ, const int nnzH,
							const double * const x,	const double * const lambda, double * const obj, double * const c,
							double * const objGrad,	double * const jac,	double * const hessian,	double * const hessVector, void * userParams);
		// optimize by Knitro
		int KTR_optimize(double* io_k, double* io_theta, double* io_t1, double* io_t2);

		// get the transformation parameter for one polygon in a local frame
		Lloyd_res optimize_one_polygon(unsigned int id, Local_frame& lf, Parameter& solution);

	private:

		static Packer *instance_;
		
		ProjectIO pio;
		std::vector<Packing_object> pack_objects;
		std::vector<Polygon_2> pgn_lib;
		TriMesh mesh;

		/** optimization **/
		Containment containment;
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
			inline Point_2 to_uv(const Point_3& p)
			{
				Vector_3 op(o, p);
				return Point_2(op*u, op*v);
			}
			inline Segment_2 to_uv(const Segment_3& s)
			{
				return Segment_2(to_uv(s.source()), to_uv(s.target()));
			}
			inline Point_3 to_xy(const Point_2& p)
			{
				return CGAL::ORIGIN + ( p.x()*u + p.y()*v );
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
		};

		// apply a 2D transformation to a 3D polygon inside the plane where the polygon is embedded
		struct Apply_transformation
		{
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
		}

	public: // for debug
		RestrictedPolygonVoronoiDiagram rpvd;


	};
}

#endif