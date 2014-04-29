#pragma once

#include <numeric>
#include "spm_cgal.h"


namespace Geex {

	class Constraint
	{
	public:
		Constraint(const Point_2& p, const Segment_2& seg, const Point_2& ref_point): p_(p), dist_(0.0), seg_(seg)
		{
			Vector_2 ref_v(seg_.source(), ref_point);
			n_ = seg_.to_vector().perpendicular(CGAL::COUNTERCLOCKWISE);
			if (n_*ref_v < 0)
				n_ = -n_;
			n_ = n_ / CGAL::sqrt(n_.squared_length());
		}
		Constraint(const Point_2& p,const Segment_2& seg ) : p_(p), dist_(0.0), seg_(seg)
		{
			Vector_2 v(seg_.source(), p_);
			n_ = seg_.to_vector().perpendicular(CGAL::COUNTERCLOCKWISE);
			if ( v*n_ < 0)
				n_ = -n_;
			n_ = n_/sqrt(n_.squared_length());  // normalize
		}
		Constraint(const Point_2& p,const Segment_2& seg, double dist ): p_(p), dist_(dist), seg_(seg)
		{
			//assume p has been already inside the region when initialization
			Vector_2 v(seg_.source(), p_);
			n_ = seg_.to_vector().perpendicular(CGAL::COUNTERCLOCKWISE);
			if ( v*n_ < 0)
				n_ = -n_;		
			n_ = n_/sqrt(n_.squared_length());  // normalize
		}
	public:
		Point_2 p_; // the transformed p should be on the positive side of the segment
		Segment_2 seg_;
		double dist_;
		Vector_2 n_; // unit vector perpendicular to the segment

		double transformed_signed_area(const Transformation_2& aff)
		{
			Point_2 pp = aff(p_);
			return (pp - seg_.source())*n_;
		}
		double transformed_signed_dist(const Transformation_2& aff, double k)
		{
			Point_2 pp = aff(p_);
			return (pp - seg_.source())*n_ - k*dist_; 
		}
		void transform(const Transformation_2& aff)
		{
			p_ = aff(p_);
			seg_ = seg_.transform(aff);
		}
	};

	class Containment
	{
		typedef CGAL::Bbox_2 Bbox;

	public:
		void load_polygons(const std::string& inner_pgn_filename,	const std::string& outer_pgn_filename);
		void load_polygons(const Polygon_2& inner_pgn, const Polygon_2& outer_pgn);
		void load_polygons(const Polygon_2& inner_pgn, const std::vector<Segment_2>& outer_pgn);// outer polygon represented as a vector of segments
		void set_constraint_list();
		void set_constraint_list_convex();           // the outer polygon is convex
		void set_constraint_list_concave();  // the outer polygon is concave
		void insert_constraint(const Point_2& p,const Segment_2& seg, double dist);  // the distance from the transformed p to the line should be greater than dist
		void centralize_constraints(const Vector_2& trans);

		void compute_constraints(const double * const x, double* c, int mc);  // interface for Knitro
		void vcompute_constraints(const double * const x, double* c, int mc);

		void compute_constraint_grads(const double * const x, double* jac);
		void vcompute_constraint_grads(const double * const x, double* jac);

		int get_constaint_list_size(){return const_list_.size();}

		// for motion direction control
		double compute_extended_object(double k, double tx, double ty);
		void compute_extended_object_grad(double k, double tx, double ty, double *grad_k, double *grad_tx, double *grad_ty);

		//void afftrans2inner_pgn(const Transformation_2& aff){inner_pgn_ = CGAL::transform(aff, inner_pgn_);}
		void afftrans2inner_pgn(const Transformation_2& aff){inner_pgn_ = CGAL::transform(aff, origin_inner_pgn_);}

		void clear();

		void get_bbox(
			double& x_min, double& y_min, double& z_min,
			double& x_max, double& y_max, double& z_max
			);
		Polygon_2& get_inner_polygon() {return inner_pgn_;}
		Polygon_2& get_outer_polygon() {return outer_pgn_;}
		Polygon_2& get_origin_inner_polygon() { return origin_inner_pgn_; }

		Vector_2 getCenteringVector() const {return vTrans_;}//get the rectification to apply to the interior polygons

	public:
		std::vector<Constraint> const_list_;

		double align_angle;

		// extension
		static double lambda;
		double region_area;
		std::vector<double> weights; // angle between triangle normal and tangent plane normal
		std::vector<Point_2> vd_vertices;
		Point_2 pgn_cent;

	private:
		Polygon_2 inner_pgn_;  // the inner polygon
		Polygon_2 outer_pgn_;  // the outer polygon
		Polygon_2 origin_inner_pgn_;

		Vector_2 vTrans_ ;  // save the translation of the polygons: the inner polygon is centralized before optimization
	};


}