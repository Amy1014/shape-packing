#pragma once

#include "spm_cgal.h"

namespace Geex {

	class Constraint
	{
	public:
		Constraint(const Point_2& p,const Segment_2& seg ) : p_(p), dist_(0.0), seg_(seg)
		{
			//if (seg_.supporting_line().has_on_negative_side(p_))
			//	seg_ = seg_.opposite();
			//Vector_2 dir = seg_.to_vector();
			//n_ = dir.perpendicular(CGAL::COUNTERCLOCKWISE);
			//Vector_2 sp(seg_.source(), p_);
			//if ( sp*n_< 0.0)
			//	n_ = -n_;
			Vector_2 v(seg_.source(), p_);
			n_ = seg_.to_vector().perpendicular(CGAL::COUNTERCLOCKWISE);
			if ( v*n_ < 0)
				n_ = -n_;
			n_ = n_/sqrt(n_.squared_length());  // normalize
		}
		Constraint(const Point_2& p,const Segment_2& seg, double dist ): p_(p), dist_(dist), seg_(seg)
		{
			//if (seg_.supporting_line().has_on_negative_side(p_))
			//	seg_ = seg_.opposite();
			//Point_2 prj = seg_.supporting_line().projection(p_);
			//assume p has been already inside the region when initialization
			Vector_2 v(seg_.source(), p_);
			n_ = seg_.to_vector().perpendicular(CGAL::COUNTERCLOCKWISE);
			if ( v*n_ < 0)
			{
				//seg_ = seg_.opposite();
				//dir = seg_.to_vector();
				n_ = -n_;
			}			
			//n_ = dir.perpendicular(CGAL::COUNTERCLOCKWISE);
			//Vector_2 sp(seg_.source(), p_);
			//if ( sp*n_ < 0.0)
			//	n_ = -n_;
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
			//Vector_2 dir = seg_.to_vector();
			//Vector_2 n = dir.perpendicular(CGAL::COUNTERCLOCKWISE);
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
		
		void set_constraint_list();
		void set_constraint_list_convex();           // the outer polygon is convex
		void set_constraint_list_concave();  // the outer polygon is concave
		void insert_constraint(const Point_2& p,const Segment_2& seg, double dist);  // the distance from the transformed p to the line should be greater than dist
		void centralize_constraints(const Vector_2& trans);

		void compute_constraints(const double * const x, double* c, int mc);  // interface for Knitro
		void compute_translation_constraint(const double *const x, double *c);
		void compute_constraint_grads(const double * const x, double* jac);
		void compute_translation_constraint_grads(const double *const x, double *jac);
		int get_constaint_list_size(){return const_list_.size();}

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

	private:
		Polygon_2 inner_pgn_;  // the inner polygon
		Polygon_2 outer_pgn_;  // the outer polygon
		Polygon_2 origin_inner_pgn_;

		Vector_2 vTrans_ ;  // save the translation of the polygons: the inner polygon is centralized before optimization
	};


}