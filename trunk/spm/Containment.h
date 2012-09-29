#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <iostream>

typedef double Coord_type;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Line_2<K> Line;
typedef CGAL::Vector_2<K> Vector2;
typedef CGAL::Segment_2<K> Segment;
typedef CGAL::Bbox_2 Bbox;
typedef CGAL::Aff_transformation_2<K>  Transformation;

namespace Geex {

	class Constraint:public Segment
	{
	public:
		Constraint(const Point& p,const Segment& seg ){
			p_ = p;
			dist_ = 0;
			seg_ = seg;

			Vector2 dir = seg.to_vector();
			n_ = dir.perpendicular(CGAL::COUNTERCLOCKWISE);
			n_ = n_/sqrt(n_.squared_length());  // normalize
		}
		Constraint(const Point& p,const Segment& seg, double dist ): Segment(seg){
			p_ = p;
			dist_= dist;
			seg_= seg;

			Vector2 dir = seg.to_vector();
			n_ = dir.perpendicular(CGAL::COUNTERCLOCKWISE);
			n_ = n_/sqrt(n_.squared_length());  // normalize
		}
	public:
		Point p_; // the transformed p should be on the positive side of the segment
		Segment seg_;
		double dist_;
		Vector2 n_; // unit vector perpendicular to the segment
		double transformed_signed_area(const Transformation& aff)
		{
			Point pp = aff(p_);
			Vector2 dir = seg_.to_vector();
			Vector2 n = dir.perpendicular(CGAL::COUNTERCLOCKWISE);
			return (pp - seg_.source())*n;
		}
		double transformed_signed_dist(const Transformation& aff, double k)
		{
			Point pp = aff(p_);
			return (pp - seg_.source())*n_ - k*dist_; 
		}
		void transform(const Transformation& aff)
		{
			p_ = aff(p_);
			seg_ = seg_.transform(aff);
		}

	};

	class Containment
	{
	public:
		Containment(void) {}
		~Containment(void) {}
	public:
		void load_polygons(const std::string& inner_pgn_filename,	const std::string& outer_pgn_filename);
		
		void set_constraint_list();
		void set_constraint_list_convex();           // the outer polygon is convex
		void set_constraint_list_concave();  // the outer polygon is concave
		void insert_constraint(const Point& p,const Segment& seg, double dist);  // the distance from the transformed p to the line should be greater than dist
		void centralize_constraints(const Vector2& trans);

		void compute_constraints(const double * const x, double* c, int mc);  // interface for Knitro
		void compute_constraint_grads(const double * const x, double* jac, int mc);
		int get_constaint_list_size(){return const_list_.size();}

		void afftrans2inner_pgn(Transformation aff){inner_pgn_ = CGAL::transform(aff, inner_pgn_);}

		void clear();

		void get_bbox(
			double& x_min, double& y_min, double& z_min,
			double& x_max, double& y_max, double& z_max
			);
		Polygon_2& get_inner_polygon() {return inner_pgn_;}
		Polygon_2& get_outer_polygon() {return outer_pgn_;}

		Vector2 getCenteringVector() const {return vTrans_;}//get the rectification to apply to the interior polygons

	public:
		std::vector<Constraint> const_list_;

	private:
		Polygon_2 inner_pgn_;  // the inner polygon
		Polygon_2 outer_pgn_;  // the outer polygon



		Vector2 vTrans_ ;  // save the translation of the polygons: the inner polygon is centralized before optimization
	};


}