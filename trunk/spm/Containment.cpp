#include "Containment.h"
#include <fstream>

namespace Geex {

void Containment::load_polygons(const std::string& inner_pgn_filename,const std::string& outer_pgn_filename)
{
	inner_pgn_.clear();
	outer_pgn_.clear();
	std::ifstream in(inner_pgn_filename.c_str()) ;
	in >> inner_pgn_;
	in.close();
	if(inner_pgn_.is_empty())
	{
		std::cout<<"Failed to load file: "<<inner_pgn_filename<<std::endl;
		return;
	}

	std::cout<<"Load file: "<<inner_pgn_filename<<std::endl;

	Bbox bbox = inner_pgn_.bbox();
	vTrans_ = Vector_2((bbox.xmin()+bbox.xmax())*0.5, (bbox.ymin()+bbox.ymax())*0.5);
	std::cout<<inner_pgn_.size()<<" vertices,   bbox center: "<<vTrans_<<std::endl ;

	in.open(outer_pgn_filename.c_str());
	in >> outer_pgn_;
	in.close();
	std::cout<<"Load file: "<<outer_pgn_filename<<std::endl;
	std::cout<<outer_pgn_.size()<<" vertices\n";
	
	if (inner_pgn_.is_clockwise_oriented())
		inner_pgn_.reverse_orientation();

	if(outer_pgn_.is_clockwise_oriented())
		outer_pgn_.reverse_orientation();

	Transformation_2 translate(CGAL::TRANSLATION, -vTrans_);
	inner_pgn_ = CGAL::transform(translate, inner_pgn_);
	outer_pgn_ = CGAL::transform(translate, outer_pgn_);

}

 void Containment::load_polygons(const Polygon_2& inner_pgn, const Polygon_2& outer_pgn)
 {	
 	inner_pgn_ = inner_pgn;
 	outer_pgn_ = outer_pgn;

	if (inner_pgn_.is_clockwise_oriented())
		inner_pgn_.reverse_orientation();

	if(outer_pgn_.is_clockwise_oriented())
		outer_pgn_.reverse_orientation();

	origin_inner_pgn_ = inner_pgn_;
 }

void Containment::get_bbox( double& x_min, double& y_min, double& z_min,
						   double& x_max, double& y_max, double& z_max )
{
	if(outer_pgn_.is_empty())
		return;
	Bbox bbox = outer_pgn_.bbox();
	x_max = bbox.xmax();
	x_min = bbox.xmin();
	y_max = bbox.ymax();
	y_min = bbox.ymin();

	z_min = z_max =0;
}

void Containment::clear()
{
	inner_pgn_.clear();
	outer_pgn_.clear();
	const_list_.clear();
}

void Containment::insert_constraint(const Point_2& p,const Segment_2& seg, double dist)
{
	const_list_.push_back(Constraint(p, seg,  dist));
}

void Containment::centralize_constraints(const Vector_2& trans)
{
	vTrans_ = trans;
	Transformation_2 translate(CGAL::TRANSLATION, vTrans_);
	for (int i=0; i<const_list_.size(); i++)
		const_list_[i].transform(translate);
}

void Containment::set_constraint_list()
{
	if (outer_pgn_.is_convex())
		set_constraint_list_convex();
	else
		set_constraint_list_concave();
}

void Containment::set_constraint_list_convex()
{
	const_list_.clear();
	for (int i=0; i<inner_pgn_.size();i++)
	{
		Point_2 p = inner_pgn_.vertex(i);
		for (unsigned int j=0; j<outer_pgn_.size(); j++)
		{
			Segment_2& seg = outer_pgn_.edge(j);
			if (!seg.is_degenerate())
				const_list_.push_back(Constraint(p, seg));
		}
	}
}

void Containment::set_constraint_list_concave()
{
	const_list_.clear();
	for (int j=0; j<outer_pgn_.size(); j++)
	{
		Point_2 q = outer_pgn_.vertex(j);
		Point_2 near_p;
		// find the nearest vertex of inner polygon
		double dist2 = std::numeric_limits<double>::max();
		for (int i=0; i<inner_pgn_.size();i++)
		{
			double temp_dist2 = (inner_pgn_.vertex(i)-q).squared_length();
			if(temp_dist2<dist2)
			{	
				dist2 = temp_dist2;
				near_p = inner_pgn_.vertex(i);
			}
		}

		Segment_2& seg_right = outer_pgn_.edge(j);
		int jj = (j==0?outer_pgn_.size()-1:j-1);
		Segment_2& seg_left = outer_pgn_.edge(jj);
		if (!seg_right.is_degenerate())
			const_list_.push_back(Constraint(near_p,seg_right ));
		if (!seg_left.is_degenerate())
			const_list_.push_back(Constraint(near_p,seg_left ));
	}

}

void Containment::compute_constraints(const double * const x, double* c, int mc)
{
	assert(mc == const_list_.size());
	// x[]: k theta t0 t1
	double fcos  = cos(x[1]); 
	double fsin  = sin(x[1]); 
	Transformation_2 aff(x[0]*fcos, -x[0]*fsin, x[2], x[0]*fsin,  x[0]*fcos, x[3]);
	for (unsigned int i=0; i<const_list_.size(); i++)
		c[i] =const_list_[i].transformed_signed_dist(aff, 1.0); // >=0	
		//c[i] = const_list_[i].transformed_signed_area(aff);
}
void Containment::compute_translation_constraint(const double *const x, double *c)
{
	Transformation_2 aff(x[0], 0.0, x[1], 0.0, x[0], x[2]);
	for (unsigned int i = 0; i < const_list_.size(); i++)
		c[i] = const_list_[i].transformed_signed_area(aff);
}
void Containment::compute_constraint_grads(const double * const x, double* jac)
{
	double fcos  = cos(x[1]); 
	double fsin  = sin(x[1]); 
	Transformation_2 rotate(fcos, -fsin, 0, fsin,  fcos, 0);
	Transformation_2 aff(-x[0]*fsin, -x[0]*fcos, 0,	x[0]*fcos,  -x[0]*fsin, 0);

	int idx =0;
	for (int i=0; i<const_list_.size(); i++)
	{
		Point_2 pp = rotate(const_list_[i].p_);
		jac[idx] = (pp-CGAL::ORIGIN)*const_list_[i].n_-const_list_[i].dist_;

		pp = aff(const_list_[i].p_);
		jac[idx+1] = (pp-CGAL::ORIGIN)*const_list_[i].n_;

		jac[idx+2] = const_list_[i].n_.x();
		jac[idx+3] = const_list_[i].n_.y();

		idx += 4;
	}
}
void Containment::compute_translation_constraint_grads(const double *const x, double *jac)
{
	int idx = 0;
	for (unsigned int i = 0; i < const_list_.size(); i++)
	{
		jac[idx] = (const_list_[i].p_-CGAL::ORIGIN)*const_list_[i].n_;
		jac[idx+1] = const_list_[i].n_.x();
		jac[idx+2] = const_list_[i].n_.y();
		idx += 3;
	}
}
} // End of namespace

