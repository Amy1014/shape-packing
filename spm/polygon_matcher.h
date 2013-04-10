#ifndef _POLYGON_MATCHER_
#define	_POLYGON_MATCHER_

/************************************************************************/
/*					Polygon Match using affine match                    */
/************************************************************************/
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cv.h>
#include <opencv2/core/types_c.h>
#include <CGAL/Dimension.h>
#include <CGAL/centroid.h>
#include "spm_cgal.h"
#include "CmAffine.h"

namespace Geex
{

template <class UserDataType>
struct Match_info_item
{
	Transformation_2 t;
	double scale;
	//double theta;
	//double tx;
	//double ty;
	double error;//match error
	UserDataType val; // index to library
};

struct Match_measure 
{
	template <class UserDataType>
	inline bool operator()(const Match_info_item<UserDataType>& i0, const Match_info_item<UserDataType>& i1) const
	{
		return i0.error > i1.error;
	}
};

class Polygon_matcher
{
	typedef std::vector<cv::Point2d> PointSet; 
public:
	Polygon_matcher(std::vector<Segment_2>& model, unsigned int _nb_sampling = 100);
	~Polygon_matcher();

	template <class UserDataType>
	Match_info_item<UserDataType> affine_match(const Polygon_2& instance, const UserDataType& data, double lambda = 0.0);

	double get_model_area() const { return model_area; }

private:
	unsigned int nb_sampling;
	PointSet modelPtsSet;
	std::vector<Vector_2> model_tangent_vec; // tangent vectors at each sampling points
	CvMat *modelMat;
	double model_area;
};



template <class UserDataType>
Match_info_item<UserDataType> Polygon_matcher::affine_match(const Polygon_2& instance_, const UserDataType& data, double lambda)
{
	Point_2 c = CGAL::centroid(instance_.vertices_begin(), instance_.vertices_end(), CGAL::Dimension_tag<0>());
	Transformation_2 to_cent(CGAL::TRANSLATION, Vector_2(c, CGAL::ORIGIN));
	Polygon_2 instance = CGAL::transform(to_cent, instance_);
	double perimeter = 0.0;
	std::vector<double> edge_lens(instance.size());
	for (unsigned int i = 0; i < instance.size(); i++)
	{
		Segment_2 s = instance.edge(i);
		edge_lens[i] = CGAL::sqrt(s.squared_length());
		perimeter += edge_lens[i];
	}
	nb_sampling = 200;
	PointSet instancePtsSet;

	std::vector<Vector_2> instance_tangent_vec;
	instance_tangent_vec.reserve(nb_sampling);

	instancePtsSet.reserve(nb_sampling);
	for (unsigned int i = 0; i < instance.size(); i++)
	{
		unsigned int n = (edge_lens[i]/perimeter) * nb_sampling;
		Segment_2 s = instance.edge(i);
		Point_2 src = s.source(), tgt = s.target();

		Vector_2 tv(src, tgt);
		cgal_vec_normalize(tv);

		for (unsigned int j = 0; j < n+1; j++)
		{
			Point_2 sp = CGAL::ORIGIN + ( (n+1-j)*(src - CGAL::ORIGIN) + j*(tgt - CGAL::ORIGIN) ) / (n+1);
			instancePtsSet.push_back(cv::Point2d(sp.x(), sp.y()));

			instance_tangent_vec.push_back(tv);
		}
	}
	CvMat *instanceMat = CmAffine::CreatePointSetMat(instancePtsSet);

	CvMat *AA = cvCreateMat(2, 2, CV_64FC1);
	CvMat *shift = cvCreateMat(2, 1, CV_64FC1);

	double theta =	CmAffine::GetAffine2D(instanceMat, modelMat, AA, shift);

	double scale = std::sqrt(model_area/std::fabs(instance.area()));

	double *r = AA->data.db;
	r[0]  = scale*std::cos(theta); r[1] = -scale*std::sin(theta);
	r[2] = -r[1];   r[3] = r[0];

	//Transformation_2 t(r[0], r[1], shift->data.db[0], r[2], r[3], shift->data.db[1]);

	CvMat *ptsMat = cvCreateMat(2, instancePtsSet.size(), CV_64FC1);
	CmAffine::FindTranslated(instanceMat, AA, shift, ptsMat);

	double *x = ptsMat->data.db;
	double *y = (double*)(ptsMat->data.ptr + ptsMat->step);
	for (int i = 0; i < ptsMat->cols; i++)
		instancePtsSet[i] = Point2d(x[i], y[i]);

	double mean_dist = 0.0;
	double mean_tangent_diff = 0.0;
	for (unsigned int i = 0; i < instancePtsSet.size(); i++)
	{
		unsigned int min_idx;
		double min_dist = std::numeric_limits<double>::max();
		for (unsigned int j = 0; j < modelPtsSet.size(); j++)
		{
			double temp = sqrDist(instancePtsSet[i], modelPtsSet[j]);
			if (min_dist > temp)
			{
				min_dist = temp;
				min_idx = j;
			}
		}
		mean_dist += std::sqrt(min_dist);
		mean_tangent_diff += std::fabs(instance_tangent_vec[i].x()*model_tangent_vec[min_dist].y() - instance_tangent_vec[i].y()*model_tangent_vec[min_dist].x());
	}
	mean_dist /= instancePtsSet.size();
	mean_tangent_diff /= instancePtsSet.size();
	mean_tangent_diff *= CGAL::sqrt(model_area);

	Match_info_item<UserDataType> res;
	//std::cout<<"( "<<mean_dist<<", "<<mean_tangent_diff<<" )\n";
	res.error = mean_dist + lambda * mean_tangent_diff; 
	res.t = Transformation_2(r[0], r[1], shift->data.db[0], r[2], r[3], shift->data.db[1])*to_cent;
	res.scale = scale;
	//res.theta = theta;
	//res.tx = shift->data.db[0];
	//res.ty = shift->data.db[1];
	res.val = data;

	cvReleaseMat(&instanceMat);
	cvReleaseMat(&shift);
	cvReleaseMat(&AA);
	cvReleaseMat(&ptsMat);

	return res;
}


}

#endif