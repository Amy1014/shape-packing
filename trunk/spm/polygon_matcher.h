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
	//Transformation_2 t;
	double k;
	double theta;
	double tx;
	double ty;
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
	Match_info_item<UserDataType> affine_match(const Polygon_2& instance, const UserDataType& data);

private:
	unsigned int nb_sampling;
	PointSet modelPtsSet;
	CvMat *modelMat;
	double model_area;
};



template <class UserDataType>
Match_info_item<UserDataType> Polygon_matcher::affine_match(const Polygon_2& instance, const UserDataType& data)
{
	double perimeter = 0.0;
	std::vector<double> edge_lens(instance.size());
	for (unsigned int i = 0; i < instance.size(); i++)
	{
		Segment_2 s = instance.edge(i);
		edge_lens[i] = CGAL::sqrt(s.squared_length());
		perimeter += edge_lens[i];
	}
	PointSet instancePtsSet;
	instancePtsSet.reserve(nb_sampling);
	for (unsigned int i = 0; i < instance.size(); i++)
	{
		unsigned int n = (edge_lens[i]/perimeter) * nb_sampling;
		Segment_2 s = instance.edge(i);
		Point_2 src = s.source(), tgt = s.target();
		for (unsigned int j = 0; j < n+1; j++)
		{
			Point_2 sp = CGAL::ORIGIN + ( (n+1-j)*(src - CGAL::ORIGIN) + j*(tgt - CGAL::ORIGIN) ) / (n+1);
			instancePtsSet.push_back(cv::Point2d(sp.x(), sp.y()));
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
	for (unsigned int i = 0; i < modelPtsSet.size(); i++)
	{
		double min_dist = std::numeric_limits<double>::max();
		for (unsigned int j = 0; j < instancePtsSet.size(); j++)
			min_dist = std::min(min_dist, sqrDist(modelPtsSet[i], instancePtsSet[j]));
		mean_dist += std::sqrt(min_dist);
	}
	mean_dist /= modelPtsSet.size();

	Match_info_item<UserDataType> res;
	res.error = mean_dist; 
	//res.t = Transformation_2(rescalor*r[0], rescalor*r[1], shift->data.db[0], rescalor*r[2], rescalor*r[3], shift->data.db[1]);
	res.k = scale;
	res.theta = theta;
	res.tx = shift->data.db[0];
	res.ty = shift->data.db[1];
	res.val = data;

	cvReleaseMat(&instanceMat);
	cvReleaseMat(&shift);
	cvReleaseMat(&AA);
	cvReleaseMat(&ptsMat);

	return res;
}


}

#endif