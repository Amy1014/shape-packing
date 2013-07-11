#include "polygon_matcher.h"

namespace Geex
{
	Polygon_matcher::Polygon_matcher(std::vector<Segment_2>& model, unsigned int _nb_samp)
	{
		// compute area of model
		model_area = 0.0;
		Point_2 p = model[0].source();
		for (unsigned int i = 0; i < model.size(); i++)
			model_area += std::fabs(CGAL::area(p, model[i].source(), model[i].target()));

		nb_sampling = _nb_samp;

		// sample the model polygon
		double perimeter = 0.0;
		vector<double> edge_lens(model.size());
		for (unsigned int i = 0; i < model.size(); i++)
		{
			edge_lens[i] = CGAL::sqrt(model[i].squared_length());
			perimeter += edge_lens[i];
		}
		modelPtsSet.reserve(nb_sampling);
		model_tangent_vec.reserve(nb_sampling);
		for (unsigned int i = 0; i < model.size(); i++)
		{
			int n = (edge_lens[i]/perimeter)*nb_sampling;
			Point_2 src = model[i].source(), tgt = model[i].target();

			Vector_2 tv(src, tgt);
			cgal_vec_normalize(tv);

			for (unsigned int j = 0; j < n+1; j++)
			{				
				Point_2 sp = CGAL::ORIGIN + ( (n+1-j)*(src - CGAL::ORIGIN) + j*(tgt - CGAL::ORIGIN) ) / (n+1);
				modelPtsSet.push_back(cv::Point2d(sp.x(), sp.y()));
				model_tangent_vec.push_back(tv);
			}
		}
		modelMat = CmAffine::CreatePointSetMat(modelPtsSet);
	}

	Polygon_matcher::Polygon_matcher(const Polygon_2& model, unsigned int _nb_sampling /* = 100 */)
	{
		model_area = std::fabs(model.area());
		nb_sampling = _nb_sampling;
		// sample the model polygon
		double perimeter = 0.0;
		vector<double> edge_lens(model.size());
		for (unsigned int i = 0; i < model.size(); i++)
		{
			edge_lens[i] = CGAL::sqrt(model.edge(i).squared_length());
			perimeter += edge_lens[i];
		}
		modelPtsSet.reserve(nb_sampling);
		model_tangent_vec.reserve(nb_sampling);
		for (unsigned int i = 0; i < model.size(); i++)
		{
			int n = (edge_lens[i]/perimeter)*nb_sampling;
			Point_2 src = model.edge(i).source(), tgt = model.edge(i).target();

			Vector_2 tv(src, tgt);
			cgal_vec_normalize(tv);

			for (unsigned int j = 0; j < n+1; j++)
			{				
				Point_2 sp = CGAL::ORIGIN + ( (n+1-j)*(src - CGAL::ORIGIN) + j*(tgt - CGAL::ORIGIN) ) / (n+1);
				modelPtsSet.push_back(cv::Point2d(sp.x(), sp.y()));
				model_tangent_vec.push_back(tv);
			}
		}
		modelMat = CmAffine::CreatePointSetMat(modelPtsSet);
	}
	Polygon_matcher::~Polygon_matcher()
	{
		modelPtsSet.clear();
		model_tangent_vec.clear();
		cvReleaseMat(&modelMat);
	}
}