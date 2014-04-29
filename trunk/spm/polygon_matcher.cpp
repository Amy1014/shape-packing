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
// 		std::vector<Point_2> ps;
// 		ps.reserve(model.size()*2);
// 		for (unsigned int i = 0; i < model.size(); i++)
// 		{
// 			ps.push_back(model[i].source());
// 			ps.push_back(model[i].target());
// 		}
// 		Point_2 p = CGAL::centroid(ps.begin(), ps.end(), CGAL::Dimension_tag<0>());
// 		for (unsigned int i = 0; i < model.size(); i++)
// 		{
// 			model_area += std::fabs(CGAL::area(p, model[i].source(), model[i].target()));
// 		}
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

	void compute_area_centroid(const Polygon_2& poly, double& area, vec2& centroid)
	{

		vec3 sarea(0,0,0);

		for (int i=1; i<poly.size()-1; i++)
		{
			vec3 p0, pi, pii;
			p0 = vec3(poly.vertex(0).x(),poly.vertex(0).y(),0);
			pi = vec3(poly.vertex(i).x(),poly.vertex(i).y(),0);
			pii = vec3(poly.vertex(i+1).x(),poly.vertex(i+1).y(),0);

			sarea += cross(pi-p0,pii-p0);
		}

		area = sarea.length();

		vec3 normal = normalize(sarea);
		vec3 cent (0,0,0);

		for (int i=1; i<poly.size()-1; i++)
		{

			vec3 p0, pi, pii;
			p0 = vec3(poly.vertex(0).x(),poly.vertex(0).y(),0);
			pi = vec3(poly.vertex(i).x(),poly.vertex(i).y(),0);
			pii = vec3(poly.vertex(i+1).x(),poly.vertex(i+1).y(),0);

			vec3 sarea_tri = cross(pi-p0,pii-p0);
			double sq_area = sarea_tri.length();
			vec3 center = 1.0/3.0*(p0+pi+pii);

			if (dot(sarea_tri,normal) > 0)
			{
				cent += sq_area*center;
			}
			else{
				cent -= sq_area*center;
			}
		}

		cent /= area;

		centroid.x = cent.x;
		centroid.y = cent.y;
	}

	void compute_area_centroid_hole(const HoleEdges& hole, double& area, vec2& centroid )
	{
		vec2 ref (0,0);

		double total_len = 0;
		for (int i=0; i<hole.size(); i++)
		{
			double len = (hole[i].first - hole[i].second).length();
			ref += 0.5*(hole[i].first + hole[i].second)*len;
			total_len += len;
		}

		ref /= total_len;

		vec3 ref_v3 (ref.x, ref.y, 0);
		area = 0;
		centroid = vec2(0,0);

		for (int i=0; i<hole.size(); i++)
		{
			vec3 pi0 (hole[i].first.x, hole[i].first.y, 0);
			vec3 pi1 (hole[i].second.x, hole[i].second.y, 0);

			double tri_area = cross(pi0-ref_v3, pi1-ref_v3).length();
			area += tri_area;
			centroid += tri_area*(ref + hole[i].first + hole[i].second)/3.0;
		}

		centroid /= area;
	}

	void sample_poly(const Polygon_2& poly, std::vector<vec2>& samples, int total_num)
	{
		//compute poly boundary length.
		double len = 0;
		for (int i=0; i<poly.size()-1; i++)
		{
			len += (to_geex_pnt_2d(poly.vertex(i)) - to_geex_pnt_2d(poly.vertex(i+1))).length();
		}

		for (int i=0; i<poly.size()-1; i++)
		{
			vec2 vi, vii;
			vi = to_geex_pnt_2d(poly.vertex(i));
			vii = to_geex_pnt_2d(poly.vertex(i+1));

			double e_len = (vii - vi).length();

			samples.push_back(vi);

			int num = int(e_len/len*total_num)-1;
			if(num > 0){
				for (int j=0; j<num; j++)
				{
					double ratio = double(j+1)/double(num+1);
					samples.push_back(ratio*vii + (1-ratio)*vi);
				}
			}
		}
	}
	void sample_hole(const HoleEdges& hole, std::vector<vec2>& samples, int total_num)
	{
		//compute poly boundary length.
		double len = 0;
		for (int i=0; i<hole.size(); i++)
		{
			len += (hole[i].first - hole[i].second).length();
		}

		for (int i=0; i<hole.size(); i++)
		{
			vec2 vi, vii;
			vi = hole[i].first;
			vii = hole[i].second;

			double e_len = (vii - vi).length();

			samples.push_back(vi);

			int num = int(e_len/len*total_num)-1;
			if(num > 0){
				for (int j=0; j<num; j++)
				{
					double ratio = double(j+1)/double(num+1);
					samples.push_back(ratio*vii + (1-ratio)*vi);
				}
			}
		}
	}


	double HausdorffDist (const std::vector<vec2>& A, const std::vector<vec2>& B)
	{
		double max_dist = 0;
		std::vector<double> Adist (A.size(), 1e10);
		std::vector<double> Bdist (B.size(), 1e10);

		for (int i=0; i<A.size(); i++)
		{
			for (int j=0; j<B.size(); j++)
			{
				double len = (A[i] - B[j]).length();
				if (len < Adist[i])
				{
					Adist[i] = len;
				}
				if (len < Bdist[j])
				{
					Bdist[j] = len;
				}
			}
		}

		for (int i=0; i<Adist.size(); i++)
		{
			max_dist = max(max_dist, Adist[i]);
		}
		for (int i=0; i<Bdist.size(); i++)
		{
			max_dist = max(max_dist, Bdist[i]);
		}

		return max_dist;
	}

	double matching_score(const Polygon_2& poly_a, const Polygon_2& poly_b, double& tx, double& ty, double& scale, double& rotate)
	{
		//compute area and centroid of two polys.
		double area_a, area_b;
		vec2 cent_a, cent_b;
		compute_area_centroid(poly_a, area_a, cent_a);
		compute_area_centroid(poly_b, area_b, cent_b);

		//generate sample points.
		const int sample_num = 100;
		std::vector<vec2> samples_a, samples_b;
		sample_poly(poly_a, samples_a, sample_num);
		sample_poly(poly_b, samples_b, sample_num);

		//scale second to make area match.
		scale = sqrt(area_a/area_b);
		for (int i=0; i<samples_a.size(); i++)
		{
			samples_a[i] = cent_a + (samples_a[i] - cent_a)/area_a;
		}
		for (int i=0; i<samples_b.size(); i++)
		{
			samples_b[i] = cent_b + (samples_b[i] - cent_b)/area_b;
		}

		//translate second sample to make centroid coincide.
		vec2 diff = cent_a - cent_b;
		tx = diff.x;
		ty = diff.y;
		for (int i=0; i<samples_b.size(); i++)
		{
			samples_b[i] += diff;
		}


		double min_dist=1e10;
		double min_dist_angle=0;
		const int slice_num = 36;//36
		for (int i=0; i<slice_num; i++)
		{
#define PI 3.14159265359
			double angle = i*2*PI/slice_num;

			vec2 rot_mat[2];
			rot_mat[0] = vec2(cos(angle), -sin(angle));
			rot_mat[1] = vec2(-rot_mat[0].y, rot_mat[0].x);

			std::vector<vec2> rot_b = samples_b;
			//rotate
			for (int j=0; j<rot_b.size(); j++)
			{
				vec2 vcp = rot_b[j]-cent_a;
				rot_b[j] = cent_a + vec2(dot(rot_mat[0],vcp), dot(rot_mat[1],vcp));
			}

			//compute Hausdorff dist.
			double hdist = HausdorffDist(samples_a, rot_b);

			if (hdist < min_dist)
			{
				min_dist = hdist;
				min_dist_angle = angle;
			}
		}

		rotate = min_dist_angle;

		return min_dist;
	}
	void transform_poly(Polygon_2& poly, double tx, double ty, double scale, double rot_angle){
		vec2 rot_mat[2];
		rot_mat[0] = vec2(cos(rot_angle), -sin(rot_angle));
		rot_mat[1] = vec2(-rot_mat[0].y, rot_mat[0].x);

		double area;
		vec2 cent;
		compute_area_centroid(poly, area, cent);

		vec2 new_cent = cent + vec2(tx,ty);

		for (Polygon_2::Vertex_iterator vitr = poly.vertices_begin(); vitr != poly.vertices_end(); vitr++)
		{
			vec2 pt = to_geex_pnt_2d(*vitr);
			vec2 diff = (pt - cent)*scale;

			diff = vec2(dot(rot_mat[0], diff), dot(rot_mat[1], diff));

			poly.set(vitr, to_cgal_pnt_2d(new_cent + diff));
		}
	}

	double matching_score(const HoleEdges& poly_a, const Polygon_2& poly_b, double& tx, double& ty, double& scale, double& rotate)
	{
		//compute area and centroid of two polys.
		double area_a, area_b;
		vec2 cent_a, cent_b;
		compute_area_centroid_hole(poly_a, area_a, cent_a);
		compute_area_centroid(poly_b, area_b, cent_b);

		//generate sample points.
		const int sample_num = 100;
		std::vector<vec2> samples_a, samples_b;
		sample_hole(poly_a, samples_a, sample_num);
		sample_poly(poly_b, samples_b, sample_num);

		//scale second to make area match.
		scale = sqrt(area_a/area_b);
		double sqrt_area_a = sqrt(area_a);
		double sqrt_area_b = sqrt(area_b);
		for (int i=0; i<samples_a.size(); i++)
		{
			samples_a[i] = cent_a + (samples_a[i] - cent_a)/sqrt_area_a;
		}
		for (int i=0; i<samples_b.size(); i++)
		{
			samples_b[i] = cent_b + (samples_b[i] - cent_b)/sqrt_area_b;
		}

		//translate second sample to make centroid coincide.
		vec2 diff = cent_a - cent_b;
		tx = diff.x;
		ty = diff.y;
		for (int i=0; i<samples_b.size(); i++)
		{
			samples_b[i] += diff;
		}


		double min_dist=1e10;
		double min_dist_angle=0;
		const int slice_num = 36;//36
		for (int i=0; i<slice_num; i++)
		{
#define PI 3.14159265359
			double angle = i*2*PI/slice_num;

			vec2 rot_mat[2];
			rot_mat[0] = vec2(cos(angle), -sin(angle));
			rot_mat[1] = vec2(-rot_mat[0].y, rot_mat[0].x);

			std::vector<vec2> rot_b = samples_b;
			//rotate
			for (int j=0; j<rot_b.size(); j++)
			{
				vec2 vcp = rot_b[j]-cent_a;
				rot_b[j] = cent_a + vec2(dot(rot_mat[0],vcp), dot(rot_mat[1],vcp));
			}

			//compute Hausdorff dist.
			double hdist = HausdorffDist(samples_a, rot_b);

			if (hdist < min_dist)
			{
				min_dist = hdist;
				min_dist_angle = angle;
			}
		}

		rotate = min_dist_angle;

		return min_dist;
	}
}