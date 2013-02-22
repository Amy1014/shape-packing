#include "SphericalBounder.h"

namespace Geex
{

#ifdef CILK
	tbb::mutex SphericalBounder::mut;
#endif

 	SphericalBounder::SphericalBounder()
 	{
 	}
	SphericalBounder::SphericalBounder(const vector<Segment_3>& aPolygon)
	{
		pgn = aPolygon;
		//msu = 0;
	}
//  	SphericalBounder::~SphericalBounder()
//  	{	}
	void SphericalBounder::setPolygon(const vector<Segment_3>& aPolygon)
	{
		pgn = aPolygon;
	}
	void SphericalBounder::genMarginBounder(const unsigned int nbSamping /* = 50 */)
	{
		if (msu.size()>0)
			return;
		//msu = new vector<Weighted_point>;
		double circumLen(0.0);
		vector<double> segLen; //length of each segment
		for (vector<Segment_3>::const_iterator it = pgn.begin(); it != pgn.end(); it++)
		{
			segLen.push_back(sqrt(it->squared_length()));
			circumLen += segLen.back();
		}
		for (unsigned int i = 0; i < pgn.size(); i++)
		{
			Point_3 s = pgn[i].source(), t = pgn[i].target();
			unsigned int n = static_cast<unsigned int>(segLen[i]/circumLen*nbSamping); //number of sampling on this segment
			n = std::max<unsigned int>(1, n);
			Vector_3 direction(s, t);
			Vector_3 step = direction/(1.0*n);
			Transformation_3 forward(CGAL::TRANSLATION, step);
			Point_3 csp(forward(s));
			Vector_3 forwardDir(csp, t);
			double squaredRadius = segLen[i]*segLen[i]/(4.0*n*n);
			for (unsigned int j = 0; j < n-1 && forwardDir*direction > 0; j++, csp = forward(csp))
			{
				Weighted_point wp(csp, squaredRadius);
				msu.push_back(wp);
				forwardDir = Vector_3(csp, t);
				assert(squared_distance(csp, t) > 1.0e-8);
			}
			msu.push_back(Weighted_point(s, squaredRadius));
		}
	}
	void SphericalBounder::genMarginBounder(double sf)
	{
		if (sf <= 0.0 || sf > 1.0)
			return;
		if (msu.size()>0)
			return;
		//shrink the polygons first
		vec3 midp(0.0, 0.0, 0.0);
		for (unsigned int i = 0; i < pgn.size(); i++)
			midp += to_geex(pgn[i].source());
		midp /= pgn.size();
		Transformation_3 move2origin(CGAL::TRANSLATION, Vector_3(-midp.x, -midp.y, -midp.z));
		Transformation_3 shrink(CGAL::SCALING, sf);
		Transformation_3 moveBack(CGAL::TRANSLATION, Vector_3(midp.x, midp.y, midp.z));
		Transformation_3 t = moveBack*(shrink*move2origin);
		vector<Segment_3> originalPolygon(pgn.size());
		for (unsigned int i = 0; i < pgn.size(); i++)
		{
			originalPolygon[i] = pgn[i];
			pgn[i] = pgn[i].transform(t);
		}
		//////////////////////// Sampling starts /////////////////////
		//msu = new vector<Weighted_point>;
		//double circumLen(0.0);
		vector<double> segLen; //length of each segment
		for (vector<Segment_3>::const_iterator it = pgn.begin(); it != pgn.end(); it++)
		{
			segLen.push_back(sqrt(it->squared_length()));
			//circumLen += segLen.back();
		}
		for (unsigned int i = 0; i < pgn.size(); i++)
		{
			Point_3 s = pgn[i].source(), t = pgn[i].target();
			double squaredRadius = CGAL::squared_distance(s, originalPolygon[i]);
			double radius = CGAL::sqrt(squaredRadius);
			unsigned int n = static_cast<unsigned int>(std::ceil(segLen[i]/(2.0*radius))); //number of sampling on this segment
			n = std::max<unsigned int>(1, n);
			Vector_3 direction(s, t);
			Vector_3 step = direction/(n-1)*0.5;
			Vector_3 step2 = direction/(n-1);
			Transformation_3 forward(CGAL::TRANSLATION, step2);
			Point_3 csp(s+step);		
			for (unsigned int j = 0; j < n-1; j++, csp = forward(csp))
			{
				Weighted_point wp(csp, squaredRadius);
				msu.push_back(wp);
				//assert(squared_distance(csp, t) > 1.0e-8);
			}
			msu.push_back(Weighted_point(s, squaredRadius));
		}
	}
	void SphericalBounder::genMarginBounder(double sf, const unsigned int nbSampling /* = 50 */)
	{
		if (msu.size() > 0)
			return;
		if (sf <=0.0 || sf >=1.0)
			return;
		//shrink the polygons first
		vec3 midp(0.0, 0.0, 0.0);
		for (unsigned int i = 0; i < pgn.size(); i++)
			midp += to_geex(pgn[i].source());
		midp /= pgn.size();
		Transformation_3 move2origin(CGAL::TRANSLATION, Vector_3(-midp.x, -midp.y, -midp.z));
		Transformation_3 shrink(CGAL::SCALING, sf);
		Transformation_3 moveBack(CGAL::TRANSLATION, Vector_3(midp.x, midp.y, midp.z));
		Transformation_3 t = moveBack*(shrink*move2origin);
		for (unsigned int i = 0; i < pgn.size(); i++)
			pgn[i] = pgn[i].transform(t);
		genMarginBounder(nbSampling);
	}
	void SphericalBounder::genInteriorMarginBounder(unsigned int nbSampling, const vec3& norm)
	{
		if (msu.size()>0)
			return;
		Vector_3 vc(0.0, 0.0, 0.0);
		for (unsigned int i = 0; i < pgn.size(); i++)
			vc = vc + (pgn[i].source() - CGAL::ORIGIN);
		Point_3 o = CGAL::ORIGIN + vc/pgn.size();
		// get the polygon normal
		Vector_3 uvz(norm.x, norm.y, norm.z);
		uvz = normalize(uvz);
		Vector_3 uvx = normalize(pgn[0].opposite().to_vector());
		Vector_3 uvy = CGAL::cross_product(uvz, uvx);
		Polygon_2 pgn_2;
		for (unsigned int i = 0; i < pgn.size(); i++)
			pgn_2.push_back(to_uv(uvx, uvy, o, pgn[i].source()));

		//compute offset distance
		double circumLen(0.0);
		vector<double> segLen(pgn.size()); //length of each segment
		for (unsigned int i = 0; i < pgn.size(); i++)
		{
			segLen[i] = sqrt(pgn[i].squared_length());
			circumLen += segLen[i];
		}
		double offset = circumLen/nbSampling*0.5;
		if (!pgn_2.is_simple())
		{
			cout<<"In spherical bounder, polygon not simple!!!\n";
			system("pause");
		}
		if (pgn_2.is_clockwise_oriented())
			pgn_2.reverse_orientation();
#ifdef CILK
		mut.lock();
#endif
		PolygonPtrVector pptr = CGAL::create_interior_skeleton_and_offset_polygons_2(offset, pgn_2); 
		//while (pptr[0]->size() != pgn_2.size())
		while (pptr.size() != 1)
		{
			offset = 0.95*offset;
			pptr = CGAL::create_interior_skeleton_and_offset_polygons_2(offset, pgn_2);
		}

#ifdef CILK
		mut.unlock();
#endif
		double squaredRadius = offset*offset;//radius of sampling spheres
		circumLen = 0.0;
		segLen.resize(pptr[0]->size());
		pgn.resize(pptr[0]->size());
		for (unsigned int i = 0; i < pptr[0]->size(); i++)
		{
			Point_2 src = pptr[0]->edge(i).source(), tgt = pptr[0]->edge(i).target();
			// preserve the sampling sphere on the original polygonal vertices
			Vector_2 v(pgn_2.vertex(i), pptr[0]->vertex(i));
			v = normalize(v);
			Point_2 vert = pgn_2.vertex(i) + v * offset;
			Point_3 pos = to_xy(uvx, uvy, o, vert);
			msu.push_back(Weighted_point(pos, squaredRadius));
			pgn[i] = Segment_3(to_xy(uvx, uvy, o, src), to_xy(uvx, uvy, o, tgt));
			segLen[i] = sqrt(pgn[i].squared_length());
			circumLen += segLen[i];
		}
		nbSampling = circumLen/offset*0.5;
		for (unsigned int i = 0; i < pgn.size(); i++)
		{
			Point_3 src = pgn[i].source(), tgt = pgn[i].target();
			//double approx_n = segLen[i]/offset*0.5;
			double approx_n = segLen[i]/circumLen*nbSampling;
			if (approx_n < 0.5)
				continue;
			unsigned int n = static_cast<unsigned int>(approx_n+0.5); //number of sampling on this segment
			n = std::max<unsigned int>(1, n);
			Vector_3 direction(src, tgt);
			Vector_3 step = direction/(1.0*n);
			Transformation_3 forward(CGAL::TRANSLATION, step);
			Point_3 csp(forward(src));
			Vector_3 forwardDir(csp, tgt);
			if (n > 1)
			{
				for (unsigned int j = 0; j < n; j++)
				{
					csp = CGAL::ORIGIN + ((n - j)*(src - CGAL::ORIGIN) + j*(tgt-CGAL::ORIGIN))/n;
					msu.push_back(Weighted_point(csp, squaredRadius));
				}
			}
			else			
				msu.push_back(Weighted_point(src+direction*0.5, squaredRadius));
		}
	}
	void SphericalBounder::genInteriorMarginBounder(unsigned int nbSampling)
	{
		if (msu.size()>0)
			return;
		//msu = new vector<Weighted_point>;
		Vector_3 vc(0.0, 0.0, 0.0);
		for (unsigned int i = 0; i < pgn.size(); i++)
			vc = vc + (pgn[i].source() - CGAL::ORIGIN);
		Point_3 o = CGAL::ORIGIN + vc/pgn.size();
		// get the polygon normal
		Vector_3 ab = pgn[0].opposite().to_vector(), ac = pgn[1].to_vector();
		Vector_3 uvz = CGAL::cross_product(ab, ac);
		uvz = normalize(uvz);
		Vector_3 uvx = normalize(ac);
		Vector_3 uvy = CGAL::cross_product(uvz, uvx);
		Polygon_2 pgn_2;
		for (unsigned int i = 0; i < pgn.size(); i++)
			pgn_2.push_back(to_uv(uvx, uvy, o, pgn[i].source()));
		
		//compute offset distance
		double circumLen(0.0);
		vector<double> segLen(pgn.size()); //length of each segment
		for (unsigned int i = 0; i < pgn.size(); i++)
		{
			segLen[i] = sqrt(pgn[i].squared_length());
			circumLen += segLen[i];
		}
		double offset = circumLen/nbSampling*0.5;
		if (!pgn_2.is_simple())
		{
			cout<<"In spherical bounder, polygon not simple!!!\n";
			system("pause");
		}
		if (pgn_2.is_clockwise_oriented())
			pgn_2.reverse_orientation();
#ifdef CILK
		mut.lock();
#endif
		PolygonPtrVector pptr = CGAL::create_interior_skeleton_and_offset_polygons_2(offset, pgn_2); 
		//while (pptr[0]->size() != pgn_2.size())
		while (pptr.size() != 1)
		{
			offset = 0.95*offset;
			pptr = CGAL::create_interior_skeleton_and_offset_polygons_2(offset, pgn_2);
		}

#ifdef CILK
		mut.unlock();
#endif
		double squaredRadius = offset*offset;//radius of sampling spheres
		circumLen = 0.0;
		segLen.resize(pptr[0]->size());
		pgn.resize(pptr[0]->size());
		for (unsigned int i = 0; i < pptr[0]->size(); i++)
		{
			//Point_3 src(pptr[0]->edge(i).source().x(), pptr[0]->edge(i).source().y(), h),
			//	tgt(pptr[0]->edge(i).target().x(), pptr[0]->edge(i).target().y(), h);
			Point_2 src = pptr[0]->edge(i).source(), tgt = pptr[0]->edge(i).target();
			// preserve the sampling sphere on the original polygonal vertices
			Vector_2 v(pgn_2.vertex(i), pptr[0]->vertex(i));
			v = normalize(v);
			Point_2 vert = pgn_2.vertex(i) + v * offset;
			//msu.push_back(Weighted_point(invt(Point_3(vert.x(), vert.y(), h)), squaredRadius));
			Point_3 pos = to_xy(uvx, uvy, o, vert);
			msu.push_back(Weighted_point(pos, squaredRadius));
			//tmp.push_back(Weighted_point(invt(Point_3(vert.x(), vert.y(), h)), squaredRadius));
			pgn[i] = Segment_3(to_xy(uvx, uvy, o, src), to_xy(uvx, uvy, o, tgt));
			segLen[i] = sqrt(pgn[i].squared_length());
			circumLen += segLen[i];
		}
		nbSampling = circumLen/offset*0.5;
		for (unsigned int i = 0; i < pgn.size(); i++)
		{
			Point_3 src = pgn[i].source(), tgt = pgn[i].target();
			//double approx_n = segLen[i]/offset*0.5;
			double approx_n = segLen[i]/circumLen*nbSampling;
			if (approx_n < 0.5)
				continue;
			unsigned int n = static_cast<unsigned int>(approx_n+0.5); //number of sampling on this segment
			n = std::max<unsigned int>(1, n);
			Vector_3 direction(src, tgt);
			Vector_3 step = direction/(1.0*n);
			Transformation_3 forward(CGAL::TRANSLATION, step);
			Point_3 csp(forward(src));
			Vector_3 forwardDir(csp, tgt);
			//Point_3 csp = s + step/2.0;
			//double squaredRadius = segLen[i]*segLen[i]/(4.0*n*n);		
			//msu->push_back(tmp[i]);
			if (n > 1)
			{
				//for (unsigned int j = 0; j < n-1 && forwardDir*direction > 0; j++, csp = forward(csp), forwardDir = Vector_3(csp, tgt))
				//{
				//	Weighted_point wp(csp, squaredRadius);
				//	msu.push_back(wp);
				//	//assert(squared_distance(csp, t) > 1.0e-8);
				//}
				//msu.push_back(Weighted_point(src, squaredRadius));
				for (unsigned int j = 0; j < n; j++)
				{
					csp = CGAL::ORIGIN + ((n - j)*(src - CGAL::ORIGIN) + j*(tgt-CGAL::ORIGIN))/n;
					msu.push_back(Weighted_point(csp, squaredRadius));
				}
			}
			else			
				msu.push_back(Weighted_point(src+direction*0.5, squaredRadius));
		}
	}

	vector<Weighted_point>* SphericalBounder::genInteriorMarginBounder(const vector<Segment_3>& plygn_, unsigned int nbSampling )
	{
		vector<Segment_3> plygn = plygn_;
		vector<Weighted_point> *bds = new vector<Weighted_point>;
		Vector_3 vc(0.0, 0.0, 0.0);
		for (unsigned int i = 0; i < plygn.size(); i++)
			vc = vc + (plygn[i].source() - CGAL::ORIGIN);
		Point_3 o = CGAL::ORIGIN + vc/plygn.size();
		// get the polygon normal
		Vector_3 ab = plygn[0].opposite().to_vector(), ac = plygn[1].to_vector();
		Vector_3 uvz = CGAL::cross_product(ab, ac);
		uvz = uvz / sqrt(uvz.squared_length());
		Vector_3 uvx = ac / sqrt(ac.squared_length());
		Vector_3 uvy = CGAL::cross_product(uvz, uvx);
		Polygon_2 pgn_2;
		for (unsigned int i = 0; i < plygn.size(); i++)
			pgn_2.push_back(to_uv(uvx, uvy, o, plygn[i].source()));

		//compute offset distance
		double circumLen(0.0);
		vector<double> segLen(plygn.size()); //length of each segment
		for (unsigned int i = 0; i < plygn.size(); i++)
		{
			segLen[i] = sqrt(plygn[i].squared_length());
			circumLen += segLen[i];
		}
		double offset = circumLen/nbSampling;
		if (!pgn_2.is_simple())
		{
			cout<<"In spherical bounder, polygon not simple!!!\n";
			system("pause");
		}
		if (pgn_2.is_clockwise_oriented())
			pgn_2.reverse_orientation();
#ifdef CILK
		mut.lock();
#endif
		//SsPrt ss = CGAL::create_interior_straight_skeleton_2(pgn_2);
		PolygonPtrVector pptr = CGAL::create_interior_skeleton_and_offset_polygons_2(offset, pgn_2); 
		//PolygonPtrVector pptr = CGAL::create_offset_polygons_2<Polygon_2>(offset, *ss);
		//Polygon_2 offsetPgn = *pptr[0];
		while (pptr[0]->size() != pgn_2.size())
		{
			offset = 0.95*offset;
			pptr = CGAL::create_interior_skeleton_and_offset_polygons_2(offset, pgn_2);
		}
#ifdef CILK
		mut.unlock();
#endif
		double squaredRadius = offset*offset;//radius of sampling spheres
		circumLen = 0.0;
		for (unsigned int i = 0; i < pptr[0]->size(); i++)
		{
			//Point_3 src(pptr[0]->edge(i).source().x(), pptr[0]->edge(i).source().y(), h),
			//	tgt(pptr[0]->edge(i).target().x(), pptr[0]->edge(i).target().y(), h);
			Point_2 src = pptr[0]->edge(i).source(), tgt = pptr[0]->edge(i).target();
			// preserve the sampling sphere on the original polygonal vertices
			Vector_2 v(pgn_2.vertex(i), pptr[0]->vertex(i));
			v = normalize(v);
			Point_2 vert = pgn_2.vertex(i) + v * offset;
			//bds.push_back(Weighted_point(invt(Point_3(vert.x(), vert.y(), h)), squaredRadius));
			Point_3 pos = to_xy(uvx, uvy, o, vert);
			bds->push_back(Weighted_point(pos, squaredRadius));
			//tmp.push_back(Weighted_point(invt(Point_3(vert.x(), vert.y(), h)), squaredRadius));
			plygn[i] = Segment_3(to_xy(uvx, uvy, o, src), to_xy(uvx, uvy, o, tgt));
			segLen[i] = sqrt(plygn[i].squared_length());
			circumLen += segLen[i];
		}
		for (unsigned int i = 0; i < plygn.size(); i++)
		{
			Point_3 src = plygn[i].source(), tgt = plygn[i].target();
			double approx_n = segLen[i]/offset*0.5;
			if (approx_n < 0.5)
				continue;
			unsigned int n = static_cast<unsigned int>(approx_n); //number of sampling on this segment
			n = std::max<unsigned int>(1, n);
			Vector_3 direction(src, tgt);
			Vector_3 step = direction/(1.0*n);
			Transformation_3 forward(CGAL::TRANSLATION, step);
			Point_3 csp(forward(src));
			Vector_3 forwardDir(csp, tgt);
			//Point_3 csp = s + step/2.0;
			//double squaredRadius = segLen[i]*segLen[i]/(4.0*n*n);		
			//bds->push_back(tmp[i]);
			if (n > 1)
			{
				for (unsigned int j = 0; j < n-1 && forwardDir*direction > 0; j++, csp = forward(csp), forwardDir = Vector_3(csp, tgt))
				{
					Weighted_point wp(csp, squaredRadius);
					bds->push_back(wp);
					//assert(squared_distance(csp, t) > 1.0e-8);
				}
				bds->push_back(Weighted_point(src, squaredRadius));
			}
			else			
				bds->push_back(Weighted_point(src+direction*0.5, squaredRadius));
		}	
		return bds;
	}

}