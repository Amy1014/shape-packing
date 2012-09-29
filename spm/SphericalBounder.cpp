#include "SphericalBounder.h"

namespace Geex
{
	SphericalBounder::SphericalBounder(const vector<Segment_3>& aPolygon)
	{
		pgn = aPolygon;
		su = 0;
		msu = 0;
	}
	SphericalBounder::~SphericalBounder()
	{
		if (su)
			su->clear();
		if (msu)
			msu->clear();
	}
	void SphericalBounder::genMarginBounder(const unsigned int nbSamping /* = 50 */)
	{
		msu = new vector<Weighted_point>;
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
			//Point_3 csp(s);//sampling start
			//for (unsigned int j = 0; j < n; j++, csp = forward(csp))
			//{
			//	Weighted_point wp(csp, segLen[i]*segLen[i]/(4.0*n*n));
			//	msu->push_back(wp);
			//}
			Point_3 csp(forward(s));
			double squaredRadius = segLen[i]*segLen[i]/(4.0*n*n);
			//double squaredRadius = 0.00005;
			//std::cout<<"========================\n";
			for (unsigned int j = 0; j < n-1; j++, csp = forward(csp))
			{
				Weighted_point wp(csp, squaredRadius);
				//std::cout<<csp<<std::endl;
				msu->push_back(wp);
				assert(squared_distance(csp, t) > 1.0e-8);
			}
			msu->push_back(Weighted_point(s, squaredRadius));
			//std::cout<<s<<std::endl;
		}
		//cout<<"Number of bounding spheres: "<< msu->size()<<endl;
		//return *msu;
	}

}