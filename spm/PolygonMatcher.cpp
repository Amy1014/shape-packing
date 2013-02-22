#include "PolygonMatcher.h"

namespace Geex
{
	// compute f(x)+s
	PiecewiseConstFunc operator+(const PiecewiseConstFunc& func, double s)
	{
		PiecewiseConstFunc newFunc = func;
		for (unsigned int i = 0; i < newFunc.size(); i++)
			newFunc[i].f += s;
		return newFunc;
	}
	PiecewiseConstFunc operator+(const PiecewiseConstFunc& f, const PiecewiseConstFunc& g)
	{
		PiecewiseConstFunc h;
		double prePoint = f[0].l, nextPoint;//assume f and g starts from the same point
		int pf = 1, pg = 1;
		// the left and right value at the current breakpoint on f and g
		double lf = f.leftValueAt(pf), lg = g.leftValueAt(pg);
		while (pf < f.size() && pg < g.size())
		{
			if (f[pf].l < g[pg].l)
			{
				nextPoint = f[pf].l;
				h.push_back(ConstPiece(lf+lg, prePoint, nextPoint));
				pf++;
				lf = f.leftValueAt(pf);
			}
			else if (f[pf].l > g[pg].l)
			{
				nextPoint = g[pg].l;
				h.push_back(ConstPiece(lf+lg, prePoint, nextPoint));
				pg++;
				lg = g.leftValueAt(pg);
			}
			else //equal
			{
				nextPoint = f[pf].l;
				h.push_back(ConstPiece(lf+lg, prePoint, nextPoint));
				pg++;
				lg = g.leftValueAt(pg);
				pf++;
				lf = f.leftValueAt(pf);
			}
			prePoint = nextPoint;
		}
		for (int i = pf; i < f.size(); i++)
		{
			nextPoint = f[i].l;
			h.push_back(ConstPiece(lf+lg, prePoint, nextPoint));
			lf = f.leftValueAt(i+1);
			prePoint = nextPoint;
		}
		for (int i = pg; i < g.size(); i++)
		{
			nextPoint = g[i].l;
			h.push_back(ConstPiece(lf+lg, prePoint, nextPoint));
			lg = g.leftValueAt(i+1);
			prePoint = nextPoint;
		}
		h.push_back(ConstPiece(f.back().f+g.back().f, prePoint, f.back().u));
		return h;
	}
	PiecewiseConstFunc operator-(const PiecewiseConstFunc& f, const PiecewiseConstFunc& g)
	{
		PiecewiseConstFunc h;
		double prePoint = f[0].l, nextPoint;//assume f and g starts from the same point
		int pf = 1, pg = 1;
		// the left and right value at the current breakpoint on f and g
		double lf = f.leftValueAt(pf), lg = g.leftValueAt(pg);
		while (pf < f.size() && pg < g.size())
		{
			if (f[pf].l < g[pg].l)
			{
				nextPoint = f[pf].l;
				h.push_back(ConstPiece(lf-lg, prePoint, nextPoint));
				pf++;
				lf = f.leftValueAt(pf);
			}
			else if (f[pf].l > g[pg].l)
			{
				nextPoint = g[pg].l;
				h.push_back(ConstPiece(lf-lg, prePoint, nextPoint));
				pg++;
				lg = g.leftValueAt(pg);
			}
			else //equal
			{
				nextPoint = f[pf].l;
				h.push_back(ConstPiece(lf-lg, prePoint, nextPoint));
				pg++;
				lg = g.leftValueAt(pg);
				pf++;
				lf = f.leftValueAt(pf);
			}
			prePoint = nextPoint;
		}
		for (int i = pf; i < f.size(); i++)
		{
			nextPoint = f[i].l;
			h.push_back(ConstPiece(lf-lg, prePoint, nextPoint));
			lf = f.leftValueAt(i+1);
			prePoint = nextPoint;
		}
		for (int i = pg; i < g.size(); i++)
		{
			nextPoint = g[i].l;
			h.push_back(ConstPiece(lf-lg, prePoint, nextPoint));
			lg = g.leftValueAt(i+1);
			prePoint = nextPoint;
		}
		h.push_back(ConstPiece(f.back().f-g.back().f, prePoint, f.back().u));
		return h;
	}
	PolygonMatcher::PolygonMatcher(const Polygon_2& model)
	{
		// scale the polygons so that the perimeter is equal to 1
		double perimeter = 0.0;
		for (Polygon_2::Edge_const_iterator eit = model.edges_begin(); eit != model.edges_end(); eit++)
			perimeter += sqrt(eit->squared_length());
		// scale transformation
		Transformation_2 scaleTrans(CGAL::SCALING, 1.0/perimeter);
		this->model = transform(scaleTrans, model);
		if (this->model.is_clockwise_oriented())
			this->model.reverse_orientation();
		constructTurningFunc(this->model, &modelTurningFunc);
		modelInt = modelTurningFunc.integral();
	}
	PolygonMatcher::PolygonMatcher(const vector<Segment_2>& model)
	{
		Vector_2 vecc(0.0, 0.0);
		unsigned int sn = 20;//sampling number
		// affine match point cloud
		for (unsigned int i = 0; i < model.size(); i++)
		{
			//modelPtsSet.push_back(cv::Point2d(model[i].source().x(), model[i].source().y()));	
			for (unsigned int j = 0; j < sn; j++)
			{
				Point_2 p = CGAL::ORIGIN+((sn-j)*(model[i].source()-CGAL::ORIGIN)+j*(model[i].target()-CGAL::ORIGIN))/sn;
				modelPtsSet.push_back(cv::Point2d(p.x(), p.y()));
			}
			vecc = vecc + (model[i].source() - CGAL::ORIGIN);
		}
		vecc = vecc / model.size();
		Point_2 c = CGAL::ORIGIN + vecc;
		// compute area of model
		modelArea = 0.0;
		for (unsigned int i = 0; i < model.size(); i++)
		{
			Triangle_2 tri(c, model[i].source(), model[i].target());
			modelArea += std::fabs(tri.area());
		}
		modelMat = CmAffine::CreatePointSetMat(modelPtsSet);
	}
	PolygonMatcher::PolygonMatcher(const vector<Point_2>& model)
	{
		for (unsigned int i = 0; i < model.size(); i++)
			modelPtsSet.push_back(cv::Point2d(model[i].x(), model[i].y()));
		// approximate the area according to convex hull
		//Polygon_2 ch;
		//CGAL::convex_hull_2(model.begin(), model.end(), std::back_inserter(ch));
		//double modelArea = std::fabs(ch.area());
		//modelMat = CmAffine::CreatePointSetMat(modelPtsSet);
		// approximate the area according to bounding box
		Bounding_box_2 bd = CGAL::bounding_box(model.begin(), model.end());
		modelArea = std::fabs(bd.area());
		modelMat = CmAffine::CreatePointSetMat(modelPtsSet);
	}
	PolygonMatcher::~PolygonMatcher()
	{
		if (!modelMat)
			cvReleaseMat(&modelMat);
	}
	inline double PolygonMatcher::radian(double x, double y)
	{
		double cs = x / std::sqrt(x*x+y*y);
		double angle = std::acos(cs);
		return (y>=0)?angle:(angle + PI);
	}
	double PolygonMatcher::insert(const Polygon_2& pgn)
	{
		// scale first
		double perimeter = 0.0;
		for (Polygon_2::Edge_const_iterator eit = pgn.edges_begin(); eit != pgn.edges_end(); eit++)
			perimeter += sqrt(eit->squared_length());
		Transformation_2 scaleTrans(CGAL::SCALING, 1.0/perimeter);
		Polygon_2 sclPgn = transform(scaleTrans, pgn);
		if (sclPgn.is_clockwise_oriented())
			sclPgn.reverse_orientation();
		//instances.push_back(sclPgn);
		///////////////////////////// Deal with the polygon without no reflection /////////////////////////////
		// get the turning function
		PiecewiseConstFunc tf;
		constructTurningFunc(sclPgn, &tf);
		double noReflectAngle = 0.0;
		double noReflectDist = computeSimilarity(tf, &noReflectAngle);
		if (noReflectDist <= 1.0e-4)
		{
			similarityDistances.push_back(noReflectDist);
			matchAngles.push_back(noReflectAngle);
			yAxisReflected.push_back(false);
			return noReflectDist;
		}
		///////////////////////////// Deal with the polygon flipped /////////////////////////////
		Transformation_2 flip(-1.0,	0.0, 0.0, 1.0);
		Polygon_2 flippedPolygon = transform(flip, sclPgn);
		if (flippedPolygon.is_clockwise_oriented())
			flippedPolygon.reverse_orientation();
		PiecewiseConstFunc rtf;
		constructTurningFunc(flippedPolygon, &rtf);
		double reflectAngle = 0.0;
		double refectDist = computeSimilarity(rtf, &reflectAngle);
		////////////////////////////// Compare /////////////////////////////////
		if (noReflectDist < refectDist)
		{
			similarityDistances.push_back(noReflectDist);
			matchAngles.push_back(noReflectAngle);
			yAxisReflected.push_back(false);
			return noReflectDist;
		}
		else
		{
			similarityDistances.push_back(refectDist);
			matchAngles.push_back(reflectAngle);
			yAxisReflected.push_back(true);
			return refectDist;
		}
	}
	double PolygonMatcher::affineMatch(const Polygon_2& instance, double *scale, double *theta, double *tx, double *ty)
	{
		// sample the input polygon
		PointSet instancePtsSet;
		const unsigned int sn = 30;//sampling number	
		// compute perimeter
		vector<double> edgePerimeter(instance.size());
		double perimeter = 0.0;
		for (unsigned int i = 0; i < instance.size(); i++)
		{
			edgePerimeter[i] = sqrt(instance.edge(i).squared_length());
			perimeter += edgePerimeter[i];
		}
		for (unsigned int i = 0; i < instance.size(); i++)
		{
			Point_2 s = instance.edge(i).source();
			Point_2 t = instance.edge(i).target();
			unsigned int n(edgePerimeter[i]/perimeter*sn);
			n += 1;//the source end point must be sampled
			for (unsigned int j = 0; j < n; j++)
			{
				Point_2 p = CGAL::ORIGIN + ( (n-j)*(s - CGAL::ORIGIN) + j*(t - CGAL::ORIGIN) )/n;
				instancePtsSet.push_back(cv::Point2d(p.x(), p.y()));
			}
		}
		CvMat *instanceMat = CmAffine::CreatePointSetMat(instancePtsSet);
		CvMat *AA = ::cvCreateMat(2, 2, CV_64FC1);
		CvMat *shift = ::cvCreateMat(2, 1, CV_64FC1);
		*theta = CmAffine::GetAffine2D(modelMat, instanceMat, AA, shift);
		double s = std::sqrt(modelArea/std::fabs(instance.area()));
		*scale = s;
		double *r = AA->data.db;
		r[0] = s*std::cos(*theta);
		r[1] = -s*std::sin(*theta);
		r[2] = -r[1];
		r[3] = r[0];
		*tx = shift->data.db[0];
		*ty = shift->data.db[1];
		CvMat *pMat = cvCreateMat(2, instancePtsSet.size(), CV_64FC1);
		CmAffine::FindTranslated(instanceMat, AA, shift, pMat);
		instancePtsSet.clear();
		double *x = pMat->data.db;
		double *y = (double*)(pMat->data.ptr + pMat->step);
		for (unsigned int i = 0; i < pMat->cols; i++)
			instancePtsSet.push_back(cv::Point2d(x[i], y[i]));
		// compute distance
		double avgDist = 0.0;
		for (unsigned int i = 0; i < instancePtsSet.size(); i++)
		{
			double minDist = DBL_MAX;
			for (unsigned int j = 0; j < modelPtsSet.size(); j++)
			{
				double tmp = sqrDist(modelPtsSet[j], instancePtsSet[i]);
				if (tmp < minDist)
					minDist = tmp;
			}
			avgDist += sqrt(minDist);
		}
		avgDist /= instancePtsSet.size();
		cvReleaseMat(&instanceMat);
		cvReleaseMat(&shift);
		cvReleaseMat(&AA);
		cvReleaseMat(&pMat);
		return avgDist;
	}
	// compute the turning angle function with x-axis, without any shift or rotation
	void PolygonMatcher::constructTurningFunc(const Polygon_2& pgn, PiecewiseConstFunc *func)
	{
		func->clear();
		double accumArcLen = 0.0, accumAngle = 0.0;
		Vector_2 preDir(1.0, 0.0);
		for (unsigned int i = 0; i < pgn.size(); i++)
		{
			double segLen = sqrt(pgn.edge(i).squared_length());
			Vector_2 curDir = pgn.edge(i).to_vector();
			double cs = curDir*preDir / (sqrt(curDir.squared_length())*sqrt(preDir.squared_length()));
			//assert(cs <=1.0 && cs>=-1.0);
			cs = std::max<double>(cs, -1.0);
			cs = std::min<double>(cs, 1.0);
			if (orientation(preDir, curDir) == CGAL::LEFT_TURN)
				accumAngle += std::acos(cs);
			else
				accumAngle -= std::acos(cs);
			ConstPiece piece(accumAngle, accumArcLen, segLen+accumArcLen);
			accumArcLen += segLen;
			func->push_back(piece);
			preDir = curDir;
		}
	}
	double PolygonMatcher::computeSimilarity(const PiecewiseConstFunc& instanceTurningFunc, double *matchAngle)
	{
		// compute the minimizer shrink t, with rotation angle equal to 0
		// construct the breakpoint of f(x) and g(x)
		vector<triplet> modelBrkPnts(modelTurningFunc.size());
		vector<triplet> insBkrPnts(instanceTurningFunc.size());
		for (unsigned int i = 0; i < modelBrkPnts.size(); i++)
			modelBrkPnts[i].val = modelTurningFunc[i].l;
		for (unsigned int i = 0; i < insBkrPnts.size(); i++)
			insBkrPnts[i].val = instanceTurningFunc[i].l;
		// sort all the breakpoints
		int pm = 0, pi = 0;//pointers to constant pieces in the two functions f(x) and g(x)
		bool preFromModel = true;//the previous one is from the model turning function
		while ( pm < modelBrkPnts.size() && pi < insBkrPnts.size())
		{
			if (modelBrkPnts[pm] <= insBkrPnts[pi])
			{
				if (preFromModel)
					modelBrkPnts[pm].left = -(pm+1); // minus the index of the previous breakpoint
				else
					modelBrkPnts[pm].left = pi-1;
				preFromModel = true;
				pm++;
			}
			else
			{
				if (!preFromModel)
					insBkrPnts[pi].left = -(pi+1); // minus the index of the previous breakpoint
				else
					insBkrPnts[pi].left = pm-1;
				preFromModel = false;
				pi++;
			}
		}
		// handle the remaining items
		for (int i = pm; i < modelBrkPnts.size(); i++)
		{
			if (!preFromModel)
				modelBrkPnts[i].left = pi-1;
			else
				modelBrkPnts[i].left = -(i+1);
			preFromModel = true;
		}
		for (int i = pi; i < insBkrPnts.size(); i++)
		{
			if (preFromModel)
				insBkrPnts[i].left = pm-1;
			else
				insBkrPnts[i].left = -(i+1);
			preFromModel = false;
		}
		// compute crucial events where breakpoints coincide
		vector<triplet> crucialEvents(modelBrkPnts.size()*insBkrPnts.size()+1);
		unsigned int idx = 0;
		crucialEvents[idx++] = triplet(0.0, -1, -1);
		for (unsigned int i = 0; i < modelBrkPnts.size(); i++)
			for (unsigned int j = 0; j < insBkrPnts.size(); j++)
			{
				if (modelBrkPnts[i].val >= insBkrPnts[j].val)
					crucialEvents[idx++] = triplet(modelBrkPnts[i].val - insBkrPnts[j].val, i, j);
				else
					crucialEvents[idx++] = triplet(1.0 + modelBrkPnts[i].val - insBkrPnts[j].val, i + modelBrkPnts.size(), j);
			}
		// determine the order in which crucial events happen
		std::sort(crucialEvents.begin(), crucialEvents.end(), cmpTriplet);
		// compute the squared height of those strips whose left and right are different breakpoints
		//correspond to strips whose left is from model polygon function and right from instance function, and vice versa
		double hmi = 0.0, him = 0.0;
		double sqhmi = 0.0, sqhim = 0.0;
		for (unsigned int i = 0; i < insBkrPnts.size(); i++)
			if (insBkrPnts[i].left >= 0)
			{
				double height = modelTurningFunc.rightValueAt(insBkrPnts[i].left) - instanceTurningFunc.leftValueAt(i);
				sqhmi += height*height;
				hmi += height;
			}
		for (unsigned int i = 0; i < modelBrkPnts.size(); i++)
			if ( modelBrkPnts[i].left >= 0)
			{
				double height = modelTurningFunc.leftValueAt(i) - instanceTurningFunc.rightValueAt(modelBrkPnts[i].left);
				sqhim += height*height;
				him += height;
			}
		// now compute the distance at each breakpoint
		double instanceInt = instanceTurningFunc.integral();
		double intDiff = modelInt - instanceInt;// int_0^1 (f(s) - g(s))ds
		PiecewiseConstFunc diffFunc = modelTurningFunc - instanceTurningFunc;
		diffFunc.square();
		double curDist = diffFunc.integral(); // int_0^1 (f(s)-g(s))^2 ds
		double minDist = curDist, minIntDiff = intDiff; // minimum distance so far
		double minShift = 0.0; // the shift minimizing the distance
		for (unsigned int i = 1; i < crucialEvents.size(); i++)
		{
			curDist += (sqhmi - sqhim)*(crucialEvents[i].val - crucialEvents[i-1].val);
			intDiff += (hmi - him)*(crucialEvents[i].val - crucialEvents[i-1].val);
			if (minDist > curDist)
			{
				minDist = curDist;
				minIntDiff = intDiff;
				minShift = crucialEvents[i].val;
			}
			// update strips and him and hmi
			int l = crucialEvents[i].left, r = crucialEvents[i].right;
			double height = modelTurningFunc.leftValueAt(l) - instanceTurningFunc.rightValueAt(r);
			sqhim -= height*height;
			him -= height;
			height = modelTurningFunc.rightValueAt(l) - instanceTurningFunc.leftValueAt(r);
			sqhmi += height*height;
			hmi += height;
			int ml = l%modelBrkPnts.size(), mr = r%insBkrPnts.size();//actual index in break point list
			int temp = insBkrPnts[mr].left;
			insBkrPnts[mr].left = l;
			modelBrkPnts[ml].left = -temp-1;
			if (modelBrkPnts[ml].left >= 0)
			{
				height = modelTurningFunc.leftValueAt(l) - instanceTurningFunc.rightValueAt(modelBrkPnts[ml].left);
				sqhim += height*height;
				him += height;
			}
			else if (modelBrkPnts[ml].left < 0)
			{
				height = modelTurningFunc.rightValueAt(-modelBrkPnts[ml].left-1) - instanceTurningFunc.leftValueAt(r);
				sqhmi -= height*height;
				hmi -= height;
			}
			if (insBkrPnts[(r+1)%insBkrPnts.size()].left == l)
			{
				height = modelTurningFunc.rightValueAt(l) - instanceTurningFunc.leftValueAt(r+1);
				sqhmi -= height*height;
				hmi -= height;
				insBkrPnts[(r+1)%insBkrPnts.size()].left = -(r+1);
			}
			else if (modelBrkPnts[(l+1)%modelBrkPnts.size()].left == -(l+1))
			{
				height = modelTurningFunc.leftValueAt(l+1) - instanceTurningFunc.rightValueAt(r);
				sqhim += height*height;
				him += height;
				modelBrkPnts[(l+1)%modelBrkPnts.size()].left = r;
			}
		}
		// compute the minimizer rotation angle
		double minRotAngle = instanceInt - modelInt - 2*PI*minShift;
		minDist = minDist + 2*minRotAngle*minIntDiff + minRotAngle*minRotAngle;
		*matchAngle = minRotAngle;
		return minDist;
	}
}