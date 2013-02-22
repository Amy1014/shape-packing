#include "Doc.h"
#include "SphericalBounder.h"
//#include <algorithm>

/* Knitro callback function that evaluates the objective and constraints */
int  callback (const int evalRequestCode,
			   const int n,
			   const int m,
			   const int nnzJ,
			   const int nnzH,
			   const double * const x,
			   const double * const lambda,
			   double * const obj,
			   double * const c,
			   double * const objGrad,
			   double * const jac,
			   double * const hessian,
			   double * const hessVector,
			   void * userParams);
namespace Geex
{
	Doc* Doc::instance_ = nil ;
	const double Doc::pi = 3.14159265358979323846;
	Doc::Doc()
	{
		gx_assert(instance_ == nil) ;
		instance_ = this ;
		minScalor = 0.3;
		maxScalor = 1.2;
		noEnlargeTimes = 0;
		areaCoverage = 0.8;
		packIterLim = 300;
		lloydUsed = &Doc::nobdLloyd;
		curvature = false;
		texture = false;
		backgroud = false;
		stopUpdate = false;
		noProgressTimes = 0;
		currentMaxAreaRatio = DBL_MIN;
	}
	Doc::~Doc() 
	{
		clear();
		instance_ = nil ; 
	}
	void Doc::loadPolygons(const string& pgnFileName)
	{
		pgnLib.clear();
		ifstream fin(pgnFileName.c_str(), std::ios_base::in);
		int nbPolygons = 0;
		fin>>nbPolygons;//total number of polygons
		for (int i = 0; i < nbPolygons; i++)
		{
			int nbEdges;
			fin>>nbEdges;
			Polygon_2 pgn;
			for (int j = 0; j < nbEdges; j++)
			{
				Polygon_2::FT x, y;
				fin>>x;
				fin>>y;
				pgn.push_back(Point_2(x, y));	
			}
			if (pgn.orientation() == CGAL::CLOCKWISE)
				pgn.reverse_orientation();
			pgnLib.push_back(pgn);
		}
		fin.close();
		// choose some number of polygons randomly
		if (meshDom.size() == 0)
		{
			std::cerr<<"Failure: No mesh loaded. Load a mesh before loading polygons!\n";
			system("pause");
			exit(1);
		}
		// compute the maximum polygon area
		maxPolygonArea = -1.0;
		double avgPolygonArea = 0.0;
		for (unsigned int i = 0; i < pgnLib.size(); i++)
		{
			double a = std::fabs(pgnLib[i].area());
			if ( a > maxPolygonArea )
				maxPolygonArea = a;
			avgPolygonArea += a;
		}
		avgPolygonArea /= pgnLib.size();
		std::cout<<"Maximum area of polygons in the library: "<<maxPolygonArea<<std::endl;
		int nbActive(meshArea/avgPolygonArea);
		//std::vector<Polygon_2> pgns;
		for (int i = 0; i < nbActive; i++)
		{
			int idx = i%pgnLib.size();
			const Polygon_2& p = pgnLib[idx];
			if (p.size()<4) //avoid triangles when initially loading
				continue;
			pgns2d.push_back(p);	
			pgnLibIdx.push_back(idx);
		}
		std::cout<<"Loaded "<<pgns2d.size()<<" polygons initially\n";
		genBoundaryPoints();
		placePolygons(pgns2d);
	}
	void Doc::loadPolygons(const string& pgnFileName, const string& restore_file)
	{
		pgnLib.clear();
		ifstream fin(pgnFileName.c_str(), std::ios_base::in);
		int nbPolygons = 0;
		fin>>nbPolygons;//total number of polygons
		pgnLib.reserve(nbPolygons);
		for (int i = 0; i < nbPolygons; i++)
		{
			int nbEdges;
			fin>>nbEdges;
			Polygon_2 pgn;
			for (int j = 0; j < nbEdges; j++)
			{
				Polygon_2::FT x, y;
				fin>>x;
				fin>>y;
				pgn.push_back(Point_2(x, y));	
			}
			if (pgn.orientation() == CGAL::CLOCKWISE)
				pgn.reverse_orientation();
			pgnLib.push_back(pgn);
		}
		fin.close();
		genBoundaryPoints();
		restore(restore_file);
		std::cout<<"Restore "<<pgns_3.size()<<" polygons initially\n";
		//std::cout<<"Area coverage: "<<areaCoverage<<std::endl;
	}
	//void Doc::loadPolygonsWithTexture(const std::vector<std::string>& filenames_)
	//{
	//	pgnLib.clear();
	//	texture = true;
	//	for (unsigned int i = 0; i < 10/*filenames_.size()*/; i++)
	//	{
	//		ifstream infile(filenames_[i].c_str(), std::ios_base::in);
	//		int nbVert;
	//		infile>>nbVert;
	//		pgnLib.push_back(Polygon_2());
	//		pgnTextCoordLib.push_back(Polygon_2());
	//		for (int j = 0; j < nbVert; j++)
	//		{
	//			double x, y, texx, texy;
	//			infile >> x;
	//			infile >> y;
	//			infile >> texx;
	//			infile >> texy;
	//			pgnLib.back().push_back(Point_2(x, y));
	//			pgnTextCoordLib.back().push_back(Point_2(texx, texy));
	//		}
	//		if (pgnLib.back().is_clockwise_oriented())
	//			pgnLib.back().reverse_orientation();
	//		infile.close();
	//	}
	//	maxPolygonArea = -1.0;
	//	double avgPolygonArea = 0.0;
	//	for (unsigned int i = 0; i < pgnLib.size(); i++)
	//	{
	//		double a = std::fabs(pgnLib[i].area());
	//		if ( a > maxPolygonArea )
	//			maxPolygonArea = a;
	//		avgPolygonArea += a;
	//	}
	//	avgPolygonArea /= pgnLib.size();
	//	std::cout<<"Maximum area of polygons in the library: "<<maxPolygonArea<<std::endl;
	//	int nbActive(meshArea/avgPolygonArea);
	//	for (int i = 0; i < nbActive; i++)
	//	{
	//		int idx = i%pgnLib.size();
	//		const Polygon_2& p = pgnLib[idx];
	//		if (p.size()<4) //avoid triangles when initially loading
	//			continue;
	//		pgns2d.push_back(p);	
	//		pgnLibIdx.push_back(idx);
	//		pgnTextCoord.push_back(pgnTextCoordLib[idx]);
	//	}
	//	std::cout<<"Loaded "<<pgns2d.size()<<" polygons initially\n";
	//	genBoundaryPoints();
	//	placePolygons(pgns2d);
	//	fbmasks.resize(pgns_3.size(), true);// treat all as foreground
	//}
	//void Doc::loadPolygonsWithTexture(const std::vector<std::string>& filenames_, const std::string& restore_file)
	//{
	//	pgnLib.clear();
	//	texture = true;
	//	for (unsigned int i = 0; i < 10/*filenames_.size()*/; i++)
	//	{
	//		ifstream infile(filenames_[i].c_str(), std::ios_base::in);
	//		int nbVert;
	//		infile>>nbVert;
	//		pgnLib.push_back(Polygon_2());
	//		pgnTextCoordLib.push_back(Polygon_2());
	//		for (int j = 0; j < nbVert; j++)
	//		{
	//			double x, y, texx, texy;
	//			infile >> x;
	//			infile >> y;
	//			infile >> texx;
	//			infile >> texy;
	//			pgnLib.back().push_back(Point_2(x, y));
	//			pgnTextCoordLib.back().push_back(Point_2(texx, texy));
	//		}
	//		if (pgnLib.back().is_clockwise_oriented())
	//			pgnLib.back().reverse_orientation();
	//		infile.close();
	//	}
	//	genBoundaryPoints();
	//	restore(restore_file);
	//	std::cout<<"Restore "<<pgns_3.size()<<" polygons\n";
	//	fbmasks.resize(pgns_3.size(), true);// treat all as foreground
	//}
	//void Doc::loadPolygonsWithBackgroud(const std::vector<string>& filenames_, const string& restore_file, 
	//									const string& bg_pgn_filename, const string& mc_pts_filename)
	//{
	//	pgnLib.clear();
	//	texture = true;
	//	double maxPgnArea = DBL_MIN;
	//	for (unsigned int i = 0; i < 10/*filenames_.size()*/; i++)
	//	{
	//		ifstream infile(filenames_[i].c_str(), std::ios_base::in);
	//		int nbVert;
	//		infile>>nbVert;
	//		pgnLib.push_back(Polygon_2());
	//		pgnTextCoordLib.push_back(Polygon_2());
	//		for (int j = 0; j < nbVert; j++)
	//		{
	//			double x, y, texx, texy;
	//			infile >> x;
	//			infile >> y;
	//			infile >> texx;
	//			infile >> texy;
	//			pgnLib.back().push_back(Point_2(x, y));
	//			pgnTextCoordLib.back().push_back(Point_2(texx, texy));
	//		}
	//		maxPolygonArea = std::max<double>(maxPolygonArea, std::fabs(pgnLib.back().area()));
	//		infile.close();
	//	}
	//	ifstream fin(bg_pgn_filename.c_str(), std::ios_base::in);
	//	int nbPolygons = 0;
	//	fin>>nbPolygons;//total number of polygons
	//	bgPgnLib.reserve(nbPolygons);
	//	bgPgnTextCoordLib.reserve(nbPolygons);
	//	for (int i = 0; i < nbPolygons; i++)
	//	{
	//		int nbVert;
	//		fin>>nbVert;
	//		Polygon_2 pgn;
	//		for (int j = 0; j < nbVert; j++)
	//		{
	//			Polygon_2::FT x, y;
	//			fin>>x;
	//			fin>>y;
	//			pgn.push_back(Point_2(x, y));	
	//		}
	//		bgPgnLib.push_back(pgn);
	//		maxPolygonArea = std::max<double>(maxPolygonArea, std::fabs(pgn.area()));
	//		// compute texture coordinates
	//		Bounding_box_2 bd = CGAL::bounding_box(pgn.vertices_begin(), pgn.vertices_end());
	//		double leftx = bd.xmin(), rightx = bd.xmax(), bottomy = bd.ymin(), topy = bd.ymax();
	//		double width = std::max<double>(rightx - leftx, topy - bottomy);
	//		Transformation_2 leftbottom2origin(CGAL::TRANSLATION, Vector_2(-leftx, -bottomy));
	//		Transformation_2 scl(CGAL::SCALING, 1.0/width);
	//		bgPgnTextCoordLib.push_back(transform(scl*leftbottom2origin, pgn));
	//	}
	//	std::cout<<"Maximum area of polygons "<<maxPolygonArea<<std::endl;
	//	fin.close();
	//	genBoundaryPoints();
	//	if (restore_file.empty())
	//		// distribute polygons at multi-class points, generate texture coordinates for background polygons
	//	{
	//		//read initial points
	//		double shrinkFactor = std::min<double>(maxPolygonArea/minFacetArea, minFacetArea/maxPolygonArea);
	//		Transformation_2 shrt(CGAL::SCALING, shrinkFactor);
	//		ifstream mcfile(mc_pts_filename.c_str(), std::ios_base::in);
	//		mcfile>>nbInitFront;
	//		mcfile>>nbInitBack;
	//		unsigned int idx = 0;
	//		pgns_3.resize(nbInitFront+nbInitBack);
	//		polygonNormals.resize(nbInitFront+nbInitBack);
	//		tangentPlanes.resize(nbInitFront+nbInitBack);
	//		projFacets.resize(nbInitFront+nbInitBack);
	//		pgnTextCoord.resize(nbInitFront+nbInitBack);
	//		pgnLibIdx.resize(nbInitFront+nbInitBack);
	//		fbmasks.resize(nbInitFront+nbInitBack, true);
	//		factors.resize(nbInitFront+nbInitBack, shrinkFactor);
	//		for (unsigned int i = 0; i < nbInitFront; i++, idx++)//distribute front polygons
	//		{
	//			double x, y, z;
	//			char dummyc;
	//			mcfile>>dummyc;
	//			mcfile>>x; mcfile>>y; mcfile>>z;
	//			vec3 location(x,y,z);
	//			meshDom.project_to_mesh(location, polygonNormals[idx], projFacets[idx]);
	//			Polygon_2 pgn = transform(shrt, pgnLib[i%pgnLib.size()]);
	//			movePolygon2Mesh(pgn, location, polygonNormals[idx], &pgns_3[idx]);
	//			//fbmasks[idx] = true;
	//			pgnTextCoord[idx] = pgnTextCoordLib[i%pgnLib.size()];
	//			pgnLibIdx[idx] = i%pgnLib.size();
	//			tangentPlanes[idx] = Plane_3(pgns_3[idx][0].source(), Vector_3(polygonNormals[idx].x,polygonNormals[idx].y,polygonNormals[idx].z));
	//		}
	//		for (unsigned int i = 0; i < nbInitBack; i++, idx++)//distribute back polygons
	//		{
	//			double x, y, z;
	//			char dummyc;
	//			mcfile>>dummyc;
	//			mcfile>>x; mcfile>>y; mcfile>>z;
	//			vec3 location(x, y, z);
	//			meshDom.project_to_mesh(location, polygonNormals[idx], projFacets[idx]);
	//			Polygon_2 pgn = transform(shrt, bgPgnLib[i%bgPgnLib.size()]);
	//			movePolygon2Mesh(pgn, location, polygonNormals[idx], &pgns_3[idx]);
	//			fbmasks[idx] = false;
	//			pgnTextCoord[idx] = bgPgnTextCoordLib[i%bgPgnLib.size()];
	//			pgnLibIdx[idx] = i%bgPgnLib.size();
	//			tangentPlanes[idx] = Plane_3(pgns_3[idx][0].source(), Vector_3(polygonNormals[idx].x,polygonNormals[idx].y,polygonNormals[idx].z));
	//		}
	//		mcfile.close();
	//	}
	//	else
	//		restore(restore_file);
	//	genSphericalBounder();
	//}
	//void Doc::loadPolygonsWithBackgroud(const std::vector<string>& filenames_, const string& restore_file, 
	//								const std::vector<string>& bg_pgn_filenames_, const string& mc_pts_filename)
	//{
	//	pgnLib.clear();
	//	texture = true;
	//	double maxPgnArea = DBL_MIN;
	//	for (unsigned int i = 0; i < 10 /*filenames_.size()*/; i++)
	//	{
	//		ifstream infile(filenames_[i].c_str(), std::ios_base::in);
	//		int nbVert;
	//		infile>>nbVert;
	//		pgnLib.push_back(Polygon_2());
	//		pgnTextCoordLib.push_back(Polygon_2());
	//		for (int j = 0; j < nbVert; j++)
	//		{
	//			double x, y, texx, texy;
	//			infile >> x;
	//			infile >> y;
	//			infile >> texx;
	//			infile >> texy;
	//			pgnLib.back().push_back(Point_2(x, y));
	//			pgnTextCoordLib.back().push_back(Point_2(texx, texy));
	//		}
	//		maxPolygonArea = std::max<double>(maxPolygonArea, std::fabs(pgnLib.back().area()));
	//		infile.close();
	//	}
	//	bgPgnLib.clear();
	//	bgPgnTextCoordLib.clear();
	//	for (unsigned int i = 0; i < 10/*bg_pgn_filenames_.size()*/; i++)
	//	{
	//		ifstream infile(bg_pgn_filenames_[i].c_str(), std::ios_base::in);
	//		int nbVert;
	//		infile>>nbVert;
	//		bgPgnLib.push_back(Polygon_2());
	//		bgPgnTextCoordLib.push_back(Polygon_2());
	//		for (int j = 0; j < nbVert; j++)
	//		{
	//			double x, y, texx, texy;
	//			infile>>x;	infile>>y;
	//			infile>>texx; infile>>texy;
	//			bgPgnLib.back().push_back(Point_2(x, y));	
	//			bgPgnTextCoordLib.back().push_back(Point_2(texx, texy));
	//		}
	//		maxPolygonArea = std::max<double>(maxPolygonArea, std::fabs(bgPgnLib.back().area()));
	//		infile.close();
	//	}
	//	std::cout<<"Maximum area of polygons "<<maxPolygonArea<<std::endl;
	//	genBoundaryPoints();
	//	if (restore_file.empty())
	//		// distribute polygons at multi-class points, generate texture coordinates for background polygons
	//	{
	//		//read initial points
	//		double shrinkFactor = std::min<double>(maxPolygonArea/minFacetArea, minFacetArea/maxPolygonArea);
	//		Transformation_2 shrt(CGAL::SCALING, shrinkFactor);
	//		ifstream mcfile(mc_pts_filename.c_str(), std::ios_base::in);
	//		mcfile>>nbInitFront;
	//		mcfile>>nbInitBack;
	//		unsigned int idx = 0;
	//		pgns_3.resize(nbInitFront+nbInitBack);
	//		polygonNormals.resize(nbInitFront+nbInitBack);
	//		tangentPlanes.resize(nbInitFront+nbInitBack);
	//		projFacets.resize(nbInitFront+nbInitBack);
	//		pgnTextCoord.resize(nbInitFront+nbInitBack);
	//		pgnLibIdx.resize(nbInitFront+nbInitBack);
	//		fbmasks.resize(nbInitFront+nbInitBack, true);
	//		factors.resize(nbInitFront+nbInitBack, shrinkFactor);
	//		for (unsigned int i = 0; i < nbInitFront; i++, idx++)//distribute front polygons
	//		{
	//			double x, y, z;
	//			char dummyc;
	//			mcfile>>dummyc;
	//			mcfile>>x; mcfile>>y; mcfile>>z;
	//			vec3 location(x,y,z);
	//			meshDom.project_to_mesh(location, polygonNormals[idx], projFacets[idx]);
	//			Polygon_2 pgn = transform(shrt, pgnLib[i%pgnLib.size()]);
	//			movePolygon2Mesh(pgn, location, polygonNormals[idx], &pgns_3[idx]);
	//			//fbmasks[idx] = true;
	//			pgnTextCoord[idx] = pgnTextCoordLib[i%pgnLib.size()];
	//			pgnLibIdx[idx] = i%pgnLib.size();
	//			tangentPlanes[idx] = Plane_3(pgns_3[idx][0].source(), Vector_3(polygonNormals[idx].x,polygonNormals[idx].y,polygonNormals[idx].z));
	//		}
	//		for (unsigned int i = 0; i < nbInitBack; i++, idx++)//distribute back polygons
	//		{
	//			double x, y, z;
	//			char dummyc;
	//			mcfile>>dummyc;
	//			mcfile>>x; mcfile>>y; mcfile>>z;
	//			vec3 location(x, y, z);
	//			meshDom.project_to_mesh(location, polygonNormals[idx], projFacets[idx]);
	//			Polygon_2 pgn = transform(shrt, bgPgnLib[i%bgPgnLib.size()]);
	//			movePolygon2Mesh(pgn, location, polygonNormals[idx], &pgns_3[idx]);
	//			fbmasks[idx] = false;
	//			pgnTextCoord[idx] = bgPgnTextCoordLib[i%bgPgnLib.size()];
	//			pgnLibIdx[idx] = i%bgPgnLib.size();
	//			tangentPlanes[idx] = Plane_3(pgns_3[idx][0].source(), Vector_3(polygonNormals[idx].x,polygonNormals[idx].y,polygonNormals[idx].z));
	//		}
	//		mcfile.close();
	//	}
	//	else
	//		restore(restore_file);
	//	genSphericalBounder();
	//}
	void Doc::loadMeshDomain(const string& bdFileName)
	{
		int len = bdFileName.length() ;
		bool is_obj = ( bdFileName.substr(len-3, 3)==string("obj"));
		if(is_obj) 
			meshDom.load(bdFileName);
		else 
			meshDom.load_tri(bdFileName) ;
		powerdiagram_ = new PowerDiagram_3(&meshDom);
		meshDom.build_kdtree();
		meshArea = 0.0;
		minFacetArea = DBL_MAX;
		for (unsigned int i = 0; i < meshDom.size(); i++)
		{
			double area = std::fabs(meshDom[i].area());
			if (area < minFacetArea)
				minFacetArea = area;
			meshArea += area;
		}
		std::cout<<"Area of underlying mesh: "<<meshArea<<std::endl;
		//get_bbox(xmin, ymin, zmin, xmax, ymax, zmax);
	}
	void Doc::genBoundaryPoints()//sampling points on boundary or feature edges
	{
		// construct boundary points
		const TriMesh::BoundaryEdgeSet& boundaryEdges = meshDom.getBoundary();
		const std::vector<pair<int, int>>& featureEdges = meshDom.getFeatureEdges();
		std::set<int> bdfeaVertIdx;//vertex indices considered as boundary or features
		if (boundaryEdges.size() > 0 || featureEdges.size() > 0)
		{
			double w = 0.0;
			int sn = 50;
			for (TriMesh::BoundaryEdgeSet::const_iterator it = boundaryEdges.begin(); it != boundaryEdges.end(); it++)
			{
				int vi0 = it->first, vi1 = it->second;
				const vec3& v0 = meshDom.vertex(vi0).point();
				const vec3& v1 = meshDom.vertex(vi1).point();
				if (bdfeaVertIdx.find(vi0) == bdfeaVertIdx.end())
				{
					bdfeaVertIdx.insert(vi0);
					boundaryPoints.push_back(Weighted_point(to_cgal(v0), w));
				}
				if (bdfeaVertIdx.find(vi1) == bdfeaVertIdx.end())
				{
					bdfeaVertIdx.insert(vi1);
					boundaryPoints.push_back(Weighted_point(to_cgal(v1), w));
				}
				for ( int i = 1; i < sn; i++ )
				{
					vec3 sv = ((sn-i)*v0 + i*v1)/sn;
					boundaryPoints.push_back(Weighted_point(to_cgal(sv), w));
				}
			}
			for (unsigned int i = 0; i < featureEdges.size(); i++)
			{
				int vi0 = featureEdges[i].first, vi1 = featureEdges[i].second;
				const vec3& v0 = meshDom.vertex(vi0).point();
				const vec3& v1 = meshDom.vertex(vi1).point();
				if (bdfeaVertIdx.find(vi0) == bdfeaVertIdx.end())
				{
					bdfeaVertIdx.insert(vi0);
					boundaryPoints.push_back(Weighted_point(to_cgal(v0), w));
				}
				if (bdfeaVertIdx.find(vi1) == bdfeaVertIdx.end())
				{
					bdfeaVertIdx.insert(vi1);
					boundaryPoints.push_back(Weighted_point(to_cgal(v1), w));
				}
				for ( int i = 1; i < sn; i++ )
				{
					vec3 sv = ((sn-i)*v0 + i*v1)/sn;
					boundaryPoints.push_back(Weighted_point(to_cgal(sv), w));
				}
			}
		}
	}
	void Doc::placePolygons(std::vector<Polygon_2> pgns)// value copy is used here on purpose
	{
		pgns_3.resize(pgns.size());
		polygonNormals.resize(pgns.size());
		tangentPlanes.resize(pgns.size());
		projFacets.resize(pgns.size());
		std::set<unsigned int> facetTranversed;
		// find the polygon with the maximum area
		const TriMesh::BoundaryEdgeSet& boundaryEdges = meshDom.getBoundary();
		const std::vector<pair<int, int>>& featureEdges = meshDom.getFeatureEdges();
		const std::vector<pair<int, int>>& featureFacets = meshDom.getFeatureFacets();
		for (unsigned int i = 0; i < meshDom.size(); i++)
		{
			const Facet& f = meshDom[i];
			int v0 = f.vertex_index[0], v1 = f.vertex_index[1], v2 = f.vertex_index[2];
			if ( boundaryEdges.find(make_pair(v0, v1)) != boundaryEdges.end() || boundaryEdges.find(make_pair(v1, v0)) != boundaryEdges.end() )
				facetTranversed.insert(i);
			else if ( boundaryEdges.find(make_pair(v0, v2)) != boundaryEdges.end() || boundaryEdges.find(make_pair(v2, v0)) != boundaryEdges.end() )
				facetTranversed.insert(i);
			else if ( boundaryEdges.find(make_pair(v1, v2)) != boundaryEdges.end() || boundaryEdges.find(make_pair(v2, v1)) != boundaryEdges.end() )
				facetTranversed.insert(i);
		} 
		for (unsigned int i = 0; i < featureFacets.size(); i++)
		{
			facetTranversed.insert(featureFacets[i].first);
			facetTranversed.insert(featureFacets[i].second);
		}
		double shrinkFactor = std::min<double>(maxPolygonArea/minFacetArea, minFacetArea/maxPolygonArea);
		shrinkFactor = std::sqrt(shrinkFactor)*0.5;
		for (unsigned int i = 0; i < pgns.size(); i++)
			shrinkPolygon(&pgns[i], shrinkFactor);
		// feature edges are given more polygons
		if ( ! curvature)
		{
			srand(time(NULL));
			for ( unsigned int i = 0; i < pgns.size(); i++ )
			{
				unsigned int fidx = ::rand() % meshDom.size();
				while ( facetTranversed.find(fidx) != facetTranversed.end() )
					fidx = ::rand() % meshDom.size();
				facetTranversed.insert(fidx);
				const Facet& f = meshDom[fidx];
				vec3 randCenter = (f.vertex[0] + f.vertex[1] + f.vertex[2])/3.0;
				vec3 n = meshDom[fidx].normal();
				//compute the shrink factor to guarantee that the polygon is totally inside a single mesh triangle
				movePolygon2Mesh(pgns[i], randCenter, n, &pgns_3[i]);
				polygonNormals[i] = n;
				// compute tangent planes at each polygons
				tangentPlanes[i] = Plane_3(pgns_3[i][0].source(), Vector_3(polygonNormals[i].x, polygonNormals[i].y, polygonNormals[i].z));
				projFacets[i] = fidx;
			}
		}
		else
		{
			double totalWeight = 0.0;
			double res = 0.0;
			int pgnIdx = 0;
			for (unsigned int i = 0; i < meshDom.nb_vertices(); i++)
				totalWeight += meshDom.vertex(i).weight();
			//////////////////////////////////////////////////////////////////////////
			const set<int>& highCurFacets = meshDom.getHighCurvatureFacets();
			set<int> vertHighCurv;
			for (set<int>::const_iterator it = highCurFacets.begin(); it != highCurFacets.end(); it++)
			{
				const Facet& f = meshDom[*it];
				for (int i = 0; i < 3; i++)
				{
					int vi = f.vertex_index[i];
					if (vertHighCurv.find(vi) != vertHighCurv.end())	
						continue;
					vertHighCurv.insert(vi);
					double nf = pgns.size()*meshDom.vertex(vi).weight()/totalWeight;
					int n(nf+res);
					if (n>=1)
					{
						res = nf + res - n;
						int nbput = std::min<int>(n, meshDom.vertex(vi).faces_.size());
						for ( int j = 0; j < nbput; j++)
						{
							int idx = meshDom.vertex(i).faces_[j];
							if (facetTranversed.find(idx) != facetTranversed.end())
								continue;
							facetTranversed.insert(idx);
							const Facet& pf = meshDom[idx];
							vec3 c = (pf.vertex[0] + pf.vertex[1] + pf.vertex[2])/3.0;
							vec3 nm = meshDom[idx].normal();
							movePolygon2Mesh(pgns[pgnIdx], c, nm, &pgns_3[pgnIdx]);
							polygonNormals[pgnIdx] = nm;
							tangentPlanes[pgnIdx] = Plane_3(pgns_3[pgnIdx][0].source(), Vector_3(nm.x, nm.y, nm.z));
							projFacets[pgnIdx] = idx;
							pgnIdx++;
						}
					}
					else
						res += nf;
				}

			}
			//////////////////////////////////////////////////////////////////////////
			//for (unsigned int i = 0; i < meshDom.nb_vertices(); i++)
			//{
			//	double nf = pgns.size()*meshDom.vertex(i).weight()/totalWeight;
			//	int n(nf + res);
			//	if (n>=1)
			//	{
			//		res = nf + res - n;
			//		int nbput = std::min<int>(n, meshDom.vertex(i).faces_.size());
			//		for ( int j = 0; j < nbput; j++)
			//		{
			//			int idx = meshDom.vertex(i).faces_[j];
			//			if (facetTranversed.find(idx) != facetTranversed.end())
			//				continue;
			//			facetTranversed.insert(idx);
			//			const Facet& f = meshDom[idx];
			//			vec3 c = (f.vertex[0] + f.vertex[1] + f.vertex[2])/3.0;
			//			vec3 nm = meshDom[idx].normal();
			//			movePolygon2Mesh(pgns[pgnIdx], c, nm, &pgns_3[pgnIdx]);
			//			polygonNormals[pgnIdx] = nm;
			//			tangentPlanes[pgnIdx] = Plane_3(pgns_3[pgnIdx][0].source(), Vector_3(nm.x, nm.y, nm.z));
			//			projFacets[pgnIdx] = idx;
			//			pgnIdx++;
			//		}
			//	}
			//	else
			//		res += nf;
			//}// pick up the lost ones
			while (pgnIdx < pgns.size())
			{
				unsigned int idx = ::rand()%meshDom.size();
				while ( facetTranversed.find(idx) != facetTranversed.end() )
					idx = ::rand()%meshDom.size();
				facetTranversed.insert(idx);
				const Facet& f = meshDom[idx];
				vec3 c = (f.vertex[0] + f.vertex[1] + f.vertex[2])/3.0;
				vec3 nm = meshDom[idx].normal();
				movePolygon2Mesh(pgns[pgnIdx], c, nm, &pgns_3[pgnIdx]);
				polygonNormals[pgnIdx] = nm;
				tangentPlanes[pgnIdx] = Plane_3(pgns_3[pgnIdx][0].source(), Vector_3(nm.x, nm.y, nm.z));
				projFacets[pgnIdx] = idx;
				pgnIdx++;
			}
		}
		factors.resize(pgns_3.size(), shrinkFactor);
		genSphericalBounder();
	}
	void Doc::setFeatureAngle(double ang)
	{
		meshDom.setFeatureAngle(ang);
		// recompute the 3d polygons and put them again
		boundaryPoints.clear();
		placePolygons(pgns2d);
	}
	void Doc::setCurvaturePercent(double percent)
	{
		if (curvature)
			meshDom.setHighCurvaturePercent(percent);
	}
	void Doc::clear()
	{
		meshDom.clear_all();
		powerdiagram_->clear();
		delete powerdiagram_;
		for (std::vector<std::vector<Segment_3>>::iterator ap = pgns_3.begin(); ap!=pgns_3.end(); ap++)
			ap->clear();
		pgns_3.clear();
	}
	void Doc::get_bbox(real& x_min, real& y_min, real& z_min, real& x_max, real& y_max, real& z_max)
	{
		x_min = y_min = z_min =  1e30;
		x_max = y_max = z_max = -1e30;
		for(unsigned int i=0; i<meshDom.size(); i++) 
		{
			for(unsigned int j=0; j<3; j++) 
			{
				const vec3& p = meshDom[i].vertex[j] ;
				x_min = gx_min(x_min, p.x) ;
				y_min = gx_min(y_min, p.y) ;
				z_min = gx_min(z_min, p.z) ;
				x_max = gx_max(x_max, p.x) ;
				y_max = gx_max(y_max, p.y) ;
				z_max = gx_max(z_max, p.z) ;
			}
		}
	}
	vec3 Doc::polygonCenter(const Polygon_2& pgn)
	{
		vec3 midp(0.0, 0.0, 0.0);
		for (Polygon_2::Vertex_const_iterator vi = pgn.vertices_begin(); vi!=pgn.vertices_end(); vi++)
		{
			vec3 v(vi->x(), vi->y(), 0.0);
			midp += v;
		}
		return midp/pgn.size();
	}
	vec3 Doc::polygonCenter(const std::vector<Segment_3>& pgn)
	{
		vec3 midp(0.0, 0.0, 0.0);
		for (std::vector<Segment_3>::const_iterator ei = pgn.begin(); ei!=pgn.end(); ei++)
			midp += to_geex(ei->source());
		return midp/pgn.size();
	}
	Point_2 Doc::polygonCentroid(const Polygon_2& p)
	{
		double polygonArea = 0.0;
		Vector_2 centroid(0.0, 0.0);
		vec3 tmp = polygonCenter(p);
		Point_2 center(tmp.x, tmp.y);
		for (unsigned int i = 0; i < p.size(); i++)
		{
			Segment_2 e = p.edge(i);
			Point_2 s = e.source(), t = e.target();
			double triArea = (s.x()-center.x())*(t.y()-center.y()) - (t.x()-center.x())*(s.y()-center.y());
			polygonArea += triArea;
			Vector_2 triCent = (center - CGAL::ORIGIN) + (s - CGAL::ORIGIN) + (t - CGAL::ORIGIN);
			centroid = centroid + triArea*triCent;
		}
		centroid = centroid / (3.0*polygonArea);
		return (CGAL::ORIGIN + centroid);
	}
	Point_3 Doc::polygonCentroid(const std::vector<Segment_3>& p)
	{
		//Vector_3 polygonArea(0.0, 0.0, 0.0);
		//Vector_3 polyCenter(0.0, 0.0, 0.0);
		//for (unsigned int i = 1; i < p.size()-1; i++)
		//{
		//	Vector_3 triArea = cross_product(Vector_3(p[0].source(), p[i].source()), Vector_3(p[0].source(),p[i+1].source()));
		//	polygonArea = polygonArea + triArea;
		//	Vector_3 triCent = (p[0].source()-CGAL::ORIGIN) + (p[i].source()-CGAL::ORIGIN) + (p[i+1].source()-CGAL::ORIGIN);
		//	polyCenter = polyCenter + CGAL::sqrt(triArea.squared_length()) * triCent;
		//}
		//polyCenter = polyCenter / (3*CGAL::sqrt(polygonArea.squared_length()));
		//return (CGAL::ORIGIN + polyCenter);
		/////////////////////////////////////////////////////////////////////////////////
		//Vector_3 polygonArea(0.0, 0.0, 0.0);
		//Vector_3 centroid(0.0, 0.0, 0.0);
		//Point_3 center = to_cgal(polygonCenter(p));
		//for (unsigned int i = 0; i < p.size(); i++)
		//{
		//	Vector_3 triArea = cross_product(Vector_3(center, p[i].source()), Vector_3(center, p[i].target()));
		//	polygonArea = polygonArea + triArea;
		//	Vector_3 triCent = (center-CGAL::ORIGIN) + (p[i].source()-CGAL::ORIGIN) + (p[i].target()-CGAL::ORIGIN);
		//	centroid = centroid + CGAL::sqrt(triArea.squared_length())*triCent;
		//}
		//centroid = centroid / (3*CGAL::sqrt(polygonArea.squared_length()));
		//return (CGAL::ORIGIN + centroid);
		////////////////////////////////////////////////////////////////////////////////
		std::vector<Point_3> verts;
		verts.reserve(p.size());
		for (unsigned int i = 0; i < p.size(); i++)
			verts.push_back(p[i].source());
		return CGAL::centroid(verts.begin(), verts.end(), CGAL::Dimension_tag<0>());
	}
	void Doc::shrinkPolygon(Polygon_2 *pgn, double f)
	{
		if (f <=0.0)
			return;
		Transformation_2 shrink(CGAL::SCALING, f);
		*pgn = CGAL::transform(shrink, *pgn);
	}
	// assume the polygon pgn is on the xy-plane, n is unit normal
	void Doc::movePolygon2Mesh(const Polygon_2& pgn, const vec3& p, const vec3& n, std::vector<Segment_3> *mvPgn)
	{
		mvPgn->clear();
		vec3 center = polygonCenter(pgn);
		Transformation_3 move_to_origin(CGAL::TRANSLATION, Vector_3(-center.x, -center.y, -center.z));
		//Point_2 cent2d = polygonCentroid(pgn);
		//Point_3 cent(cent2d.x(), cent2d.y(), 0.0);
		//Transformation_3 move_to_origin(CGAL::TRANSLATION, Vector_3(cent, CGAL::ORIGIN));
		// assume the normal vector n form angle theta with z axis, its projection on xOy form angle phi with x axis
		//sin and cos for transformation
		double cos_theta = n.z;
		double sin_theta = std::sqrt(n.x*n.x + n.y*n.y);
		// construct the rotation matrix
		Transformation_3 align_normal(CGAL::IDENTITY);
		if (n.x != 0.0 || n.y != 0.0) // if the facet has already had the same normal as z axis
		{
			double cos_phi = n.x / std::sqrt(n.x*n.x + n.y*n.y);
			double sin_phi = n.y / std::sqrt(n.x*n.x + n.y*n.y);
			align_normal = Transformation_3(cos_phi*cos_theta, -sin_phi, cos_phi*sin_theta, sin_phi*cos_theta, 
				cos_phi, sin_phi*sin_theta, -sin_theta, 0.0, cos_theta);
		}

		//Vector_3 translation(rotated_center_3d, to_cgal(p));
		Transformation_3 move_to_p(CGAL::TRANSLATION, Vector_3(p.x, p.y, p.z));
		Transformation_3 t = move_to_p*(align_normal*move_to_origin);
		//for (Polygon_2::Edge_const_iterator ei = pgn.edges_begin(); ei != pgn.edges_end(); ei++)
		for (unsigned int i = 0; i < pgn.size(); i++)
		{
			Point_2 src = pgn.vertex(i), tgt = pgn.vertex((i+1)%pgn.size());
			Segment_3 edge_3d(Point_3(src.x(), src.y(), 0.0), Point_3(tgt.x(), tgt.y(), 0.0));
			mvPgn->push_back(edge_3d.transform(t));
		}
	}
	/*
	** n: normal of the point p on mesh
	** m: normal of polygon plane
	*/
	void Doc::movePolygon2Mesh(const std::vector<Segment_3>& pgn, const vec3& p, vec3& n, vec3& m, std::vector<Segment_3> *mvPng)
	{
		mvPng->clear();
		// normalize
		n = n/n.length();
		m = m/m.length();
		// first, compute the rotation matrix so the two normals are aligned
		if (m.x*n.x+m.y*n.y+m.z*n.z < 0.0)
			m = -m;
		double cos_theta_m = m.z, sin_theta_m = sqrt(m.x*m.x+m.y*m.y);
		double cos_theta_n = n.z, sin_theta_n = std::sqrt(n.x*n.x + n.y*n.y);
		Transformation_3 rotzm(CGAL::IDENTITY);//rotate m on to the xOz plane
		if (m.x!=0 || m.y!=0)
		{
			double cos_phi_m = m.x/sin_theta_m, sin_phi_m = m.y/sin_theta_m;
			rotzm = Transformation_3(cos_phi_m, sin_phi_m, 0.0, -sin_phi_m, cos_phi_m, 0.0, 0.0, 0.0, 1.0);
		}
		double m00 = cos_theta_n*cos_theta_m + sin_theta_n*sin_theta_m,
			   m02 = sin_theta_n*cos_theta_m - cos_theta_n*sin_theta_m;
		Transformation_3 rotymn(m00, 0.0, m02, 0.0, 1.0, 0.0, -m02, 0.0, m00);//rotate m's mirror on xOz to n's mirror on xOz
		Transformation_3 rotzn(CGAL::IDENTITY);//rotate n's mirror on xOz to n
		if (n.x!=0 || n.y!=0)
		{
			double cos_phi_n = n.x/sin_theta_n, sin_phi_n = n.y/sin_theta_n;
			rotzn = Transformation_3(cos_phi_n, -sin_phi_n, 0.0, sin_phi_n, cos_phi_n, 0.0, 0.0, 0.0, 1.0);
		}
		Transformation_3 align_mn = rotzn*(rotymn*rotzm);
		//vec3 center = polygonCenter(pgn);
		vec3 center = to_geex(polygonCentroid(pgn));
		Transformation_3 center2origin(CGAL::TRANSLATION, Vector_3(-center.x, -center.y, -center.z));
		Transformation_3 move2p(CGAL::TRANSLATION, Vector_3(p.x, p.y, p.z));
		Transformation_3 t = move2p*(align_mn*center2origin);
		for (unsigned int i = 0; i < pgn.size(); i++)
			mvPng->push_back(pgn[i].transform(t));
	}
	void Doc::adjustPolygons(std::vector<Segment_3>& pgn, double f)
	{
		Point_3 pc = polygonCentroid(pgn);
		Transformation_3 move_to_origin(CGAL::TRANSLATION, Vector_3(pc, CGAL::ORIGIN));
		Transformation_3 move_back(CGAL::TRANSLATION, Vector_3(CGAL::ORIGIN, pc));
		Transformation_3 scale(CGAL::SCALING, f);
		Transformation_3 t = move_back*(scale*move_to_origin);
		for (std::vector<Segment_3>::iterator it = pgn.begin(); it!=pgn.end(); it++)
			*it = it->transform(t);
	}
	inline void Doc::adjustPolygons(unsigned int nb)
	{
#ifdef CILK
		cilk_for( unsigned int i = 0; i < nb; i++)
#else
		for ( unsigned int i = 0; i < nb; i++ )
#endif
		{
			adjustPolygons(pgns_3[i], 0.9);
			factors[i] *= 0.9;
		}
	}
	void Doc::adjustPolygons( double f )
	{
		if (f == 1.0 || f <= 1.0e-3)
			return;
		for ( unsigned int i = 0; i < pgns_3.size(); i++)
			adjustPolygons(pgns_3[i], f);
		for (unsigned int i = 0; i < factors.size(); i++)
			factors[i] *= f;
		genSphericalBounder();
	}
	void Doc::genSphericalBounder(unsigned int sampleGrain)
	{
		if (pgns_3.size() == 0)
			return;
		if (stopUpdate)
			return;
		std::vector<SphericalBounder> bList(pgns_3.size());
		//std::cout<<"Start generating spherical bounders\n";
#ifdef CILK
		cilk_for( unsigned int i = 0; i < bList.size(); i++)
#else
//#pragma omp parallel for 
		for (unsigned int i = 0; i < bList.size(); i++)
#endif
		{
			bList[i].setPolygon(pgns_3[i]);
			//bList[i].genInteriorMarginBounder(24);
			bList[i].genInteriorMarginBounder(sampleGrain, polygonNormals[i]);
		}
		// form the power diagram of these sphere8
		powerdiagram_->clear();
		powerdiagram_->begin_insert();
		for (unsigned int i = 0; i < bList.size(); i++)
		{
			powerdiagram_->insert_grouped_points(bList[i].getMarginBounder(), tangentPlanes[i]);
		}
		if (boundaryPoints.size() > 0)
		{
			//std::cout<<"Inserting boundary points\n";
			for (unsigned int i = 0; i < boundaryPoints.size(); i++)
				powerdiagram_->insert(boundaryPoints[i]);
		}
		powerdiagram_->end_insert();
		//std::cerr<<"Finish inserting points.\n";
	}
	Doc::LloydStatus Doc::nobdLloyd()
	{
		std::vector<Parameter> optiSolutions(pgns_3.size());
		polygonsReserved.clear();
		std::vector<bool> preserveMask(pgns_3.size(), false), failMask(pgns_3.size(), false);
		double kmin = 1.0e30;
		// optimize every group	
		//std::cout<<"Knitro optimization starts..."<<std::endl;
		const set<int>& fhc = meshDom.getHighCurvatureFacets();
		std::vector<Point_3> rotCents(pgns_3.size());
		std::vector<Vector_3> uvx(pgns_3.size()), uvy(pgns_3.size()), uvz(pgns_3.size());	
		bool noProgressFlag = false;
#ifdef CILK
		containments.resize(pgns_3.size());
		cilk::reducer_opor<bool> failFlag(false);
		cilk::reducer_opand<bool> stopFlag(true);
		cilk_for (unsigned int i = 0; i < pgns_3.size(); i++)
#else
		bool failFlag = false, stopFlag = true;
		for (unsigned int i = 0; i < pgns_3.size(); i++)
#endif
		{
#ifdef CILK
			containments[i].clear();
#else
			contaiment_.clear();
#endif
			uvx[i] = pgns_3[i][0].to_vector();
			//uvx[i] = Vector_3(pgns_3[i][0].source(), pgns_3[i][pgns_3[i].size()/2].source());
			uvx[i] = uvx[i]/sqrt(uvx[i].squared_length());
			uvz[i] = Vector_3(polygonNormals[i].x, polygonNormals[i].y, polygonNormals[i].z);
			uvz[i] = uvz[i]/sqrt(uvz[i].squared_length());
			uvy[i] = CGAL::cross_product(uvz[i], uvx[i]);
			Parameter solution(1.0, 0.0, 0.0, 0.0);
			int tryTimes = 1, knitroRes;
#ifdef CILK
			if (!stopUpdate)
				powerdiagram_->setConstraints(i, &containments[i].const_list_, uvx[i], uvy[i], &rotCents[i]);			
			else
				powerdiagram_->setConstraints(i, &containments[i].const_list_, uvx[i], uvy[i], &rotCents[i], pgns_3[i]);
			while ( (knitroRes=optimize_by_Knitro(&solution.k, &solution.theta, &solution.t1, &solution.t2, i)) && tryTimes <= 10 )
#else
			if (!stopUpdate)
				powerdiagram_->setConstraints(i, &contaiment_.const_list_, uvx[i], uvy[i], &rotCents[i]);
			else
				powerdiagram_->setConstraints(i, &contaiment_.const_list_, uvx[i], uvy[i], &rotCents[i], pgns_3[i]);
			while ( ((knitroRes=optimize_by_Knitro(&solution.k, &solution.theta, &solution.t1, &solution.t2)) ) && tryTimes <= 10 )
#endif
			{
				solution = Parameter(1.0, (Numeric::float64() - 0.5)/3.0*pi, 0.0, 0.0);
				tryTimes++;
			}
			if (!knitroRes)
			{
				optiSolutions[i] = solution;
				if (solution.k < kmin)
					kmin = solution.k;
				preserveMask[i] = true;
			}
			else
			{
				std::cerr<< "Failed to optimize group "<<i<<". return value: "<<knitroRes<<", k = "<<solution.k<<std::endl;
				//double m = (factors[i] + minScalor)/2;
				if ( factors[i] >= (minScalor+maxScalor)/2.0 )
					//polygonsReserved.push_back(i);
					preserveMask[i] = true;
				if (curvature && fhc.find(projFacets[i]) != fhc.end())//polygons at high curvature will not removed
					preserveMask[i] = true;
				double m = std::min<double>((minScalor+factors[i])/2.0, factors[i]);
				adjustPolygons(pgns_3[i], m/factors[i]);
				factors[i] = m;
				failMask[i] = true;
				optiSolutions[i].k = 1.0;
				optiSolutions[i].theta = optiSolutions[i].t1 = optiSolutions[i].t2 = 0.0;
				failFlag |= true;
			}
		}
		for (unsigned int i = 0; i < preserveMask.size(); i++)
			if (preserveMask[i])
				polygonsReserved.push_back(i);

		double k = 1.0;//actual enlargement scale factor
		if (kmin >= 1.0005) // polygons can be enlarged
		{
			if ( kmin >= 1.05)
				//k = (kmin-1.0)*0.8+1;
				k = kmin;
			else
				k = kmin;
			if (noEnlargeTimes >= 2)
			{
				std::cout<<"Enlarging...\n";
				for (unsigned int j = 0; j < optiSolutions.size(); j++)
				{
					if (!failMask[j])
						optiSolutions[j].k = k;
					if ( curvature && fhc.find(projFacets[j]) != fhc.end())
					{
						double mappedScale = slope*meshDom[projFacets[j]].weight() + y_intercept;
						optiSolutions[j].k = std::min<double>(optiSolutions[j].k, mappedScale/factors[j]);
					}
				}
			}
			else
			{
				for (unsigned int j = 0; j < optiSolutions.size(); j++)
					optiSolutions[j].k = 1.0;
				noEnlargeTimes++;
			}
			stopUpdate = false;
		}
		else
		{
			for (unsigned int j = 0; j < optiSolutions.size(); j++)
			{
				if (optiSolutions[j].k <= 1.0)
				{
					optiSolutions[j].t1 = optiSolutions[j].t2 = optiSolutions[j].theta = 0.0;
					optiSolutions[j].k = 1.0;
				}
				else
				{				
					if (curvature && fhc.find(projFacets[j]) != fhc.end())
					{
						double mappedScale = slope*meshDom[projFacets[j]].weight() + y_intercept;
						// facets at high curvature will not be enlarged too much
						optiSolutions[j].k = std::min<double>(mappedScale/factors[j], optiSolutions[j].k);
					}
					else
						optiSolutions[j].k = std::min<double>(optiSolutions[j].k, maxScalor/factors[j]);
				}
			}
			noProgressFlag = true;
			stopUpdate = false;
		}
#ifdef CILK
		cilk_for(unsigned int i = 0; i < pgns_3.size(); i++)
#else
		for (unsigned int i = 0; i < pgns_3.size(); i++)
#endif
		{
			if (factors[i] >= maxScalor)	
			{
				stopFlag &= true;
				optiSolutions[i].k = 1.0;
			}
			else
				stopFlag &= false;
			std::vector<Segment_3>& pgn = pgns_3[i];
			double fcos = cos(optiSolutions[i].theta);
			double fsin = sin(optiSolutions[i].theta);
			Polygon_2 planarPolygon;
			for (unsigned int j = 0; j < pgn.size(); j++)
				planarPolygon.push_back(to_uv(uvx[i], uvy[i], rotCents[i], pgn[j].source()));
			// coordinates of rotation center in the local system
			double ik = optiSolutions[i].k;
			Transformation_2 enlarge(ik*fcos, -ik*fsin, optiSolutions[i].t1, ik*fsin, ik*fcos, optiSolutions[i].t2);
			planarPolygon = transform(enlarge, planarPolygon);
			std::vector<Segment_3> floatPgn;
			for (unsigned int j = 0; j < planarPolygon.size(); j++)
			{
				Point_2 s = planarPolygon.edge(j).source(), t = planarPolygon.edge(j).target();
				Point_3 src = rotCents[i] + s.x()*uvx[i] + s.y()*uvy[i];
				Point_3 tgt = rotCents[i] + t.x()*uvx[i] + t.y()*uvy[i];
				floatPgn.push_back(Segment_3(src, tgt));
			}
			if (!stopUpdate)
			{
				vec3 center = to_geex(polygonCentroid(floatPgn));
				vec3 prjPnt =  meshDom.project_to_mesh(center, polygonNormals[i], projFacets[i]);		
				movePolygon2Mesh(floatPgn, prjPnt, polygonNormals[i], vec3(uvz[i].x(), uvz[i].y(), uvz[i].z()), &pgns_3[i]);
				Vector_3 n(polygonNormals[i].x, polygonNormals[i].y, polygonNormals[i].z);
				tangentPlanes[i] = Plane_3(pgns_3[i][0].source(), n);
			}
			else
				pgns_3[i] = floatPgn;
			factors[i] *= ik;
		}

		double polygonAreaSum = 0.0;
		for (unsigned int i = 0; i < pgns_3.size(); i++)
		{
			polygonAreaSum += computePolygonArea(pgns_3[i]);
		}
		areaRatio = polygonAreaSum/meshArea;
		std::cout<<"Area Ratio: "<<areaRatio<<std::endl;

		if (!stopUpdate)
			genSphericalBounder();
		else
			std::cout<<"No update\n";
#ifdef CILK
		if (stopFlag.get_value())
#else
		if (stopFlag)
#endif
			return LLOYD_STOP;

#ifdef CILK
		else if (noProgressFlag || failFlag.get_value() )
#else
		else if (noProgressFlag || failFlag )
#endif	
		{

			if ( areaRatio > areaCoverage)
				return LLOYD_CONVERGE;
#ifdef CILK
			else if (failFlag.get_value())
#else
			else if (failFlag)
#endif
				return LLOYD_FAILURE;
			else
				return LLOYD_NO_PROGRESS;
		}
		else
			return LLOYD_SUCCESS;
	}
	void Doc::postLloyd(unsigned int nbBeforeFilling)
	{
		std::vector<Parameter> optiSolutions(pgns_3.size());
		polygonsReserved.clear();
		std::vector<bool> masks(pgns_3.size(), false);
		double kmin = DBL_MAX;
		std::vector<Point_3> rotCents(pgns_3.size());
		std::vector<Vector_3> uvx(pgns_3.size()), uvy(pgns_3.size()), uvz(pgns_3.size());
		const set<int>& hcf = meshDom.getHighCurvatureFacets();
#ifdef CILK
		containments.resize(pgns_3.size());
		cilk::reducer_opor<bool> failFlag(false);
		cilk_for (unsigned int i = 0; i < pgns_3.size(); i++)
#else
		bool failFlag = false;
		for (unsigned int i = 0; i < pgns_3.size(); i++)
#endif
		{
#ifdef CILK
			containments[i].clear();
#else
			contaiment_.clear();
#endif
			uvx[i] = pgns_3[i][0].to_vector();
			uvx[i] = uvx[i]/sqrt(uvx[i].squared_length());
			uvz[i] = Vector_3(polygonNormals[i].x, polygonNormals[i].y, polygonNormals[i].z);
			uvz[i] = uvz[i]/sqrt(uvz[i].squared_length());
			uvy[i] = CGAL::cross_product(uvz[i], uvx[i]);
			Parameter solution(1.0, 0.0, 0.0, 0.0);
			int tryTimes = 1, knitroRes;
#ifdef CILK
			powerdiagram_->setConstraints(i, &containments[i].const_list_, uvx[i], uvy[i], &rotCents[i]);
			//powerdiagram_->setConstraints(i, &containments[i].const_list_, uvx[i], uvy[i], &rotCents[i], pgns_3[i]);
			while ( ((knitroRes=optimize_by_Knitro(&solution.k, &solution.theta, &solution.t1, &solution.t2, i)) ) && tryTimes <= 10 )
#else
			powerdiagram_->setConstraints(i, &contaiment_.const_list_, uvx[i], uvy[i], &rotCents[i]);
			//powerdiagram_->setConstraints(i, &contaiment_.const_list_, uvx[i], uvy[i], &rotCents[i], pgns_3[i]);
			while ( ((knitroRes=optimize_by_Knitro(&solution.k, &solution.theta, &solution.t1, &solution.t2)) ) && tryTimes <= 10 )
#endif			
			{
				solution = Parameter(1.0, (Numeric::float64() - 0.5)/3.0*pi, 0.0, 0.0);
				tryTimes++;
			}
			if (tryTimes <= 10)
			{
				optiSolutions[i] = solution;
				if (solution.k < kmin)
					kmin = solution.k;
				masks[i] = true;
			}
			else 
			{
				std::cerr<< "Failed to optimize group "<<i<<". Now shrink it "<<std::endl;
				//if ( i < nbBeforeFilling )
				//	//polygonsReserved.push_back(i);
				//	masks[i] = true;
				//else
				//{
				//	double m = std::min<double>(minScalor, (minScalor+factors[i])/2.0);
				//	adjustPolygons(pgns_3[i], m/factors[i]);
				//	factors[i] = m;
				//}
				optiSolutions[i].k = 1.0;
				optiSolutions[i].theta = optiSolutions[i].t1 = optiSolutions[i].t2;
				failFlag |= true;
			}
		}
		for (unsigned int i = 0; i < masks.size(); i++)
			if (masks[i])
				polygonsReserved.push_back(i);
		for (unsigned int j = 0; j < optiSolutions.size(); j++)
		{
			if (optiSolutions[j].k < 1.0 || !masks[j])
			{
				optiSolutions[j].t1 = optiSolutions[j].t2 = optiSolutions[j].theta = 0.0;
				optiSolutions[j].k = 1.0;
			}
			//if (curvature && hcf.find(projFacets[j]) != hcf.end())
			//{
			//	double mappedScale = slope*meshDom[projFacets[j]].weight() + y_intercept;
			//	optiSolutions[j].k = std::min<double>(mappedScale/factors[j], optiSolutions[j].k);
			//}
		}
#ifdef CILK
		cilk_for (unsigned int i = 0; i < pgns_3.size(); i++)
#else
		for (unsigned int i = 0; i < pgns_3.size(); i++)
#endif
		{
			std::vector<Segment_3>& pgn = pgns_3[i];
			double fcos = cos(optiSolutions[i].theta);
			double fsin = sin(optiSolutions[i].theta);
			Polygon_2 planarPolygon;
			for (unsigned int j = 0; j < pgn.size(); j++)
				planarPolygon.push_back(to_uv(uvx[i], uvy[i], rotCents[i], pgn[j].source()));
			// coordinates of rotation center in the local system
			double ik = optiSolutions[i].k;
			Transformation_2 enlarge(ik*fcos, -ik*fsin, optiSolutions[i].t1, ik*fsin, ik*fcos, optiSolutions[i].t2);
			planarPolygon = transform(enlarge, planarPolygon);
			std::vector<Segment_3> floatPgn;
			for (unsigned int j = 0; j < planarPolygon.size(); j++)
			{
				Point_2 s = planarPolygon.edge(j).source(), t = planarPolygon.edge(j).target();
				Point_3 src = rotCents[i] + s.x()*uvx[i] + s.y()*uvy[i];
				Point_3 tgt = rotCents[i] + t.x()*uvx[i] + t.y()*uvy[i];
				floatPgn.push_back(Segment_3(src, tgt));
			}
			//vec3 center = polygonCenter(floatPgn);
			vec3 center = to_geex(polygonCentroid(floatPgn));
			vec3 prjPnt =  meshDom.project_to_mesh(center, polygonNormals[i], projFacets[i]);
			movePolygon2Mesh(floatPgn, prjPnt, polygonNormals[i], vec3(uvz[i].x(), uvz[i].y(), uvz[i].z()), &pgns_3[i]);
			//Vector_3 n = cross_product(pgns_3[i][0].to_vector(), pgns_3[i][1].opposite().to_vector());
			//polygonNormals[i] = to_geex(CGAL::ORIGIN + n);
			Vector_3 n(polygonNormals[i].x, polygonNormals[i].y, polygonNormals[i].z);
			tangentPlanes[i] = Plane_3(pgns_3[i][0].source(), n);
			factors[i] *= ik;
		}
	}
	void Doc::postProcess(unsigned int nbBeforeFilling)
	{
		std::cout<<"Post processing starts...\n";
		//find the overlap pieces
			for (int i = 0; i < 2; i++)
			{
				postLloyd(nbBeforeFilling);
				genSphericalBounder();
				glut_viewer_redraw();
			}
/*			if (polygonsReserved.size() != pgns_3.size())
			{
				std::cout<<"Removing the recently added holes that overlap with others\n";
				remove();
				genSphericalBounder();
				glut_viewer_redraw();
				std::cout<<"Removing ends.\n";
			}*/		
/*			{
				std::cout<<"Final filling\n";
				detectHoles();
				mergeHoles();
				nbBeforeFilling = pgns_3.size();
				int nbFilled = fill(0.8);
				std::cout<<"Filling hole ends. "<<nbFilled<<" holes are added\n";
				genSphericalBounder();
			}*/	
		std::cout<<"Post processing ends.\n";
	}
	int Doc::replace()
	{
#ifdef CILK
		cilk::reducer_opadd<int> nbReplaced;
#else
		int nbReplaced = 0;
#endif
		// search the polygon library and find the best matcher
		// compute the 2D mirror of bisectors on xOy plane
		const std::vector<std::vector<Bisector>>& clippedPolygons = powerdiagram_->getClippedPolygons();
		const set<int>& hcf = meshDom.getHighCurvatureFacets();
#ifdef CILK
		cilk_for(unsigned int j = 0; j < pgns_3.size(); j++)
#else
		for (unsigned int j = 0; j < pgns_3.size(); j++)
#endif
		{	
			if (curvature && hcf.find(projFacets[j]) != hcf.end())//polygons at high curvature remains
				continue;
			const std::vector<Bisector>& clippedPolygon = clippedPolygons[j];
			Vector_3 uvx = pgns_3[j][0].to_vector();
			uvx = uvx/sqrt(uvx.squared_length());
			Vector_3 uvz(polygonNormals[j].x, polygonNormals[j].y, polygonNormals[j].z);
			uvz = uvz/sqrt(uvz.squared_length());
			Vector_3 uvy = CGAL::cross_product(uvz, uvx);
			Point_3 o = to_cgal(polygonCenter(pgns_3[j]));
			Polygon_2 base;
			for (unsigned int k = 0; k < clippedPolygon.size(); k++)
				base.push_back(to_uv(uvx, uvy, o, clippedPolygon[k].source()));
			if (!base.is_simple())
			{
				std::cout<<"Non-simple! Ignore it and continue...\n";
				continue;
			}
			PolygonMatcher pm(base);
			std::vector<SimilarityInfo> simDist;
			int cnt = 0;
			for ( int k = 0; k < pgnLib.size(); k++)
			{
				double dist = pm.insert(pgnLib[k]);
				simDist.push_back(SimilarityInfo(dist, k, cnt));
				cnt++;
			}
			//search for the best matcher
			std::sort(simDist.begin(), simDist.end(), similarityCmp);
			int tryTimes = std::min<int>(cnt, 50);
			Point_2 cbase = polygonCentroid(base), ccandidate;
			Polygon_2 candidate;
			int m;
			double optiScalor = 1.0;	
			bool flipped = false;
#ifdef CILK
			bool found = false;//the best match is found
			m = 0;
			while ( m < tryTimes && !found)
#else
			for (m = 0; m < tryTimes; m++)
#endif
			{
				const Polygon_2& p = pgnLib[simDist[m].libIndex];
				double minAngle = -pm.getMatchAngle(simDist[m].insertionOrder);
				Transformation_2 rot(CGAL::ROTATION, std::sin(minAngle), std::cos(minAngle));
				if (pm.getyAxisReflection(simDist[m].insertionOrder))
				{
					candidate = transform(Transformation_2(-1.0, 0.0, 0.0, 1.0), p);
					candidate = transform(rot, candidate);
					candidate.reverse_orientation();
					flipped = true;
				}
				else
				{
					candidate = transform(rot, p);
					flipped = false;
				}
				ccandidate = polygonCentroid(candidate);
				Transformation_2 align(CGAL::TRANSLATION, Vector_2(ccandidate, cbase));
				candidate = transform(align, candidate);
				ccandidate = ccandidate.transform(align);
#ifdef CILK
				found = softFill(base, ccandidate, candidate, optiScalor);
				m++;
#else
				if ( softFill(base, ccandidate, candidate, optiScalor))
					break;
#endif
			}
			if ( m == tryTimes )//no good match
				continue;
			std::vector<Segment_3> tmppgn;
			for (unsigned int k = 0; k < candidate.size(); k++)
			{
				Point_2 s = candidate.edge(k).source(), t = candidate.edge(k).target();
				Point_3 src = o + s.x()*uvx + s.y()*uvy;
				Point_3 tgt = o + t.x()*uvx + t.y()*uvy;
				tmppgn.push_back(Segment_3(src, tgt));
			}
			factors[j] = optiScalor;
			vec3 prjp = meshDom.project_to_mesh(to_geex(to_xy(uvx, uvy, o, ccandidate)), polygonNormals[j], projFacets[j]);
			movePolygon2Mesh(tmppgn, prjp, polygonNormals[j], vec3(uvz.x(), uvz.y(), uvz.z()), &pgns_3[j]);
			tangentPlanes[j] = Plane_3(pgns_3[j][0].source(), Vector_3(polygonNormals[j].x, polygonNormals[j].y, polygonNormals[j].z));			
			if (texture)
			{
				pgnLibIdx[j] = simDist[m-1].libIndex;
				if (flipped)
				{
					pgnTextCoord[j].clear();
					const Polygon_2& tp = pgnTextCoordLib[simDist[m-1].libIndex];
					for (unsigned int k = 0; k < tp.size(); k++)
						pgnTextCoord[j].push_back(tp.vertex((tp.size()-k)%tp.size()));
				}
				else
					pgnTextCoord[j] = pgnTextCoordLib[simDist[m-1].libIndex];
			}		
			nbReplaced += 1;	
		}
		//genSphericalBounder();
#if CILK
		return nbReplaced.get_value();
#else
		return nbReplaced;
#endif
	}
	int Doc::noDupReplace()
	{
		std::cout<<"start replace with less duplication\n";
		std::vector<std::vector<int>> candidatesIndex(pgns_3.size());//candidates' indices to polygon library
		std::vector<std::vector<Polygon_2>> candidatePolygons(pgns_3.size());
		std::vector<std::vector<bool>> flipMasks(pgns_3.size());
		std::vector<set<int>> neighborIndex(pgns_3.size());
		std::vector<std::vector<double>> candidateFactors(pgns_3.size());
		// form the neighborhood relationship
#ifdef CILK
		cilk::reducer_opadd<int> nbReplaced(0);
		cilk_for (unsigned int i = 0; i < pgns_3.size(); i++)
#else
		int nbReplaced = 0;
		for (unsigned int i = 0; i < pgns_3.size(); i++)			
#endif
		{
			const PowerDiagram_3::VGroup& vg = powerdiagram_->getGroup(i);
			for (unsigned int j = 0; j < vg.size(); j++)
			{
				std::vector<int>& dualVertGroupId = vg[j]->dualVert;
				for (unsigned int k = 0; k < dualVertGroupId.size(); k++)
					if (dualVertGroupId[k] >= 0)
						neighborIndex[i].insert(dualVertGroupId[k]);
			}
		}
		const std::vector<std::vector<Bisector>>& clippedPolygons = powerdiagram_->getClippedPolygons();
		const set<int>& hcf = meshDom.getHighCurvatureFacets();
		// compute all candidates for each cell
#ifdef CILK
		cilk_for(unsigned int j = 0; j < pgns_3.size(); j++)
#else
		for (unsigned int j = 0; j < pgns_3.size(); j++)
#endif
		{
			if (curvature && hcf.find(projFacets[j]) != hcf.end())//polygons at high curvature remains
				continue;
			const std::vector<Bisector>& clippedPolygon = clippedPolygons[j];
			Vector_3 uvx = pgns_3[j][0].to_vector();
			uvx = uvx/sqrt(uvx.squared_length());
			Vector_3 uvz(polygonNormals[j].x, polygonNormals[j].y, polygonNormals[j].z);
			uvz = uvz/sqrt(uvz.squared_length());
			Vector_3 uvy = CGAL::cross_product(uvz, uvx);
			Point_3 o = to_cgal(polygonCenter(pgns_3[j]));
			Polygon_2 base;
			for (unsigned int k = 0; k < clippedPolygon.size(); k++)
				base.push_back(to_uv(uvx, uvy, o, clippedPolygon[k].source()));
			if (!base.is_simple())
			{
				std::cout<<"Non-simple! Ignore it and continue...\n";
				continue;
			}
			PolygonMatcher pm(base);
			std::vector<SimilarityInfo> simDist;
			std::vector<Polygon_2> *pl;
			if (fbmasks[j])		pl = &pgnLib;
			else				pl = &bgPgnLib;	
			for ( int k = 0; k < pl->size(); k++)
			{
				double dist = pm.insert(pl->at(k));
				simDist.push_back(SimilarityInfo(dist, k, k));
			}
			//search for the best matcher
			std::sort(simDist.begin(), simDist.end(), similarityCmp);
			int tryTimes = std::min<int>(pl->size(), 50);
			Point_2 cbase = polygonCentroid(base), ccandidate;
			Polygon_2 candidate;
			double optiScalor = 1.0;		
			for (int m = 0; m < tryTimes; m++)
			{
				bool flipped = false;
				const Polygon_2& p = pl->at(simDist[m].libIndex);
				double minAngle = -pm.getMatchAngle(simDist[m].insertionOrder);
				Transformation_2 rot(CGAL::ROTATION, std::sin(minAngle), std::cos(minAngle));
				if (pm.getyAxisReflection(simDist[m].insertionOrder))
				{
					candidate = transform(Transformation_2(-1.0, 0.0, 0.0, 1.0), p);
					candidate = transform(rot, candidate);
					candidate.reverse_orientation();
					flipped = true;
				}
				else
				{
					candidate = transform(rot, p);
					flipped = false;
				}
				ccandidate = polygonCentroid(candidate);
				Transformation_2 align(CGAL::TRANSLATION, Vector_2(ccandidate, cbase));
				candidate = transform(align, candidate);
				ccandidate = ccandidate.transform(align);
				if (softFill(base, ccandidate, candidate, optiScalor))
				{
					candidatePolygons[j].push_back(candidate);
					candidatesIndex[j].push_back(simDist[m].libIndex);
					candidateFactors[j].push_back(optiScalor);
					flipMasks[j].push_back(flipped);
				}
			}
		}
		// determine which to replace
		std::vector<int> rr(pgns_3.size(), -1);//resulting replace, candidates' indices to polygon library
		std::vector<int> rri(pgns_3.size(), -1);//candidates' indices to candidateIndex
		for (unsigned int i = 0; i < pgns_3.size(); i++)
		{
			if (candidatesIndex[i].size() == 0)
				continue;//no replace candidate
			const set<int>& neighbor = neighborIndex[i];
			unsigned int j;
			for ( j = 0; j < candidatesIndex[i].size(); j++)
			{
				set<int>::const_iterator it;
				for ( it = neighbor.begin(); it != neighbor.end(); it++)
					if (rr[*it] == candidatesIndex[i][j])//it already exists in its neighbors
						break;
				if (it == neighbor.end())
				{
					rr[i] = candidatesIndex[i][j];
					rri[i] = j;
					break;
				}
			}
			if (j == candidatesIndex[i].size())//if every candidates occurs in neighbor, choose the first candidate
			{
				rr[i] = candidatesIndex[i][0];
				rri[i] = 0;
			}
		}
		// do replacement
#ifdef CILK
		cilk_for (unsigned int i = 0; i < pgns_3.size(); i++)
#else
		for (unsigned int i = 0; i < pgns_3.size(); i++)
#endif
		{
			if (rri[i] < 0)
				continue;
			Vector_3 uvx = pgns_3[i][0].to_vector();
			uvx = uvx/sqrt(uvx.squared_length());
			Vector_3 uvz(polygonNormals[i].x, polygonNormals[i].y, polygonNormals[i].z);
			uvz = uvz/sqrt(uvz.squared_length());
			Vector_3 uvy = CGAL::cross_product(uvz, uvx);
			Point_3 o = to_cgal(polygonCenter(pgns_3[i]));
			std::vector<Segment_3> tmppgn;
			Polygon_2& candidate = candidatePolygons[i][rri[i]];
			for (unsigned int k = 0; k < candidate.size(); k++)
			{
				Point_2 s = candidate.edge(k).source(), t = candidate.edge(k).target();
				Point_3 src = o + s.x()*uvx + s.y()*uvy;
				Point_3 tgt = o + t.x()*uvx + t.y()*uvy;
				tmppgn.push_back(Segment_3(src, tgt));
			}
			factors[i] = candidateFactors[i][rri[i]];
			Point_2 ccandidate = polygonCentroid(candidate);
			vec3 prjp = meshDom.project_to_mesh(to_geex(to_xy(uvx, uvy, o, ccandidate)), polygonNormals[i], projFacets[i]);
			movePolygon2Mesh(tmppgn, prjp, polygonNormals[i], vec3(uvz.x(), uvz.y(), uvz.z()), &pgns_3[i]);
			tangentPlanes[i] = Plane_3(pgns_3[i][0].source(), Vector_3(polygonNormals[i].x, polygonNormals[i].y, polygonNormals[i].z));			
			if (texture)
			{
				if (rr[i] != pgnLibIdx[i])
					nbReplaced += 1;
				pgnLibIdx[i] = rr[i];
				Polygon_2 texCoord;
				if (fbmasks[i])
					texCoord = pgnTextCoordLib[rr[i]];
				else
					texCoord = bgPgnTextCoordLib[rr[i]];
				if (flipMasks[i][rri[i]])
				{
					pgnTextCoord[i].clear();
					for (unsigned int k = 0; k < texCoord.size(); k++)
						pgnTextCoord[i].push_back(texCoord.vertex((texCoord.size()-k)%texCoord.size()));
				}
				else 
					pgnTextCoord[i] = texCoord;
			}	
		}
#ifdef CILK
		return nbReplaced.get_value();
#else
		return nbReplaced;
#endif
	}
	int Doc::bgreplace()
	{
		std::cout<<"Start background replacing\n";
#ifdef CILK
		cilk::reducer_opadd<int> nbReplaced;
#else
		int nbReplaced = 0;
#endif
		const std::vector<std::vector<Bisector>>& clippedPolygons = powerdiagram_->getClippedPolygons();
		const set<int>& hcf = meshDom.getHighCurvatureFacets();
#ifdef CILK
		cilk_for(unsigned int j = 0; j < pgns_3.size(); j++)
#else
		for (unsigned int j = 0; j < pgns_3.size(); j++)
#endif
		{	
			if (curvature && hcf.find(projFacets[j]) != hcf.end())//polygons at high curvature remains
				continue;
			const std::vector<Bisector>& clippedPolygon = clippedPolygons[j];
			Vector_3 uvx = pgns_3[j][0].to_vector();
			uvx = uvx/sqrt(uvx.squared_length());
			Vector_3 uvz(polygonNormals[j].x, polygonNormals[j].y, polygonNormals[j].z);
			uvz = uvz/sqrt(uvz.squared_length());
			Vector_3 uvy = CGAL::cross_product(uvz, uvx);
			Point_3 o = to_cgal(polygonCenter(pgns_3[j]));
			Polygon_2 base;
			for (unsigned int k = 0; k < clippedPolygon.size(); k++)
				base.push_back(to_uv(uvx, uvy, o, clippedPolygon[k].source()));
			if (!base.is_simple())
			{
				std::cout<<"Non-simple! Ignore it and continue...\n";
				continue;
			}
			PolygonMatcher pm(base);
			std::vector<SimilarityInfo> simDist;
			unsigned int size = 0;
			if (fbmasks[j])
			{
				size = pgnLib.size();
				for (int k = 0; k < pgnLib.size(); k++)
				{
					double dist = pm.insert(pgnLib[k]);
					simDist.push_back(SimilarityInfo(dist, k, k));
				}
			}
			else
			{
				size = bgPgnLib.size();
				for (int k = 0; k < bgPgnLib.size(); k++)
				{
					double dist = pm.insert(bgPgnLib[k]);
					simDist.push_back(SimilarityInfo(dist, k, k));
				}
			}
			//search for the best matcher
			std::sort(simDist.begin(), simDist.end(), similarityCmp);
			int tryTimes = std::min<int>(size, 50);
			Point_2 cbase = polygonCentroid(base), ccandidate;
			Polygon_2 candidate;
			int m;
			double optiScalor = 1.0;	
			bool flipped = false;
#ifdef CILK
			bool found = false;//the best match is found
			m = 0;
			while ( m < tryTimes && !found)
#else
			for (m = 0; m < tryTimes; m++)
#endif
			{
				Polygon_2 p;
				if (fbmasks[j])
					p = pgnLib[simDist[m].libIndex];
				else
					p = bgPgnLib[simDist[m].libIndex];
				double minAngle = -pm.getMatchAngle(simDist[m].insertionOrder);
				Transformation_2 rot(CGAL::ROTATION, std::sin(minAngle), std::cos(minAngle));
				if (pm.getyAxisReflection(simDist[m].insertionOrder))
				{
					candidate = transform(Transformation_2(-1.0, 0.0, 0.0, 1.0), p);
					candidate = transform(rot, candidate);
					candidate.reverse_orientation();
					flipped = true;
				}
				else
				{
					candidate = transform(rot, p);
					flipped = false;
				}
				ccandidate = polygonCentroid(candidate);
				Transformation_2 align(CGAL::TRANSLATION, Vector_2(ccandidate, cbase));
				candidate = transform(align, candidate);
				ccandidate = ccandidate.transform(align);
#ifdef CILK
				found = softFill(base, ccandidate, candidate, optiScalor);
				m++;
#else
				if ( softFill(base, ccandidate, candidate, optiScalor))
					break;
#endif
			}
			if ( m == tryTimes )//no good match
				continue;
			std::vector<Segment_3> tmppgn;
			for (unsigned int k = 0; k < candidate.size(); k++)
			{
				Point_2 s = candidate.edge(k).source(), t = candidate.edge(k).target();
				Point_3 src = o + s.x()*uvx + s.y()*uvy;
				Point_3 tgt = o + t.x()*uvx + t.y()*uvy;
				tmppgn.push_back(Segment_3(src, tgt));
			}
			factors[j] = optiScalor;
			vec3 prjp = meshDom.project_to_mesh(to_geex(to_xy(uvx, uvy, o, ccandidate)), polygonNormals[j], projFacets[j]);
			movePolygon2Mesh(tmppgn, prjp, polygonNormals[j], vec3(uvz.x(), uvz.y(), uvz.z()), &pgns_3[j]);
			tangentPlanes[j] = Plane_3(pgns_3[j][0].source(), Vector_3(polygonNormals[j].x, polygonNormals[j].y, polygonNormals[j].z));			
			if (texture)
			{
				pgnLibIdx[j] = simDist[m-1].libIndex;
				Polygon_2 texCoord;
				if (fbmasks[j])
					texCoord = pgnTextCoordLib[simDist[m-1].libIndex];
				else
					texCoord = bgPgnTextCoordLib[simDist[m-1].libIndex];
				if (flipped)
				{
					pgnTextCoord[j].clear();
					for (unsigned int k = 0; k < texCoord.size(); k++)
						pgnTextCoord[j].push_back(texCoord.vertex((texCoord.size()-k)%texCoord.size()));
				}
				else 
					pgnTextCoord[j] = texCoord;
			}		
			nbReplaced += 1;	
		}
#if CILK
		return nbReplaced.get_value();
#else
		return nbReplaced;
#endif
	}
	int Doc::affineReplace()
	{
		
		const std::vector<std::vector<Bisector>>& clippedPolygons = powerdiagram_->getClippedPolygons();
#ifdef CILK
		cilk::reducer_opadd<int> nbReplaced(0);
		cilk_for ( unsigned int i = 0; i < clippedPolygons.size(); i++ )
#else
		int nbReplaced = 0;
		for (unsigned int i = 0; i < clippedPolygons.size(); i++)
#endif
		{
			const std::vector<Bisector>& clippedPolygon = clippedPolygons[i];
			Vector_3 uvx = pgns_3[i][0].to_vector();
			uvx = uvx / sqrt(uvx.squared_length());
			Vector_3 uvz(polygonNormals[i].x, polygonNormals[i].y, polygonNormals[i].z);
			uvz = uvz / sqrt(uvz.squared_length());
			Vector_3 uvy = CGAL::cross_product(uvz, uvx);
			Point_3 o = to_cgal(polygonCenter(pgns_3[i]));
			std::vector<Segment_2> base;
			Polygon_2 region;
			for (unsigned int j = 0; j < clippedPolygon.size(); j++)
			{
				region.push_back(to_uv(uvx, uvy, o, clippedPolygon[j].source()));
				base.push_back(Segment_2(to_uv(uvx, uvy, o, clippedPolygon[j].source()), to_uv(uvx, uvy, o, clippedPolygon[j].target())));
			}
			if (!region.is_simple())
			{
				std::cout<<"Non-simple! Ignore it and continue...\n";
				continue;
			}
			PolygonMatcher pm(base);
			std::vector<SimilarityInfo> simDist;
			std::vector<Parameter> paras(pgnLib.size());
			for (int k = 0; k < pgnLib.size(); k++)
			{
				double dist = pm.affineMatch(pgnLib[k], &paras[k].k, &paras[k].theta, &paras[k].t1, &paras[k].t2);
				simDist.push_back(SimilarityInfo(dist, k, k));
			}
			std::sort(simDist.begin(), simDist.end(), similarityCmp);
			Polygon_2 candidate;
			unsigned int m;
			double optiScalor = 1.0;
			bool found = false;
			for (m=0; m<pgnLib.size() && !found; m++)
			{
				int idx = simDist[m].insertionOrder;
				candidate.clear();
				candidate = pgnLib[idx];				
				found = softFill(region, candidate, paras[idx], optiScalor);
			}
			if (m == pgnLib.size())
				continue;
			pgns_3[i].clear();
			for (unsigned int j = 0; j < candidate.size(); j++)
			{
				Point_2 s = candidate.edge(j).source(), t = candidate.edge(j).target();
				pgns_3[i].push_back(Segment_3(to_xy(uvx, uvy, o, s), to_xy(uvx, uvy, o, t)));
			}
			factors[i] = optiScalor;
			nbReplaced += 1;
		}
#ifdef CILK
		return nbReplaced.get_value();
#else
		return nbReplaced;
#endif
		//genSphericalBounder();
	}
	void Doc::pack()
	{
		const std::vector<pair<int,int>>& featureFacets = meshDom.getFeatureFacets();
		if (featureFacets.size() == 0)
			Doc::lloydUsed = &Doc::nobdLloyd;
		else
			Doc::lloydUsed = &Doc::nobdLloyd;
		if (curvature)
			setCurvScaleMap();
		rpack();
	}
	void Doc::remove()
	{
		for (unsigned int i = 0; i < polygonsReserved.size(); i++)
		{
			if (i != polygonsReserved[i])
			{
				pgns_3[i] = pgns_3[polygonsReserved[i]];
				polygonNormals[i] = polygonNormals[polygonsReserved[i]];
				tangentPlanes[i] = tangentPlanes[polygonsReserved[i]];
				factors[i] = factors[polygonsReserved[i]];
				projFacets[i] = projFacets[polygonsReserved[i]];
				if (texture)
				{
					pgnLibIdx[i] = pgnLibIdx[polygonsReserved[i]];
					pgnTextCoord[i] = pgnTextCoord[polygonsReserved[i]];
				}
				if (backgroud)
					fbmasks[i] = fbmasks[polygonsReserved[i]];
			}
		}
		pgns_3.resize(polygonsReserved.size());
		polygonNormals.resize(polygonsReserved.size());
		tangentPlanes.resize(polygonsReserved.size());
		factors.resize(polygonsReserved.size());
		projFacets.resize(polygonsReserved.size());
		if (texture)
		{
			pgnLibIdx.resize(polygonsReserved.size());
			pgnTextCoord.resize(polygonsReserved.size());
		}
		if (backgroud)
			fbmasks.resize(polygonsReserved.size());

	}
	void Doc::rpack()
	{
		pgns2d.clear();
		static int lloydIterTimes = 0;
		static int ft = 0; //fail times
		std::cout<<"Packing starts..."<<std::endl;
		static unsigned int nbBeforeFilling;
		static unsigned int replaceTimes = 0;
		std::cout<<pgns_3.size()<<" polygons initially."<<std::endl;
		static int noProgressTimes = 0;
		static LloydStatus ls, prels;
		static int rmTimes = 0;
		int n;
		for (n = 0; n < packIterLim; n++)
		{
			std::cout<<":::::: Lloyd Iteration "<<lloydIterTimes++<<" ::::::"<<std::endl;
			//while ( (ls = lloyd()) == LLOYD_SUCCESS)//magnify to the original size, then stop
			//{
			ls = lloyd();
			if (ls == LLOYD_SUCCESS)
			{
				std::cout<<"...Lloyd succeeded...\n";
				glut_viewer_redraw();
				//std::cout<<":::::: Lloyd Iteration "<<lloydIterTimes++<<" ::::::"<<std::endl;
				ft = 0;
				prels = ls;
			}
			//}
			else if (ls == LLOYD_FAILURE)//overlap occurs
			{
				std::cout<<"...Lloyd failed...\n";
				ft++;
				if (ft == 2)//remove those polygons that cannot grow
				{
					std::cout<<"Failed twice. Remove polygons that cannot grow...\n";
					remove();
					genSphericalBounder();
					std::cout<<"Removing polygons ends.\n";
					rmTimes++;
					ft = 0;
				}
				if (rmTimes >= 3)
				{
					std::cout<<"Removal happens more than 3 times.\n";
					if (replaceTimes < 3)
					{
						std::cout<<"Start replacing...\n";
						replaceTimes++;
						int nbReplaced = 0;
						if (!texture)
							nbReplaced = replace();
						else
							//nbReplaced = bgreplace();
							nbReplaced = noDupReplace();
						std::cout<<"Replacing ends. "<<nbReplaced<<" polygons are replaced\n";
					}
					else
					{
						std::cout<<"Filling holes starts..."<<std::endl;
						detectHoles();
						mergeHoles();
						nbBeforeFilling = pgns_3.size();
						int nbFilled = 0;
						if (!backgroud)
							nbFilled = fill(0.3);
						else
							nbFilled = bgfill(0.3);
						std::cout<<"Filling holes ends. "<<nbFilled<<" holes are added."<<std::endl;
						adjustPolygons(nbBeforeFilling);
						noEnlargeTimes = 0;
						replaceTimes = 0;
					}
					rmTimes = 0;
					ft = 0;
					genSphericalBounder();
				}
				glut_viewer_redraw();
			}
			else if (ls == LLOYD_NO_PROGRESS)
			{
				std::cout<<"Lloyd made no progress. \n";
				noProgressTimes++;
				if (noProgressTimes <= 3 /*|| prels == LLOYD_FAILURE*/ /*&& ft <= 2 */)
				{
					std::cout<<"Shrink some polygons...\n";
					adjustPolygons(nbBeforeFilling);
					for (unsigned int i = 0; i < pgns_3.size(); i++)
						if (factors[i]>maxScalor)
						{
							adjustPolygons(pgns_3[i], maxScalor/factors[i]);
							factors[i] = maxScalor;
						}	
					noEnlargeTimes = 0;
				}
				else
				{
					if (replaceTimes < 3)
					{
						std::cout<<"Start replacing...\n";
						replaceTimes++;
						int nbReplaced = 0;
						if (!texture)
							nbReplaced = replace();
						else
							//nbReplaced = bgreplace();
							nbReplaced = noDupReplace();
						std::cout<<"Replacing ends. "<<nbReplaced<<" polygons are replaced\n";
					}
					else
					{
						std::cout<<"Filling holes starts..."<<std::endl;
						detectHoles();
						mergeHoles();
						nbBeforeFilling = pgns_3.size();
						int nbFilled = 0;
						if (!backgroud)
							nbFilled = fill(0.3);
						else
							nbFilled = bgfill(0.3);
						std::cout<<"Filling holes ends. "<<nbFilled<<" holes are added."<<std::endl;
						adjustPolygons(nbBeforeFilling);
						noEnlargeTimes = 0;
						replaceTimes = 0;
					}
					noProgressTimes = 0;
					ft = 0;
				}
				genSphericalBounder();
				glut_viewer_redraw();
			}
			else if ( ls == LLOYD_CONVERGE)
				break;
			else// (ls == LLOYD_STOP)
			{
				std::cout<<"...Lloyd stopped...\n";
				// shrink those polygons which are enlarged over maxScale
				for ( unsigned int i = 0; i < pgns_3.size(); i++ )
					if (factors[i] > maxScalor)
					{
						adjustPolygons(pgns_3[i], 1.0/factors[i]);
						factors[i] = maxScalor;
					}
					genSphericalBounder();
					if (replaceTimes < 3)
					{
						std::cout<<"Replacing the current polygons...\n";
						int nbReplaced = 0;
						if (!backgroud)
							nbReplaced = replace();
						else
							//nbReplaced = bgreplace();
							nbReplaced = noDupReplace();
						std::cout<<"Replacing finished. "<<nbReplaced<<" polygons are replaced\n";
						replaceTimes++;
					}
					else
					{
						detectHoles();
						mergeHoles();
						nbBeforeFilling = pgns_3.size();
						int nbFilled = 0;
						if (!backgroud)
							nbFilled = fill(0.3);
						else
							nbFilled = bgfill(0.3);
						std::cout<<"Filling holes ends. "<<nbFilled<<" holes are added."<<std::endl;
						adjustPolygons(nbBeforeFilling);
						noEnlargeTimes=0;
						replaceTimes = 0;
					}
					genSphericalBounder();
					glut_viewer_redraw();
			}
			stringstream ss;
			ss << n;
			dump(immResDir + ss.str() +".dat");
			prels = ls;
			if (areaRatio <= currentMaxAreaRatio)
				noProgressTimes++;
			else
			{
				currentMaxAreaRatio = areaRatio;
				noProgressTimes = 0;
			}
			if (noProgressTimes > 10)
				break;
		}
		if (n == packIterLim)
			return;
		postProcess(nbBeforeFilling);
		double minf, maxf, avgf;
		getFactorStatistic(&maxf, &minf, &avgf);
		double polygonAreaSum = 0.0;
		for (unsigned int i = 0; i < pgns_3.size(); i++)
			polygonAreaSum += computePolygonArea(pgns_3[i]);
		double areaRatio = polygonAreaSum/meshArea;
		std::cout<<"Packing ends."<<std::endl;
		std::cout<<"Report: \n";
		std::cout<<"\t-"<<pgns_3.size()<<" polygons are packed."<<std::endl;
		std::cout<<"\t-"<<"Max: "<<maxf*100.0<<'%'<<std::endl;
		std::cout<<"\t-"<<"Min: "<<minf*100.0<<'%'<<std::endl;
		std::cout<<"\t-"<<"Average: "<<avgf*100.0<<'%'<<std::endl;
		std::cout<<"\t-"<<"Area ratio: "<<areaRatio*100.0<<'%'<<std::endl;
		glut_viewer_redraw();
	}
	bool Doc::inside(const Polygon_2& outer, const Polygon_2& inner)
	{
		unsigned int i;
		for ( i = 0; i < inner.size(); i++)
			if ( outer.has_on_unbounded_side(inner.vertex(i)) )
				break;
		if ( i == inner.size() )
		{
			// further check
			for (i = 0; i < inner.size(); i++)
				for (unsigned int j = 0; j < outer.size(); j++)
					if (CGAL::do_intersect(inner.edge(i), outer.edge(j)))
						return false;
			return true;
		}
		else
			return false;
	}
	bool Doc::inside(const std::vector<Point_2>& outer, const Polygon_2& inner)
	{
		unsigned int i;
		for ( i = 0; i < outer.size(); i++)
			if (inner.has_on_bounded_side(outer[i]))
				break;
		if (i==outer.size())
			return false;
		else
			return true;
	}
	bool Doc::softFill(const Polygon_2& region, Polygon_2& filler, const Parameter& para, double &optiScalor)
	{
		double theta = para.theta;
		double c = std::cos(theta), s = std::sin(theta);
		Transformation_2 t(maxScalor*c, -maxScalor*s, para.t1, maxScalor*s, maxScalor*c, para.t2);
		Polygon_2 sfiller = transform(t, filler);
		if (inside(region, sfiller))
		{
			filler = transform(t, filler);
			optiScalor = maxScalor;
			return true;
		}
		t = Transformation_2(minScalor*c, -minScalor*s, para.t1, minScalor*s, minScalor*c, para.t2);
		sfiller.clear();
		sfiller = transform(t, filler);
		if (!inside(region, sfiller))
			return false;
		double u = maxScalor, l = minScalor, m;
		optiScalor = minScalor;
		int seartchTimes = 20, i = 0;
		while (i < seartchTimes)
		{
			m = (u+l)/2.0;
			sfiller.clear();
			t = Transformation_2(m*c, -m*s, para.t1, m*s, m*c, para.t2);
			sfiller = transform(t, filler);
			if (inside(region, sfiller))
			{
				optiScalor = m;
				l = m;
			}
			else
				u = m;
			i++;
		}
		t = Transformation_2(optiScalor*c, -optiScalor*s, para.t1, optiScalor*s, optiScalor*c, para.t2);
		filler = transform(t, filler);
		return true;
	}
	bool Doc::softFill(const Polygon_2 hole, const Point_2& c, Polygon_2& filler, double& optiScalor)
	{
		Transformation_2 move2origin(CGAL::TRANSLATION, Vector_2(c, CGAL::ORIGIN));
		Transformation_2 shrink(CGAL::SCALING, maxScalor);
		Transformation_2 moveBack(CGAL::TRANSLATION, c - CGAL::ORIGIN);
		Transformation_2 t = moveBack*(shrink*move2origin);
		Polygon_2 sfiller = transform(t, filler);
		if (inside(hole, sfiller))
		{
			filler = transform(t, filler);
			optiScalor = maxScalor;
			return true;
		}
		shrink = Transformation_2(CGAL::SCALING, minScalor);
		t = moveBack*(shrink*move2origin);
		sfiller = transform(t, filler);
		if (!inside(hole, sfiller))
			return false;
		double u = maxScalor, l = minScalor, m;
		optiScalor = minScalor;
		int searchTimes = 20, i = 0;
		while ( i < searchTimes)
		{
			m = (u+l)/2.0;
			sfiller.clear();
			shrink = Transformation_2(CGAL::SCALING, m);
			t = moveBack*(shrink*move2origin);
			sfiller = transform(t,filler);
			if (inside(hole, sfiller))
			{
				l = m;
				optiScalor = m;
			}
			else
				u = m;
			i++;
		}
		shrink = Transformation_2(CGAL::SCALING, optiScalor);
		t = moveBack*(shrink*move2origin);
		filler = transform(t, filler);
		return true;
	}
	int Doc::fill(double percent)
	{
		std::cerr<<"Start filling holes..."<<std::endl;
		const std::vector<std::vector<Point_3>>& holes = powerdiagram_->getHoles();
		int nbFilled = 0;
		unsigned int nbLargeHoles = unsigned int(percent*holes.size());
		for (unsigned int i = 0; i < nbLargeHoles; i++)
		{
			const std::vector<Point_3>& hole = holes[i];
			Point_3 o = hole[1];
			Vector_3 uvx(hole[1],hole[0]);
			uvx = uvx/sqrt(uvx.squared_length());
			Vector_3 uvz = CGAL::cross_product(Vector_3(hole[1], hole[0]), Vector_3(hole[1], hole[2]));
			uvz = uvz/sqrt(uvz.squared_length());
			Vector_3 uvy = CGAL::cross_product(uvz, uvx);
			Polygon_2 uvHole;
			for (unsigned int k = 0; k < hole.size(); k++)
				uvHole.push_back(to_uv(uvx, uvy, o, hole[k]));
			if (!uvHole.is_simple())
			{
				std::cout<<"uvHole not simple!!!\n";
				system("pause");
			}
			PolygonMatcher pm(uvHole);
			std::vector<SimilarityInfo> simDist;//similarity distance and index
			int cnt = 0;
			for ( int k = 0; k < pgnLib.size(); k++)
			{
				if (minScalor * minScalor * std::abs(pgnLib[k].area())>= std::abs(uvHole.area()))
					continue;
				double dist = pm.insert(pgnLib[k]);
				simDist.push_back(SimilarityInfo(dist, k, cnt));
				cnt++;
			}
			std::sort(simDist.begin(), simDist.end(), similarityCmp);
			int tryTimes = std::min<int>(cnt, 50);
			// find the best matcher polygon in the polygon library
			// in the case of shape similarity, then check containment
			// there is a maximum trial time, avoid severe time overhead
			Point_2 chole = polygonCentroid(uvHole), cfiller;
			int j;
			Polygon_2 filler;
			double optiScalor = 1.0;
			bool flipped = false;
			for (j = 0; j < tryTimes; j++)
			{
				const Polygon_2& p = pgnLib[simDist[j].libIndex];
				double minAngle = -pm.getMatchAngle(simDist[j].insertionOrder);
				Transformation_2 rot_2(CGAL::ROTATION, std::sin(minAngle), std::cos(minAngle));
				if (pm.getyAxisReflection(simDist[j].insertionOrder))
				{
					filler = transform(Transformation_2(-1.0, 0.0, 0.0, 1.0), p);
					filler = transform(rot_2, filler);
					filler.reverse_orientation();
					flipped = true;
				}
				else
				{
					filler = transform(rot_2, p);
					flipped = false;
				}
				// check if filler is contained inside the current hole
				// if contained, break
				cfiller = polygonCentroid(filler);
				Transformation_2 align(CGAL::TRANSLATION, Vector_2(cfiller, chole));
				filler = transform(align, filler);
				cfiller = cfiller.transform(align);
				if ( softFill(uvHole, cfiller, filler, optiScalor) )
					break;
			}
			if ( j == tryTimes )//no good match
				continue;
			//pgns.push_back(filler);
			//if (filler.is_clockwise_oriented())
			//	filler.reverse_orientation();
			std::vector<Segment_3> pgn_3, tmppgn;
			for ( unsigned int k = 0; k < filler.size(); k++ )
			{
				Point_2 s = filler.vertex(k), t = filler.vertex((k+1)%filler.size());
				Point_3 src = o + s.x()*uvx + s.y()*uvy;
				Point_3 tgt = o + t.x()*uvx + t.y()*uvy;
				tmppgn.push_back(Segment_3(src,tgt));
			}
			factors.push_back(optiScalor);
			vec3 n;
			int pf;
			vec3 prjp = meshDom.project_to_mesh(to_geex(to_xy(uvx, uvy, o, cfiller)), n, pf);
			movePolygon2Mesh(tmppgn, prjp, n, vec3(uvz.x(), uvz.y(), uvz.z()), &pgn_3);
			polygonNormals.push_back(n);
			tangentPlanes.push_back(Plane_3(pgn_3[0].source(), Vector_3(n.x, n.y, n.z)));
			pgns_3.push_back(pgn_3);
			projFacets.push_back(pf);	
			if (texture)
			{
				pgnLibIdx.push_back(simDist[j].libIndex);
				if (flipped)
				{
					const Polygon_2& tp = pgnTextCoordLib[simDist[j].libIndex];
					pgnTextCoord.push_back(Polygon_2());
					for (unsigned int k = 0; k < tp.size(); k++)
						pgnTextCoord.back().push_back(tp.vertex((tp.size()-k)%tp.size()));
				}
				else
					pgnTextCoord.push_back(pgnTextCoordLib[simDist[j].libIndex]);
			}

			nbFilled++;
		}
		return nbFilled;
	}
	int Doc::bgfill(double percent)
	{
		std::cerr<<"Start background filling holes..."<<std::endl;
		const std::vector<std::vector<Point_3>>& holes = powerdiagram_->getHoles();
		int nbFilled = 0;
		unsigned int nbLargeHoles = unsigned int(percent*holes.size());
		for (unsigned int i = 0; i < nbLargeHoles; i++)
		{
			const std::vector<Point_3>& hole = holes[i];
			Point_3 o = hole[1];
			Vector_3 uvx(hole[1],hole[0]);
			uvx = uvx/sqrt(uvx.squared_length());
			Vector_3 uvz = CGAL::cross_product(Vector_3(hole[1], hole[0]), Vector_3(hole[1], hole[2]));
			uvz = uvz/sqrt(uvz.squared_length());
			Vector_3 uvy = CGAL::cross_product(uvz, uvx);
			Polygon_2 uvHole;
			for (unsigned int k = 0; k < hole.size(); k++)
				uvHole.push_back(to_uv(uvx, uvy, o, hole[k]));
			if (!uvHole.is_simple())
			{
				std::cout<<"uvHole not simple!!!\n";
				system("pause");
			}
			PolygonMatcher pm(uvHole);
			std::vector<SimilarityInfo> simDist;//similarity distance and index
			int cnt = 0;
			for ( int k = 0; k < bgPgnLib.size(); k++)
			{
				if (minScalor * minScalor * std::abs(bgPgnLib[k].area())>= std::abs(uvHole.area()))
					continue;
				double dist = pm.insert(bgPgnLib[k]);
				simDist.push_back(SimilarityInfo(dist, k, cnt));
				cnt++;
			}
			std::sort(simDist.begin(), simDist.end(), similarityCmp);
			int tryTimes = std::min<int>(cnt, 50);
			Point_2 chole = polygonCentroid(uvHole), cfiller;
			int j;
			Polygon_2 filler;
			double optiScalor = 1.0;
			bool flipped = false;
			for (j = 0; j < tryTimes; j++)
			{
				const Polygon_2& p = bgPgnLib[simDist[j].libIndex];
				double minAngle = -pm.getMatchAngle(simDist[j].insertionOrder);
				Transformation_2 rot_2(CGAL::ROTATION, std::sin(minAngle), std::cos(minAngle));
				if (pm.getyAxisReflection(simDist[j].insertionOrder))
				{
					filler = transform(Transformation_2(-1.0, 0.0, 0.0, 1.0), p);
					filler = transform(rot_2, filler);
					filler.reverse_orientation();
					flipped = true;
				}
				else
				{
					filler = transform(rot_2, p);
					flipped = false;
				}
				cfiller = polygonCentroid(filler);
				Transformation_2 align(CGAL::TRANSLATION, Vector_2(cfiller, chole));
				filler = transform(align, filler);
				cfiller = cfiller.transform(align);
				if ( softFill(uvHole, cfiller, filler, optiScalor) )
					break;
			}
			if ( j == tryTimes)//no good match
				continue;
			std::vector<Segment_3> pgn_3, tmppgn;
			for ( unsigned int k = 0; k < filler.size(); k++ )
			{
				Point_2 s = filler.vertex(k), t = filler.vertex((k+1)%filler.size());
				Point_3 src = o + s.x()*uvx + s.y()*uvy;
				Point_3 tgt = o + t.x()*uvx + t.y()*uvy;
				tmppgn.push_back(Segment_3(src,tgt));
			}
			factors.push_back(optiScalor);
			vec3 n;
			int pf;
			vec3 prjp = meshDom.project_to_mesh(to_geex(to_xy(uvx, uvy, o, cfiller)), n, pf);
			movePolygon2Mesh(tmppgn, prjp, n, vec3(uvz.x(), uvz.y(), uvz.z()), &pgn_3);
			polygonNormals.push_back(n);
			tangentPlanes.push_back(Plane_3(pgn_3[0].source(), Vector_3(n.x, n.y, n.z)));
			pgns_3.push_back(pgn_3);
			projFacets.push_back(pf);	
			fbmasks.push_back(false);
			if (texture)
			{
				pgnLibIdx.push_back(simDist[j].libIndex);
				if (flipped)
				{
					const Polygon_2& tp = bgPgnTextCoordLib[simDist[j].libIndex];
					pgnTextCoord.push_back(Polygon_2());
					for (unsigned int k = 0; k < tp.size(); k++)
						pgnTextCoord.back().push_back(tp.vertex((tp.size()-k)%tp.size()));
					//pgnTextCoord.push_back(bgPgnTextCoordLib[simDist[j].libIndex]);
				}
				else
					pgnTextCoord.push_back(bgPgnTextCoordLib[simDist[j].libIndex]);
			}

			nbFilled++;
		}
		return nbFilled;
	}
#ifdef CILK
	int Doc::optimize_by_Knitro(double* io_k, double* io_theta, double* io_t1, double* io_t2, unsigned int groupId)
#else
	int Doc::optimize_by_Knitro(double* io_k, double* io_theta, double* io_t1, double* io_t2)
#endif
	{
		int  nStatus;
		/* variables that are passed to KNITRO */
		KTR_context *kc;
		int n, m, nnzJ, nnzH, objGoal, objType;
		int *cType;
		int *jacIndexVars, *jacIndexCons;
		double obj, *x, *lambda;
		double *xLoBnds, *xUpBnds, *xInitial, *cLoBnds, *cUpBnds;
		int i, j, k; // convenience variables

		/*problem size and mem allocation */
		n = 4;
#ifdef CILK
		m = containments[groupId].get_constaint_list_size();
#else
		m = contaiment_.get_constaint_list_size();
#endif
		nnzJ = n*m;
		nnzH = 0;
		x = (double *)malloc(n * sizeof(double));
		lambda = (double *)malloc((m+n) * sizeof(double));

		xLoBnds      = (double *)malloc(n * sizeof(double));
		xUpBnds      = (double *)malloc(n * sizeof(double));
		xInitial     = (double *)malloc(n * sizeof(double));
		cType        = (int    *)malloc(m * sizeof(int));
		cLoBnds      = (double *)malloc(m * sizeof(double));
		cUpBnds      = (double *)malloc(m * sizeof(double));
		jacIndexVars = (int    *)malloc(nnzJ * sizeof(int));
		jacIndexCons = (int    *)malloc(nnzJ * sizeof(int));

		/* objective type */
		//objType = KTR_OBJTYPE_GENERAL;
		objType = KTR_OBJTYPE_LINEAR;
		objGoal = KTR_OBJGOAL_MAXIMIZE;

		/* bounds and constraints type */
		xLoBnds[0] = 0.0;
		xUpBnds[0] = KTR_INFBOUND;
		xLoBnds[1] = -Doc::pi/6.0;
		xUpBnds[1] = Doc::pi/6.0;
		for (i = 2; i < n; i++) {
			xLoBnds[i] = -KTR_INFBOUND;
			xUpBnds[i] = KTR_INFBOUND;
		}

		for (j = 0; j < m; j++) {
			cType[j] = KTR_CONTYPE_GENERAL;
			cLoBnds[j] = 0.0;
			cUpBnds[j] = KTR_INFBOUND;
		}

		/* initial point */
		xInitial[0] = *io_k;
		xInitial[1] = *io_theta;
		xInitial[2] = *io_t1;
		xInitial[3] = *io_t2;

		/* sparsity pattern (here, of a full matrix) */
		k = 0;
		for (j = 0; j < m; j++) 
			for (i = 0; i < n; i++){
				jacIndexCons[k] = j;
				jacIndexVars[k] = i;
				k++;
			}

		/* create a KNITRO instance */
		kc = KTR_new();
		if (kc == NULL)
			return -1 ; // probably a license issue

		/* set options: automatic gradient and hessian matrix */
		if (KTR_set_int_param_by_name(kc, "gradopt", KTR_GRADOPT_EXACT) != 0)
			return -1 ;
		if (KTR_set_int_param_by_name(kc, "hessopt", KTR_HESSOPT_LBFGS) != 0)
			return -1 ;
		if (KTR_set_int_param_by_name(kc, "outlev", 0) != 0)
			return -1 ;
		if (KTR_set_int_param_by_name(kc, "maxit", 1000) != 0)
			return -1 ;

		/* register the callback function */
		if (KTR_set_func_callback(kc, &callback) != 0)
			return -1 ;
		if (KTR_set_grad_callback(kc, &callback) != 0)
			return -1 ;

		/* pass the problem definition to KNITRO */
		nStatus = KTR_init_problem(kc, n, objGoal, objType, xLoBnds, xUpBnds, m, cType, cLoBnds, cUpBnds,
									nnzJ, jacIndexVars, jacIndexCons, nnzH, NULL, NULL, xInitial, NULL);

		//		nStatus = KTR_check_first_ders(kc, xInitial, 2, 1e-6, 1e-6, 0, 0.0, NULL, NULL, NULL, NULL);
		/* free memory (KNITRO maintains its own copy) */
		free(xLoBnds);
		free(xUpBnds);
		free(xInitial);
		free(cType);
		free(cLoBnds);
		free(cUpBnds);
		free(jacIndexVars);
		free(jacIndexCons);

		/* solver call */
#ifdef CILK
		nStatus = KTR_solve(kc, x, lambda, 0, &obj, NULL, NULL, NULL, NULL, NULL, &groupId);
#else
		nStatus = KTR_solve(kc, x, lambda, 0, &obj, NULL, NULL, NULL, NULL, NULL, NULL);
#endif

		if(nStatus ==0)
		{	
			*io_k =  x[0];
			*io_theta = x[1];
			*io_t1 = x[2];
			*io_t2 = x[3];
		}

		/* delete the KNITRO instance and primal/dual solution */
		KTR_free(&kc);
		free(x);
		free(lambda);
		return nStatus;
	}

	void Doc::saveCurrentPolygons(const string& fn) const
	{
		ofstream ofile(fn.c_str(), std::ios_base::out);
		ofile<< pgns_3.size()<<std::endl;
		for (int i = 0; i < pgns_3.size(); i++)
		{
			ofile<<pgns_3[i].size()<<std::endl;
			for (int j = 0; j < pgns_3[i].size(); j++)
			{
				Point_3 p = pgns_3[i][j].source();
				ofile<<p.x()<<" "<<p.y()<<" "<<p.z()<<std::endl;
			}
		}
		ofile.close();
	}

	double Doc::computePolygonArea( const std::vector<Segment_3>& p )
	{
		Vector_3 polygonArea(0.0, 0.0, 0.0);
		for (unsigned int i = 1; i < p.size()-1; i++)
		{
			Vector_3 triArea = cross_product(Vector_3(p[0].source(), p[i].source()), Vector_3(p[0].source(),p[i+1].source()));
			polygonArea = polygonArea + triArea;
		}
		return 0.5*sqrt(polygonArea.squared_length());
	}

	void Doc::getFactorStatistic(double *maxFactor, double *minFactor, double *avgFactor) const
	{
		*maxFactor = DBL_MIN;
		*minFactor = DBL_MAX;
		double sum = 0.0;
		for (unsigned int i = 0; i < factors.size(); i++)
		{
			if (*maxFactor < factors[i])
				*maxFactor = factors[i];
			if (*minFactor > factors[i])
				*minFactor = factors[i];
			sum += factors[i];
		}
		*avgFactor = sum/factors.size();
	}
	void Doc::setCurvScaleMap()
	{
		// linear mapping
		// highest curvature facet -> minScalor
		// lowest curvature facet -> (minScalor+maxScalor)/2
		double wmax = meshDom.getMaxFacetWeight();
		double wmin = meshDom.getMinFacetWeight();
		if (wmin == wmax)
		{
			slope = 0.0;
			y_intercept = (minScalor + maxScalor)/2.0;
		}
		else
		{
			slope = (minScalor - maxScalor) / (2.0*(wmax - wmin));
			y_intercept = minScalor - slope*wmax;
		}
	}
	void Doc::flipMeshNormals()
	{
		meshDom.flip_normals();
		for (unsigned int i = 0; i < polygonNormals.size(); i++)
			polygonNormals[i] = -polygonNormals[i];
	}

	void Doc::dump(const string& fn)
	{
		Archive arch(fn, Archive::STORE);
		//arch.append_var(areaCoverage);
		arch.append_vector(pgns_3);
		arch.append_vector(polygonNormals);
		arch.append_vector(projFacets);
		arch.append_vector(factors);
		if (texture)
		{
			arch.append_vector(pgnLibIdx);
			arch.append_vector(pgnTextCoord);
		}
		if (backgroud)
			arch.write(fbmasks);
	}
	void Doc::restore(const string& fn)
	{
		Archive arch(fn, Archive::LOAD);
		//arch.restore_var(areaCoverage);
		arch.restore_vector(pgns_3);
		arch.restore_vector(polygonNormals);
		arch.restore_vector(projFacets);
		arch.restore_vector(factors);
		if (texture)
		{
			arch.restore_vector(pgnLibIdx);
			arch.restore_vector(pgnTextCoord);
		}
		if (backgroud)
			arch.read(fbmasks);
		// restore tangent plane
		for (unsigned int i = 0; i < pgns_3.size(); i++)
			tangentPlanes.push_back(Plane_3(pgns_3[i][0].source(), Vector_3(polygonNormals[i].x, polygonNormals[i].y, polygonNormals[i].z)));
		genSphericalBounder();
	}
	void Doc::makeClipPlane(unsigned int pgnId)
	{
		// get the neighbor spheres in triangulation
		const PowerDiagram_3::VGroup& vg = powerdiagram_->getGroup(pgnId);
		std::vector<Point_3> fitPnts;
		for (unsigned int i = 0; i < vg.size(); i++)		
		{
			// get all adjacent vertices with vg[i]
			std::vector<PowerDiagram_3::Vertex_handle> adjv;
			powerdiagram_->adjacent_vertices(vg[i], back_inserter(adjv));
			// find the nearest one
			std::map<double, PowerDiagram_3::Vertex_handle> dists;
			for (unsigned int j = 0; j < adjv.size(); j++)
			{
				if (adjv[j]->group_id == vg[i]->group_id)
					continue;
				double d = CGAL::squared_distance(vg[i]->point(), adjv[j]->point());
				dists[d] = adjv[j];
			}
			fitPnts.push_back(CGAL::midpoint(dists.begin()->second->point(), vg[i]->point()));
		}
		std::vector<vec3> prjPnts(fitPnts.size());
		for (unsigned int i = 0; i < fitPnts.size(); i++)
		{
			vec3 n;
			prjPnts[i] = meshDom.project_to_mesh(to_geex(fitPnts[i]), n);
			fitPnts[i] = to_cgal(prjPnts[i]);
		}
		fpts.push_back(fitPnts);
		Point_3 centroid;
		double quality = CGAL::linear_least_squares_fitting_3(fitPnts.begin(), fitPnts.end(), tangentPlanes[pgnId], centroid, CGAL::Dimension_tag<0>());
		// move polygon
		// compute new normal vector
		vec3 newnormal;
		meshDom.project_to_mesh(to_geex(centroid), newnormal, projFacets[pgnId]);
		Vector_3 pnormal = tangentPlanes[pgnId].orthogonal_vector();			
		std::vector<Segment_3> tmp;
		movePolygon2Mesh(pgns_3[pgnId], to_geex(centroid), 
					vec3(pnormal.x(), pnormal.y(), pnormal.z()), polygonNormals[pgnId], &tmp);
		pgns_3[pgnId] = tmp;
		if (newnormal.x*pnormal.x() + newnormal.y*pnormal.y()+ newnormal.z*pnormal.z() <= 0.0)
			polygonNormals[pgnId] = vec3(-pnormal.x(), -pnormal.y(), -pnormal.z());
		else
			polygonNormals[pgnId] = vec3(pnormal.x(), pnormal.y(), pnormal.z());
		polygonNormals[pgnId] = polygonNormals[pgnId]/polygonNormals[pgnId].length();
	}
	void Doc::fitPlanes()
	{
//#ifdef CILK
//		cilk_for (unsigned int i = 0; i < pgns_3.size(); i++)
//#else
		for (unsigned int i = 0; i < pgns_3.size(); i++) 
//#endif
			makeClipPlane(i);
		genSphericalBounder();
		glut_viewer_redraw();
	}
	int Doc::affineFill()
	{
		int nbFilled = 0;
		const std::vector<std::vector<PowerDiagram_3::Vertex_handle>>& pholes = powerdiagram_->holeBorders;
		std::vector<Parameter> paras(pgnLib.size());
		std::vector<SimilarityInfo> simInfo(pgnLib.size());
		for (unsigned int i = 0; i < pholes.size(); i++)
		{
			const std::vector<PowerDiagram_3::Vertex_handle>& ahole = pholes[i];
			// transform to 2d
			std::vector<Point_3> pos;
			//std::vector<pair<Point_3, double>> pos;
			for (unsigned int j = 0; j < ahole.size(); j++)
				pos.push_back(ahole[j]->point());
				//pos.push_back(make_pair(ahole[j]->point(), 1.0/ahole.size()));
			Point_3 cent = CGAL::centroid(pos.begin(), pos.end(), CGAL::Dimension_tag<0>());
			//Point_3 cent = CGAL::barycenter(pos.begin(), pos.end());
			vec3 n; 
			int fid;
			vec3 prjcent = meshDom.project_to_mesh(to_geex(cent), n, fid);
			Vector_3 uvz(n.x, n.y, n.z);
			Point_3 o(to_cgal(prjcent));
			Plane_3 prjPlane(o, uvz);
			Vector_3 uvx = prjPlane.base1();
			uvx = uvx / sqrt(uvx.squared_length());
			Vector_3 uvy = cross_product(uvz, uvx);
			// project onto 2d plane
			std::vector<Point_2> hole2d;
			for (unsigned int j = 0; j < ahole.size(); j++)
				hole2d.push_back( to_uv(uvx, uvy, o, prjPlane.projection(ahole[j]->point())) );
			PolygonMatcher pm(hole2d);
			for (unsigned int j = 0; j < pgnLib.size(); j++)
			{
				simInfo[j].dist = pm.affineMatch(pgnLib[j], &paras[j].k, &paras[j].theta, &paras[j].t1, &paras[j].t2);
				simInfo[j].libIndex = simInfo[j].insertionOrder = j;
				
			}
			std::sort(simInfo.begin(), simInfo.end(), similarityCmp);
			int searchtimes = 10;
			double maxs = 1.0, mins = 0.6, optiScalor;
			unsigned int j;
			Polygon_2 testPgn;
			std::vector<Polygon_2> candidates;
			std::vector<unsigned int> candidatesIndex;
			for ( j = 0; j < simInfo.size(); j++)
			{
				//const Polygon_2& candidate = pgnLib[simInfo[j].libIndex];
				const Parameter& para = paras[simInfo[j].insertionOrder];
				double cs = para.k*std::cos(para.theta), ss = para.k*std::sin(para.theta);
				Transformation_2 t(cs, -ss, para.t1, ss, cs, para.t2);
				Polygon_2 candidate = CGAL::transform(t, pgnLib[simInfo[j].libIndex]);
				// binary search
				Point_2 pgnCentroid = CGAL::centroid(candidate.vertices_begin(), candidate.vertices_end(), CGAL::Dimension_tag<0>());
				Transformation_2 move2origin(CGAL::TRANSLATION, Vector_2(pgnCentroid, CGAL::ORIGIN));
				Transformation_2 scale(CGAL::SCALING, mins);
				testPgn = CGAL::transform(scale*move2origin, candidate);
				if (inside(hole2d, testPgn))
					continue;
				optiScalor = mins;
				int k = 0;
				double mid, u = maxs, l = mins;
				while (k < searchtimes)
				{
					mid = (u + l)/2;
					scale = Transformation_2(CGAL::SCALING, mid);
					testPgn = CGAL::transform(scale*move2origin, candidate);
					if (!inside(hole2d, testPgn))
					{
						l = mid;
						optiScalor = mid;
					}
					else
						u = mid;
					k++;
				}
				//break;
				candidates.push_back(testPgn);
				candidatesIndex.push_back(j);
			}
			if (candidates.size() == 0)
				continue;
			// containment check
			// check neighbor polygons
			set<int> neighborPolygonId;
			for (unsigned int k = 0; k < ahole.size(); k++)
				neighborPolygonId.insert(ahole[k]->group_id);
			map<int, Polygon_2> neighborPolygon2d;
			for (set<int>::iterator it = neighborPolygonId.begin(); it != neighborPolygonId.end(); it++)
			{
				const std::vector<Segment_3>& np = pgns_3[*it];
				Polygon_2 np2d;
				for (unsigned int k = 0; k < np.size(); k++)
				{
					Point_3 vert = prjPlane.projection(np[k].source());
					np2d.push_back(to_uv(uvx, uvy, o, vert));
				}
				neighborPolygon2d[*it] = np2d;
			}
			unsigned int l;
			for (l = 0; l < candidates.size(); l++)
			{
				set<int>::const_iterator it;
				const Polygon_2& cand = candidates[l];
				for ( it = neighborPolygonId.begin(); it != neighborPolygonId.end(); it++)
				{
					const Polygon_2& np2d = neighborPolygon2d[*it];
					int insideVertNb = 0;
					for (unsigned int k = 0; k < cand.size(); k++)
						if (np2d.has_on_bounded_side(cand.vertex(k)))
							insideVertNb++;
					if (insideVertNb > 1)//too much overlap(one vertex overlap is tolerated)
						break;
				}
				if (it == neighborPolygonId.end())
					break;
			}
			if ( l == candidates.size())//no good filler, even in many candidates
				continue;
			testPgn = candidates[l];
			j = candidatesIndex[l];
			std::vector<Segment_3> tmppgn;
			// transform to global system
			for (unsigned int m = 0; m < testPgn.size(); m++)
			{
				Point_2 src = testPgn.vertex(m), tgt = testPgn.vertex((m+1)%testPgn.size());
				tmppgn.push_back(Segment_3(to_xy(uvx, uvy, o, src), to_xy(uvx, uvy, o, tgt)));
			}
			cent = polygonCentroid(tmppgn);
			prjcent = meshDom.project_to_mesh(to_geex(cent), n, fid);
			pgns_3.push_back(std::vector<Segment_3>());
			movePolygon2Mesh(tmppgn, prjcent, n, vec3(uvz.x(), uvz.y(), uvz.z()), &pgns_3.back());
			factors.push_back(paras[simInfo[j].insertionOrder].k*optiScalor);
			polygonNormals.push_back(n);
			projFacets.push_back(fid);
			tangentPlanes.push_back(prjPlane);
			if (texture)
			{
				pgnLibIdx.push_back(simInfo[j].libIndex);
				pgnTextCoord.push_back(pgnTextCoordLib[simInfo[j].libIndex]);
			}	
			nbFilled++;
		}
		//std::cout<<"Start gen\n";
		genSphericalBounder();
		return nbFilled;
	}

	void Doc::postOptimization(unsigned int nbHoleFilled)
	{
		//shrink recently added polygons
		//std::cout<<"Number of holes: "<<nbHoleFilled<<std::endl;
		unsigned int base = pgns_3.size()-nbHoleFilled;
#ifdef CILK
		cilk_for (unsigned int i = 0; i < nbHoleFilled; i++)
#else
		for (unsigned int i = 0; i < nbHoleFilled; i++)
#endif
		{
			adjustPolygons(pgns_3[i+base], 0.2);
			factors[i+base] *= 0.2;
		}
		std::cout<<"End scaling\n";
		genSphericalBounder();
		std::cout<<"End sphere generating\n";
		glut_viewer_redraw();
		std::vector<Parameter> optiSolutions(nbHoleFilled);
		std::vector<Point_3> rotCents(nbHoleFilled);
		std::vector<Vector_3> uvx(nbHoleFilled), uvy(nbHoleFilled), uvz(nbHoleFilled);
		
		for (unsigned int iter = 0; iter < 6; iter++)
		{

		
#ifdef CILK
		containments.resize(nbHoleFilled);
		cilk_for (unsigned int i = 0; i < nbHoleFilled; i++)
#else
		for (unsigned int i = 0; i < nbHoleFilled; i++)
#endif
		{
#ifdef CILK
			containments[i].clear();
#else
			contaiment_.clear();
#endif
			uvx[i] = pgns_3[i+base][0].to_vector();
			uvx[i] = uvx[i] / sqrt(uvx[i].squared_length());
			uvz[i] = Vector_3(polygonNormals[i+base].x, polygonNormals[i+base].y, polygonNormals[i+base].z);
			uvy[i] = cross_product(uvz[i], uvx[i]);
			Parameter solution(1.0, 0.0, 0.0, 0.0);
			int tryTimes = 1, knitroRes;
#ifdef CILK
			powerdiagram_->setConstraints(i+base, &containments[i].const_list_, uvx[i], uvy[i], &rotCents[i]);
			while ( (knitroRes = optimize_by_Knitro(&solution.k, &solution.theta, &solution.t1, &solution.t2, i)) && tryTimes <= 10 )
#else
			powerdiagram_->setConstraints(i, &contaiment_.const_list_, uvx[i], uvy[i], &rotCents[i]);
			while ( (knitroRes=optimize_by_Knitro(&solution.k, &solution.theta, &solution.t1, &solution.t2)) && tryTimes <= 10 )
#endif
			{
				solution = Parameter(1.0, (Numeric::float64()-0.5)/3.0*pi, 0.0, 0.0);
				tryTimes++;
			}
			if (!knitroRes && solution.k > 1.0)
			{
				if (iter < 3)
					optiSolutions[i].k = 1.0;//only rotation and translation at the first half steps
				else
					optiSolutions[i] = solution;
			}
			else
			{
				optiSolutions[i].k = 1.0;
				optiSolutions[i].theta = optiSolutions[i].t1 = optiSolutions[i].t2 = 0.0;
			}
		}
#ifdef CILK
		cilk_for (unsigned int i = 0; i < nbHoleFilled; i++)
#else
		for (unsigned int i = 0; i < nbHoleFilled; i++)
#endif
		{
			if (optiSolutions[i].k <= 1.0)
				continue;
			std::vector<Segment_3>& pgn = pgns_3[i+base];
			double fcos = optiSolutions[i].k*cos(optiSolutions[i].theta), fsin = optiSolutions[i].k*sin(optiSolutions[i].theta);
			Polygon_2 planarPolygon;
			for (unsigned int j = 0; j < pgn.size(); j++)
				planarPolygon.push_back(to_uv(uvx[i], uvy[i], rotCents[i], pgn[j].source()));
			Transformation_2 t(fcos, -fsin, optiSolutions[i].t1, fsin, fcos, optiSolutions[i].t2);
			planarPolygon = transform(t, planarPolygon);
			std::vector<Segment_3> floatPgn;
			for (unsigned int j = 0; j < planarPolygon.size(); j++)
			{
				Point_2 src = planarPolygon.edge(j).source(), tgt = planarPolygon.edge(j).target();
				floatPgn.push_back(Segment_3(to_xy(uvx[i], uvy[i], rotCents[i], src), to_xy(uvx[i], uvy[i], rotCents[i], tgt)));
			}
			vec3 cent = to_geex(polygonCentroid(floatPgn));
			vec3 prjcent = meshDom.project_to_mesh(cent, polygonNormals[i+base], projFacets[i+base]);
			movePolygon2Mesh(floatPgn, prjcent, polygonNormals[i+base], vec3(uvz[i].x(), uvz[i].y(), uvz[i].z()), &pgns_3[i+base]);
			tangentPlanes[i+base] = Plane_3(pgns_3[i+base][0].source(), Vector_3(polygonNormals[i+base].x, polygonNormals[i+base].y, polygonNormals[i+base].z));
			factors[i+base] *= optiSolutions[i].k;
		}
		genSphericalBounder();
		glut_viewer_redraw();
		}
	}

}//namespace Geex

int  callback (const int evalRequestCode,
			   const int n,
			   const int m,
			   const int nnzJ,
			   const int nnzH,
			   const double * const x,
			   const double * const lambda,
			   double * const obj,
			   double * const c,
			   double * const objGrad,
			   double * const jac,
			   double * const hessian,
			   double * const hessVector,
			   void * userParams) 
{

   Geex::Doc* doc = Geex::Doc::instance() ;
   if (evalRequestCode == KTR_RC_EVALFC)
   {
	   /* objective function */
	   // x[0]- x[4]: k cos sin t1 t2
	   *obj    = x[0];
#ifdef CILK
		doc->containments[*(unsigned int*)userParams].compute_constraints(x,c,m);
#else
	   doc->contaiment_.compute_constraints(x,c,m);
#endif
	   return(0);
   }
   else if (evalRequestCode == KTR_RC_EVALGA)
   {
	   objGrad[0] = 1.0;
	   objGrad[1] = 0.0;
	   objGrad[2] = 0.0;
	   objGrad[3] = 0.0;
#ifdef CILK
	   doc->containments[*(unsigned int*)userParams].compute_constraint_grads(x, jac);
#else
	   doc->contaiment_.compute_constraint_grads(x, jac);
#endif
	   return(0);
   }
   else 
   {
	   printf ("Wrong evalRequestCode in callback function.\n");
	   return(-1);
   }
}

