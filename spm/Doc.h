#ifndef _DOC_H_
#define  _DOC_H_

#include <fstream>
#include <cmath>
#include <algorithm>
//#include <vector>
#include <iterator>
//#include <string>
#include <cfloat>
//#include <boost/shared_ptr.hpp>
#include <CGAL/create_offset_polygons_2.h>
#include <glut_viewer/glut_viewer.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Dimension.h>
#include <CGAL/barycenter.h>
#include "knitro.h"
//#include "spm_cgal.h"

//#include "meshes.h"


#include "archive.h"
#include "power_diagram.h"
#include "PolygonMatcher.h"

#ifdef CILK
#include <cilk/cilk.h>
#include <cilk/reducer_opor.h>
#include <cilk/reducer_opand.h>
#include <cilk/reducer_opadd.h>
#endif

namespace Geex {

	class Parameter
	{
	public:
		Parameter(){
			k = 1.0;
			theta = 0.0;
			t1 = 0;
			t2 = 0;
		}
		Parameter(double k, double theta, double t1, double t2)
		{
			this->k = k;
			this->theta = theta;
			this->t1 = t1;
			this->t2 = t2;
		}
		double k ;
		double theta;
		double t1;
		double t2 ;
	};
	// ===================================================================================
	struct SimilarityInfo
	{
		double dist;
		int libIndex;
		int insertionOrder;
		SimilarityInfo(double d, int li, int insertOrder): dist(d), libIndex(li), insertionOrder(insertOrder){}
		SimilarityInfo() : dist(0.0), libIndex(-1), insertionOrder(-1) {}
	};
	inline bool similarityCmp(const SimilarityInfo& s1, const SimilarityInfo& s2)
	{
		return (s1.dist <= s2.dist);
	}	
	class Doc 
	{
		//typedef boost::shared_ptr<Geex::Polygon_2> PolygonPrt;
		//typedef std::vector<PolygonPrt> PolygonPtrVector;	
	public:
		typedef enum {LLOYD_SUCCESS, LLOYD_STOP, LLOYD_FAILURE, LLOYD_NO_PROGRESS, LLOYD_CONVERGE} LloydStatus;
		Doc() ;
		~Doc();
		/* set methods */
		void loadPolygons( const string& pgnFileName);
		//void loadPolygonsWithTexture(const vector<std::string>& filenames);
		void loadPolygons(const string& pgnFileName, const string& restore_file);
		//void loadPolygonsWithTexture(const vector<std::string>& filenames, const std::string& restore_file);
		//void loadPolygonsWithBackgroud(const vector<string>& filenames, const string& restore_file,
		//								const string& bg_pgn_filename, const string& mc_pts_filename);
		//void loadPolygonsWithBackgroud(const vector<string>& filenames, const string& restore_file,
		//								const vector<string>& bg_pgn_filenames, const string& mc_pts_filename);
		void loadMeshDomain(const string& bdFileName);
		void flipMeshNormals();
		void clear();
		void adjustPolygons(double f);
		double& setMinScalor() { return minScalor; }
		double& setMaxScalor()	{ return maxScalor; }
		double& setAreaCoverage() { return areaCoverage; }
		int& setPackIterLim() { return packIterLim; }
		void setFeatureAngle(double ang);// { meshDom.setFeatureAngle(ang); }
		void setCurvaturePercent(double percent);
		/* access methods */
		static Doc* instance() { return instance_ ; }
		const PowerDiagram_3& getPowerDiagram() const { return *powerdiagram_; }
		const vector<vector<Segment_3>>& getPolygons() const {return pgns_3;}
		const vector<int>& getPolygonLibIndex() const { return pgnLibIdx; }
		const vector<Polygon_2>& getTextCoords() const { return pgnTextCoord; }
		void saveCurrentPolygons(const string& fn) const;
		const vector<vec3>& getPolygonNormals() const {return polygonNormals;}
		const vector<Weighted_point>& getBoundaryPoints() const { return boundaryPoints; }
		inline double getFeatureAngle() const { return meshDom.getFeatureAngle(); }
		inline double getCurvaturePercet() const { return meshDom.getHighCurvaturePercent(); }
		inline bool hasCurvature() const { return curvature; }
		inline const vector<bool>& getFBMasks() const { return fbmasks; }
		inline const bool hasBackground() const { return backgroud; }
		//debug
		const vector<int>& getProjectFacets() const { return projFacets; }
		TriMesh& getMeshDomain() {return meshDom;}
		void fitPlanes();
		const vector<vector<Point_3>>& fitPoints() { return fpts; }
/*		double getShrinkFactor() const { return shrinkFactor; }*/
		void getFactorStatistic(double *maxFactor, double *minFactor, double *avgFactor) const;
		void setCurvature() { curvature = true; }
		void setBackground() { backgroud = true; }
		//void correctNormals(); // Laplacian-like smoothness of normals
		//void smoothNormals();
		void dump(const string& fn);
		void restore(const string& fn);
		/* computation */
#ifdef CILK
		int optimize_by_Knitro(double* io_k, double* io_theta, double* io_t1, double* io_t2, unsigned int groupId);
#else
		int optimize_by_Knitro(double* io_k, double* io_theta, double* io_t1, double* io_t2);
#endif
		void get_bbox(real& x_min, real& y_min, real& z_min, real& x_max, real& y_max, real& z_max);
		inline Doc::LloydStatus lloyd() { return (this->*lloydUsed)(); }
		//LloydStatus bdLloyd();//lloyd with bounding box
		LloydStatus nobdLloyd();
		void pack();	
		void rpack();
		int fill(double percent);
		int bgfill(double percent);
		int noDupReplace();//replace with less duplication with 1-ring neighbors
		int replace();
		int bgreplace();//replace with the same kind of polygons
		int affineReplace();
		int affineFill();
		void postOptimization(unsigned int nbHoleFilled);
		void detectHoles() { powerdiagram_->detectHoles();  }
		void clusterDetectHoles() { genSphericalBounder(28); powerdiagram_->clusterHoleDetect(); }
		void mergeHoles() { powerdiagram_->mergeHoles(); glut_viewer_redraw();}
	private:
		void postProcess(unsigned int nbBeforeFilling);
		void postLloyd(unsigned int nbBeforeFilling);
		double computePolygonArea(const vector<Segment_3>& pgn);
		void adjustPolygons(vector<Segment_3>& pgn, double f);//adjust a single polygon
		inline void adjustPolygons(unsigned int nb);//adjust the only front nb polygons
		inline bool inside(const Polygon_2& outer, const Polygon_2& inner);
		inline bool inside(const vector<Point_2>& outer, const Polygon_2& inner);
		bool softFill(const Polygon_2& region, Polygon_2& filler, const Parameter& para, double &optiScalor);
		bool softFill(const Polygon_2 hole, const Point_2& c, Polygon_2& filler, double& optiScalor);
		void placePolygons(vector<Polygon_2> pgns);
		void mcPlacePolygons(const vector<Polygon_2>& pgns);//place polygons at multi-class points
		vec3 polygonCenter(const Polygon_2& pgn);
		vec3 polygonCenter(const vector<Segment_3>&);
		void shrinkPolygon(Polygon_2 *pgn, double f);
		//move a polygon to point p on mesh, with its normal aligned with n, which is p's normal on mesh
		void movePolygon2Mesh(const Polygon_2& pgn, const vec3& p, const vec3& n, vector<Segment_3> *mvPgn);
		//void movePolygon2Mesh(const Polygon_2& pgn, const Facet& f, vector<Segment_3> *mvPgn);
		//move an arbitrary polygon in 3D to a point p on mesh with normal n
		void movePolygon2Mesh(const vector<Segment_3>& pgn, const vec3& p, vec3& n, vec3& m, vector<Segment_3> *mvPng);
		void genSphericalBounder(unsigned int sampleGrain = 24);
		void genBoundaryPoints();
		Point_3 polygonCentroid(const vector<Segment_3>& p);
		//Point_3 polygonCentroid(const vector<Point_3>& p);
		Point_2 polygonCentroid(const Polygon_2& p);
		void remove();//remove polygons
		void setCurvScaleMap();// compute scale factor based on curvature at high curvature
		void makeClipPlane(unsigned int pgnId);
	public:
#ifdef CILK
		vector<Containment> containments;
		string immResDir;//immediate results
#else
		Containment contaiment_;
		string immResDir;//immediate results
#endif
	private:
		TriMesh meshDom;
		double meshArea;
		double minFacetArea;
		double maxPolygonArea;
		PowerDiagram_3 *powerdiagram_;
		vector<Polygon_2> pgnLib;
		vector<Polygon_2> pgnTextCoordLib;
		vector<Polygon_2> bgPgnLib;
		vector<Polygon_2> bgPgnTextCoordLib;
		vector<bool> fbmasks;
		vector<Polygon_2> pgns2d;
		vector<vector<Segment_3>> pgns_3;//3d polygons
		vector<Polygon_2> pgnTextCoord;
		vector<int> pgnLibIdx;//indices to the polygon library
		vector<int> projFacets;//indices of facets onto which polygons are projected
		vector<Weighted_point> boundaryPoints;
		vector<vec3> polygonNormals;
		vector<Plane_3> tangentPlanes;
		vector<double> factors;
	
		vector<unsigned int> polygonsReserved;
		static Doc* instance_ ;
		static const double pi /*= 3.14159265358979323846*/;
		double minScalor;
		double maxScalor;
		bool curvature;//has curvature input
		bool texture;//has texture input
		bool backgroud;//has back and front ground differences
		unsigned int nbInitBack;
		unsigned int nbInitFront;
		bool stopUpdate;//stop updating the power diagram
		//debug	
		vector<Point_3> bdboxVerts;
		double xmin, xmax, ymin, ymax, zmin, zmax;
		double areaCoverage;
		double areaRatio;
		double currentMaxAreaRatio;
		int noProgressTimes;
		int noEnlargeTimes;
		int packIterLim;
		string archFile;//dump and restore file name
		vector<vector<Point_3>> fpts;
		

	private:
		typedef LloydStatus (Doc::*ptr_memfunction)();
		//LloydStatus (Doc::*lloydUsed)();// the lloyd used actually
		ptr_memfunction lloydUsed;
		//for curvature to scale mapping
		double slope;
		double y_intercept;
	} ;


}

#endif