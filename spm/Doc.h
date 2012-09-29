#pragma  once
#include <fstream>
#include <cmath>
#include "spm_cgal.h"
#include "SphericalBounder.h"
#include "power_diagram.h"
//#include "knitro.h"
#include "meshes.h"
using Geex::Numeric::is_nan;
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
		double k ;
		double theta;
		double t1;
		double t2 ;
	};
	// ===================================================================================
	class Doc {
	public:
		typedef enum {LLOYD_SUCCESS, LLOYD_STOP, LLOYD_FAILURE} LloydStatus;
		Doc() ;
		~Doc();
		/* set methods */
		void loadPolygons( const string& pgnFileName );
		void loadMeshDomain(const string& bdFileName);
		void setOptimal();
		bool& setRotate() {return rotateAllowed;}
		void clear();
		void adjustPolygons(double f);
		/* access methods */
		static Doc* instance() { return instance_ ; }
		const PowerDiagram_3& getPowerDiagram() { return *powerdiagram_; }
		const vector<vector<Segment_3>>& getPolygons() const {return pgns_3;}
		const vector<vec3>& getPolygonNormals() const {return polygonNormals;}
		const vector<SphericalBounder*>& getMA() const {return bList;}//get the skeleton of MA
		const TriMesh& getMeshDomain() {return meshDom;}
		const vector<SphericalBounder*>& getSphericalBounder() {return bList;}
		const vector<int>& getSelectedTriangles()	{return facetTranversed;}
		/* computation */
		int optimize_by_Knitro(double* io_k, double* io_theta, double* io_t1, double* io_t2);
		LloydStatus lloyd(int nb_it);
		void shrink( double shrinkFactor = 0.9);
		void transformPolygons(int group_id, const vec2& midp, const Transformation& trans, const double k);
		void permute();
		void get_bbox(real& x_min, real& y_min, real& z_min,
						real& x_max, real& y_max, real& z_max);
	private:
		void genBounder();  //generate the medial axis of the interior polygons
		vec3 polygonCenter(const Polygon_2& pgn);
		void shrinkPolygon(Polygon_2 *pgn, double f);
		//move a polygon to point p on mesh, with its normal aligned with n, which is p's normal on mesh
		//vector<Segment_3> movePolygon2Mesh(const Polygon_2& pgn, const vec3& p, const vec3& n);
		//vector<Segment_3> movePolygon2Mesh(const vector<Segment_3>& pgn, const vec3& p, const vec3& n);
		void movePolygon2Mesh(const Polygon_2& pgn, const vec3& p, const vec3& n, vector<Segment_3> *mvPgn);
		void movePolygon2Mesh(const Polygon_2& pgn, const Facet& f, vector<Segment_3> *mvPgn);
		void genSphericalBounder();
		inline void clearBounder();//clear all the bounding spheres for update use
	private:
		TriMesh meshDom;
		PowerDiagram_3 *powerdiagram_;
		Containment contaiment_;
		vector<Polygon_2> pgns;
		vector<unsigned int> multiplicity;
		vector<vector<Segment_3>> pgns_3;//3d polygons
		vector<vec3> polygonNormals;
		vector<Polygon_2> optimalPolygons;
		vector<SphericalBounder*> bList;
		vector<Parameter> optimalConfig;//the optimal configuration so far
		double optimalScaling; //the optimal scaling factor k so far
		double totalEnlargement;
		bool rotateAllowed;
		static Doc* instance_ ;

		//debug
		vector<int> facetTranversed;
	} ;


}