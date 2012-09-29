#include "Doc.h"
#include <glut_viewer/glut_viewer.h>
#include <algorithm>
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

	Doc::Doc()
	{
		gx_assert(instance_ == nil) ;
		instance_ = this ;
		totalEnlargement = 1.0;
		optimalScaling = 1.0;
		rotateAllowed = true;
	}
	Doc::~Doc() 
	{
		clear();
		instance_ = nil ; 
	}
	void Doc::loadPolygons( const string& pgnFileName )
	{
		pgns.clear();
		for (vector<Polygon_2>::iterator it = pgns.begin(); it!= pgns.end(); it++)
			it->clear();
		ifstream fin(pgnFileName.c_str(), std::ios_base::in);
		int nbPolygons = 0;
		fin>>nbPolygons;
		for (int i = 0; i < nbPolygons; i++)
		{
			int nbEdges;
			double percent;
			fin>>percent;
			fin>>nbEdges;
			Polygon_2 pgn;
			for (int j = 0; j < nbEdges; j++)
			{
				Polygon_2::FT x, y;
				fin>>x;
				fin>>y;
				pgn.push_back(Point(x, y));	
			}
			if (pgn.orientation() == CGAL::CLOCKWISE)
				pgn.reverse_orientation();
			pgns.push_back(pgn);
			multiplicity.push_back(std::max<unsigned int>(1, static_cast<unsigned int>(percent*128)));
		}
		fin.close();
		//size_t size = pgns.size();
		//for (unsigned int i = 0; i < size; i++ )
		//	pgns.push_back(pgns[i]);
		// distribute all the polygons over the mesh, this step transforms polygons in 2D to 3D
		// distribute the same number of points over the mesh
		unsigned int totalSum = 0;
		for (unsigned int n = 0; n < multiplicity.size(); n++)
			totalSum += multiplicity[n];
		pgns_3.resize(totalSum);
		polygonNormals.resize(totalSum);
		bList.resize(totalSum);
		srand(time(NULL));
		unsigned int idx = 0;
		for ( unsigned int i = 0; i < pgns.size(); i++)
		{
			double a = std::fabs(pgns[i].area());
			double fa = std::fabs(meshDom[0].area());// need to be more robust
			double factor = std::min<double>(a/fa, fa/a);
			shrinkPolygon(&pgns[i], std::sqrt(factor)*0.5);
			for (unsigned int j = 0; j < multiplicity[i]; j++)
			{
				int fidx = ::rand() % meshDom.size();
				while ( std::find(facetTranversed.begin(), facetTranversed.end(), fidx) != facetTranversed.end() ) fidx = ::rand() % meshDom.size();
				facetTranversed.push_back(fidx);
				const Facet& f = meshDom[fidx];
				vec3 randCenter = (f.vertex[0] + f.vertex[1] + f.vertex[2])/3.0;
				vec3 n = -meshDom[fidx].normal();
				//compute the shrink factor to guarantee that the polyon is totally inside a single mesh triangle
				movePolygon2Mesh(pgns[i], randCenter, n, &pgns_3[idx]);
				for (unsigned int k = 0; k < pgns_3[idx].size(); k++)
				{
					const Segment_3& edge = pgns_3[idx][k];
					const Point_3& s = edge.source(), t = edge.target();
					if (is_nan(s.x()) || is_nan(s.y()) || is_nan(s.z()) || is_nan(t.x()) || is_nan(t.y()) || is_nan(t.z()))
					{
						cout << "============ " << i << " ============\n";
						system("pause");
					}
				}
				//movePolygon2Mesh(pgns[i], f, &pgns_3[i]);
				polygonNormals[idx] = n;
				idx++;
			}
			
		}
		genSphericalBounder();
	}
	void Doc::loadMeshDomain(const string& bdFileName)
	{
		vec3 g ;
		//double R ;
		int len = bdFileName.length() ;
		bool is_obj = ( bdFileName.substr(len-3, 3)==string("obj"));

		if(is_obj) 
		{
			//boundary_.load(filename, &boundary_planes_) ;
			meshDom.load(bdFileName);
		} 
		else 
		{
			meshDom.load_tri(bdFileName) ;
		}
		powerdiagram_ = new PowerDiagram_3(&meshDom);
	}
	void Doc::clear()
	{
		meshDom.clear_all();
		powerdiagram_->clear();
		delete powerdiagram_;
		pgns.clear();
		for (vector<vector<Segment_3>>::iterator ap = pgns_3.begin(); ap!=pgns_3.end(); ap++)
			ap->clear();
		pgns_3.clear();
		optimalPolygons.clear();
		clearBounder();
		bList.clear();
		optimalConfig.clear();

		glut_viewer_redraw();
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
	void Doc::shrinkPolygon(Polygon_2 *pgn, double f)
	{
		if (f <=0.0)
			return;
		vec3 center = polygonCenter(*pgn);
		Vector_2 translation(center.x, center.y);
		Transformation_2 move_to_origin(CGAL::TRANSLATION, -translation);
		Transformation_2 move_back(CGAL::TRANSLATION, translation);
		Transformation_2 shrink(CGAL::SCALING, f);
		Transformation_2 t = move_back*(shrink*move_to_origin);
		*pgn = CGAL::transform(t, *pgn);
	}
	// assume the polygon pgn is on the xy-plane, n is unit normal
	void Doc::movePolygon2Mesh(const Polygon_2& pgn, const vec3& p, const vec3& n, vector<Segment_3> *mvPgn)
	{
		mvPgn->clear();
		vec3 center = polygonCenter(pgn);
		Transformation_3 move_to_origin(CGAL::TRANSLATION, Vector_3(-center.x, -center.y, -center.z));
		// assume the normal vector n form angle theta with z axis, its projection on xOy form angle phi with x axis
		//sin and cos for transformation
		double cos_theta = n.z;
		double sin_theta = std::sqrt(n.x*n.x + n.y*n.y);
		double cos_phi = n.x / std::sqrt(n.x*n.x + n.y*n.y);
		double sin_phi = std::abs(n.y) / std::sqrt(n.x*n.x + n.y*n.y);
		// construct the rotation matrix
		Transformation_3 align_normal(CGAL::IDENTITY);
		if (n.x != 0.0 || n.y != 0.0) // if the facet has already had the same normal as z axis
		{
			if (n.y > 0.0)
				align_normal = Transformation_3(cos_phi*cos_theta, -sin_phi, cos_phi*sin_theta, sin_phi*cos_theta, 
				cos_phi, sin_phi*sin_theta, -sin_theta, 0.0, cos_theta);
			else
				align_normal = Transformation_3(cos_phi*cos_theta, sin_phi, cos_phi*sin_theta, -sin_phi*cos_theta, 
				cos_phi, -sin_phi*sin_theta, -sin_theta, 0.0, cos_theta);
		}

		//Vector_3 translation(rotated_center_3d, to_cgal(p));
		Transformation_3 move_to_p(CGAL::TRANSLATION, Vector_3(p.x, p.y, p.z));
		Transformation_3 t = move_to_p*(align_normal*move_to_origin);
		//Point_3 rotated_center_3d = t(to_cgal(center));
		//vec3 vr = to_geex(rotated_center_3d);
		for (Polygon_2::Edge_const_iterator ei = pgn.edges_begin(); ei != pgn.edges_end(); ei++)
		{
			Segment_3 edge_3d(Point_3(ei->source().x(), ei->source().y(), 0.0), Point_3(ei->target().x(), ei->target().y(), 0.0));
			mvPgn->push_back(edge_3d.transform(t));
		}
	}

	void Doc::movePolygon2Mesh( const Polygon_2& pgn, const Facet& f, vector<Segment_3> *mvPgn )
	{
		// move facet to the origin
		vec3 fCenter = f.center();
		Facet fo(f.vertex[0]-fCenter, f.vertex[1]-fCenter, f.vertex[2]-fCenter);
		vec3 n = fo.normal();
		movePolygon2Mesh(pgn, fCenter, n, mvPgn);
	}

	void Doc::adjustPolygons( double f )
	{
		if (f == 1.0 || f <= 1.0e-3)
			return;
		for ( unsigned int i = 0; i < pgns_3.size(); i++)
		{
			vector<Segment_3>& ap = pgns_3[i];
			vec3 midp(0.0, 0.0, 0.0);
			for (vector<Segment_3>::const_iterator it = ap.begin(); it!=ap.end(); it++)
				midp += to_geex(it->source());
			midp /= ap.size();
			Transformation_3 move_to_origin(CGAL::TRANSLATION, Vector_3(-midp.x, -midp.y, -midp.z));
			Transformation_3 move_back(CGAL::TRANSLATION, Vector_3(midp.x, midp.y, midp.z));
			Transformation_3 scale(CGAL::SCALING, f);
			Transformation_3 t = move_back*(scale*move_to_origin);
			for (vector<Segment_3>::iterator it = ap.begin(); it!=ap.end(); it++)
				*it = it->transform(t);
		}
		clearBounder();
		genSphericalBounder();
	}
	void Doc::genSphericalBounder()
	{
		if (pgns_3.size() == 0)
			return;
		for (unsigned int i = 0; i < pgns_3.size(); i++)
		{
			SphericalBounder *bounder = new SphericalBounder(pgns_3[i]);
			bounder->genMarginBounder(24);//10, 14
			bList[i] = bounder;
		}
		// form the power diagram of these sphere8
		std::cerr<< "Start inserting points into power diagram......\n";
		powerdiagram_->begin_insert();
		for (unsigned int i = 0; i < bList.size(); i++)
		{
			SphericalBounder* aBounder = bList[i];
			//powerdiagram_->insert_grouped_points(aBounder->getMarginBounder());
			Vector_3 cgalVec(polygonNormals[i].x, polygonNormals[i].y, polygonNormals[i].z);
			//if (is_nan(polygonNormals[i].x) || is_nan(polygonNormals[i].y) || is_nan(polygonNormals[i].z))
			//{
			//	std::cout << "========== " << i << " =========="<<std::endl;
			//	system("pause");
			//}
			Plane_3 tangentPlane(pgns_3[i][0].source(), cgalVec);
			double sx = pgns_3[i][0].source().x(), sy = pgns_3[i][0].source().y(), sz = pgns_3[i][0].source().z();
			using Geex::Numeric::is_nan;
			if (is_nan(tangentPlane.a()) || is_nan(tangentPlane.b()) || is_nan(tangentPlane.c()) || is_nan(tangentPlane.d()))
			{
				std::cout << "========== " << i << " =========="<<std::endl;
				system("pause");
			}
			powerdiagram_->insert_grouped_points(aBounder->getMarginBounder(), tangentPlane/*Plane_3(pgns_3[i][0].source(), cgalVec)*/);
		}
		powerdiagram_->end_insert();
		std::cerr<<"Finish inserting points.\n";
	}
	inline void Doc::clearBounder()
	{
		for (vector<SphericalBounder*>::iterator it = bList.begin(); it!=bList.end(); it++)
			delete *it;
	}
}//namespace Geex
