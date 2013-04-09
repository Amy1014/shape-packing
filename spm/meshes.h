/*
*  _____ _____ ________  _
* /  __//  __//  __/\  \//
* | |  _|  \  |  \   \  / 
* | |_//|  /_ |  /_  /  \ 
* \____\\____\\____\/__/\\
*
* Graphics Environment for EXperimentations.
*  Copyright (C) 2006 INRIA - Project ALICE
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program; if not, write to the Free Software
*  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*
*  If you modify this software, you should include a notice giving the
*  name of the person performing the modification, the date of modification,
*  and the reason for such modification.
*
*  Contact: 
*
*     ALICE Project - INRIA
*     INRIA Lorraine, 
*     Campus Scientifique, BP 239
*     54506 VANDOEUVRE LES NANCY CEDEX 
*     FRANCE
*
*  Note that the GNU General Public License does not permit incorporating
*  the Software into proprietary programs. 
*/

#ifndef GEOM_TYPES
#define GEOM_TYPES

#include <set>
#include <map>
#include <cctype>
#include <CGAL/MP_Float.h> // by Wenchao Hu, for more accurate computation of cross product
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Vector_3.h>
#include <Geex/mathematics/glsl_linear.h>
#include <Geex/CVT/geometry.h>
//#include "oriented_line.h"
#include <Geex/CVT/ann_kdtree.h>
//#include "spm_cgal.h"
using CGAL::normal;
using CGAL::MP_Float;
using CGAL::to_double;


namespace Geex {


	class Convex : public std::vector< Plane<real> > {
	public:
		inline bool contains(const vec3& p) {
			for(unsigned int i=0; i<size(); i++) {
				if((*this)[i].side(p) < 0) {
					return false ;
				}
			}
			return true ;
		}
	} ;

	class MeshVertex {
	public:

		MeshVertex(const vec3& p) {
			pos_ = p;
			tex_coord_ = 0. ;
		}

		MeshVertex() {
			pos_ = vec3(0, 0, 0) ;
			tex_coord_ = 0. ;
		}

		vec3& point()    { return pos_ ; } 
		const vec3& point() const { return pos_;}
		double& weight() { return weight_ ; }
		//MeshVertex& operator = (const MeshVertex& rhs) {
		//    pos_ = rhs.pos_ ;
		//    weight_ = rhs.weight_ ;
		//    faces_ = rhs.faces_ ;
		//    return *this ;
		//}

	public:
		vec3             pos_;
		double           weight_ ;
		std::vector<int> faces_;
		real             tex_coord_ ;
		double			curvature;
		//	int              gidx_;  // index of nearest delaunay vertex
	};

	//class MeshEdge {
	//public:
	//	int halfedge_[2];
	//	int vertex_[2];
	//	bool flag_; // 
	//};

	//class MeshHalfEdge {
	//public:
	//	int fidx_; // equals -1 if no face
	//	int vidx_; // the vertex which halfedge points to
	//	int eidx_; // the edge idx which halfedge belongs to
	//	int next_; // index of next halfedge
	//	int prev_; // previous halfedge
	//	int oppo_; // opposite halfedge
	//	bool flag_; //
	//};

	class LineSeg {
	public:
		LineSeg(const vec3& p0, const vec3& p1, real w0=1, real w1=1, int n0=-1, int n1=-1, int pidx=-1) {
			v0_ = p0 ;
			v1_ = p1 ;
			w0_ = w0 ;
			w1_ = w1 ;
			nclip_[0] = n0 ;
			nclip_[1] = n1 ;
			pidx_ = pidx ;
		}

		double length() { return (v1_-v0_).length() ; }
		vec3   midpoint() { return 0.5*(v0_+v1_) ; }
		int    parentidx() { return pidx_ ; }

	public:
		vec3 v0_, v1_ ;
		real w0_, w1_ ;  // weigth of each vertex
		int  nclip_[2] ; // neighbor clipping seed
		int  pidx_ ;     // parent feature line index
	} ;

	class Facet {
	public:
		enum VertexType {V_VERTEX, V_EDGE, V_FACE, V_AUX} ; // three types of vertex: 
		enum EdgeType   {E_EDGE, E_CLIP, E_AUX, E_FEA} ;   // original edge, edge by plane/face intersect, aux edge, feature edge

		Facet(
			const vec3& p1, const vec3& p2, const vec3& p3,
			int e1, int e2, int e3,
			int i1, int i2, int i3
			) {
				vertex[0] = p1 ;
				vertex[1] = p2 ;
				vertex[2] = p3 ;
				edge_flag[0] = e1 ;
				edge_flag[1] = e2 ;
				edge_flag[2] = e3 ;
				vertex_index[0] = i1 ;
				vertex_index[1] = i2 ;
				vertex_index[2] = i3 ; // fix bug. dmyan
				face_index[0] = -1 ;
				face_index[1] = -1 ;
				face_index[2] = -1 ;
				//	    set_vertex_type(V_VERTEX, V_VERTEX, V_VERTEX) ;
				//	    set_edge_type(E_EDGE, E_EDGE, E_EDGE) ;
				ftype_ = -1 ; // F_AUX
				//	    clear_clip_planes() ;
				set_clip_neighbor(-1, -1, -1) ;
				set_vertex_weight(1, 1, 1) ;
		}

		Facet(
			const vec3& p1, const vec3& p2, const vec3& p3,
			int e1, int e2, int e3
			) {
				vertex[0] = p1 ;
				vertex[1] = p2 ;
				vertex[2] = p3 ;
				edge_flag[0] = e1 ;
				edge_flag[1] = e2 ;
				edge_flag[2] = e3 ;
				face_index[0] = -1 ;
				face_index[1] = -1 ;
				face_index[2] = -1 ;
				//	    flag = false ;
				ftype_ = -1 ; // F_AUX
				//	    clear_clip_planes() ;
				set_clip_neighbor(-1, -1, -1) ;
				set_vertex_weight(1, 1, 1) ;
		}

		Facet(const vec3& p1, const vec3& p2, const vec3& p3) {
			vertex[0] = p1 ;
			vertex[1] = p2 ;
			vertex[2] = p3 ;		
			face_index[0] = -1 ;
			face_index[1] = -1 ;
			face_index[2] = -1 ;
			ftype_ = -1 ; // F_AUX
			clear_edge_flags() ;
			set_clip_neighbor(-1, -1, -1) ;
			set_vertex_weight(1, 1, 1) ;
		}

		Facet() { 
			face_index[0] = -1 ;
			face_index[1] = -1 ;
			face_index[2] = -1 ;
			ftype_ = -1 ; // F_AUX
			clear_edge_flags() ;
			set_clip_neighbor(-1, -1, -1) ;
			set_vertex_weight(1, 1, 1) ;
		} 

		void clear_edge_flags() {
			edge_flag[0] = 0 ;
			edge_flag[1] = 0 ;
			edge_flag[2] = 0 ;
		}

		void set_clip_neighbor(int n0, int n1, int n2) {
			nclip_[0] = n0 ;
			nclip_[1] = n1 ;
			nclip_[2] = n2 ;
		}

		void set_vertex_weight(real w0, real w1, real w2) {
			vertex_weight[0] = w0 ;
			vertex_weight[1] = w1 ;
			vertex_weight[2] = w2 ;
		}

		void set_face_neighbor(int f0, int f1, int f2) {
			face_index[0] = f0 ; 
			face_index[1] = f1 ;
			face_index[2] = f2 ;
		}

		VertexType intersect_type(EdgeType et) {
			switch(et) {
		case E_EDGE:
			return V_EDGE ;
		case E_CLIP:
			return V_FACE ;
		case E_AUX:
			return V_AUX ;
			}
		}

		void get_edge(int i, vec3 seg[2]) {
			seg[0] = vertex[(i+1)%3] ;
			seg[1] = vertex[(i+2)%3] ;
		}
		vec3 normal() const { /*why?? dmyan*/
			//return normalize( cross(vertex[2] - vertex[0], vertex[1] - vertex[0]) ) ; 
			return norm;
		}
		
		// by Wenchao Hu
		vec3 MP_Float_normal() const
		{
			MP_Float p0x(vertex[0].x), p0y(vertex[0].y), p0z(vertex[0].z),
					 p1x(vertex[1].x), p1y(vertex[1].y), p1z(vertex[1].z), 
					 p2x(vertex[2].x), p2y(vertex[2].y), p2z(vertex[2].z);
			MP_Float p01x = p1x - p0x, p01y = p1y - p0y, p01z = p1z - p0z, 
					 p02x = p2x - p0x, p02y = p2y - p0y, p02z = p2z - p0z;
			MP_Float cpx =  p01y*p02z - p01z*p02y, cpy = p01z*p02x - p01x*p02z, cpz = p01x*p02y - p01y*p02x;
			MP_Float l2 = cpx*cpx + cpy*cpy + cpz*cpz;
			real x2 = to_double(cpx*cpx)/to_double(l2), y2 = to_double(cpy*cpy)/to_double(l2), z2 = to_double(cpz*cpz)/to_double(l2);
			return vec3(std::sqrt(x2), std::sqrt(y2), std::sqrt(z2));

		}

		//vec3 CGAL_Normal() const
		//{
		//	typedef CGAL::V
		//	Vector_3 n = CGAL::normal(to_cgal(vertex[0]), to_cgal(vertex[1]), to_cgal(vertex[2]));
		//	double sl = n.x()*n.x() + n.y()*n.y() + n.z()*n.z();
		//	return vec3(n.x(), n.y(), n.z());
		//}

		vec3 xdir() const {
			return normalize( vertex[1] - vertex[0] ) ;
		}

		vec3 ydir() const {
			return cross(xdir() ,normal() ) ;
		}

		vec3 org() const {
			return vertex[0] ;
		}

		vec3 center() const {
			return 1.0/3.0*(vertex[0]+vertex[1]+vertex[2]) ;
		}

		Plane<real> plane() const { 
			return Plane<real>(vertex[0], normal()) ; 
		}

		real weight() {
			//return (vertex_weight[0]+vertex_weight[1]+vertex_weight[2])/3.0 ;
			double avg_cur = (vertex_weight[0]+vertex_weight[1]+vertex_weight[2])/3.0;
			return avg_cur*avg_cur;
		}

		real weight(int v0, int v1, vec3& p) {
			real u = distance(p, vertex[v0])/distance(vertex[v0], vertex[v1]) ;
			return (1.-u)*vertex_weight[v0] + u*vertex_weight[v1] ;
		}

		//inline const vec3& vertex(int i) const { return vertex[i%3]; }

		bool has_vertex(int vidx) { return vertex_index[0]==vidx || vertex_index[1]==vidx || vertex_index[2]==vidx ; }
		int  find_edge(int v0, int v1) {
			for(int i=0; i<3; ++i) {
				if(vertex_index[i]==v0 && vertex_index[(i+1)%3]==v1 
					|| vertex_index[i]==v1 && vertex_index[(i+1)%3]==v0) 
					return (i+2)%3 ;
			}
			return -1 ;
		}

		real area() const
		{
			vec3 ba = vertex[1] - vertex[0];
			vec3 ca = vertex[2] - vertex[0];
			vec3 cp = cross(ba, ca);
			return 0.5 * cp.length();
		}

		Facet& operator = (const Facet& rhs) {
			ftype_ = rhs.ftype_ ;
			for(int i=0; i<3; ++i) {
				vertex_index[i] = rhs.vertex_index[i] ;
				vertex[i] = rhs.vertex[i] ;
				vertex_weight[i] = rhs.vertex_weight[i] ;
				face_index[i] = rhs.face_index[i] ;
				edge_flag[i] = rhs.edge_flag[i] ;
				nclip_[i] = rhs.nclip_[i] ;
			}

			return *this ;
		} ;
		int  vertex_index[3] ; // only used for boundary triangle
		vec3 vertex[3] ;
		real vertex_weight[3] ;
		int  face_index[3] ; // indices of three neighbor faces
		int  edge_flag[3] ;
		int  nclip_[3] ;
		int  ftype_ ; // used in volume clipping. face type could be F_BDY, F_AUX and F_CLIP
		vec3 norm;
	} ;

	// FaceVertex is used to merge the same vertices in clipped mesh
	class FaceVertex {
	public:
		vec3 pos_ ;
		int  fidx_ ;  // face index it belongs to
		int  fvidx_ ; // the index in the face  0, 1, or 2

	public:				
		FaceVertex() {} ;

		FaceVertex(vec3& pos, int fidx, int fvidx) {
			pos_  = pos ;
			fidx_ = fidx ;
			fvidx_  = fvidx ; 
		}

		bool operator < (const FaceVertex& rhs) const {
			if(pos_.x < rhs.pos_.x) 
				return true ;
			else if (pos_.x > rhs.pos_.x) 
				return false ;
			else {
				if(pos_.y < rhs.pos_.y) 
					return true ;
				else if (pos_.y > rhs.pos_.y) 
					return false ;
				else {
					if(pos_.z < rhs.pos_.z) 
						return true ;
					else if (pos_.z > rhs.pos_.z) 
						return false ;
					else 
						return false ;
				}
			}
		}
	} ; 

	class Grid { // used to merge same vertices
	public:
		void add_face_vertex(FaceVertex fv) {
			fvs_.push_back(fv) ;
		}
	public:
		std::vector<FaceVertex> fvs_ ;
	} ;

	using std::pair;
	using std::set;
	class TriMesh : public std::vector<Facet> {
	public:
		enum ClippingMode { CLIP_BNDRY=0, CLIP_CONVEX=1} ;

		TriMesh() ;
		~TriMesh() ;

		// reads Alias|Wavefront .obj file format.
		void load(const std::string& filename, Convex* planes = nil) ;
		void load(std::istream& in, Convex* planes = nil) ;
		bool load_feature(const std::string& filename) ;

		void save(const std::string& filename) ;
		void save(std::ostream& out) ;
		void save_obj(const std::string& filename) ;

		void convex_clip(TriMesh& to, const Plane<real>& P, ClippingMode mode = CLIP_CONVEX) const {
			switch(mode) {
			case CLIP_BNDRY:
				convex_clip(to, P, false) ;
				break ;
			case CLIP_CONVEX:
				convex_clip(to, P, true) ;
				break ;
			}
		}

		void flip_normals();

		void convex_clip(TriMesh& to, const Plane<real>& P, bool close = true) const ;
		bool contains(const vec3& p) const ;

		void normalize(vec3& g, double& R, bool disable_normalize_ = false) ;
		void compute_bounding_box() ;

		int nb_vertices() const { return nb_vertices_ ; }
		MeshVertex& vertex(int idx) { return vertices_[idx]; }
		inline bool is_on_boundary(int idx) const { return vert_on_boundary[idx]; }
		inline bool is_on_feature(int idx) const { return vert_on_feature[idx]; }
		inline const MeshVertex& vertex(int idx) const { return vertices_[idx]; }
		void add_vertex(vec3 p) { vertices_.push_back(MeshVertex(p)); }
		void clear_all() {
			clear();
			vertices_.clear();
			nb_vertices_ = 0;
		}

		void get_bounding_box(vec3& minb, vec3& maxb) { 
			minb = bbox_[0] ;
			maxb = bbox_[1] ;
		}
		real xmin() const { return bbox_[0][0] ; } 
		real xmax() const { return bbox_[1][0] ; } 
		real ymin() const { return bbox_[0][1] ; } 
		real ymax() const { return bbox_[1][1] ; } 
		real zmin() const { return bbox_[0][2] ; } 
		real zmax() const { return bbox_[1][2] ; } 

		void build_kdtree() ;
		vec3 project_to_mesh(const vec3& pt, vec3& norm) ;
		vec3 project_to_mesh(const vec3& pt, vec3& norm, int &nearestFacet);

		bool load_density(const std::string& filename, double gamma=1.) ;
		void set_uniform_density() ;
		void merge_same_vertices(const int res, const real toler=1e-8) ; 
		void set_facet_neighbor() ;
		void set_feature_edgeflag() ; // 

		void measure_quality() ;
		void perturb_facet_vertices(double perturb_ratio = 0.001) ;
		void save_tri(const std::string& filename) ;
		void load_tri(const std::string& filename) ;

		// for feature extracting, by Wenchao Hu
		void setFeatureAngle(double ang);
		inline double getFeatureAngle() const { return featureAngleCriterion; }
		inline double getHighCurvaturePercent() const { return highCurvatureCriterion; }
		void setHighCurvaturePercent(double percent);
		const set<int>& getHighCurvatureFacets() const {return facetHighCurvature;}
		const vector<pair<int, int>>& getFeatureFacets() const { return featureFacets; }
		const vector<pair<int, int>>& getFeatureEdges() const { return featureEdges; }
		inline double getMaxFacetWeight() const { return maxFacetWeight; }
		inline double getMinFacetWeight() const { return minFacetWeight; }

		inline double curvature_at_vertex(int idx) const { return vertices_[idx].curvature; }
		//inline double get_max_curvature() const { return max_vert_curvature; }
		//inline double get_min_curvature() const { return min_vert_curvature; }
 
	protected:
		int                       nb_vertices_ ;
		std::vector<MeshVertex>   vertices_ ;
		vec3                      bbox_[2] ; // bounding box of boundary mesh
	private:
		class PairIntIntCmp
		{
		public:
			inline bool operator()(const std::pair<int, int>& p0, const std::pair<int, int>& p1) const
			{
				if (p0.first < p1.first)
					return true;
				else if ( p0.first > p1.first)
					return false;
				else if ( p0.second < p1.second)
					return true;
				else 
					return false;			
			}
		};
		class PairDoubleIntCmp
		{
		public:
			inline bool operator()(const std::pair<double, int>& p0, const std::pair<double, int>& p1) const
			{
				return p0.first < p1.first;
			}
		};
		struct BihedralAngle
		{
			//friend class TriMesh;
			int findex0;
			int findex1;
			int vindex0;
			int vindex1;
			double cosAngle;
			BihedralAngle( int faceIndex0, int faceIndex1, double angleCos ) : findex0(faceIndex0), findex1(faceIndex1), 
																vindex0(-1), vindex1(-1), cosAngle(angleCos) {}
			BihedralAngle( int faceIndex0, int faceIndex1, int vertIndex0, int vertIndex1, double angleCos )
					:findex0(faceIndex0), findex1(faceIndex1), vindex0(vertIndex0), vindex1(vertIndex1), cosAngle(angleCos) {}
		};
		std::set<pair<int, int>, PairIntIntCmp> boundaryEdges;//record the boundary edge, by Wenchao Hu
		//std::set<pair<int, int>, PairIntIntCmp> bihedralAngles;//dihedral angles to extract sharp features
		vector<BihedralAngle> bihedralAngles;
		double featureAngleCriterion;//the bihedral angle criterion determining sharp feature edges, in degree [0, 180]
		double highCurvatureCriterion;
		//std::multimap<std::pair<int, int>, int, Pa
		vector<pair<int, int>> featureFacets;//updated whenever featureAngleCriterion is set by setFeatureAngle
		vector<pair<int, int>> featureEdges;
		set<int> facetHighCurvature;
		double maxFacetWeight;
		double minFacetWeight;
		//double max_vert_curvature;
		//double min_vert_curvature;
		vector<pair<double, int>> facetWeight;
		static double PI;

		// record whether a vertex is on boundary or feature
		vector<bool> vert_on_boundary;
		vector<bool> vert_on_feature;

	public:
		typedef std::set<std::pair<int, int>, PairIntIntCmp> BoundaryEdgeSet;//Wenchao
		const BoundaryEdgeSet& getBoundary() const { return boundaryEdges; } //Wenchao	

	public:
		std::vector<int>          corners_ ;  // constrained vertex indices
		std::vector<std::vector<int> > fealines_ ;
		std::vector<LineSeg>      features_ ;    // feature edges
		ANNKdTree_3               *kdtree_ ;
		int                       samplerate_ ;
		
		
	} ;

	//class PluckerSegment {
	//public:
	//	PluckerSegment(
	//		const vec3& p1, const vec3& p2
	//		) : L(p1, p2){
	//			P[0] = p1 ;
	//			P[1] = p2 ;
	//	}
	//	vec3 P[2] ;
	//	OrientedLine<real> L ;
	//} ;


	//class PluckerFacet {
	//public:
	//	PluckerFacet(const vec3& v1, const vec3& v2, const vec3& v3) {
	//		P = Plane<real>(v1, cross(v2-v1, v3-v1)) ;
	//		E[0] = OrientedLine<real>(v2,v3) ;
	//		E[1] = OrientedLine<real>(v3,v1) ;
	//		E[2] = OrientedLine<real>(v1,v2) ;
	//	}

	//	double seg_side(const PluckerSegment& e, const OrientedLine<real>& l) const {
	//		return side(l,e.L) ;
	//	}

	//	double plane_side(const vec3& p, const Plane<real>& P) const {
	//		return P.side(p) ;
	//	}

	//	//bool intersects(const PluckerSegment& e) const {
	//	//	double s1 = plane_side(e.P[0], P) ;
	//	//	double s2 = plane_side(e.P[1], P) ;
	//	//	if(s1 > 0 && s2 > 0 || s1 < 0 && s2 < 0) { return false ; }
	//	//	double sp1 = seg_side(e, E[0]) ;
	//	//	double sp2 = seg_side(e, E[1]) ;
	//	//	double sp3 = seg_side(e, E[2]) ;
	//	//	return (
	//	//		(sp1 > 0 && sp2 > 0 && sp3 > 0) ||
	//	//		(sp1 < 0 && sp2 < 0 && sp3 < 0)
	//	//		) ;
	//	//}

	//	Plane<real> P ;
	//	//OrientedLine<real> E[3] ;
	//} ;

	//class PluckerMesh : public std::vector<PluckerFacet> {
	//public:
	//	bool contains(const vec3& P) const ;
	//} ;

}

#endif
