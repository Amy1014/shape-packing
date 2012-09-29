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

#include <Geex/mathematics/glsl_linear.h>
#include <Geex/CVT/geometry.h>
//#include "oriented_line.h"
#include <Geex/CVT/ann_kdtree.h>


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
			return normalize( cross(vertex[2] - vertex[0], vertex[1] - vertex[0]) ) ; 
		}

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
			return (vertex_weight[0]+vertex_weight[1]+vertex_weight[2])/3.0 ;
		}

		real weight(int v0, int v1, vec3& p) {
			real u = distance(p, vertex[v0])/distance(vertex[v0], vertex[v1]) ;
			return (1.-u)*vertex_weight[v0] + u*vertex_weight[v1] ;
		}

		bool has_vertex(int vidx) { return vertex_index[0]==vidx || vertex_index[1]==vidx || vertex_index[2]==vidx ; }
		int  find_edge(int v0, int v1) {
			for(int i=0; i<3; ++i) {
				if(vertex_index[i]==v0 && vertex_index[(i+1)%3]==v1 
					|| vertex_index[i]==v1 && vertex_index[(i+1)%3]==v0) 
					return (i+2)%3 ;
			}
			return -1 ;
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

		void convex_clip(TriMesh& to, const Plane<real>& P, bool close = true) const ;
		bool contains(const vec3& p) const ;

		void normalize(vec3& g, double& R, bool disable_normalize_ = false) ;
		void compute_bounding_box() ;

		int nb_vertices() const { return nb_vertices_ ; }
		MeshVertex& vertex(int idx) { return vertices_[idx]; }
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
		real xmin() { return bbox_[0][0] ; } 
		real xmax() { return bbox_[1][0] ; } 
		real ymin() { return bbox_[0][1] ; } 
		real ymax() { return bbox_[1][1] ; } 
		real zmin() { return bbox_[0][2] ; } 
		real zmax() { return bbox_[1][2] ; } 

		void build_kdtree() ;
		vec3 project_to_mesh(const vec3& pt, vec3& norm) ;

		bool load_density(const std::string& filename, double gamma=1.) ;
		void set_uniform_density() ;
		void merge_same_vertices(const int res, const real toler=1e-8) ; 
		void set_facet_neighbor() ;
		void set_feature_edgeflag() ; // 

		void measure_quality() ;
		void perturb_facet_vertices(double perturb_ratio = 0.001) ;
		void save_tri(const std::string& filename) ;
		void load_tri(const std::string& filename) ;


	protected:
		int                       nb_vertices_ ;
		std::vector<MeshVertex>   vertices_ ;
		vec3                      bbox_[2] ; // bounding box of boundary mesh

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
