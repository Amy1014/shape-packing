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

#ifndef __RVD_CCVT_H__
#define __RVD_CCVT_H__

#include "meshes.h"
#include "rvc_ccvt.h"
#include "geometry_ccvt.h"
#include <Geex/graphics/opengl.h>

namespace Geex {

	class Delaunay ;
	class TopoPolyMesh ;
	class DelaunaySkeleton ;

	enum RVDMode {TRI_KD_INCIDENT, POLY_KD_INCIDENT, TRI_LINEAR, POLY_LINEAR, TRI_KD_SPLIT, PARA_TRI, PARA_POLY} ;
#define ClipModeNames "Tri_Incident" | "Poly_Incident" | "Tri_Linear" | "Poly_Linear" | "Tri_Split" | "Para_Tri" | "Para_Poly" 

	class RestrictedVoronoiDiagram {
	public:

		RestrictedVoronoiDiagram() : dt_(nil), mesh_(nil), kdtree_(nil), rvd_mode_(TRI_KD_INCIDENT) {}
		RestrictedVoronoiDiagram(MyDelaunay* dt, Geex::TriMesh* mesh, RVDMode mode=TRI_KD_INCIDENT) ;
		~RestrictedVoronoiDiagram() ;

		GLenum& rvd_mode()                            { return rvd_mode_ ; } 
		std::vector<RestrictedVoronoiCell>& rvcells() { return rvcells_ ; }
		Geex::TriMesh& dual_mesh()                    { return remesh_ ; }

		void compute(Array1d<Geex::TopoPolyMesh>& parts, MyDelaunay* dt) ;
		void compute(std::vector<Geex::TriMesh>& parts, MyDelaunay* dt) ;
		void compute(Geex::TriMesh* mesh, MyDelaunay* dt) ;
		void extract_mesh() ;
		void extract_dual() ;
        void voronoi_filtering(double input_theta_degree); // 0 <= input_theta_degree <= 90 degree
        void alpha_shape();

		// -------------------------------------------------------------------------------------
		// ---------------------------- Intefaces for CVT iteration ----------------------------
		// -------------------------------------------------------------------------------------

		inline void get_fg(int vidx, double& f, vec3& grad_f)  { rvcells_[vidx].get_fg(f, grad_f) ; }
		inline void get_centroid(int vidx, vec3& g, double& V) { rvcells_[vidx].get_centroid(g, V) ; }
		inline void get_fgw(int vidx, double& f, vec3& grad_f)  { rvcells_[vidx].get_fgw(f, grad_f) ; }
		inline void get_centroidw(int vidx, vec3& g, double& V) { rvcells_[vidx].get_centroidw(g, V) ; }
        inline vec3 project_to_mesh(const vec3& pt, vec3& norm) { return mesh_->project_to_mesh(pt, norm) ; } 

		//interface for cvt
		void get_centroids(std::vector<vec3>& gs, std::vector<double>& Vs, bool weight=false) ;
		void get_constrained_centroids(std::vector<vec3>& gs, std::vector<double>& Vs, bool weight=false) ;
//		void get_gfs(double* f, std::vector<vec3>& grad_fs, bool constrained=false) ;
		void get_fgs(double* f, double* grad_fs, bool weight=false) ;
		void get_constrained_fgs(int N, double* x, double* f, double* g, bool weight=false) ;

	protected:
        void init_vcells() ; // initialize Voronoi cells
		void build_kdtree() ;
		void collect_crease_segments() ;

		// -------------------------------------------------------------------------------------
		// --------------------- surface RVD computation ---------------------------------------
		// -------------------------------------------------------------------------------------

		// linear rvd computation by polygon clipping
		void clip_poly_linear(TopoPolyMesh* poly, MyDelaunay* dt, DelaunaySkeleton* skel) ;
		// linear rvd computation by polygon clipping
		void clip_poly_kdtree(TopoPolyMesh* poly, MyDelaunay* dt, DelaunaySkeleton* skel, ANNKdTree_3* kdtree) ;


        // linear voronoi clip by propagation
        void voronoi_clip_linear() ;
        void clip_surface_linear() ;
        void clip_tri_by_cell(int fidx, int cell_idx, std::vector<Geex::Facet>& from) ;
        void clip_tri_by_plane(Geex::Facet& F, const Geex::Plane<real>& P, int oppvidx, std::vector<Geex::Facet>& inlist, std::vector<Geex::Facet>& outlist) ;

        // voronoi clip by kd-tree query and using incident cells
        void clip_surface_mlogn() ;

        // voronoi clip by kd-tree query and triangle splitting
        void voronoi_clip_splitting() ;
        void clip_surface_splitting() ;
        void subdivide_tri(const Geex::Facet& F, int fidx) ;
        void subdivide_mesh(int vid, Geex::TriMesh& mesh, int fidx) ;

		//----------feature & corner clipping algorithm has no use now  -----------------
		//----------keeping for reference -----------------------------------------------

		//------------------------ feature edge clipping algorithm ----------------------
        void clip_feature() ;
        void clip_edge_by_cell(int eidx, int vidx,  std::set<int> edge_cells) ;
        void clip_edge_by_plane(Geex::LineSeg& E, const Geex::Plane<real>& P, int oppvidx, std::vector<Geex::LineSeg>& inlist) ;
        void refine_feature() ; // insert new seeds if one seed occupied more than 1 feature lineseg
        //------------------------ corner clipping algorithm ----------------------
        void clip_corner() ; // used for identify corners after load data

	protected:
		MyDelaunay             *dt_ ;
		Geex::TriMesh          *mesh_ ;
		ANNKdTree_3            *kdtree_ ; // only used for kd-tree based clipping
		Geex::TriMesh          remesh_ ;
		GLenum                 rvd_mode_ ;
		std::vector<RestrictedVoronoiCell>  rvcells_ ;

        static int configs[8][3] ;
        static int configv[16][4] ;

		friend class RestrictedVoronoiCell ; 
		friend class MyDelaunay ;
	} ;

}


#endif
