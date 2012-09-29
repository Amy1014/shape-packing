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

#ifndef __RVC_H__
#define __RVC_H__

#include "delaunay_cgal.h" 
#include "meshes.h"
#include <Geex/CVT/topo_poly_mesh.h>

namespace Geex {

	class MyDelaunay ;

	class RestrictedVoronoiCell {
	public:
		RestrictedVoronoiCell() {}
		RestrictedVoronoiCell(Vertex_handle vh) : vh_(vh) {}

		bool insert_clipped_tet(int tid)  { return clip_tet_list_.insert(tid).second ; }
		void remove_incident_tet(int tid) { clip_tet_list_.erase(tid) ; } 
		bool insert_clipped_seg(int sid)  { return clip_seg_list_.insert(sid).second ; }
		bool insert_clipped_facet(int i)  { return clip_facet_list_.insert(i).second ; }
		void insert_facet( const Geex::Facet& F) { mesh_.push_back(F) ; } 
		void insert_lineseg(Geex::LineSeg& L)    { linesegs_.push_back(L) ; } 
		void add_mesh(Geex::TriMesh& m)          { vol_meshes_.push_back(m) ; }
		std::set<int>& incident_facets()         { return clip_facet_list_ ; }

		vec3 point()                      { return to_geex(vh_->point()) ; }
		Vertex_handle& vh()               { return vh_ ; }
		Geex::TriMesh& mesh()             { return mesh_ ; }
		std::vector<Geex::TriMesh>& vol_meshes() { return vol_meshes_ ; }      
		void append(const RestrictedVoronoiCell& rhs) ;

		bool is_boundary() ;
		bool is_crease() ;
		bool is_corner() ;

		void get_fg(double& f, vec3& grad) ; // interface for lbfgs
		void get_fgw(double& f, vec3& grad) ;
		void get_centroid(vec3& g, double& V) ;
		void get_centroidw(vec3& g, double& V) ;
		void get_constrained_centroid(TriMesh& boundary, vec3& g, vec3& n) ;
		void get_constrained_centroidw(TriMesh& boundary, vec3& g, vec3& n) ;
		vec3 project_to_crease(vec3& pt, vec3& tang) ;
		void collect_crease_segments() ;

		// interface for volume CVT
		void get_volume_centroid(MyDelaunay* dt, vec3& g, double& V) ;
		void get_volume_centroidw(MyDelaunay* dt, vec3& g, double& V) ;

	protected:
		// centroid computation for surface
		void get_boundary_centroid(vec3& g, double& V) ;
		void get_boundary_centroidw(vec3&g, double& V) ;
		void get_smooth_centroid(vec3& g, double& V) ;
		void get_smooth_centroidw(vec3& g, double& V) ;
		vec3 get_crease_tangent() ;
		void get_crease_centroid(vec3& g, double& L) ;
		void get_crease_centroidw(vec3& g, double& L) ;
		void get_crease_centroids(std::vector<vec3>& gs) ; // used for feature splitting
		void get_perturbed_crease_centroid(vec3& g, double& L) ;
		void get_perturbed_crease_centroidw(vec3& g, double& L) ;

		// centroid computation for volume
		void get_inner_centroid(MyDelaunay* dt, vec3& g, double& V) ;
		void get_inner_centroidw(MyDelaunay* dt, vec3& g, double& V) ;
		void get_inner_fg(MyDelaunay* dt, double& f, vec3& grad) ;
		void get_boundary_fg(MyDelaunay* dt, double& f, vec3& grad, bool clip_boundary_only=false) ;
		void get_inner_fgw(MyDelaunay* dt, double& f, vec3& grad) ;
		void get_boundary_fgw(MyDelaunay* dt, double& f, vec3& grad, bool clip_boundary_only=false) ;
		
		// ODT centroid computation for volume
		void get_odt_inner_centroid(MyDelaunay* dt, vec3& g, double& V) ;
		void get_odt_inner_centroidw(MyDelaunay* dt, vec3& g, double& V) ;

	protected:
		Vertex_handle                   vh_ ;
		Geex::TriMesh                   mesh_ ;            // for boundary cell
		std::vector<Geex::TriMesh>      vol_meshes_ ;      // for boundary cell
		std::vector<Geex::LineSeg>      linesegs_ ;        // for feature cell
		std::set<int>                   clip_facet_list_ ; // clipped facets by this seed 
		std::set<int>                   clip_seg_list_ ;   // clipped feature segment index  
		std::set<int>                   clip_tet_list_ ;   //
//		simple_vector<int>              polymesh_ ;

//		friend class TopoPolyMesh ;
		friend class MyDelaunay ;
		friend class RestrictedVoronoiDiagram ;
		friend class RestrictedVoronoiDiagram3D ;
	} ;

} // end of namespace Geex
#endif // __RVC_H__
