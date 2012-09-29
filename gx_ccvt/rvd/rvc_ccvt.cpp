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

#include "rvc_ccvt.h"
#include "delaunay.h"
#include "lloyd_energy.h"
#include "geometry_ccvt.h"

namespace Geex {

	bool RestrictedVoronoiCell::is_boundary() {
		return mesh_.size() > 0 ;
	}

	bool RestrictedVoronoiCell::is_crease() {
		return linesegs_.size() > 0 && !is_corner() ;
	}

	bool RestrictedVoronoiCell::is_corner() {
		return vh_->type() == V_CORNER ;
	} 

	void RestrictedVoronoiCell::append(const RestrictedVoronoiCell& rhs) {
		gx_assert(vh_ == rhs.vh_) ;
		mesh_.insert(mesh_.end(), rhs.mesh_.begin(), rhs.mesh_.end()) ;
		vol_meshes_.insert(vol_meshes_.end(), rhs.vol_meshes_.begin(), rhs.vol_meshes_.end()) ;
		linesegs_.insert(linesegs_.end(), rhs.linesegs_.begin(), rhs.linesegs_.end()) ; 
		clip_seg_list_.insert(rhs.clip_seg_list_.begin(), rhs.clip_seg_list_.end()) ;
		clip_facet_list_.insert(rhs.clip_facet_list_.begin(), rhs.clip_facet_list_.end()) ;
		clip_tet_list_.insert(rhs.clip_tet_list_.begin(), rhs.clip_tet_list_.end()) ;
	}

	void RestrictedVoronoiCell::collect_crease_segments() {
		int nf = (int) mesh_.size() ;

		linesegs_.clear() ;
		for(int i=0; i<nf; ++i) {
			Geex::Facet& F = mesh_[i] ;
			for(int j=0; j<3; ++j) {
				if(F.edge_flag[j] == Geex::Facet::E_FEA) 
					linesegs_.push_back(Geex::LineSeg(F.vertex[(j+2)%3], F.vertex[(j+1)%3], F.vertex_weight[(j+2)%3], F.vertex_weight[(j+1)%3])) ;
			}			
		}
	}

	// ---------------------------------------------------------------------------------------
	// ---------------- functions for centroid and energy computation. -----------------------
	// ---------------------------------------------------------------------------------------

	void RestrictedVoronoiCell::get_fg(double& f, vec3& grad) {
		vec3 p0 = to_geex(vh_->point()) ;
		double V = 0.0 ;
		f = 0.0 ;		

		vec3 gradi;
		double Vi = 0;
		grad.x = 0; grad.y = 0; grad.z = 0;
		// compute energy
		for(unsigned int i=0; i<mesh_.size(); ++i) {

			f += Lloyd_energy(p0, mesh_[i].vertex[0], mesh_[i].vertex[1], mesh_[i].vertex[2], gradi, Vi);
			
			V += Vi;
			grad += gradi;

		}

		// compute gradient
		if(is_corner()) {
			grad = vec3(0, 0, 0) ;
		}
		else if(is_crease()) {
			vec3 tangent;
			project_to_crease(p0, tangent);
			grad = dot(grad, tangent)*tangent;

			//vec3 gc ;
			//double Lc ;
			//get_crease_centroid(gc, Lc) ;
			//grad = 2.0 * V*(p0 - gc) ; // gradient for features
		}
	}

	// current implementation without feature segments, do it later
	void RestrictedVoronoiCell::get_fgw(double& f, vec3& grad) {
		vec3 p0 = to_geex(vh_->point()) ;
		double V = 0.0 ;
		//bool  non = false ;
		f = 0.0 ;		
		
		vec3 gradi;
		double Vi;
		
		grad.x = 0; grad.y = 0; grad.z = 0;
		for(unsigned int i=0; i<mesh_.size(); ++i) {
			Geex::Facet& F = mesh_[i] ;

			f+= Lloyd_energy(p0, F.vertex[0], F.vertex[1], F.vertex[2], 
				F.vertex_weight[0], F.vertex_weight[1], F.vertex_weight[2], gradi, Vi) ; 
			grad += gradi;
			V += Vi;
		}
		if(is_corner()) {
			grad = vec3(0, 0, 0) ;
		} 
		else if(is_crease()) {
			//vec3 tangent;
			//project_to_crease(p0, tangent);
			//grad = dot(grad, tangent)*tangent;

			vec3 gc ;
			double Lc ;
			get_crease_centroidw(gc, Lc) ;
			grad = 2.0 * V*(p0 - gc) ; // gradient for features
		}

	}

	void RestrictedVoronoiCell::get_centroid(vec3& g, double& V) {
		if(is_corner())
			g = to_geex(vh_->point()) ;
		else if(is_crease()) 
			get_crease_centroid(g, V) ;
		else
			get_boundary_centroid(g, V) ;
	}

	void RestrictedVoronoiCell::get_centroidw(vec3& g, double& V) {
		if(is_corner())
			g = to_geex(vh_->point()) ;
		else if(is_crease()) 
			get_crease_centroidw(g, V) ;
		else
			get_boundary_centroidw(g, V) ;
	}

	void RestrictedVoronoiCell::get_boundary_centroid(vec3& g, double& V) {
		g.x = 0.0 ; g.y = 0.0 ; g.z = 0.0 ;
		V = 0.0 ;
		vec3 gi ;
		double Vi ;
		for(unsigned int i=0; i<mesh_.size(); i++) {
			Geex::Facet& F = mesh_[i] ;
			tri_centroid(F.vertex[0], F.vertex[1], F.vertex[2],	gi, Vi) ;
			g += Vi * gi ;
			V += Vi ;
		}
		if (V > 0)
		{
			g /= V;
		}
		else
		{
			g = to_geex(vh_->point());
		}
		
	}

	void RestrictedVoronoiCell::get_boundary_centroidw(vec3&g, double& V) {
		g.x = 0.0 ; g.y = 0.0 ; g.z = 0.0 ;
		V = 0.0 ;
		vec3 gi ; 
		double Vi ; 
		for(unsigned int i=0; i<mesh_.size(); i++) {
			Geex::Facet& F = mesh_[i] ;
			tri_centroid(F.vertex[0], F.vertex[1], F.vertex[2], 
				F.vertex_weight[0], F.vertex_weight[1], F.vertex_weight[2],
				gi, Vi) ;
			g += Vi * gi ;
			V += Vi ;
		}
		if (V > 0)
		{
			g /= V;
		}
		else
		{
			g = to_geex(vh_->point());
		}
	}

	void RestrictedVoronoiCell::get_crease_centroid(vec3& g, double& L) {
		g = vec3(0, 0, 0) ;
		L = 0 ;

		for(unsigned int i=0; i<linesegs_.size(); ++i) {
			g = g+linesegs_[i].length()*linesegs_[i].midpoint() ;
			L += linesegs_[i].length() ;
		}
		g /= L; //g = 1.0/L * g ;
		// perturb later outside.
		//if(dt_->perturb_feature_)
		//	g = g + 1e-6*vec3(Numeric::random_float64(), Numeric::random_float64(), Numeric::random_float64()) ; // pertubation
	}

	void RestrictedVoronoiCell::get_crease_centroidw(vec3& g, double& L) { // implement later
		g = vec3(0, 0, 0) ;
		L = 0 ;

		// use the weight of middle point, implement exact integration later.
		for(unsigned int i=0; i<linesegs_.size(); ++i) {
			g = g+0.5*(linesegs_[i].w0_+linesegs_[i].w1_)*linesegs_[i].length()*linesegs_[i].midpoint() ;
			L += 0.5*(linesegs_[i].w0_+linesegs_[i].w1_)*linesegs_[i].length() ;
		}
		g /= L; //g = 1.0/L * g ;
	}

	vec3 RestrictedVoronoiCell::get_crease_tangent() {
		vec3 T(0, 0, 0) ;
		for(unsigned int i=0; i<linesegs_.size(); ++i) {
			vec3   T = T+(linesegs_[i].v1_-linesegs_[i].v0_) ;
		}
		Geex::normalize(T) ;
		return T ;
	}

	vec3 RestrictedVoronoiCell::project_to_crease(vec3& pt, vec3& tan) {
		int idx ;
		double mindist2 = 1e10 , dist2;
		vec3   fp ; // foot point 
		vec3 p ;
		for(unsigned int i=0; i<linesegs_.size(); ++i) {
			project_to_linesegment(pt, linesegs_[i].v0_, linesegs_[i].v1_, p, dist2) ;
			if(mindist2 > dist2) {
				mindist2 = dist2 ;
				idx = i ;
				fp = p ;
			}
		}

		tan = linesegs_[idx].v1_-linesegs_[idx].v0_ ;
		Geex::normalize(tan) ;
		//g = fp ;
		return fp ;
	}

	//void RestrictedVoronoiCell::get_crease_centroids(std::vector<vec3>& gs) {
	//	std::map<int, std::vector<LineSeg> > lines ;
	//	typedef std::map<int, std::vector<LineSeg> >::iterator LineIter ;
	//	typedef std::pair<int, std::vector<LineSeg> > Fea_Pair ;

	//	for(unsigned int i=0; i<linesegs_.size(); ++i) {
	//		if(lines.find(linesegs_[i].parentidx())==lines.end()) {
	//			std::vector<LineSeg> line ;
	//			line.push_back(linesegs_[i]) ;
	//			lines.insert(Fea_Pair(linesegs_[i].parentidx(), line)) ;
	//		} else {
	//			lines[linesegs_[i].parentidx()].push_back(linesegs_[i]) ;
	//		}
	//	}

	//	gs.clear() ;
	//	for(LineIter lit=lines.begin(); lit!=lines.end(); ++lit) {
	//		std::vector<LineSeg>& line = lit->second ;
	//		vec3 g = vec3(0, 0, 0) ;
	//		double L = 0 ;

	//		for(unsigned int i=0; i<line.size(); ++i) {
	//			g = g+line[i].length()*line[i].midpoint() ;
	//			L += line[i].length() ;
	//		}
	//		g = 1.0/L * g ;

	//		gs.push_back(g) ;
	//	}
	//}

	void RestrictedVoronoiCell::get_constrained_centroid(TriMesh& boundary, vec3& g, vec3& n) {
		double V ;
		get_centroid(g, V) ;
		if(is_corner()) //do nothing
			{}
		else if(is_crease()) 
			project_to_crease(g, n) ;
		else
			boundary.project_to_mesh(g, n) ;
	}

	void RestrictedVoronoiCell::get_constrained_centroidw(TriMesh& boundary, vec3& g, vec3& n) {
		double V ;
		get_centroidw(g, V) ;
		if(is_corner()) //do nothing
			{}
		else if(is_crease()) 
			project_to_crease(g, n) ;
		else
			boundary.project_to_mesh(g, n) ;
	}

	// ------------------------------------------------------------------------------------
	// ---------------------------- Interface of Volume CVT -------------------------------
	// ------------------------------------------------------------------------------------
#ifndef F_BDY 
#define F_BDY -2  // boundary face
#endif 
#ifndef F_AUX 
#define F_AUX -1  // aux face
#endif 

	void RestrictedVoronoiCell::get_volume_centroid(MyDelaunay* dt, vec3& g, double& V) {
		if(is_boundary()) {
			vec3 p0 = to_geex(vh_->point()) ;
			vec3 Vg(0.0, 0.0, 0.0) ;
			V = 0.0 ;

//			gx_assert(vol_meshes_.size()>0) ;

			for(unsigned int i=0; i<vol_meshes_.size(); ++i) {
				Geex::TriMesh& mi = vol_meshes_[i] ;
				for(unsigned int j=0; j<mi.size(); ++j) {
					Geex::Facet& Fij = mi[j] ;
					if(Fij.ftype_!=F_AUX)  {
						double Vij = tetra_signed_volume(p0, Fij.vertex[0], Fij.vertex[1], Fij.vertex[2]) ;
						vec3 gij = tetra_centroid(p0, Fij.vertex[0], Fij.vertex[1], Fij.vertex[2]) ;
						V += Vij ;
						Vg += Vij * gij ;
					}
				}
			}
			g = 1./V * Vg ;
		} else {
			get_inner_centroid(dt, g, V) ;
		}
	}

	void RestrictedVoronoiCell::get_inner_centroid(MyDelaunay* dt, vec3& g, double& V) {
		std::vector<MyDelaunay::Edge>& edges = dt->stars()[vh_->id] ;
		vec3 p0 = to_geex(vh_->point()) ;
		g.x = 0.0 ; g.y = 0.0 ; g.z = 0.0 ;
        V = 0.0 ;
        for(unsigned int i=0; i<edges.size(); i++) {
            MyDelaunay::Edge e = edges[i] ;
            vec3 pts[2] ;
            int index = 0 ;
            MyDelaunay::Cell_circulator ccir = dt->incident_cells(e) ;
            MyDelaunay::Cell_circulator cdone = ccir ;
            do {
                vec3 p = ccir->dual ;
                if(index<2) {
                    pts[index] = p ; 
                    index++ ;
                } else {
                    double Vi = fabs(tetra_signed_volume(p0, pts[0], pts[1], p)) ;
                    vec3 gi = tetra_centroid(p0, pts[0], pts[1], p) ;
                    V += Vi ;
                    g += Vi * gi ;
                    pts[1] = p ;
                }
                ccir++ ;
            } while(ccir != cdone) ;
        }
        double s = (1.0 / V) ;
        g.x = s*g.x ;
        g.y = s*g.y ;
        g.z = s*g.z ;	
	}

	void RestrictedVoronoiCell::get_volume_centroidw(MyDelaunay* dt, vec3& g, double& V) {
		if(is_boundary()) {
			vec3 p0 = to_geex(vh_->point()) ;
			vec3 Vg(0.0, 0.0, 0.0) ;
			V = 0.0 ;

//			gx_assert(vol_meshes_.size()>0) ;

			for(unsigned int i=0; i<vol_meshes_.size(); ++i) {
				Geex::TriMesh& mi = vol_meshes_[i] ;
				for(unsigned int j=0; j<mi.size(); ++j) {
					Geex::Facet& Fij = mi[j] ;
					if(Fij.ftype_!=F_AUX)  {
						double w = vh_->weight() ; // test
						double Vij = tetra_signed_volume(p0, Fij.vertex[0], Fij.vertex[1], Fij.vertex[2])*w ;
						vec3 gij = tetra_centroid(p0, Fij.vertex[0], Fij.vertex[1], Fij.vertex[2]) ;
						V += Vij ;
						Vg += Vij * gij ;
					}
				}
			}
			g = 1./V * Vg ;
		} else {
			get_inner_centroidw(dt, g, V) ;
		}
	}

	void RestrictedVoronoiCell::get_inner_centroidw(MyDelaunay* dt, vec3& g, double& V) {
		std::vector<MyDelaunay::Edge>& edges = dt->stars()[vh_->id] ;
		vec3 p0 = to_geex(vh_->point()) ;
		real w0 = vh_->weight() ;
		g.x = 0.0 ; g.y = 0.0 ; g.z = 0.0 ;
        V = 0.0 ;
        for(unsigned int i=0; i<edges.size(); i++) {
            MyDelaunay::Edge e = edges[i] ;
            vec3 pts[2] ;
			real ws[2] ;
            int index = 0 ;
            MyDelaunay::Cell_circulator ccir = dt->incident_cells(e) ;
            MyDelaunay::Cell_circulator cdone = ccir ;
            do {
                vec3 p = ccir->dual ;
				real w = ccir->weight() ;
                if(index<2) {
                    pts[index] = p ; 
					ws[index] = w ;
                    index++ ;
                } else {
                    double Vi = fabs(tetra_signed_volume(p0, pts[0], pts[1], p))*(w0+w+ws[0]+ws[1])*0.25 ;
                    vec3 gi = tetra_centroid(p0, pts[0], pts[1], p) ;
                    V += Vi ;
                    g += Vi * gi ;
                    pts[1] = p ;
					ws[1] = w ;
                }
                ccir++ ;
            } while(ccir != cdone) ;
        }
        double s = (1.0 / V) ;
        g.x = s*g.x ;
        g.y = s*g.y ;
        g.z = s*g.z ;	
	}

	void RestrictedVoronoiCell::get_boundary_fg(MyDelaunay* dt, double& f, vec3& grad, bool clip_boundary_only) {
	    vec3 p0 = to_geex(vh_->point()) ;
        double V = 0.0 ;
        f = 0.0 ;
        grad.x = grad.y = grad.z = 0;
		for(unsigned int i=0; i<vol_meshes_.size(); ++i) {
			Geex::TriMesh& mi = vol_meshes_[i] ;
			for(unsigned int j=0; j<mi.size(); ++j) {
				Geex::Facet& Fij = mi[j] ;
				double Vij = 0 ;
				vec3   gij(0,0,0) ;
				if(Fij.ftype_!=F_AUX)  {
                    f += Lloyd_energy_tet(p0, p0, Fij.vertex[0], Fij.vertex[1], Fij.vertex[2], gij, Vij);
                    grad += gij;
					V += Vij ;
				}
			}
		}
		if(clip_boundary_only) {
			vec3 g ;
			double A ;
			this->get_boundary_centroid(g, A) ;
            //g -= p0;
            //normalize(g);
            //grad = dot(grad, g) * g;
			grad = 2.0 * V * (p0 - g) ;
		}
	}

	void RestrictedVoronoiCell::get_inner_fg(MyDelaunay* dt, double& f, vec3& grad) {
		std::vector<MyDelaunay::Edge>& edges = dt->stars()[vh_->id] ;
        vec3 p0 = to_geex(vh_->point()) ;
        vec3 g(0.0, 0.0, 0.0) ;
		f = 0.0 ;
        double V;
        grad.x = grad.y = grad.z = 0;
        double fi;
        for(unsigned int i=0; i<edges.size(); i++) {
            MyDelaunay::Edge e = edges[i] ;
            int index = 0 ;
            MyDelaunay::Cell_circulator ccir = dt->incident_cells(e) ;
            MyDelaunay::Cell_circulator cdone = ccir ;
            
            vec3 p = ccir->dual ;
            ccir++;
            MyDelaunay::Cell_circulator ccir2 = ccir;
            ccir2++;
            
            do {
                fi = Lloyd_energy_tet(p0, p0, ccir->dual, ccir2->dual, p, g, V) ;
                if (fi < 0)
                {
                    f -= fi;
                    grad -= g;
                }
                else
                {
                    f += fi;
                    grad += g;
                }
                ccir++;
                ccir2++;
            } while(ccir2 != cdone) ;
        }

	}

	void RestrictedVoronoiCell::get_boundary_fgw(MyDelaunay* dt, double& f, vec3& grad, bool clip_boundary_only) {
	    vec3 p0 = to_geex(vh_->point()) ;
        vec3 g(0.0, 0.0, 0.0) ;
        double V = 0.0 ;
        f = 0.0 ;
        grad.x = grad.y = grad.z = 0;
		for(unsigned int i=0; i<vol_meshes_.size(); ++i) {
			Geex::TriMesh& mi = vol_meshes_[i] ;
			for(unsigned int j=0; j<mi.size(); ++j) {
				Geex::Facet& Fij = mi[j] ;
				if(Fij.ftype_!=F_AUX)  {
                    f += Lloyd_energy_tet(p0, p0, Fij.vertex[0], Fij.vertex[1], Fij.vertex[2], g, V)*vh_->weight();
                    grad += vh_->weight()*g;
				}
			}
		}

		if(clip_boundary_only) {
			vec3 g ;
			double A ;
			this->get_boundary_centroidw(g, A) ;
            //g -= p0;
            //normalize(g);
            //grad = dot(grad, g) * g;
			grad = 2.0 * vh_->weight() * V * (p0 - g) ;
		}
	}

	void RestrictedVoronoiCell::get_inner_fgw(MyDelaunay* dt, double& f, vec3& grad) {
		std::vector<MyDelaunay::Edge>& edges = dt->stars()[vh_->id] ;
        vec3 p0 = to_geex(vh_->point()) ;
		real w0 = vh_->weight() ;
        vec3 g(0.0, 0.0, 0.0) ;
		f = 0.0 ;
        double V;
        grad.x = grad.y = grad.z = 0;
        double fi;
        for(unsigned int i=0; i<edges.size(); i++) {
            MyDelaunay::Edge e = edges[i] ;
            int index = 0 ;
            MyDelaunay::Cell_circulator ccir = dt->incident_cells(e) ;
            MyDelaunay::Cell_circulator cdone = ccir ;
            
            vec3 p = ccir->dual ;
			real w = ccir->weight() ;
            ccir++;
            MyDelaunay::Cell_circulator ccir2 = ccir;
            ccir2++;
            
            do {
				real w1= ccir->weight() ;
				real w2= ccir2->weight() ;
                fi = Lloyd_energy_tet(p0, p0, ccir->dual, ccir2->dual, p, g, V)*(w0+w1+w2+w)*0.25 ;
                if (fi >= 0)
                {
                    f += fi;
                    grad += (w0+w1+w2+w)*0.25*g;
                }
                else
                {
                    f -= fi;
                    grad -= (w0+w1+w2+w)*0.25*g;
                }
                
                ccir++;
                ccir2++;
            } while(ccir2 != cdone) ;
        }
	}

	// ODT Inter face
	void RestrictedVoronoiCell::get_odt_inner_centroid(MyDelaunay* dt, vec3& g, double& V) {
		std::vector<MyDelaunay::Edge>& edges = dt->stars()[vh_->id] ;
		vec3 p0 = to_geex(vh_->point()) ;
		g.x = 0.0 ; g.y = 0.0 ; g.z = 0.0 ;
        V = 0.0 ;
        for(unsigned int i=0; i<edges.size(); i++) {
            MyDelaunay::Edge e = edges[i] ;
            vec3 pts[2] ;
            int index = 0 ;
            MyDelaunay::Cell_circulator ccir = dt->incident_cells(e) ;
            MyDelaunay::Cell_circulator cdone = ccir ;
            do {
                vec3 p = ccir->dual ;
				double Vi = fabs(tetra_signed_volume(to_geex(ccir->vertex(0)->point()), 
					                                 to_geex(ccir->vertex(1)->point()), 
													 to_geex(ccir->vertex(2)->point()),
													 to_geex(ccir->vertex(2)->point()))) ;
                V += Vi ;
                g += Vi * p ;
                ccir++ ;
            } while(ccir != cdone) ;
        }
        double s = (1.0 / V) ;
        g.x = s*g.x ;
        g.y = s*g.y ;
        g.z = s*g.z ;	
	}

	void RestrictedVoronoiCell::get_odt_inner_centroidw(MyDelaunay* dt, vec3& g, double& V) {
		std::vector<MyDelaunay::Edge>& edges = dt->stars()[vh_->id] ;
		vec3 p0 = to_geex(vh_->point()) ;
		real w0 = vh_->weight() ;
		g.x = 0.0 ; g.y = 0.0 ; g.z = 0.0 ;
        V = 0.0 ;
        for(unsigned int i=0; i<edges.size(); i++) {
            MyDelaunay::Edge e = edges[i] ;
            vec3 pts[2] ;
//			real ws[2] ;
            int index = 0 ;
            MyDelaunay::Cell_circulator ccir = dt->incident_cells(e) ;
            MyDelaunay::Cell_circulator cdone = ccir ;
            do {
                vec3 p = ccir->dual ;
				real w = ccir->weight() ;
				double Vi = fabs(tetra_signed_volume(to_geex(ccir->vertex(0)->point()), 
					                                 to_geex(ccir->vertex(1)->point()), 
													 to_geex(ccir->vertex(2)->point()),
													 to_geex(ccir->vertex(2)->point()))) * w ;
				V += Vi ;
                g += Vi * p ;
                ccir++ ;
            } while(ccir != cdone) ;
        }
        double s = (1.0 / V) ;
        g.x = s*g.x ;
        g.y = s*g.y ;
        g.z = s*g.z ;	
	}

} // end of namespace Geex

