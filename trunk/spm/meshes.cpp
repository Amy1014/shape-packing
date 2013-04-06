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

#include "meshes.h"
#include "geometry_ccvt.h"
#include <Geex/basics/line_stream.h>
#include <fstream>
#include <algorithm>
#include <queue>
#include <assert.h>

namespace Geex {

	double TriMesh::PI = 3.14159265358979323846;

    TriMesh::TriMesh() {
        nb_vertices_ = 0 ;
		samplerate_  = 36 ; 
		kdtree_      = nil ;
		featureAngleCriterion = 0.0;
		highCurvatureCriterion = 0.2;
		maxFacetWeight = minFacetWeight = 0.0;

		//min_vert_curvature = max_vert_curvature = 0.0;
    }
	TriMesh::~TriMesh() {
		if(kdtree_) delete kdtree_ ;
	}

    void TriMesh::load(const std::string& filename, Convex* planes) {
        std::cerr << "Loading " << filename << std::endl ;
        clear() ;
		nb_vertices_ = 0;  // dmyan
		vertices_.clear(); // dmyan
        std::ifstream in(filename.c_str()) ;
        if(!in) {
            std::cerr << "Could not open " << filename << std::endl ;
            return ;
        }
        load(in, planes) ;
    }

    void TriMesh::load(std::istream& input, Convex* planes) {

        if(planes != nil) { planes->clear() ; }

        LineInputStream in(input) ;
        std::vector<vec3> vertices ;
//		vertices_.clear();
		std::map<pair<int, int>, pair<int, int>, PairIntIntCmp> edgeValence;//to compute boundary and feature edges, by Wenchao
        while(!in.eof()) {
            in.get_line() ;
            std::string keyword ;
            
            in >> keyword ;
            
            if(keyword == "v") {
                vec3 p ;
                in >> p ;
                vertices.push_back(p) ;
				vertices_.push_back(MeshVertex(p)) ;
				nb_vertices_++ ;				
            } else if(keyword == "f") {
                std::vector<int> cur_facet ;
                while(!in.eol()) {
                    std::string s ;
                    in >> s ;
                    if(s.length() > 0) {
                        std::istringstream v_input(s) ;
                        int index ;
                        v_input >> index ;
                        if(index < 1 || index > (int)vertices_.size()) {
                            std::cerr << "Out of bounds vertex index" << std::endl ;
                        } else {
                            cur_facet.push_back(index - 1) ;
                        }
                        char c ;
                        v_input >> c ;
                        // tex vertex, ignored
                        if(c == '/') {
                            v_input >> index ;
                        }
                    }
                }
                if(cur_facet.size() < 3) {
                    std::cerr << "facet with less than 3 vertices, ignoring" << std::endl ;
                } else {
                    for(unsigned int i=1; i<cur_facet.size() - 1; i++) {
                        int id0 = cur_facet[0] ;
                        int id1 = cur_facet[i] ;
                        int id2 = cur_facet[i+1] ;
                        //const vec3& p0 = vertices_[ id0 ].pos_ ;
                        //const vec3& p1 = vertices_[ id1 ].pos_ ;
                        //const vec3& p2 = vertices_[ id2 ].pos_ ;
                        const vec3& p0 = vertices[ id0 ] ;
                        const vec3& p1 = vertices[ id1 ] ;
                        const vec3& p2 = vertices[ id2 ] ;
						// compute the number of adjacent faces for each edge, by Wenchao
						using std::make_pair;
						std::map<pair<int, int>, pair<int, int>, PairIntIntCmp>::iterator findIter0, findIter1;
						findIter0 = edgeValence.find(make_pair(id0, id1));
						findIter1 = edgeValence.find(make_pair(id1, id0));
						if (findIter0 == edgeValence.end() && findIter1 == edgeValence.end())
						{
							edgeValence[make_pair(id0, id1)] = make_pair(int(size()), -1);
							//edgeValence[make_pair(id0, id1)]++;
						}
						else if ( findIter0 == edgeValence.end() && findIter1 != edgeValence.end())
						{
							//findIter1->second++;
							findIter1->second.second = int(size());
						}
						else //if ( findIter0 != edgeValence.end() && findIter1 == edgeValence.end() )
						{
							//findIter0->second++;
							findIter0->second.second = int(size());
						}
						findIter0 = edgeValence.find(make_pair(id0, id2));
						findIter1 = edgeValence.find(make_pair(id2, id0));
						if (findIter0 == edgeValence.end() && findIter1 == edgeValence.end())
						{
							//edgeValence[make_pair(id0, id2)]++;
							edgeValence[make_pair(id0, id2)] = make_pair(int(size()), -1);
						}
						else if ( findIter0 == edgeValence.end() && findIter1 != edgeValence.end())
						{
							//findIter1->second++;
							findIter1->second.second = int(size());
						}
						else //if ( findIter0 != edgeValence.end() && findIter1 == edgeValence.end() )
						{
							//findIter0->second++;
							findIter0->second.second = int(size());
						}
						findIter0 = edgeValence.find(make_pair(id1, id2));
						findIter1 = edgeValence.find(make_pair(id2, id1));
						if (findIter0 == edgeValence.end() && findIter1 == edgeValence.end())
						{
							//edgeValence[make_pair(id1, id2)]++;
							edgeValence[make_pair(id1, id2)] = make_pair(int(size()), -1);
						}
						else if ( findIter0 == edgeValence.end() && findIter1 != edgeValence.end())
						{
							//findIter1->second++;
							findIter1->second.second = int(size());
						}
						else //if ( findIter0 != edgeValence.end() && findIter1 == edgeValence.end() )
						{
							//findIter0->second++;
							findIter0->second.second = int(size());
						}

//						std::cout <<"M:" <<pt0 << " " << pt1 << " " << pt2 << "\nR:" << p0 <<" " <<p1 << " " << p2 << std::endl;
                        bool e0 = true ;
                        bool e1 = (i == cur_facet.size() - 2) ;
                        bool e2 = (i == 1) ;
                        push_back(Facet(p0, p1, p2, e0, e1, e2, id0, id1, id2)) ;
						back().norm = Geex::normalize(cross(p1-p0, p2-p0));
                        if(i == 1 && planes != nil) {
                            planes->push_back(rbegin()->plane()) ;
                        }

						vertices_[ id0 ].faces_.push_back((int)size()-1);
						vertices_[ id1 ].faces_.push_back((int)size()-1);
						vertices_[ id2 ].faces_.push_back((int)size()-1);
                    }
                } 
            }
        }
		// collect the boundary edges and compute bihedral angles
		//for ( std::map<std::pair<int, int>, int, PairIntIntCmp>::const_iterator it = edgeValence.begin(); it!=edgeValence.end(); it++)
		//	if (it->second == 1)
		//		boundaryEdges.insert(it->first);
		vert_on_boundary.resize(nb_vertices(), false);
		vert_on_feature.resize(nb_vertices(), false);
		for (std::map<pair<int, int>, pair<int, int>, PairIntIntCmp>::const_iterator it=edgeValence.begin(); it!=edgeValence.end(); it++)
		{
			if (it->second.second < 0)
			{
				boundaryEdges.insert(it->first);
				vert_on_boundary[it->first.second] = true;
				vert_on_boundary[it->first.first] = true;
			}
			else
			{
				int fi0 = it->second.first, fi1 = it->second.second;
				const Facet& f0 = (*this)[fi0];
				const Facet& f1 = (*this)[fi1];
				// common edge between the two triangular facets
				int vi0 = it->first.first, vi1 = it->first.second;
				const vec3& v0 = vertices_[vi0].point();
				const vec3& v1 = vertices_[vi1].point();
				// the different vertices
				int pi0, pi1;
				if ( f0.vertex_index[0] != vi0 && f0.vertex_index[0] != vi1)		pi0 = f0.vertex_index[0];
				else if ( f0.vertex_index[1] != vi0 && f0.vertex_index[1] != vi1)	pi0 = f0.vertex_index[1];
				else																pi0 = f0.vertex_index[2];
				if (f1.vertex_index[0] != vi0 && f1.vertex_index[0] != vi1)			pi1 = f1.vertex_index[0];
				else if (f1.vertex_index[1] != vi0 && f1.vertex_index[1] != vi1)	pi1 = f1.vertex_index[1];
				else																pi1 = f1.vertex_index[2];
				const vec3& p0 = vertices_[pi0].point();
				const vec3& p1 = vertices_[pi1].point();
				// this more complicate formula is used to avoid error introduced by subtracting two small vectors
				vec3 n0 = cross(p0, v1) - (cross(v0, v1) + cross(p0, v0));
				n0 /= n0.length();
				vec3 n1 = cross(v1, p1) - (cross(v0, p1) + cross(v1, v0));
				n1 /= n1.length();
				//bihedralAngles[make_pair(fi0, fi1)] = -dot(n0, n1);
				bihedralAngles.push_back(BihedralAngle(fi0, fi1, vi0, vi1, -dot(n0, n1)));
			}
		}
		setFeatureAngle(featureAngleCriterion);
		set_facet_neighbor() ; // set three neighbors of each face

		std::cout << "nb_vertices: " <<nb_vertices_ << "; nb_faces: " <<size()<< std::endl;
		std::cout << "number of boundary vertices: "<<boundaryEdges.size()<<std::endl;
    }

	void TriMesh::setFeatureAngle(double ang)
	{
		if (ang < 0.0 || ang > 180.0)
			return;	
		vert_on_feature.assign(nb_vertices(), false);
		featureFacets.clear();
		featureEdges.clear();
		featureAngleCriterion = ang;
		double cosCriterion = std::cos(featureAngleCriterion*PI/180.0);
		using std::make_pair;
		for (unsigned int i = 0; i < bihedralAngles.size(); i++)
		{
			if ( bihedralAngles[i].cosAngle >= cosCriterion )
			{
				featureFacets.push_back(make_pair(bihedralAngles[i].findex0, bihedralAngles[i].findex1));
				featureEdges.push_back(make_pair(bihedralAngles[i].vindex0, bihedralAngles[i].vindex1));
				vert_on_feature[bihedralAngles[i].vindex0] = vert_on_feature[bihedralAngles[i].vindex1] = true;
			}
		}
	}
	void TriMesh::flip_normals()
	{
		for (unsigned int i = 0; i < size(); i++)
			operator[](i).norm = -operator[](i).norm;
	}
	void TriMesh::setHighCurvaturePercent(double percent /* = 0.2 */)
	{
		facetHighCurvature.clear();
		if (percent <= 0.0 || percent > 1.0 || facetWeight.size() == 0)
		{
			highCurvatureCriterion = 0.0;
			maxFacetWeight = minFacetWeight = 0.0;
			return;
		}
		int n(percent*size());
		for (int i = size()-1; i > size()-n-1; i--)
			facetHighCurvature.insert(facetWeight[i].second);
		highCurvatureCriterion = percent;
		maxFacetWeight = facetWeight[size()-1].first;
		minFacetWeight = facetWeight[size()-n].first;

	}
	void TriMesh::set_facet_neighbor() {
		for(int i=0; i<(int)size(); ++i) {
			Geex::Facet& F = (*this)[i] ;
			for(int j=0; j<3; ++j) {
				int vidx = F.vertex_index[j] ;
				std::vector<int>& vfacets = vertex(vidx).faces_ ;
				for(int k=0; k<(int)vfacets.size(); ++k) {
					if((*this)[vfacets[k]].has_vertex(F.vertex_index[(j+1)%3]) && vfacets[k]!=i) {
						F.face_index[(j+2)%3] = vfacets[k] ;
						break ;
					}
				}
			}
		}
	}

	void TriMesh::set_feature_edgeflag() {
		for(int i=0; i<(int)fealines_.size(); ++i) {
			std::vector<int>& line = fealines_[i] ;
			for(int j=0; j<(int)line.size()-1; ++j) {
				int vidx0 = line[j] ;
				int vidx1 = line[j+1] ;
				MeshVertex& v0 = vertex(vidx0) ;
				for(int k=0; k<(int)v0.faces_.size(); ++k) {
					Facet& F = (*this)[v0.faces_[k]];
					if(F.has_vertex(vidx1)) {
						int eidx = F.find_edge(vidx0, vidx1) ;
						if(eidx>=0) F.edge_flag[eidx] = Facet::E_FEA ; 
					}
				}
			}
		}
	}

	bool TriMesh::load_feature(const std::string& filename) {
		std::ifstream in(filename.c_str()) ;
		int nv, /*ne,*/ nl ;

		if( !in.is_open() )
			return false ;
		
		corners_.clear() ;
		fealines_.clear() ;
		features_.clear() ;

		in >> nv >> nl ;
		
		// load constrained vertex indices
		for(int i=0; i<nv; ++i) {
			int vidx ;
			in >> vidx ;
			corners_.push_back(vidx) ;
		}

		// load feature lines
		for(int i=0; i<nl; ++i) {
			int ni ;
			std::vector<int> line ;
			in >> ni ;
			for(int j=0; j<ni; ++j) {
				int vidx ;
				in >> vidx ;
				line.push_back(vidx) ;
			}
			fealines_.push_back(line) ;
		}
		
		in.close();

		// create feature edges
		for(int i=0; i<(int)fealines_.size(); ++i) {
			for(int j=0; j<(int)fealines_[i].size()-1; ++j) {
				int p0 = fealines_[i][j] ;
				int p1 = fealines_[i][j+1] ;
				features_.push_back(LineSeg(vertices_[p0].point(), vertices_[p1].point(), 
					                 vertices_[p0].weight(), vertices_[p1].weight(), -1, -1, i)) ;
			}
		}

		// mark feature edges for each facet
		set_feature_edgeflag() ;

		return true ;
	}

	void TriMesh::set_uniform_density() {
		for(unsigned int i=0; i<vertices_.size(); ++i) 
			vertices_[i].weight() = 1.f ;
		for(unsigned int i=0; i<size(); ++i) 
			for(int j=0; j<3; ++j) 
			(*this)[i].vertex_weight[j] = 1.f ;
	}

	bool TriMesh::load_density(const std::string& filename, double gamma)
	{
		std::ifstream in(filename.c_str()) ;
		if(!in.is_open()) {
			std::cerr << "could not open " << filename << std::endl ;
			return false ;
		}

		char c;
		in >> c;
		while (!std::isdigit(c) && in)
			in >> c;

		if (!in)
		{
			std::cerr<<"Invalid curvature file: no data input.\n";
			return false;
		}
		
		// begin read curvature
		std::vector<double> min_cur, max_cur;
		min_cur.reserve(nb_vertices());
		max_cur.reserve(nb_vertices());

		while (in)
		{
			double val;
			// min curvature
			in >> val;
			min_cur.push_back(val);
			// ignore direction
			in >> val; in >> val; in >> val;
			// max curvature
			in >> val;
			max_cur.push_back(val);
			// ignore direction
			in >> val; in >> val; in >> val;
			// go to the next vertex
			int nxt_vert_idx;
			in >> nxt_vert_idx;
		}
		if (min_cur.size() != nb_vertices() || max_cur.size() != nb_vertices())
		{
			std::cerr<<"Caution: number of vertices in curvature and mesh file does not match. Curvature will not be used!\n";
			return false;
		}

		// take the larger curvature value in terms of absolute value
		std::vector<double> curs(nb_vertices());
		for (unsigned int i = 0; i < nb_vertices(); i++)
			curs[i] = std::max(std::fabs(min_cur[i]), std::fabs(max_cur[i]));
		std::vector<double> avgcurs(nb_vertices(), 0.0);
		std::vector<int> nb_neighbors(nb_vertices(), 0);
		for (unsigned int i = 0; i < size(); i++)
		{
			const Facet& f = operator[](i);
			int vi0 = f.vertex_index[0], vi1 = f.vertex_index[1], vi2 = f.vertex_index[2];
			avgcurs[vi0] += (curs[vi1]+curs[vi2]);
			nb_neighbors[vi0] += 2;
			avgcurs[vi1] += (curs[vi0]+curs[vi2]);
			nb_neighbors[vi1] += 2;
			avgcurs[vi2] += (curs[vi0]+curs[vi1]);
			nb_neighbors[vi2] += 2;
		}
		// for the boundary, sum once more
		for (BoundaryEdgeSet::const_iterator it = boundaryEdges.begin(); it != boundaryEdges.end(); it++)
		{
			int vi0 = it->first, vi1 = it->second;
			avgcurs[vi0] += curs[vi1];
			nb_neighbors[vi0]++;
			avgcurs[vi1] += curs[vi0];
			nb_neighbors[vi1]++;
		}
		//std::ofstream test_file("test.txt");
		for (unsigned int i = 0; i < nb_vertices(); i++)
		{
			avgcurs[i] /= nb_neighbors[i];
			vertices_[i].curvature = avgcurs[i]; // load curvature at this vertex
			//if (min_vert_curvature > avgcurs[i])
			//	min_vert_curvature = avgcurs[i];
			//if (max_vert_curvature < avgcurs[i])
			//	max_vert_curvature = avgcurs[i];
			//test_file<<vertices_[i].curvature<<std::endl;
		}
		for (unsigned int i = 0; i < nb_vertices(); i++)
			vertices_[i].weight() = pow(avgcurs[i]*2, gamma);

		// compute face's weight by averaging weight of vertices
		for(unsigned int i=0; i < size(); ++i)
		{
			Facet& F = (*this)[i] ;
			for(int j=0; j<3; ++j) 
				F.vertex_weight[j] = vertices_[F.vertex_index[j]].weight() ;
		}
		// compute facet weight and sort
		facetWeight.resize(size());
		for (int i = 0; i < size(); i++)
		{
			const Facet& f = operator[](i);
			facetWeight[i].first = (f.vertex_weight[0] + f.vertex_weight[1] + f.vertex_weight[2])/3.0;
			facetWeight[i].second = i;
		}
		std::sort(facetWeight.begin(), facetWeight.end(), PairDoubleIntCmp());
		setHighCurvaturePercent(highCurvatureCriterion);
		return true ;
	}

#if 0
    bool TriMesh::load_density(const std::string& filename, double gamma) {
		std::ifstream in(filename.c_str()) ;
		if(!in.is_open()) {
			std::cerr << "could not open " << filename << std::endl ;
			return false ;
		}

		// each vertex has a density value
		int nv ;
		in >> nv ;
		if(nv!=nb_vertices()) {
			std::cerr << "number of vertices doesn't match. " << std::endl ;
			return false ;
		}

		std::cout << "loading density for boundary mesh..." << std::endl ;

		double wmax=0, wmin=1e20 ;
		//smooth the curvature of each vertex by averaging its neighbors
		vector<double> curs(nv);//curvatures of vertices, by Wenchao Hu
		for (unsigned int i = 0; i < curs.size(); i++)
			in>>curs[i];
		vector<double> avgcurs(nv, 0.0);
		vector<int> nbNeighbors(nv, 0);///number of neighbors
		for (unsigned int i = 0; i < this->size(); i++)
		{
			// compute curvature sum of neighbors
			const Facet& f = operator[](i);
			int v0 = f.vertex_index[0], v1 = f.vertex_index[1], v2 = f.vertex_index[2];
			avgcurs[v0] += (curs[v1]+curs[v2]);
			nbNeighbors[v0] += 2;
			avgcurs[v1] += (curs[v0]+curs[v2]);
			nbNeighbors[v1] += 2;
			avgcurs[v2] += (curs[v0]+curs[v1]);
			nbNeighbors[v2] += 2;
		}
		// for the boundary, sum once more
		for (BoundaryEdgeSet::const_iterator it = boundaryEdges.begin(); it != boundaryEdges.end(); it++)
		{
			int v0 = it->first, v1 = it->second;
			avgcurs[v0] += curs[v1];
			nbNeighbors[v0]++;
			avgcurs[v1] += curs[v0];
			nbNeighbors[v1]++;
		}
		for (unsigned int i = 0; i < avgcurs.size(); i++)
		{
			avgcurs[i] /= nbNeighbors[i];
			vertices_[i].curvature = avgcurs[i]; // load curvature at this vertex
		}
//		for(int i=0; i<nv; ++i) {
//		    double w ; 
//		    in >> w ;
////		    if(wmin > w) wmin = w;
////		    if(wmax < w) wmax = w ;
////		    vertices_[i].weight() = w; 
//		    vertices_[i].weight() = pow(w*10, gamma) ;
//		}
		for (unsigned int i = 0; i< avgcurs.size(); i++)
			vertices_[i].weight() = pow(avgcurs[i]*5, gamma);

		// compute face's weight by averaging weight of vertices
		for(unsigned int i=0; i<size(); ++i) {
			Geex::Facet& F = (*this)[i] ;
//			double w = 0 ;
			for(int j=0; j<3; ++j) 
				F.vertex_weight[j] = vertices_[F.vertex_index[j]].weight() ;
//			w /= 3.0 ;

//			F.weight() = w ;
		}
		// compute facet weight and sort
		facetWeight.resize(size());
		for (int i = 0; i < size(); i++)
		{
			const Facet& f = operator[](i);
			//double facetMass = tri_mass(f.vertex[0], f.vertex[1], f.vertex[2], 
			//							f.vertex_weight[0], f.vertex_weight[1], f.vertex_weight[2]);
			facetWeight[i].first = (f.vertex_weight[0] + f.vertex_weight[1] + f.vertex_weight[2])/3.0;
			facetWeight[i].second = i;
		}
		std::sort(facetWeight.begin(), facetWeight.end(), PairDoubleIntCmp());
		setHighCurvaturePercent(highCurvatureCriterion);
		return true ;
	}
#endif

    void TriMesh::save(const std::string& filename) {
        std::ofstream out(filename.c_str()) ;
        if(!out) {
            std::cerr << "could not open " << filename << std::endl ;
            return ;
        }
        save(out) ;
    }
    
    void TriMesh::save(std::ostream& out) {
        int cur_v = 1 ;
        for(unsigned int i=0; i<size(); i++) {
            const Facet& F = (*this)[i] ;
            for(unsigned int j=0; j<3; j++) {
                out << "v " << F.vertex[j] << std::endl ;
            }
            out << "f " << cur_v << " " << cur_v+1 << " " << cur_v+2 << std::endl ;
            cur_v += 3 ;
        }
    }

	void TriMesh::save_obj(const std::string& filename) {
		std::cout << "saving mesh: " << filename << std::endl ;
		std::ofstream ofile(filename.c_str()) ;
		for(unsigned int i=0; i<vertices_.size(); ++i) 
			ofile <<"v " << vertices_[i].pos_ << std::endl ;
		for(unsigned int i=0; i<size(); ++i) 
			ofile <<"f " << (*this)[i].vertex_index[0]+1 << " " << (*this)[i].vertex_index[1]+1 << " " << (*this)[i].vertex_index[2]+1 << std::endl ;
		ofile.close() ;
		std::cout << "end save mesh" << std::endl ;
	}
	void TriMesh::build_kdtree() {
		if(kdtree_ != nil) delete kdtree_;
		int curridx = 0;
		double * data = new double[(int)size()*(samplerate_+1)*(samplerate_+2)*3/2];
		for(int i=0; i<(int)size(); ++i) {
			const Geex::Facet& F = (*this)[i] ;
			Geex::vec3 centro;
			centro = F.center() ; //./3. * (boundary_[i].vertex[0] + boundary_[i].vertex[1] + boundary_[i].vertex[2]);

			for(int j = 0; j <= samplerate_; j ++)
				for(int k = 0; k <= samplerate_ - j; k ++)
				{
					Geex::vec3 samplepoint;
					samplepoint = 
						0.99 * j / ((double)samplerate_) * (F.vertex[0] - centro) +
						0.99 * k / ((double)samplerate_) * (F.vertex[1] - centro) +
						0.99 * (1 - j / ((double)samplerate_) - k / ((double)samplerate_)) * (F.vertex[2] - centro);
					
					samplepoint += centro;
					data[curridx*3+0] = samplepoint[0];
					data[curridx*3+1] = samplepoint[1];
					data[curridx*3+2] = samplepoint[2];
					curridx ++;
				}		
		}

		kdtree_ = new ANNKdTree_3(data, (int)size()*(samplerate_+1)*(samplerate_+2)/2, 3);
		delete [] data;	
	}

	vec3 TriMesh::project_to_mesh(const vec3& pt, vec3& norm) {
			Geex::vec3 cent = pt;
			int fid;
			fid  = kdtree_->queryNearestNeighbor(cent[0],cent[1],cent[2]) / 
				((samplerate_ + 1)*(samplerate_ + 2)/2);

			double dist = 1e20;
			Geex::vec3 fp;
			const Geex::Facet& F = (*this)[fid];
			Geex::vec3 v0,v1,v2;
			v0 = F.vertex[0];
			v1 = F.vertex[1];
			v2 = F.vertex[2];
			Geex::vec3 v01, v02;
			v01 = F.vertex[1] - F.vertex[0];
			v02 = F.vertex[2] - F.vertex[0];
			//Geex::vec3 n = cross(v1-v0, v2-v0);
			//n = Geex::normalize(n);
			vec3 n = F.normal();
			n = n / n.length();
			Geex::vec3 tempfp;
			tempfp = cent - dot(cent-F.vertex[0],n)*n;	
			fp = tempfp;
			double a0 = dot(cross(v0-tempfp, v1-tempfp), n); 
			double a1 = dot(cross(v1-tempfp, v2-tempfp), n); 
			double a2 = dot(cross(v2-tempfp, v0-tempfp), n);
			double threshold = -1.0e-14;
			if(a0<threshold || a1<threshold || a2<threshold) 
			{
				Geex::vec3 fp0,fp1,fp2;
				double dist0,dist1,dist2;
				project_to_linesegment(tempfp,v0,v1,fp0,dist0);
				project_to_linesegment(tempfp,v1,v2,fp1,dist1);
				project_to_linesegment(tempfp,v0,v2,fp2,dist2);
				if( (dist0 < dist) || (dist1 < dist) || (dist2 < dist))
				{
					if( (dist0 <= dist1) && (dist0 <= dist2))
					{
						dist = dist0;
						fp = fp0;
						norm = n;
					}
					else if( (dist1 <= dist0) && (dist1 <= dist2))
					{
						dist = dist1;
						fp = fp1;
						norm = n;
					}
					else
					{
						dist = dist2;
						fp = fp2;
						norm = n;
					}
				}
			}
			else 
			{
				double dist0;
				dist0 = (tempfp-cent).length2();
				if(dist0 < dist)
				{
					dist = dist0;
					fp = tempfp;
					norm = n;
				}
			}	
			//std::cerr << norm[0] << ' ' << norm[1] << ' ' << norm[2] << std::endl;
			return fp;
	}
	vec3 TriMesh::project_to_mesh(const vec3& pt, vec3& norm, int &nearestFacet)
	{
		Geex::vec3 cent = pt;
		int fid;
		fid  = kdtree_->queryNearestNeighbor(cent[0],cent[1],cent[2]) / 
			((samplerate_ + 1)*(samplerate_ + 2)/2);

		double dist = 1e20;
		Geex::vec3 fp;
		const Geex::Facet& F = (*this)[fid];
		Geex::vec3 v0,v1,v2;
		v0 = F.vertex[0];
		v1 = F.vertex[1];
		v2 = F.vertex[2];
		Geex::vec3 v01, v02;
		v01 = F.vertex[1] - F.vertex[0];
		v02 = F.vertex[2] - F.vertex[0];
		//Geex::vec3 n = cross(v1-v0, v2-v0);
		//n = Geex::normalize(n);
		vec3 n = F.normal();
		n = n / n.length();
		Geex::vec3 tempfp;
		tempfp = cent - dot(cent-F.vertex[0],n)*n;	
		fp = tempfp;
		double a0 = dot(cross(v0-tempfp, v1-tempfp), n); 
		double a1 = dot(cross(v1-tempfp, v2-tempfp), n); 
		double a2 = dot(cross(v2-tempfp, v0-tempfp), n);
		double threshold = -1.0e-14;
		if(a0<threshold || a1<threshold || a2<threshold) 
		{
			Geex::vec3 fp0,fp1,fp2;
			double dist0,dist1,dist2;
			project_to_linesegment(tempfp,v0,v1,fp0,dist0);
			project_to_linesegment(tempfp,v1,v2,fp1,dist1);
			project_to_linesegment(tempfp,v0,v2,fp2,dist2);
			if( (dist0 < dist) || (dist1 < dist) || (dist2 < dist))
			{
				if( (dist0 <= dist1) && (dist0 <= dist2))
				{
					dist = dist0;
					fp = fp0;
					norm = n;
				}
				else if( (dist1 <= dist0) && (dist1 <= dist2))
				{
					dist = dist1;
					fp = fp1;
					norm = n;
				}
				else
				{
					dist = dist2;
					fp = fp2;
					norm = n;
				}
			}
		}
		else 
		{
			double dist0;
			dist0 = (tempfp-cent).length2();
			if(dist0 < dist)
			{
				dist = dist0;
				fp = tempfp;
				norm = n;
			}
		}	
		//std::cerr << norm[0] << ' ' << norm[1] << ' ' << norm[2] << std::endl;
		nearestFacet = fid;
		return fp;
	}

    void TriMesh::convex_clip(TriMesh& to, const Plane<real>& P, bool close) const {
        to.clear() ;

        bool first = true ;
        vec3 firstP ;

        for(unsigned int i=0; i<size(); i++) {
            const Facet& F = (*this)[i] ;
            int config = 
                ((P.side(F.vertex[0]) > 0.0) ? 1 : 0) |
                ((P.side(F.vertex[1]) > 0.0) ? 2 : 0) |
                ((P.side(F.vertex[2]) > 0.0) ? 4 : 0) ;
            
            switch(config) {

            case 0: {
                // Completely outside
            } break ;

            case 1: {
                vec3 P01 = seg_plane_intersection(F.vertex[0], F.vertex[1], P) ;
                vec3 P20 = seg_plane_intersection(F.vertex[2], F.vertex[0], P) ;
                to.push_back(Facet(F.vertex[0], P01, P20, 1, F.edge_flag[1], F.edge_flag[2])) ;
                if(close) {
                    if(first) { 
                        first = false ; firstP = P01 ; 
                    } else {
                        to.push_back(Facet(firstP, P20, P01, 0, 0, 0)) ;
                    }
                }
            } break ;

            case 2: {
                vec3 P01 = seg_plane_intersection(F.vertex[0], F.vertex[1], P) ;
                vec3 P12 = seg_plane_intersection(F.vertex[1], F.vertex[2], P) ;
                to.push_back(Facet(P01, F.vertex[1], P12, F.edge_flag[0], 1, F.edge_flag[2])) ;
                if(close) {
                    if(first) { 
                        first = false ; firstP = P12 ; 
                    } else {
                        to.push_back(Facet(firstP, P01, P12, 0, 0, 0)) ;
                    }
                }
            } break ;

            case 3: {
                vec3 P12 = seg_plane_intersection(F.vertex[1], F.vertex[2], P) ;
                vec3 P20 = seg_plane_intersection(F.vertex[2], F.vertex[0], P) ;
                to.push_back(Facet(F.vertex[0], P12, P20, 1, F.edge_flag[1], 0)) ;
                to.push_back(Facet(F.vertex[0], F.vertex[1], P12,  F.edge_flag[0], 0, F.edge_flag[2])) ;
                if(close) {
                    if(first) { 
                        first = false ; firstP = P12 ; 
                    } else {
                        to.push_back(Facet(firstP, P20, P12, 0, 0, 0)) ;
                    }
                }
            } break ;

            case 4: {
                vec3 P12 = seg_plane_intersection(F.vertex[1], F.vertex[2], P) ;
                vec3 P20 = seg_plane_intersection(F.vertex[2], F.vertex[0], P) ;
                to.push_back(Facet(P20, P12, F.vertex[2], F.edge_flag[0], F.edge_flag[1], 1)) ;                
                if(close) {
                    if(first) { 
                        first = false ; firstP = P20 ; 
                    } else {
                        to.push_back(Facet(firstP, P12, P20, 0, 0, 0)) ;
                    }
                }
            } break ;

            case 5: {
                vec3 P01 = seg_plane_intersection(F.vertex[0], F.vertex[1], P) ;
                vec3 P12 = seg_plane_intersection(F.vertex[1], F.vertex[2], P) ;
                to.push_back(Facet(F.vertex[0], P01, P12, 1, 0, F.edge_flag[2])) ;                
                to.push_back(Facet(F.vertex[0], P12, F.vertex[2], F.edge_flag[0], F.edge_flag[1], 0)) ;                
                if(close) {
                    if(first) { 
                        first = false ; firstP = P01 ; 
                    } else {
                        to.push_back(Facet(firstP, P12, P01, 0, 0, 0)) ;
                    }
                }
            } break ;

            case 6: {
                vec3 P01 = seg_plane_intersection(F.vertex[0], F.vertex[1], P) ;
                vec3 P20 = seg_plane_intersection(F.vertex[2], F.vertex[0], P) ;
                to.push_back(Facet(P01, F.vertex[1], F.vertex[2], F.edge_flag[0], 0, F.edge_flag[2])) ;  
                to.push_back(Facet(P01, F.vertex[2], P20, F.edge_flag[1], 1, false)) ;                  
                if(close) {
                    if(first) { 
                        first = false ; firstP = P20 ; 
                    } else {
                        to.push_back(Facet(firstP, P01, P20, 0, 0, 0)) ;
                    }
                }
            } break ;

            case 7: {
                // Completely inside
                to.push_back(F) ;
            } break ;
            }
        }
    }

    //inline bool line_intersects_triangle(
    //    const OrientedLine<real>& L, const vec3& q1, const vec3& q2,
    //    const vec3& p1, const vec3& p2, const vec3& p3
    //) {
    //    vec3 N = cross(p2-p1, p3-p1) ;
    //    double side1 = dot(q1-p1,N) ;
    //    double side2 = dot(q2-p1,N) ;
    //    if(
    //        side1 > 0 && side2 > 0 ||
    //        side1 < 0 && side2 < 0
    //    ) { 
    //        return false ; 
    //    }

    //    OrientedLine<real> e1(p1,p2) ;
    //    OrientedLine<real> e2(p2,p3) ;
    //    OrientedLine<real> e3(p3,p1) ;
    //    double s1 = side(e1, L) ;
    //    double s2 = side(e2, L) ;
    //    double s3 = side(e3, L) ;
    //    return (
    //        (s1 > 0 && s2 > 0 && s3 > 0) ||
    //        (s1 < 0 && s2 < 0 && s3 < 0) 
    //    ) ;
    //}

    //bool TriMesh::contains(const vec3& P) const {
    //    int nb_isect = 0 ;
    //    vec3 Q = P ;
    //    Q.z += 1e3 ;
    //    OrientedLine<real> L(P,Q) ;
    //    for(unsigned int i=0; i<size(); i++) {
    //        const Facet& F = (*this)[i] ;
    //        vec3 I ;
    //        if(
    //            line_intersects_triangle(
    //                L, P, Q, F.vertex[0], F.vertex[1], F.vertex[2]
    //            )
    //        ) { nb_isect++ ; }
    //    }
    //    return ((nb_isect & 1) != 0) ;
    //}

    void TriMesh::normalize(vec3& g, double& R, bool disable_normalize_) {
        if (disable_normalize_)
        {
            compute_bounding_box() ;
            return;
        }

        g = vec3(0.0, 0.0, 0.0) ;
        double cnt = 0.0 ;
        for(unsigned int i=0; i<size(); i++) {
            for(unsigned int j=0; j<3; j++) {
                g += (*this)[i].vertex[j] ;
                cnt += 1.0 ;
            }
        }
        if(cnt != 0.0) {
            g = (1.0 / cnt) * g ;
        }
        R = 0.0 ;
        for(unsigned int i=0; i<size(); i++) {
            for(unsigned int j=0; j<3; j++) {
                (*this)[i].vertex[j] -= g ;
                R = gx_max(R,((*this)[i].vertex[j]).length()) ;
            }
        }
        double scal = (R != 0.0) ? (1.0 / R) : 0.0 ;
        for(unsigned int i=0; i<size(); i++) {
            for(unsigned int j=0; j<3; j++) {
                (*this)[i].vertex[j] *= scal ;
                (*this)[i].vertex[j] += vec3(1.0, 1.0, 1.0) ;
                (*this)[i].vertex[j] *= 0.5 ;                

				if(vertices_.size()>0)
					vertex((*this)[i].vertex_index[j]).pos_ = (*this)[i].vertex[j]; // dmyan
            }
        }        

		compute_bounding_box() ;
    }

	void TriMesh::compute_bounding_box() {
		bbox_[0] = vec3( 1e10,  1e10,  1e10) ;
		bbox_[1] = vec3(-1e10, -1e10, -1e10) ;

		for(unsigned int i=0; i<size(); ++i) 
			for(int j=0; j<3; ++j) {
#ifdef min
				bbox_[0][0] = min(bbox_[0][0], (*this)[i].vertex[j][0]) ;
				bbox_[0][1] = min(bbox_[0][1], (*this)[i].vertex[j][1]) ;
				bbox_[0][2] = min(bbox_[0][2], (*this)[i].vertex[j][2]) ;
				bbox_[1][0] = max(bbox_[1][0], (*this)[i].vertex[j][0]) ;
				bbox_[1][1] = max(bbox_[1][1], (*this)[i].vertex[j][1]) ;
				bbox_[1][2] = max(bbox_[1][2], (*this)[i].vertex[j][2]) ;
#else
				bbox_[0][0] = std::min(bbox_[0][0], (*this)[i].vertex[j][0]) ;
				bbox_[0][1] = std::min(bbox_[0][1], (*this)[i].vertex[j][1]) ;
				bbox_[0][2] = std::min(bbox_[0][2], (*this)[i].vertex[j][2]) ;
				bbox_[1][0] = std::max(bbox_[1][0], (*this)[i].vertex[j][0]) ;
				bbox_[1][1] = std::max(bbox_[1][1], (*this)[i].vertex[j][1]) ;
				bbox_[1][2] = std::max(bbox_[1][2], (*this)[i].vertex[j][2]) ;
#endif
			}

//		for(unsigned int i=0; i<vertices_.size(); ++i) {
//#ifdef min
//		    bbox_[0][0] = min(bbox_[0][0], vertices_[i].pos_[0]) ;
//			bbox_[0][1] = min(bbox_[0][1], vertices_[i].pos_[1]) ;
//			bbox_[0][2] = min(bbox_[0][2], vertices_[i].pos_[2]) ;
//			bbox_[1][0] = max(bbox_[1][0], vertices_[i].pos_[0]) ;
//			bbox_[1][1] = max(bbox_[1][1], vertices_[i].pos_[1]) ;
//			bbox_[1][2] = max(bbox_[1][2], vertices_[i].pos_[2]) ;
//#else
//		    bbox_[0][0] = std::min(bbox_[0][0], vertices_[i].pos_[0]) ;
//			bbox_[0][1] = std::min(bbox_[0][1], vertices_[i].pos_[1]) ;
//			bbox_[0][2] = std::min(bbox_[0][2], vertices_[i].pos_[2]) ;
//			bbox_[1][0] = std::max(bbox_[1][0], vertices_[i].pos_[0]) ;
//			bbox_[1][1] = std::max(bbox_[1][1], vertices_[i].pos_[1]) ;
//			bbox_[1][2] = std::max(bbox_[1][2], vertices_[i].pos_[2]) ;
//#endif
//		}
	}

	void TriMesh::merge_same_vertices(const int res, const real toler) {
		Grid   ***grid ; //[res][res][res] ;
		
		grid = new Grid**[res] ;
		for(int i=0; i<res; ++i) {
			grid[i] = new Grid*[res] ;
			for(int j=0; j<res; ++j) {
				grid[i][j] = new Grid[res] ;
			}
		}

		// compute bounding box of mesh
		bbox_[0] = vec3( 1e10,  1e10,  1e10) ;
		bbox_[1] = vec3(-1e10, -1e10, -1e10) ;

		for(unsigned int i=0; i<size(); ++i) {
			for(int j=0; j<3; ++j) {
	#ifdef min
				bbox_[0][0] = min(bbox_[0][0], (*this)[i].vertex[j][0]) ;
				bbox_[0][1] = min(bbox_[0][1], (*this)[i].vertex[j][1]) ;
				bbox_[0][2] = min(bbox_[0][2], (*this)[i].vertex[j][2]) ;
				bbox_[1][0] = max(bbox_[1][0], (*this)[i].vertex[j][0]) ;
				bbox_[1][1] = max(bbox_[1][1], (*this)[i].vertex[j][1]) ;
				bbox_[1][2] = max(bbox_[1][2], (*this)[i].vertex[j][2]) ;
	#else
				bbox_[0][0] = std::min(bbox_[0][0], (*this)[i].vertex[j][0]) ;
				bbox_[0][1] = std::min(bbox_[0][1], (*this)[i].vertex[j][1]) ;
				bbox_[0][2] = std::min(bbox_[0][2], (*this)[i].vertex[j][2]) ;
				bbox_[1][0] = std::max(bbox_[1][0], (*this)[i].vertex[j][0]) ;
				bbox_[1][1] = std::max(bbox_[1][1], (*this)[i].vertex[j][1]) ;
				bbox_[1][2] = std::max(bbox_[1][2], (*this)[i].vertex[j][2]) ;
	#endif			
			}
		}

		double xstep = (bbox_[1][0]-bbox_[0][0])/res ; 
		double ystep = (bbox_[1][1]-bbox_[0][1])/res ; 
		double zstep = (bbox_[1][2]-bbox_[0][2])/res ; 

//		double toler = 1e-7 ; 
		std::priority_queue<FaceVertex> all_verts ;
		std::vector<FaceVertex> face_verts ;

		vertices_.clear() ;

		// partition point data
		for(unsigned int i=0; i<size(); ++i) {
			Geex::Facet& F = (*this)[i] ;
			for(int j=0; j<3; ++j) {
				int xidx = int((F.vertex[j][0]-bbox_[0][0])/xstep) ;
				int yidx = int((F.vertex[j][1]-bbox_[0][1])/ystep) ;
				int zidx = int((F.vertex[j][2]-bbox_[0][2])/zstep) ;
				if(xidx >= res) xidx = res-1 ;
				if(yidx >= res) yidx = res-1 ;
				if(zidx >= res) zidx = res-1 ;
				grid[xidx][yidx][zidx].add_face_vertex(FaceVertex(F.vertex[j], i, j)) ;
			}
		}

		// merge same point with grid
		int count = 0 ;
		for(int i=0; i<res; ++i) 
			for(int j=0; j<res; ++j) 
				for(int k=0; k<res; ++k) {
					std::vector<FaceVertex>& fvs = grid[i][j][k].fvs_ ;
					count += (int)fvs.size(); 
		
					if(fvs.size() > 0) {
						unsigned int start = (int)vertices_.size() ;
						vertices_.push_back(fvs[0].pos_) ;

						for(unsigned int gi=0; gi<fvs.size(); ++gi) {
							FaceVertex& fv = fvs[gi] ; 
							bool bsame = false ;
							for(unsigned int it=start; it<vertices_.size(); ++it) {
								if(distance(fvs[gi].pos_, vertices_[it].pos_) < toler) {
									(*this)[fvs[gi].fidx_].vertex_index[fvs[gi].fvidx_] = it ;
									(*this)[fvs[gi].fidx_].vertex[fvs[gi].fvidx_] = fvs[gi].pos_ ;
									bsame = true ;
									break ;
								}
							}
							if(!bsame) {
								(*this)[fvs[gi].fidx_].vertex_index[fvs[gi].fvidx_] = (int)vertices_.size() ;
								(*this)[fvs[gi].fidx_].vertex[fvs[gi].fvidx_] = fvs[gi].pos_ ;
								vertices_.push_back(fvs[gi].pos_) ;
							}
						}
					}
				}

				assert(count == size()*3) ;

	/*	for(unsigned int i=0; i<size(); ++i) {
			Geex::Facet& F = (*this)[i] ;
			for(int j=0; j<3; ++j) 
				face_verts.push_back(FaceVertex(F.vertex[j], i, j)) ;
		}

		add_vertex(face_verts[0].pos_) ;
		for(unsigned int i=0; i<face_verts.size(); ++i) {
			bool bfind = false ;
			for(unsigned int j=0; j<vertices_.size(); ++j) {
				if(distance(face_verts[i].pos_, vertices_[j].pos_) < toler) {
					(*this)[face_verts[i].fidx_].vertex_index[face_verts[i].fvidx_] = j ;
					(*this)[face_verts[i].fidx_].vertex[face_verts[i].fvidx_] = face_verts[i].pos_ ;
					bfind = true ;
					break ;
				}
			}
			if(!bfind) {
				vertices_.push_back(face_verts[i].pos_) ;
				(*this)[face_verts[i].fidx_].vertex_index[face_verts[i].fvidx_] = vertices_.size()-1 ;
				(*this)[face_verts[i].fidx_].vertex[face_verts[i].fvidx_] = face_verts[i].pos_ ;
			}
		}*/

		


		
		//for(int i=0; i<size(); ++i) {
		//	Geex::Facet& F = (*this)[i] ;
		//	for(int j=0; j<3; ++j) 
		//		all_verts.push(FaceVertex(F.vertex[j], i, j)) ;
		//}

		//FaceVertex pre_v = all_verts.top() ;
		//all_verts.pop() ;
		//add_vertex(pre_v.pos_) ;
		//(*this)[pre_v.fidx_].vertex_index[pre_v.fvidx_] = vertices_.size()-1 ;

		//while(!all_verts.empty()) {
		//	FaceVertex cur_v = all_verts.top() ;
		//	all_verts.pop() ;
		//	
		//	if(!distance(pre_v.pos_, cur_v.pos_)<toler)  {
		//		add_vertex(cur_v.pos_) ;
		//		(*this)[cur_v.fidx_].vertex_index[cur_v.fvidx_] = vertices_.size()-1 ;
		//		(*this)[cur_v.fidx_].vertex[cur_v.fvidx_] = cur_v.pos_ ;

		//		pre_v = cur_v ;
		//	}
		//	else {
		//		(*this)[cur_v.fidx_].vertex_index[cur_v.fvidx_] = vertices_.size()-1 ;
		//		(*this)[cur_v.fidx_].vertex[cur_v.fvidx_] = cur_v.pos_ ;
		//	}
		//}

		for(unsigned int i=0; i<size(); ++i) {
			Geex::Facet& F = (*this)[i] ; 
			for(int j=0; j<3; ++j) {
				vertices_[ F.vertex_index[j] ].faces_.push_back(i);
				F.edge_flag[j] = 1 ;
			}
		}

		nb_vertices_ = (int)vertices_.size() ;

		for(int i=0; i<res; ++i) {
			for(int j=0; j<res; ++j) {
				delete [] grid[i][j] ;
			}
			delete [] grid[i] ;
		}
		delete [] grid ;
	}

    void TriMesh::measure_quality(){
	real qmin=1e10, qav=0, anglemin=1e10, anglemax=0, minangle_ave=0 ;
	int  nless30=0 ;
	real alpha = 4*sqrt(3.0) ;
	real angle_coef = 180. / M_PI ;

	for(int i=0; i<(int)size(); ++i) {
	    Geex::Facet& F = (*this)[i] ;
	    vec3 ab = F.vertex[1] - F.vertex[0] ;
	    vec3 bc = F.vertex[2] - F.vertex[1] ;
	    vec3 ca = F.vertex[0] - F.vertex[2] ;
	    real lab = length(ab) ;
	    real lbc = length(bc) ;
	    real lca = length(ca) ;
	    real anglea = acos(-dot(ca, ab)/(lca*lab))*angle_coef ;
	    real angleb = acos(-dot(ab, bc)/(lab*lbc))*angle_coef ;
	    real anglec = acos(-dot(bc, ca)/(lbc*lca))*angle_coef ;

	    if(anglea < 30) nless30 ++ ;
	    if(angleb < 30) nless30 ++ ;
	    if(anglec < 30) nless30 ++ ;
	    real minangle = gx_min(gx_min(anglea, angleb), anglec) ;
	    anglemin = gx_min(anglemin, minangle) ;
	
		anglemax = gx_max(gx_max(gx_max(anglea, angleb), anglec), anglemax);

	    minangle_ave += minangle ;
	
	    real area = tri_area(F.vertex[0], F.vertex[1], F.vertex[2]) ;
	    real Qi = alpha*area/((lab+lbc+lca)*gx_max(lab, gx_max(lbc, lca))) ;
	    qmin = gx_min(qmin, Qi) ;
	    qav += Qi ;
	}
	qav /= size() ;
	minangle_ave /= size() ;

	std::cout << "Mesh geometry quality: " << std::endl ;
	std::cout << "  AngleMax = " << anglemax << std::endl ;
	std::cout << "  AngleMin = " << anglemin << std::endl ;
	std::cout << "  AverageMinAngle = " << minangle_ave << std::endl ;
	std::cout << "  Qmin = " << qmin << std::endl ;
	std::cout << "  Qav = " << qav << std::endl ;
	std::cout << "  P30 = " << (double)nless30 / (size()*3.) << std::endl << std::endl ;
    }

	void TriMesh::perturb_facet_vertices(double perturb_ratio) {
		compute_bounding_box() ;
		double dist = perturb_ratio*Geex::distance(bbox_[0], bbox_[1]) ;
		for(unsigned int i=0; i<size(); ++i) {
			for(int j=0; j<3; ++j) {
				vec3 delta = dist*vec3(Numeric::random_float64(), Numeric::random_float64(), Numeric::random_float64()) ;
				(*this)[i].vertex[j] += delta ;
			}
		}		
	}

	void TriMesh::save_tri(const std::string& filename) {
		std::ofstream out(filename.c_str()) ;
		out << size() << std::endl ;
		for(unsigned int i=0; i<size(); ++i) {
			const Geex::Facet& F = (*this)[i] ;
			out.precision(16) ;
			out << F.vertex[0] << " " << F.vertex[1] << " " << F.vertex[2] << std::endl ;
		}
			
		out.close() ;
	}

	void TriMesh::load_tri(const std::string& filename) {
		std::ifstream in(filename.c_str()) ;
		int nf ;
		in >> nf ;
		for(int i=0; i<nf; ++i) {
			vec3 p0, p1, p2 ;
			in >> p0 >> p1 >> p2 ;
			push_back(Geex::Facet(p0, p1, p2, 1, 1, 1)) ;
		}
		in.close() ;
	}
    // Note: PluckerMesh could be made a bit faster, by sharing plucker
    // coordinates for each segment shared by two edges.

    //bool PluckerMesh::contains(const vec3& P) const {
    //    vec3 Q = P ;
    //    Q.z += 1e3 ;
    //    PluckerSegment L(P,Q) ;
    //    int nb_isect = 0 ;
    //    for(unsigned int i=0; i<size(); i++) {
    //        if((*this)[i].intersects(L)) {
    //            nb_isect++ ;
    //        }
    //    }
    //    return ((nb_isect & 1) != 0) ;
    //}

}
