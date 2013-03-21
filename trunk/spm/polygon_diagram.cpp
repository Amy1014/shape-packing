#include "polygon_diagram.h"

namespace Geex
{
	
	RestrictedPolygonVoronoiDiagram::RestrictedPolygonVoronoiDiagram() 
	{
		mesh = 0;
		trimesh = 0;
		nb_groups = 0;
		open = false;
	}

	RestrictedPolygonVoronoiDiagram::~RestrictedPolygonVoronoiDiagram()
	{
		mesh->clear();
		std::for_each(samp_pnts.begin(), samp_pnts.end(), std::mem_fun_ref(&VertGroup::clear));
		samp_pnts.clear();
		bounding_pnts.clear();
		edge_face_adjacency.clear();
		rdt_ds.clear();
	}

	void RestrictedPolygonVoronoiDiagram::set_mesh(const std::string& mesh_file_name)
	{
		assert(!mesh);
		mesh = new TopoPolyMesh;
		mesh->load(mesh_file_name);
		//rvd = new RestrictedVoronoiDiagram(del, mesh);
	}
	
	unsigned int RestrictedPolygonVoronoiDiagram::insert_polygons(const Polygon_3& polygon, unsigned int group_id, unsigned int samp_nb)
	{	
		samp_pnts.push_back(VertGroup());
		samp_pnts.back().reserve(samp_nb+1);
		unsigned int nb_pnts = 0; // for debug
		// sample on polygon edges
		double perimeter = 0.0;
		std::vector<double> edge_len(polygon.size());
		for ( unsigned int i = 0; i < polygon.size(); i++ )
		{
			Segment_3 e = polygon.edge(i);
			edge_len[i] = CGAL::sqrt(e.squared_length());
			perimeter += edge_len[i];
		}

		for ( unsigned int i = 0; i < polygon.size(); i++ )
		{
			Segment_3 e = polygon.edge(i);
			unsigned int n = edge_len[i] / perimeter * samp_nb;
			Point_3 src = e.source(), tgt = e.target();
			for ( unsigned int j = 0; j < n+1; j++ )
			{
				Point_3 p = CGAL::ORIGIN + ( (n+1-j)*(src - CGAL::ORIGIN) + j*(tgt - CGAL::ORIGIN) )/(n+1);
				RDT_data_structure::Vertex_handle vh = rdt_ds.create_vertex(RDT_data_structure::Vertex(p, group_id));
				samp_pnts.back().push_back(vh);
				vec3 n;
				vec3 pp = trimesh->project_to_mesh(to_geex_pnt(p), n);
				vh->mp = to_cgal_pnt(pp);
				vh->n = to_cgal_vec(n);
				nb_pnts++;
			}
		}
		return nb_pnts;
	}
	
	void RestrictedPolygonVoronoiDiagram::insert_bounding_points(unsigned int samp_nb)
	{
		const TriMesh::BoundaryEdgeSet& bes = trimesh->getBoundary(); // boundary edges
		const std::vector<std::pair<int,int>>& fes = trimesh->getFeatureEdges(); // feature edges
		unsigned int nbe = insert_bounding_edges(bes.begin(), bes.end(), samp_nb);
		//std::cout<<nbe<<" bounding edges\n";
		unsigned int nfe = insert_bounding_edges(fes.begin(), fes.end(), samp_nb);
		//std::cout<<nfe<<" feature edges\n";
	}
	void RestrictedPolygonVoronoiDiagram::begin_insert()
	{
		assert(!open);
		open = true;
		std::for_each(samp_pnts.begin(), samp_pnts.end(), std::mem_fun_ref(&VertGroup::clear));
		samp_pnts.clear();
		bounding_pnts.clear();
		edge_face_adjacency.clear();
		edge_valance.clear();
		//for (Vertex_iterator vit = rdt_ds.vertices_begin(); vit != rdt_ds.vertices_end(); ++vit)
		//	vit->vd_vertices.clear();
		rdt_ds.clear();
	}
	
	void RestrictedPolygonVoronoiDiagram::end_insert()
	{
		assert(trimesh); 
		unsigned int nb_pnts = rdt_ds.number_of_vertices();
		std::cout<<"Number of sample points: "<<nb_pnts<<std::endl;
		double *pnt_coord = new double[nb_pnts*3];
		double *ptr = pnt_coord;
		std::vector<Vertex_handle> all_vertices;
		all_vertices.reserve(nb_pnts);
		for ( unsigned int i = 0; i < samp_pnts.size(); i++ )
			for (unsigned int j = 0; j < samp_pnts[i].size(); j++)
			{
				//Point_3 p = samp_pnts[i][j]->point_3();
				Point_3 p = samp_pnts[i][j]->mp;
				//vec3 dummyn;
				//vec3 prjp = trimesh->project_to_mesh(to_geex_pnt(p), dummyn);
				*ptr++ = p.x(); *ptr++ = p.y(); *ptr++ = p.z();
				//*ptr++ = prjp.x; *ptr++ = prjp.y; *ptr++ = prjp.z;
				all_vertices.push_back(samp_pnts[i][j]);
			}
		for ( unsigned int i = 0; i < bounding_pnts.size(); i++ )
		{
			Point_3 p = bounding_pnts[i]->point_3();
			*ptr++ = p.x();	*ptr++ = p.y();	*ptr++ = p.z();
			all_vertices.push_back(bounding_pnts[i]);
		}
		Delaunay *del = Delaunay::create("CGAL");
		del->set_vertices(nb_pnts, pnt_coord);
		rdt_ds.set_dimension(2);
		RestrictedVoronoiDiagram *rvd = new RestrictedVoronoiDiagram(del, mesh);
		rvd->for_each_primal_triangle(Construct_RDT_structure(rdt_ds, all_vertices, edge_face_adjacency, edge_valance));
		for (std::map<Vertex_pair, int>::iterator it = edge_valance.begin(); it != edge_valance.end(); ++it)
			if (it->second == 1)
			{
				Vertex_handle v0 = it->first.first, v1 = it->first.second;
				Face_handle infinite_face = rdt_ds.create_face(v0, v1, rdt_ds.infinite_vertex());
				Face_handle f = edge_face_adjacency[it->first];
				int i0 = f->index(v0), i1 = f->index(v1);
				bool mask[] = {true, true, true};
				mask[i0] = mask[i1] = false;
				if (mask[0])	
					rdt_ds.set_adjacency(f, 0, infinite_face, infinite_face->index(rdt_ds.infinite_vertex()));
				if (mask[1])	
					rdt_ds.set_adjacency(f, 1, infinite_face, infinite_face->index(rdt_ds.infinite_vertex()));
				if (mask[2])	
					rdt_ds.set_adjacency(f, 2, infinite_face, infinite_face->index(rdt_ds.infinite_vertex()));
			}
		delete rvd;
		delete del;
		delete pnt_coord;
		open = false;
		edge_valance.clear();
	}
	void RestrictedPolygonVoronoiDiagram::compute_clipped_VD(std::vector<Plane_3>& clipping_planes)
	{
		if (nb_groups == 0)
			return;
		assert(clipping_planes.size() == nb_groups);
		for (Vertex_iterator vit = rdt_ds.vertices_begin(); vit != rdt_ds.vertices_end(); ++vit)
			vit->vd_vertices.clear();
		for (unsigned int i = 0; i < nb_groups; i++)
		{
			const VertGroup& vg = samp_pnts[i];
			const Plane_3& pln = clipping_planes[i];
			for (unsigned int j = 0; j < vg.size(); j++)
			{
				//RDT_data_structure::Edge_circulator cur_edge, nxt_edge, end, start_edge, end_edge;
				//Plane_3 curpln, nxtpln;
				//bool start_found = false, end_found = false;
				//end = cur_edge = rdt_ds.incident_edges(vg[j]);
				//nxt_edge = cur_edge;
				//++nxt_edge;
				//do 
				//{
				//	Face_handle cur_fc = cur_edge->first, nxt_fc = nxt_edge->first;
				//	int cur_idx = cur_edge->second, nxt_idx = nxt_edge->second;
				//	Vertex_handle cur_other_vert = cur_fc->vertex(cur_fc->cw(cur_idx))==vg[j]?
				//						cur_fc->vertex(cur_fc->ccw(cur_idx)):cur_fc->vertex(cur_fc->cw(cur_idx));
				//	Vertex_handle nxt_other_vert = nxt_fc->vertex(nxt_fc->cw(nxt_idx))==vg[j]?
				//						nxt_fc->vertex(nxt_fc->ccw(nxt_idx)):nxt_fc->vertex(nxt_fc->cw(nxt_idx));
				//	if ( cur_other_vert->group_id == vg[j]->group_id 
				//			&& nxt_other_vert->group_id != vg[j]->group_id )
				//	{
				//		start_edge = cur_edge;
				//		start_found = true;
				//		if (end_found)
				//			break;					
				//	}
				//	else if ( cur_other_vert->group_id != vg[j]->group_id
				//			&& nxt_other_vert->group_id == vg[j]->group_id )
				//	{
				//		end_edge = cur_edge;
				//		end_found = true;
				//		if (start_found)
				//			break;					
				//	}
				//	++cur_edge; ++nxt_edge;

				//} while (cur_edge != end );
				//if (!start_found && !end_found )
				//	continue;
				//nxt_edge = cur_edge = start_edge;
				//++nxt_edge;
				//++end_edge;
				//do 
				//{
				//	// check whether inside the same group
				//	Face_handle cur_fc = cur_edge->first, nxt_fc = nxt_edge->first;
				//	int cur_idx = cur_edge->second, nxt_idx = nxt_edge->second;
				//	Vertex_handle cur_other_vert = cur_fc->vertex(cur_fc->cw(cur_idx))==vg[j]?
				//								cur_fc->vertex(cur_fc->ccw(cur_idx)):cur_fc->vertex(cur_fc->cw(cur_idx));
				//	Vertex_handle nxt_other_vert = nxt_fc->vertex(nxt_fc->cw(nxt_idx))==vg[j]?
				//								nxt_fc->vertex(nxt_fc->ccw(nxt_idx)):nxt_fc->vertex(nxt_fc->cw(nxt_idx));
				//	assert(cur_other_vert != nxt_other_vert);
				//	if (vg[j]->group_id != cur_other_vert->group_id || vg[j]->group_id != nxt_other_vert->group_id)
				//	{
				//		// get two adjacent bisector planes
				//		curpln = get_bisector(vg[j], cur_other_vert), nxtpln = get_bisector(vg[j], nxt_other_vert);
				//		CGAL::Object o = CGAL::intersection(pln, curpln, nxtpln);
				//		if (const Point_3 *interp = CGAL::object_cast<Point_3>(&o))
				//			vg[j]->vd_vertices.push_back(*interp);
				//	}
				//	++cur_edge;	++nxt_edge;
				//} while (cur_edge != end_edge);
				int gid = vg[j]->group_id;
				RDT_data_structure::Face_circulator fc = rdt_ds.incident_faces(vg[j]), end, nxt_fc;
				nxt_fc = end = fc;
				++nxt_fc;
				do 
				{
					//Vertex_handle v0 = fc->vertex(0), v1 = fc->vertex(1), v2 = fc->vertex(2);
					//if ((v0->group_id == gid && v1->group_id == gid && v2->group_id == gid))
					//	break;
					Vertex_handle v = fc->vertex(fc->ccw(fc->index(vg[j])));
					if (v->group_id != vg[j]->group_id)
						break;
					++fc; ++nxt_fc;
				} while (fc != end);
				//if (fc == end)
				//	continue;
				//end = fc = nxt_fc;
				end = fc;
				do 
				{
					//Vertex_handle v0 = fc->vertex(0), v1 = fc->vertex(1), v2 = fc->vertex(2);
					Vertex_handle v = fc->vertex(fc->ccw(fc->index(vg[j])));
					//if (v0->group_id != gid || v1->group_id != gid || v2->group_id != gid)
					if (v->group_id != vg[j]->group_id)
					{
						Vertex_handle v0 = fc->vertex(0), v1 = fc->vertex(1), v2 = fc->vertex(2);
						Point_3 c = CGAL::circumcenter(v0->point_3(), v1->point_3(), v2->point_3());
						vg[j]->vd_vertices.push_back(c);
					}
					++fc;
				} while (fc != end);
				std::vector<Point_3>& vd_vertices = vg[j]->vd_vertices;
				for (unsigned int k = 0; k <  vd_vertices.size(); k++)
				{
					vd_vertices[k] = pln.projection(vd_vertices[k]);
				}
			}
		}
	}

	bool RestrictedPolygonVoronoiDiagram::is_delaunay_edge(const Edge& e)
	{
		Face_handle fi = e.first;
		int i = e.second;
		Face_handle fj = fi->neighbor(i);
		int j = fj->index(fi);
		Vertex_handle vi = fi->vertex(i), vi_cw = fi->vertex(fi->cw(i)), vi_ccw = fi->vertex(fi->ccw(i));
		Vertex_handle vj = fj->vertex(j), vj_cw = fj->vertex(fj->cw(j)), vj_ccw = fj->vertex(fj->ccw(j));
		//double a2 = CGAL::squared_distance(vi->point_3(), vi_cw->point_3()),
		//	b2 = CGAL::squared_distance(vi->point_3(), vi_ccw->point_3()),
		//	c2 = CGAL::squared_distance(vi_cw->point_3(), vi_ccw->point_3());
		double a2 = CGAL::squared_distance(vi->mp, vi_cw->mp),
			b2 = CGAL::squared_distance(vi->mp, vi_ccw->mp),
			c2 = CGAL::squared_distance(vi_cw->mp, vi_ccw->mp);
		double cos_ang_i = (a2 + b2 - c2) / (2*std::sqrt(a2*b2));
		double sin_ang_i = std::sqrt(1-cos_ang_i*cos_ang_i);
		//a2 = CGAL::squared_distance(vj->point_3(), vj_cw->point_3());
		//b2 = CGAL::squared_distance(vj->point_3(), vj_ccw->point_3());
		//c2 = CGAL::squared_distance(vj_cw->point_3(), vj_ccw->point_3());
		a2 = CGAL::squared_distance(vj->mp, vj_cw->mp);
		b2 = CGAL::squared_distance(vj->mp, vj_ccw->mp);
		c2 = CGAL::squared_distance(vj_cw->mp, vj_ccw->mp);
		double cos_ang_j = (a2 + b2 - c2)/(2*std::sqrt(a2*b2));
		double sin_ang_j = std::sqrt(1 - cos_ang_j*cos_ang_j);
		double sin_ij = sin_ang_i*cos_ang_j + sin_ang_j*cos_ang_i;
		if (sin_ij >= 0.0) // sum of two opposite angles not larger than pi
			return true;
		else
			return false;
	}

	inline bool RestrictedPolygonVoronoiDiagram::is_interior_edge(const Edge& e)
	{
		Face_handle f = e.first;
		int i = e.second;
		Vertex_handle v0 = f->vertex(f->cw(i)), v1 = f->vertex(f->ccw(i));
		return (v0->group_id >= 0) && (v1->group_id >= 0);
	}
	void RestrictedPolygonVoronoiDiagram::add_quadrilateral_edge(Face_handle f, int i, std::queue<Edge>& q)
	{
		Face_handle tf = f->neighbor(i);
		assert(tf != Face_handle());
		int ti = tf->index(f);
		Edge te(tf, ti);
		if (is_interior_edge(te) && !tf->visited[ti] && !rdt_ds.is_infinite(tf) && !rdt_ds.is_infinite(f))
		{
			q.push(te);
			assert(is_interior_edge(te));
		}
		else
			f->visited[i] = tf->visited[ti] = true;		
	}
	void RestrictedPolygonVoronoiDiagram::add_quadrilateral_edge(Vertex_handle v0, Vertex_handle v1, 
												std::queue<Vertex_pair>&q, std::set<Vertex_pair>& s)
	{
		std::map<Vertex_pair, Face_handle>::iterator it;
		if ( (it = edge_face_adjacency.find(Vertex_pair(v0, v1))) == edge_face_adjacency.end() )
			it = edge_face_adjacency.find(Vertex_pair(v1, v0));
		assert(it != edge_face_adjacency.end());

		if (s.find(Vertex_pair(v0, v1)) != s.end())
			return;
		if (s.find(Vertex_pair(v1, v0)) != s.end())
			return;
		if (v0->group_id < 0 && v1->group_id < 0)
			return;
		q.push(Vertex_pair(v0, v1));
		s.insert(Vertex_pair(v0, v1));
	}
	void RestrictedPolygonVoronoiDiagram::iDT_update()
	{
		std::set<Vertex_pair> edges_visited;
		//for (Edge_iterator eit = rdt_ds.edges_begin(); eit != rdt_ds.edges_end(); ++eit)
		//{
		//	Face_handle f = eit->first;
		//	int i = eit->second;
		//	Vertex_handle v0 = f->vertex(f->cw(i)), v1 = f->vertex(f->ccw(i));
		//	Vertex_pair e(v0, v1);
		//	std::map<Vertex_pair, Face_handle>::iterator it;
		//	if ( (it = edge_face_adjacency.find(e)) == edge_face_adjacency.end() )
		//		it = edge_face_adjacency.find(Vertex_pair(v1, v0));
		//	assert(it != edge_face_adjacency.end());
		//}
		RDT_data_structure::Edge_iterator eit = rdt_ds.edges_begin();
		while (true)
		{
			if (is_interior_edge(*eit) && !rdt_ds.is_infinite(*eit))
				break;
			++eit;
		}
		std::queue<Vertex_pair> q;
		Face_handle f = eit->first;
		int i = eit->second;
		Vertex_handle v0 = f->vertex(f->cw(i)), v1 = f->vertex(f->ccw(i));
		q.push(Vertex_pair(v0, v1));
		edges_visited.insert(Vertex_pair(v0, v1));
		while (!q.empty())
		{
			Vertex_pair e = q.front();
			Vertex_handle v0 = e.first, v1 = e.second;
			std::map<Vertex_pair, Face_handle>::iterator it = edge_face_adjacency.find(Vertex_pair(v0, v1));
			if ( it == edge_face_adjacency.end() )
				it = edge_face_adjacency.find(Vertex_pair(v1, v0));
			assert(it != edge_face_adjacency.end());
			Face_handle f = it->second;
			bool mask[] = {true, true, true};
			mask[f->index(v0)] = mask[f->index(v1)] = false;
			int i;
			if (mask[0]) i = 0;
			if (mask[1]) i = 1;
			if (mask[2]) i = 2;
			Face_handle nf = f->neighbor(i);
			int ni = nf->index(f);
			Vertex_handle P = f->vertex(i), Q = nf->vertex(ni);
			//Vertex_handle r = nf->vertex(nf->cw(ni)), s = nf->vertex(nf->ccw(ni));
			if (!is_delaunay_edge(Edge(f, i)))
			{
				rdt_ds.flip(f, i);
				edge_face_adjacency.erase(it);
				edge_face_adjacency[Vertex_pair(P,Q)] = f;
				//edge_face_adjacency[Vertex_pair(v0, P)] = edge_face_adjacency[Vertex_pair(v0, Q)] = f;
				if (f->has_vertex(v0))
				{
					edge_face_adjacency[Vertex_pair(v0, P)] = edge_face_adjacency[Vertex_pair(v0, Q)] = f;
					assert(f->has_vertex(v0) && f->has_vertex(P) && f->has_vertex(Q));
					edge_face_adjacency[Vertex_pair(v1, P)] = edge_face_adjacency[Vertex_pair(v1, Q)] = nf;
					assert(nf->has_vertex(v1) && nf->has_vertex(P) && nf->has_vertex(Q));
				}
				else if (f->has_vertex(v1))
				{
					edge_face_adjacency[Vertex_pair(v1, P)] = edge_face_adjacency[Vertex_pair(v1, Q)] = f;
					assert(f->has_vertex(v1) && f->has_vertex(P) && f->has_vertex(Q));
					edge_face_adjacency[Vertex_pair(v0, P)] = edge_face_adjacency[Vertex_pair(v0, Q)] = nf;
					assert(nf->has_vertex(v0) && nf->has_vertex(P) && nf->has_vertex(Q));
				}
				edges_visited.insert(Vertex_pair(P, Q));
			}
			//edges_visited.insert(e);
			add_quadrilateral_edge(v0, P, q, edges_visited);
			add_quadrilateral_edge(v0, Q, q, edges_visited);
			add_quadrilateral_edge(v1, P, q, edges_visited);
			add_quadrilateral_edge(v1, Q, q, edges_visited);
			
			q.pop();
		}
	}
	void RestrictedPolygonVoronoiDiagram::Construct_RDT_structure::find_adjacency(const Vertex_handle& v0, const Vertex_handle& v1, int v2, const Face_handle& f) const
	{
		Vertex_pair vp01(v0, v1), vp10(v1, v0);
		std::map<Vertex_pair, Face_handle>::iterator it;
		if ( (it = edge_face_adjacency.find(vp01)) != edge_face_adjacency.end() )
		{
			int idx0 = it->second->index(v0), idx1 = it->second->index(v1);
			bool mask[] = {true, true, true};
			mask[idx0] = false;
			mask[idx1] = false;
			if (mask[0])	rdt_ds.set_adjacency(it->second, 0, f, v2);
			if (mask[1])	rdt_ds.set_adjacency(it->second, 1, f, v2);
			if (mask[2])	rdt_ds.set_adjacency(it->second, 2, f, v2);
			edge_valance[vp01]++;
		}
		else if ( (it = edge_face_adjacency.find(vp10)) != edge_face_adjacency.end() )
		{
			int idx0 = it->second->index(v0), idx1 = it->second->index(v1);
			bool mask[] = {true, true, true};
			mask[idx0] = false;
			mask[idx1] = false;
			if (mask[0])	rdt_ds.set_adjacency(it->second, 0, f, v2);
			if (mask[1])	rdt_ds.set_adjacency(it->second, 1, f, v2);
			if (mask[2])	rdt_ds.set_adjacency(it->second, 2, f, v2);
			edge_valance[vp10]++;
		}
		else
		{
			edge_face_adjacency[vp01] = f;
			edge_valance[vp01] = 1;
			//F_face = ace_handle infiniterdt_ds.create_face(v0, v1, rdt_ds.infinite_vertex());
			//rdt_ds.set_adjacency(infinite_face, infinite_face->index(rdt_ds.infinite_vertex()), f, v2);
		}
	}
	void RestrictedPolygonVoronoiDiagram::Construct_RDT_structure::operator()(unsigned int i, unsigned int j, unsigned int k) const
	{
		Face_handle f = rdt_ds.create_face(samp_pnts[i], samp_pnts[j], samp_pnts[k]);
		//Vector_3 n = samp_pnts[i]->n + samp_pnts[j]->n + samp_pnts[k]->n;
		//Vector_3 nij(samp_pnts[i]->point_3(), samp_pnts[j]->point_3());
		//Vector_3 njk(samp_pnts[j]->point_3(), samp_pnts[k]->point_3());
		//Vector_3 det = CGAL::cross_product(nij, njk);
		//if (det*n > 0.0)
		//{
		//	f->set_vertex(0, samp_pnts[i]);
		//	f->set_vertex(f->ccw(0), samp_pnts[j]);
		//	f->set_vertex(f->cw(0), samp_pnts[k]);
		//}
		//else
		//{
		//	f->set_vertex(0, samp_pnts[i]);
		//	f->set_vertex(f->cw(0), samp_pnts[j]);
		//	f->set_vertex(f->ccw(0), samp_pnts[k]);
		//}
		find_adjacency(samp_pnts[i], samp_pnts[j], f->index(samp_pnts[k]), f);
		find_adjacency(samp_pnts[j], samp_pnts[k], f->index(samp_pnts[i]), f);
		find_adjacency(samp_pnts[k], samp_pnts[i], f->index(samp_pnts[j]), f);
		samp_pnts[i]->set_face(f);
		samp_pnts[j]->set_face(f);
		samp_pnts[k]->set_face(f);
	}

}