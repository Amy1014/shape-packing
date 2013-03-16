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
		//std::for_each(clipped_VD.begin(), clipped_VD.end(), std::mem_fun_ref(&std::vector<Point_3>::clear));
		//clipped_VD.clear();
		glDeleteLists(RDT_disp_list, 1);
		glDeleteLists(clipped_VD_disp_list, 1);
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
				nb_pnts++;
			}
		}

		cents.push_back(polygon.centroid());
		return nb_pnts;
	}
	
	void RestrictedPolygonVoronoiDiagram::insert_bounding_points(unsigned int samp_nb)
	{
		const TriMesh::BoundaryEdgeSet& bes = trimesh->getBoundary(); // boundary edges
		const std::vector<std::pair<int,int>>& fes = trimesh->getFeatureEdges(); // feature edges
		unsigned int nbe = insert_bounding_edges(bes.begin(), bes.end(), samp_nb);
		std::cout<<nbe<<" bounding edges\n";
		unsigned int nfe = insert_bounding_edges(fes.begin(), fes.end(), samp_nb);
		std::cout<<nfe<<" feature edges\n";
	}
	void RestrictedPolygonVoronoiDiagram::begin_insert()
	{
		assert(!open);
		open = true;
		std::for_each(samp_pnts.begin(), samp_pnts.end(), std::mem_fun_ref(&VertGroup::clear));
		samp_pnts.clear();
		bounding_pnts.clear();
		edge_face_adjacency.clear();
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
				Point_3 p = samp_pnts[i][j]->point_3();
				vec3 dummyn;
				vec3 prjp = trimesh->project_to_mesh(to_geex_pnt(p), dummyn);
				*ptr++ = prjp.x; *ptr++ = prjp.y; *ptr++ = prjp.z;
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
		//rvd->for_each_primal_triangle(Construct_RDT_structure(rdt_ds, all_vertices, edge_face_adjacency));
		RestrictedVoronoiDiagram *rvd = new RestrictedVoronoiDiagram(del, mesh);
		rvd->for_each_primal_triangle(Construct_RDT_structure(rdt_ds, all_vertices, edge_face_adjacency));
		delete rvd;
		delete del;
		delete pnt_coord;
		open = false;

	}
	void RestrictedPolygonVoronoiDiagram::compute_clipped_VD(std::vector<Plane_3>& clipping_planes)
	{
		if (nb_groups == 0)
			return;
		assert(clipping_planes.size() == nb_groups);
		//clipped_VD.reserve(nb_groups);
		for (unsigned int i = 0; i < nb_groups; i++)
		{
			const VertGroup& vg = samp_pnts[i];
			const Plane_3& pln = clipping_planes[i];
			//clipped_VD.push_back(std::vector<Point_3>());
			for (unsigned int j = 0; j < vg.size(); j++)
			{
				RDT_data_structure::Edge_circulator cur_edge, nxt_edge, end;
				Plane_3 curpln, nxtpln;
				//assert(vg[j]->face() != Face_handle());
				end = cur_edge = rdt_ds.incident_edges(vg[j]);
				nxt_edge = cur_edge;
				++nxt_edge; 			
				do 
				{
					// check whether inside the same group
					Face_handle cur_fc = cur_edge->first, nxt_fc = nxt_edge->first;
					int cur_idx = cur_edge->second, nxt_idx = nxt_edge->second;
					Vertex_handle cur_other_vert = cur_fc->vertex(cur_fc->cw(cur_idx))==vg[j]?
												cur_fc->vertex(cur_fc->ccw(cur_idx)):cur_fc->vertex(cur_fc->cw(cur_idx));
					Vertex_handle nxt_other_vert = nxt_fc->vertex(nxt_fc->cw(nxt_idx))==vg[j]?
												nxt_fc->vertex(nxt_fc->ccw(nxt_idx)):nxt_fc->vertex(nxt_fc->cw(nxt_idx));
					assert(cur_other_vert != nxt_other_vert);
					if (vg[j]->group_id != cur_other_vert->group_id || vg[j]->group_id != nxt_other_vert->group_id)
					{
						// get two adjacent bisector planes
						curpln = get_bisector(vg[j], cur_other_vert), nxtpln = get_bisector(vg[j], nxt_other_vert);
						CGAL::Object o = CGAL::intersection(pln, curpln, nxtpln);
						if (const Point_3 *interp = CGAL::object_cast<Point_3>(&o))
						{
							//clipped_VD.back().push_back(*interp);
							vg[j]->vd_vertices.push_back(*interp);
						}
					}
					++cur_edge;	++nxt_edge;
				} while (cur_edge != end);
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
		double a2 = CGAL::squared_distance(vi->point_3(), vi_cw->point_3()),
			b2 = CGAL::squared_distance(vi->point_3(), vi_ccw->point_3()),
			c2 = CGAL::squared_distance(vi_cw->point_3(), vi_ccw->point_3());
		double cos_ang_i = (a2 + b2 - c2) / (2*std::sqrt(a2*b2));
		double sin_ang_i = std::sqrt(1-cos_ang_i*cos_ang_i);
		a2 = CGAL::squared_distance(vj->point_3(), vj_cw->point_3());
		b2 = CGAL::squared_distance(vj->point_3(), vj_ccw->point_3());
		c2 = CGAL::squared_distance(vj_cw->point_3(), vj_ccw->point_3());
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
		return (v0->group_id >= 0) || (v1->group_id >= 0);
	}
	inline void RestrictedPolygonVoronoiDiagram::add_quadrilateral_edge(Face_handle f, int i, std::queue<Edge>& q)
	{
		Face_handle t = f->neighbor(i);
		int ti = t->index(f);
		Edge te(t, ti);
		if (!t->visited[ti] && is_interior_edge(te))
			q.push(te);
		else
			f->visited[i] = true;		
	}
	void RestrictedPolygonVoronoiDiagram::iDT_update()
	{
		// mark all edges as unvisited
		for (RDT_data_structure::Face_iterator fit = rdt_ds.faces_begin(); fit != rdt_ds.faces_end(); ++fit)
			fit->visited[0] = fit->visited[1] = fit->visited[2] = false;
		//choose an arbitrary interior edge
		RDT_data_structure::Edge_iterator eit = rdt_ds.edges_begin();
		while (true)
		{
			if (is_interior_edge(*eit))
				break;
			++eit;
		}
		std::queue<Edge> q;
		q.push(*eit);
		while (!q.empty())
		{
			Edge e = q.front();
			q.pop();
			Face_handle f = e.first;
			int i = e.second;
			Face_handle nf = f->neighbor(i);
			int ni = nf->index(f);
			if (!is_delaunay_edge(e))
			{
				rdt_ds.flip(e.first, e.second);
				// mark this diagonal as visisted, in both faces, respectively
				f->visited[f->ccw(i)] = nf->visited[nf->ccw(ni)] = true;
				// add the four edges of the quadrilateral
				add_quadrilateral_edge(f, i, q);
				add_quadrilateral_edge(f, f->cw(i), q);
				add_quadrilateral_edge(nf, ni, q);
				add_quadrilateral_edge(nf, nf->cw(ni), q);
			}
			else
			{
				f->visited[i] = nf->visited[ni] = true;
				add_quadrilateral_edge(f, f->cw(i), q);
				add_quadrilateral_edge(f, f->ccw(i), q);
				add_quadrilateral_edge(nf, nf->cw(ni), q);
				add_quadrilateral_edge(nf, nf->ccw(ni), q);
			}
		}
	}

//#if _DEBUG
	void RestrictedPolygonVoronoiDiagram::draw_clipped_VD()
	{
		if (mesh)
		{
			if (!glIsList(clipped_VD_disp_list))
			{
				unsigned int max_uint = std::numeric_limits<unsigned int>::max();
				clipped_VD_disp_list = glGenLists(1);
				glNewList(clipped_VD_disp_list, GL_COMPILE);
				glDisable(GL_LIGHTING);
				for (unsigned int i = 0; i < nb_groups; i++)
				{
					setGroupColor(i, nb_groups);
					const VertGroup& vg = samp_pnts[i];
					glBegin(GL_TRIANGLES);
					for (unsigned int j = 0; j < vg.size(); j++)
					{	
						const std::vector<Point_3>& verts = vg[j]->vd_vertices;
						unsigned int sz = verts.size();
						if (sz > 1)
							for (unsigned int k = 0; k < sz; k++)
							{
								glVertex3d(cents[i].x(), cents[i].y(), cents[i].z());
								glVertex3d(verts[k].x(), verts[k].y(), verts[k].z());
								glVertex3d(verts[(k+1)%sz].x(), verts[(k+1)%sz].y(), verts[(k+1)%sz].z());
							}
					}
					glEnd();
				}
				glEnable(GL_LIGHTING);
				glEndList();
			}
			glCallList(clipped_VD_disp_list);
		}
	}
//#endif


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
		}
		else
			edge_face_adjacency[vp01] = f;
	}
	void RestrictedPolygonVoronoiDiagram::Construct_RDT_structure::operator()(unsigned int i, unsigned int j, unsigned int k) const
	{
		Face_handle f = rdt_ds.create_face(samp_pnts[i], samp_pnts[j], samp_pnts[k]);
		assert( f != Face_handle() );
		find_adjacency(samp_pnts[i], samp_pnts[j], f->index(samp_pnts[k]), f);
		find_adjacency(samp_pnts[j], samp_pnts[k], f->index(samp_pnts[i]), f);
		find_adjacency(samp_pnts[k], samp_pnts[i], f->index(samp_pnts[j]), f);
		samp_pnts[i]->set_face(f);
		samp_pnts[j]->set_face(f);
		samp_pnts[k]->set_face(f);
	}

}