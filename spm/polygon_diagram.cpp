#include "polygon_diagram.h"

namespace Geex
{
	
	RestrictedPolygonVoronoiDiagram::RestrictedPolygonVoronoiDiagram() 
	{
		mesh = 0;
		trimesh = 0;
		//rvd = 0;
		//pnt_coord = 0;
		nb_groups = 0;
		open = false;
		//del = Delaunay::create("CGAL");
	}

	RestrictedPolygonVoronoiDiagram::~RestrictedPolygonVoronoiDiagram()
	{
		//delete rvd;
		//delete del;
		//delete pnt_coord;
		mesh->clear();
		std::for_each(samp_pnts.begin(), samp_pnts.end(), std::mem_fun_ref(&VertGroup::clear));
		samp_pnts.clear();
		std::for_each(clipped_VD.begin(), clipped_VD.end(), std::mem_fun_ref(&std::vector<Point_3>::clear));
		clipped_VD.clear();
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

	void RestrictedPolygonVoronoiDiagram::begin_insert()
	{
		assert(!open);
		open = true;
		std::for_each(samp_pnts.begin(), samp_pnts.end(), std::mem_fun_ref(&VertGroup::clear));
		samp_pnts.clear();
		std::for_each(clipped_VD.begin(), clipped_VD.end(), std::mem_fun_ref(&std::vector<Point_3>::clear));
		clipped_VD.clear();
		//if (!pnt_coord) delete pnt_coord;
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
		clipped_VD.reserve(nb_groups);
		for (unsigned int i = 0; i < nb_groups; i++)
		{
			const VertGroup& vg = samp_pnts[i];
			const Plane_3& pln = clipping_planes[i];
			clipped_VD.push_back(std::vector<Point_3>());
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
							clipped_VD.back().push_back(*interp);
							vg[j]->vd_vertices.push_back(*interp);
						}
					}
					++cur_edge;	++nxt_edge;
				} while (cur_edge != end);
			}
		}
	}

//#if _DEBUG
	void RestrictedPolygonVoronoiDiagram::draw_DT()
	{
		if (mesh)
		{
			if (!glIsList(RDT_disp_list))
			{
				RDT_disp_list = glGenLists(1);
				glNewList(RDT_disp_list, GL_COMPILE);
				glLineWidth(1.5f);
				glColor3f(0.8f, 0.0f, 0.0f);
				glDisable(GL_LIGHTING);
				glBegin(GL_LINES);
				for (RDT_data_structure::Edge_iterator eit = rdt_ds.edges_begin(); eit != rdt_ds.edges_end(); ++eit)
				{
					RDT_data_structure::Face_handle f = eit->first;
					int i = eit->second;
					RDT_data_structure::Vertex_handle vh0 = f->vertex(f->cw(i)), vh1 = f->vertex(f->ccw(i));
					Point_3 p0 = vh0->point_3(), p1 = vh1->point_3();
					glVertex3d(p0.x(), p0.y(), p0.z());
					glVertex3d(p1.x(), p1.y(), p1.z());				
				}
				glEnd();
				glEnable(GL_LIGHTING);
				glEndList();
			}
			glCallList(RDT_disp_list);
		}
	}

	void RestrictedPolygonVoronoiDiagram::draw_clipped_VD()
	{
		if (mesh && clipped_VD.size() > 0)
		{
			if (!glIsList(clipped_VD_disp_list))
			{
				unsigned int max_uint = std::numeric_limits<unsigned int>::max();
				clipped_VD_disp_list = glGenLists(1);
				glNewList(clipped_VD_disp_list, GL_COMPILE);
				glDisable(GL_LIGHTING);
				//for (unsigned int i = 0; i < clipped_VD.size(); i++)
				//{
				//	setGroupColor(i, clipped_VD.size());
				//	Point_3 c = CGAL::centroid(clipped_VD[i].begin(), clipped_VD[i].end(), CGAL::Dimension_tag<0>());
				//	glBegin(GL_TRIANGLES);
				//	unsigned int sz = clipped_VD[i].size();
				//	for (unsigned int j = 0; j < sz; j++)
				//	{
				//		glVertex3d(c.x(), c.y(), c.z());
				//		glVertex3d(clipped_VD[i][j].x(), clipped_VD[i][j].y(), clipped_VD[i][j].z());
				//		glVertex3d(clipped_VD[i][(j+1)%sz].x(), clipped_VD[i][(j+1)%sz].y(), clipped_VD[i][(j+1)%sz].z());
				//	}
				//	glEnd();
				//}
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