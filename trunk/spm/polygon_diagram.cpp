#include "polygon_diagram.h"

namespace Geex
{
	int nb_invalid_edges = 0;

	double RestrictedPolygonVoronoiDiagram::pi = 3.141592653589793;

	//std::ofstream of("test_rdt.obj");

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
		//samp_pnts.push_back(VertGroup());
		//samp_pnts.back().reserve(samp_nb+1);
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
				//all_points.push_back(p);
				//p->group_id = group_id;
				p = e.supporting_line().projection(p);
				vec3 dv;
				vec3 pp = trimesh->project_to_mesh(to_geex_pnt(p), dv);
				all_points.push_back(MyPoint(p, to_cgal_pnt(pp), group_id));
				//all_points.back().mp = to_cgal_pnt(pp);
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
		unsigned int nfe = insert_bounding_edges(fes.begin(), fes.end(), samp_nb);
	}
	void RestrictedPolygonVoronoiDiagram::begin_insert()
	{
		assert(!open);
		open = true;
		std::for_each(samp_pnts.begin(), samp_pnts.end(), std::mem_fun_ref(&VertGroup::clear));
		samp_pnts.clear();
		bounding_pnts.clear();
		all_points.clear();
		rdt_ds.clear();
	}
	
	void RestrictedPolygonVoronoiDiagram::end_insert()
	{
		assert(trimesh); 
		//unsigned int nb_pnts = rdt_ds.size_of_vertices();
		unsigned int nb_pnts = all_points.size();
		std::cout<<"Number of sample points: "<<nb_pnts<<std::endl;
		double *pnt_coord = new double[nb_pnts*3];
		double *ptr = pnt_coord;
		for (unsigned int i = 0; i < all_points.size(); i++)
		{
			Point_3 p = all_points[i].prj_pnt;
			*ptr++ = p.x(); *ptr++ = p.y(); *ptr++ = p.z();
		}
		Delaunay *del = Delaunay::create("CGAL");
		del->set_vertices(nb_pnts, pnt_coord);
		//rdt_ds.set_dimension(2);
		RestrictedVoronoiDiagram *rvd = new RestrictedVoronoiDiagram(del, mesh);
		rdt_ds.delegate(Builder(*rvd, all_points, *trimesh));
		assert(nb_pnts == rdt_ds.size_of_vertices());
		assert(rdt_ds.is_pure_triangle());
		assert(rdt_ds.is_valid(false, 0));
		for (Vertex_iterator vit = rdt_ds.vertices_begin(); vit != rdt_ds.vertices_end(); ++vit)
		{
			if (vit->group_id >= 0)
				samp_pnts[vit->group_id].push_back(vit);
			else
				bounding_pnts.push_back(vit);
		}
		delete rvd;
		delete del;
		delete pnt_coord;
		all_points.clear();
		open = false;
	}

	void RestrictedPolygonVoronoiDiagram::compute_clipped_VD(std::vector<Plane_3>& clipping_planes)
	{
		if (nb_groups == 0)
			return;
		assert(clipping_planes.size() == nb_groups);
		for (Vertex_iterator vit = rdt_ds.vertices_begin(); vit != rdt_ds.vertices_end(); ++vit)
		{
			vit->vd_vertices.clear();
			vit->contain_non_delaunay_facet = false;
			vit->penetration = true;
			//vit->weights.clear();
		}
		int nb_lost_facets = 0;
		typedef RDT_data_structure::Halfedge_around_vertex_circulator Edge_circulator;
		for (unsigned int i = 0; i < nb_groups; i++)
		{
			const VertGroup& vg = samp_pnts[i];
			const Plane_3& pln = clipping_planes[i];

			// for curvature correction optimization
			Vector_3 pln_nm = pln.orthogonal_vector();
			cgal_vec_normalize(pln_nm);
			for (unsigned int j = 0; j < vg.size(); j++)
			{
				int gid = vg[j]->group_id;
				Edge_circulator start_edge = vg[j]->vertex_begin();
				Edge_circulator current_edge = start_edge;
				do 
				{
					Vertex_handle v_pre = current_edge->prev()->vertex();
					//Vertex_handle v_pre = current_edge->opposite()->vertex();
					if ( v_pre->group_id == vg[j]->group_id)
					{
						vg[j]->penetration = false;
						break;
					}
					++current_edge;
				} while (current_edge != start_edge);
				Edge_circulator end = current_edge;
				do 
				{
					Vertex_handle v_pre = current_edge->prev()->vertex();
					Vertex_handle v_nxt = current_edge->next()->vertex();
					if ( v_pre->group_id != vg[j]->group_id || v_nxt->group_id != vg[j]->group_id)
					{	
						Point_3 c;
						if (current_edge->facet() == Facet_handle())
						{
							//std::cout<<"Caution: invalid surface mesh. poor sampling!.\n";
							nb_lost_facets++;
						}
						else if (current_edge->facet()->is_delaunay)
						{
							c = CGAL::circumcenter(v_pre->mp, v_nxt->mp, vg[j]->mp);
							vg[j]->vd_vertices.push_back(c);

							// for curvature correction optimization
							//Vector_3 tri_nm = CGAL::orthogonal_vector(v_pre->mp, v_nxt->mp, vg[j]->mp);
							//cgal_vec_normalize(tri_nm);
							//double w = std::fabs(pln_nm * tri_nm);
							//vg[j]->weights.push_back(w);
						}
						else
						{
							c = CGAL::centroid(v_pre->mp, v_nxt->mp, vg[j]->mp);
							vg[j]->vd_vertices.push_back(c);
							vg[j]->contain_non_delaunay_facet = true;
							// for curvature corretion optimization
							//Vector_3 tri_nm = CGAL::orthogonal_vector(v_pre->mp, v_nxt->mp, vg[j]->mp);
							//cgal_vec_normalize(tri_nm);
							//double w = std::fabs(pln_nm * tri_nm);
							//vg[j]->weights.push_back(w);
						}				
					}
					++current_edge;
				} while (current_edge != end);
				std::vector<Point_3>& vd_vertices = vg[j]->vd_vertices;
				for (unsigned int k = 0; k <  vd_vertices.size(); k++)
				{
					vd_vertices[k] = pln.projection(vd_vertices[k]);
				}
			}
		}
		std::cout<<"Number of lost facets (due to wrong RDT): "<<nb_lost_facets<<std::endl;
	}
	
	bool RestrictedPolygonVoronoiDiagram::is_delaunay_edge(Halfedge_handle e)
	{
		Vertex_handle vi = e->next()->vertex(), vi_cw = e->vertex(), vi_ccw = e->prev()->vertex();
		Halfedge_handle oe = e->opposite();
		Vertex_handle vj = oe->next()->vertex(), vj_cw = oe->vertex(), vj_ccw = oe->prev()->vertex();
		Point_3 P = vi->mp, Q = vj->mp, R = vi_cw->mp, S = vi_ccw->mp;
		Vector_3 PR(P, R), PS(P, S), QR(Q, R), QS(Q, S);
		PR = PR / CGAL::sqrt(PR.squared_length());
		PS = PS / CGAL::sqrt(PS.squared_length());
		QR = QR / CGAL::sqrt(QR.squared_length());
		QS = QS / CGAL::sqrt(QS.squared_length());
		double cos_ang_P = PR*PS, sin_ang_P = std::fabs(CGAL::sqrt(CGAL::cross_product(PR, PS).squared_length()));
		double cos_ang_Q = QS*QR, sin_ang_Q = std::fabs(CGAL::sqrt(CGAL::cross_product(QS, QR).squared_length()));
		double sin_PQ = sin_ang_P*cos_ang_Q + sin_ang_Q*cos_ang_P;
		if (sin_PQ >= 0.0) // sum of two opposite angles not larger than pi
			return true;
		else
			return false;
	}
	inline void RestrictedPolygonVoronoiDiagram::add_quadrilateral_edge(Halfedge_handle e, std::stack<Halfedge_handle>& q, std::set<Halfedge_handle>& s)
	{
		if (!e->is_border() && !e->opposite()->is_border())
		{
			std::set<Halfedge_handle>::iterator it = s.find(e);
			std::set<Halfedge_handle>::iterator oit = s.find(e->opposite());
			if (it != s.end() && oit != s.end())
			{
				q.push(e);
				q.push(e->opposite());
				s.erase(it);
				s.erase(oit);
			}
		}
	}
	void RestrictedPolygonVoronoiDiagram::iDT_update()
	{
		typedef RDT_data_structure::Halfedge_around_vertex_circulator Edge_circulator;
		int nb_unflippable_edges = 0;
		std::set<Halfedge_handle> visited_edges;
		std::stack<Halfedge_handle> q;
		std::set<Halfedge_handle> added_set;
		typedef std::pair<Vertex_handle, Vertex_handle> Edge;
		std::set<Edge> vpedges;
		for (Halfedge_iterator eit = rdt_ds.halfedges_begin(); eit != rdt_ds.halfedges_end(); ++eit)
		{
			if (!eit->is_border() && !eit->opposite()->is_border())
			{
				std::set<Halfedge_handle>::iterator it = added_set.find(eit);
				std::set<Halfedge_handle>::iterator oit = added_set.find(eit->opposite());
				if (it == added_set.end() && oit == added_set.end())
				{
					q.push(eit);
					q.push(eit->opposite());
					vpedges.insert(Edge(eit->vertex(), eit->opposite()->vertex()));
					vpedges.insert(Edge(eit->opposite()->vertex(), eit->vertex()));
					added_set.insert(eit);
					added_set.insert(eit->opposite());
				}
			}
		}
		for (Facet_iterator fit = rdt_ds.facets_begin(); fit != rdt_ds.facets_end(); ++fit)
		{
			fit->is_delaunay = true;
		}
		while (!q.empty())
		{
			Halfedge_handle e = q.top();
			visited_edges.insert(e);
			q.pop();
			Halfedge_handle oe = q.top();
			visited_edges.insert(oe);
			q.pop();
			if (!is_delaunay_edge(e))
			{
				// check whether flippable
				Vertex_handle v = e->next()->vertex();
				Vertex_handle u = e->opposite()->next()->vertex();
				//Edge_circulator surround_edge = v->vertex_begin();
				//Edge_circulator end_edge = surround_edge;
				//do 
				//{
				//	if (surround_edge->is_border() || surround_edge->opposite()->is_border())
				//		system("pause");
				//	if (surround_edge->opposite()->vertex() == u && surround_edge->vertex() == v
				//		|| surround_edge->vertex() == u && surround_edge->opposite()->vertex() == v)
				//	{
				//		std::cout<<"Unflippable edge!\n";
				//		e->facet()->is_delaunay = oe->facet()->is_delaunay = false;
				//	}
				//	++surround_edge;
				//}while (surround_edge != end_edge);
				if (vpedges.find(Edge(u, v)) != vpedges.end())
				{
					//std::cout<<"Unflippable edge!\n";
					nb_unflippable_edges++;
					e->facet()->is_delaunay = oe->facet()->is_delaunay = false;
				}

				if (!e->facet()->is_delaunay && !oe->facet()->is_delaunay)
					continue;

				/** test validity of pointers after split and joint **/
				Halfedge_handle pre_e2 = oe->prev(), pre_e3 = e->next();
				Halfedge_handle pre_e0 = e->prev(), pre_e1 = oe->next();
				//std::cout<<&(*pre_e0)<<", "<<&(*pre_e1)<<", "<<&(*pre_e2)<<", "<<&(*pre_e3)<<std::endl;
				vpedges.erase(Edge(e->vertex(), e->opposite()->vertex()));
				vpedges.erase(Edge(e->opposite()->vertex(), e->vertex()));
				Halfedge_handle fe = rdt_ds.join_facet(e);
				fe = rdt_ds.split_facet(fe->next(), fe->prev());
				visited_edges.erase(e);
				visited_edges.erase(oe);

				visited_edges.insert(fe);
				visited_edges.insert(fe->opposite());

				vpedges.insert(Edge(fe->vertex(), fe->opposite()->vertex()));
				vpedges.insert(Edge(fe->opposite()->vertex(), fe->vertex()));
				//assert(!fe->is_border() && !fe->opposite()->is_border());
				Halfedge_handle e0 = fe->next(), e1 = fe->prev();
				Halfedge_handle e2 = fe->opposite()->next(), e3 = fe->opposite()->prev();
				//std::cout<<&(*e0)<<", "<<&(*e1)<<", "<<&(*e2)<<", "<<&(*e3)<<std::endl;
				add_quadrilateral_edge(e0, q, visited_edges);
				add_quadrilateral_edge(e1, q, visited_edges);
				add_quadrilateral_edge(e2, q, visited_edges);
				add_quadrilateral_edge(e3, q, visited_edges);
				if (e0 != pre_e0 || e1 != pre_e1 || e2 != pre_e2 || e3 != pre_e3)
					system("pause");
				fe->facet()->is_delaunay = fe->opposite()->facet()->is_delaunay = true;
			}
			e->facet()->is_delaunay = true;
			oe->facet()->is_delaunay = true;
		}
		std::cout<<"Number of unflippable edges: "<<nb_unflippable_edges<<std::endl;
	}

	void RestrictedPolygonVoronoiDiagram::save_triangulation(const std::string fn)
	{
		std::ofstream os(fn.c_str());
		std::vector<Point_3> pnts(rdt_ds.size_of_vertices());
		for (Vertex_iterator vit = rdt_ds.vertices_begin(); vit != vertices_end(); ++vit)
			pnts[vit->idx] = vit->mp;
		for (unsigned int i = 0; i < pnts.size() ; i++)
			os<<"v "<<pnts[i].x()<<' '<<pnts[i].y()<<' '<<pnts[i].z()<<'\n';
		for (Facet_iterator fit = rdt_ds.facets_begin(); fit != rdt_ds.facets_end(); ++fit)
		{
			Halfedge_handle eh = fit->halfedge();
			int i0 = eh->vertex()->idx;
			eh = eh->next();
			int i1 = eh->vertex()->idx;
			eh = eh->next();
			int i2 = eh->vertex()->idx;		
			os<<"f "<<i0+1<<' '<<i1+1<<' '<<i2+1<<'\n';			
		}
	}
}