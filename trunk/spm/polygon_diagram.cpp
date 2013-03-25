#include "polygon_diagram.h"

namespace Geex
{
	int nb_invalid_edges = 0;

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
		//edge_face_adjacency.clear();
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
				vec3 n;
				vec3 pp = trimesh->project_to_mesh(to_geex_pnt(p), n);
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
		open = false;
	}

	void RestrictedPolygonVoronoiDiagram::compute_clipped_VD(std::vector<Plane_3>& clipping_planes)
	{
		if (nb_groups == 0)
			return;
		assert(clipping_planes.size() == nb_groups);
		for (Vertex_iterator vit = rdt_ds.vertices_begin(); vit != rdt_ds.vertices_end(); ++vit)
			vit->vd_vertices.clear();
		typedef RDT_data_structure::Halfedge_around_vertex_circulator Edge_circulator;
		for (unsigned int i = 0; i < nb_groups; i++)
		{
			const VertGroup& vg = samp_pnts[i];
			const Plane_3& pln = clipping_planes[i];
			for (unsigned int j = 0; j < vg.size(); j++)
			{
				int gid = vg[j]->group_id;
				Edge_circulator start_edge = vg[j]->vertex_begin();
				Edge_circulator current_edge = start_edge;
				do 
				{
					Vertex_handle v_pre = current_edge->prev()->vertex();
					if ( v_pre->group_id == vg[j]->group_id)
						break;
					++current_edge;
				} while (current_edge != start_edge);
				Edge_circulator end = current_edge;
				do 
				{
					Vertex_handle v_pre = current_edge->prev()->vertex();
					Vertex_handle v_nxt = current_edge->next()->vertex();
					if ( v_pre->group_id != vg[j]->group_id || v_nxt->group_id != vg[j]->group_id)
					{					
						Point_3 c = CGAL::circumcenter(v_pre->mp, v_nxt->mp, vg[j]->mp);
						vg[j]->vd_vertices.push_back(c);
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
	}

	bool RestrictedPolygonVoronoiDiagram::is_delaunay_edge(HalfEdge_handle e)
	{
		Vertex_handle vi = e->next()->vertex(), vi_cw = e->vertex(), vi_ccw = e->prev()->vertex();
		HalfEdge_handle oe = e->opposite();
		Vertex_handle vj = oe->next()->vertex(), vj_cw = oe->vertex(), vj_ccw = oe->prev()->vertex();
		double a2 = CGAL::squared_distance(vi->mp, vi_cw->mp),
			b2 = CGAL::squared_distance(vi->mp, vi_ccw->mp),
			c2 = CGAL::squared_distance(vi_cw->mp, vi_ccw->mp);
		double cos_ang_i = (a2 + b2 - c2) / (2*std::sqrt(a2*b2));
		double sin_ang_i = std::sqrt(1-cos_ang_i*cos_ang_i);
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

	//inline void RestrictedPolygonVoronoiDiagram::add_quadrilateral_edge(HalfEdge_handle e, std::queue<HalfEdge_handle>& q, std::set<HalfEdge_handle>& s)
	inline void RestrictedPolygonVoronoiDiagram::add_quadrilateral_edge(HalfEdge_handle e, std::stack<HalfEdge_handle>& q, std::set<HalfEdge_handle>& s)
	{
		if (!e->is_border() && !e->opposite()->is_border())
		{
			std::set<HalfEdge_handle>::iterator it = s.find(e);
			std::set<HalfEdge_handle>::iterator oit = s.find(e->opposite());
			if (it != s.end() && oit != s.end())
			{
				q.push(e);
				q.push(e->opposite());
				s.erase(it);
				s.erase(oit);
			}
		}
	}
#if 0
	void RestrictedPolygonVoronoiDiagram::add_quadrilateral_edge(std::pair<Vertex_handle, Vertex_handle> e, 
		std::stack<std::pair<Vertex_handle, Vertex_handle>>& q, 
		std::set<std::pair<Vertex_handle, Vertex_handle>>& s)
	{
		typedef std::pair<Vertex_handle, Vertex_handle> Edge;
		using std::make_pair;

		Vertex_handle v = e.first;
		HalfEdge_handle h = v->halfedge();
		HalfEdge_handle oh = h->opposite();
		while ( oh->vertex() != e.second)
		{
			if (h->is_border() || oh->is_border())
				break;
			h = h->next()->opposite();
			oh = h->opposite();
		}
		if (oh->vertex() != e.second)
		{
			h = h->opposite()->prev();
			oh = h->opposite();
			while ( oh->vertex() != e.second)
			{
				h = h->opposite()->prev();
				oh = h->opposite();
			}
		}
		if (!h->is_border() && !oh->is_border())
		{
			std::set<Edge>::iterator it = s.find(e);
			if (it == s.end())
				it = s.find(make_pair(e.second, e.first));
			if (it != s.end())
			{
				q.push(e);
				s.erase(it);
			}
		}
	}
	void RestrictedPolygonVoronoiDiagram::iDT_update()
	{
		typedef std::pair<Vertex_handle, Vertex_handle> Edge;
		using std::make_pair;
		std::stack<Edge> q;
		std::set<Edge> marked_edges;
		std::set<Edge> added_set;
		for (Halfedge_iterator eit = rdt_ds.halfedges_begin(); eit != rdt_ds.halfedges_end(); ++eit)
		{
			if (!eit->is_border() && !eit->opposite()->is_border())
			{
				Vertex_handle v0 = eit->vertex(), v1 = eit->opposite()->vertex();
				std::set<Edge>::iterator it, oit;
				it = added_set.find(make_pair(v0, v1));
				oit = added_set.find(make_pair(v1, v0));
				if (it == added_set.end() && oit == added_set.end())
				{
					q.push(make_pair(v0, v1));
					added_set.insert(make_pair(v0, v1));
				}
			}
		}
		while (!q.empty())
		{
			Edge e = q.top();
			q.pop();
			Vertex_handle v0 = e.first, v1 = e.second;
			HalfEdge_handle eh = v0->halfedge();
			HalfEdge_handle oeh = eh->opposite();
			while ( oeh->vertex() != v1)
			{
				if (eh->is_border() || oeh->is_border())
					break;
				eh = eh->next()->opposite();
				oeh = eh->opposite();
			}
			if (oeh->vertex() != v1)
			{
				eh = eh->opposite()->prev();
				oeh = eh->opposite();
				while ( oeh->vertex() != v1)
				{
					eh = eh->opposite()->prev();
					oeh = eh->opposite();
				}
			}
			Vertex_handle v2 = eh->next()->vertex(), v3 = oeh->next()->vertex();
			marked_edges.insert(e);
			if (!is_delaunay_edge(eh))
			{
				rdt_ds.flip_edge(eh);
				marked_edges.erase(e);
				marked_edges.insert(Edge(v2, v3));
				add_quadrilateral_edge(Edge(v0, v2), q, marked_edges);
				add_quadrilateral_edge(Edge(v0, v3), q, marked_edges);
				add_quadrilateral_edge(Edge(v1, v2), q, marked_edges);
				add_quadrilateral_edge(Edge(v1, v3), q, marked_edges);
			}
		}
	}
#endif

	void RestrictedPolygonVoronoiDiagram::iDT_update()
	{
		std::set<HalfEdge_handle> visited_edges;
		std::stack<HalfEdge_handle> q;
		std::set<HalfEdge_handle> added_set;
		for (Halfedge_iterator eit = rdt_ds.halfedges_begin(); eit != rdt_ds.halfedges_end(); ++eit)
		{
			if (!eit->is_border() && !eit->opposite()->is_border())
			{
				std::set<HalfEdge_handle>::iterator it = added_set.find(eit);
				std::set<HalfEdge_handle>::iterator oit = added_set.find(eit->opposite());
				if (it == added_set.end() && oit == added_set.end())
				{
					q.push(eit);
					q.push(eit->opposite());
					added_set.insert(eit);
					added_set.insert(eit->opposite());
				}
			}
		}	
		while (!q.empty())
		{
			//HalfEdge_handle e = q.front();
			HalfEdge_handle e = q.top();
			visited_edges.insert(e);
			q.pop();
			//HalfEdge_handle oe = q.front();
			HalfEdge_handle oe = q.top();
			q.pop();
			visited_edges.insert(oe);
			//HalfEdge_handle pre_e0 = e->prev(), pre_e1 = e->next();
			//HalfEdge_handle pre_e2 = oe->prev(), pre_e3 = oe->next();
			if (!is_delaunay_edge(e))
			{
				//HalfEdge_handle fe = rdt_ds.flip_edge(e);
				HalfEdge_handle fe = rdt_ds.join_facet(e);
				fe = rdt_ds.split_facet(fe->next(), fe->prev());
				visited_edges.erase(e);
				visited_edges.erase(oe);
				visited_edges.insert(fe);
				visited_edges.insert(fe->opposite());
				assert(!fe->is_border() && !fe->opposite()->is_border());
				HalfEdge_handle e0 = fe->next(), e1 = fe->prev();
				HalfEdge_handle e2 = fe->opposite()->next(), e3 = fe->opposite()->prev();
				//if (!(e0 == pre_e0 && e3 == pre_e1 && e2 == pre_e2 && e1 == pre_e3))
				//	std::cout<<"HalfEdge incosistency!\n";
				add_quadrilateral_edge(e0, q, visited_edges);
				add_quadrilateral_edge(e1, q, visited_edges);
				add_quadrilateral_edge(e2, q, visited_edges);
				add_quadrilateral_edge(e3, q, visited_edges);
			}
		}
	}

	void RestrictedPolygonVoronoiDiagram::save_triangulation(const std::string fn)
	{
		std::ofstream os(fn.c_str());
		std::vector<Point_3> pnts(rdt_ds.size_of_vertices());
		for (Vertex_iterator vit = rdt_ds.vertices_begin(); vit != vertices_end(); ++vit)
			pnts[vit->idx] = vit->point();
		for (unsigned int i = 0; i < pnts.size() ; i++)
			os<<"v "<<pnts[i].x()<<' '<<pnts[i].y()<<' '<<pnts[i].z()<<'\n';
		for (Facet_iterator fit = rdt_ds.facets_begin(); fit != rdt_ds.facets_end(); ++fit)
		{
			HalfEdge_handle eh = fit->halfedge();
			int i0 = eh->vertex()->idx;
			eh = eh->next();
			int i1 = eh->vertex()->idx;
			eh = eh->next();
			int i2 = eh->vertex()->idx;		
			os<<"f "<<i0+1<<' '<<i1+1<<' '<<i2+1<<'\n';			
		}
	}
}