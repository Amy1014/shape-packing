#include "polygon_diagram.h"

namespace Geex
{

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
		//mesh->clear();
		delete mesh;
		std::for_each(samp_pnts.begin(), samp_pnts.end(), std::mem_fun_ref(&VertGroup::clear));
		samp_pnts.clear();
		bounding_pnts.clear();
		rdt_ds.clear();
	}

	void RestrictedPolygonVoronoiDiagram::set_mesh(const std::string& mesh_file_name)
	{
		//assert(!mesh);
		delete mesh;
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
				assert(!Numeric::is_nan(p.x()) && !Numeric::is_nan(p.y()) && !Numeric::is_nan(p.z()));
				//all_points.push_back(p);
				//p->group_id = group_id;
				p = e.supporting_line().projection(p);
				assert(!Numeric::is_nan(p.x()) && !Numeric::is_nan(p.y()) && !Numeric::is_nan(p.z()));
				vec3 dv;
				vec3 pp = trimesh->project_to_mesh(to_geex_pnt(p), dv);
				assert(!Numeric::is_nan(pp.x) && !Numeric::is_nan(pp.y) && !Numeric::is_nan(pp.z));
				all_points.push_back(MyPoint(p, to_cgal_pnt(pp), group_id));
				//all_points.back().mp = to_cgal_pnt(pp);
				nb_pnts++;
			}
		}
		return nb_pnts;
	}

	unsigned int RestrictedPolygonVoronoiDiagram::insert_polygons(const Polygon_3& polygon, unsigned int group_id)
	{
		unsigned int nb_pnts = 0; 
		for ( unsigned int i = 0; i < polygon.size(); i++ )
		{
			Point_3 p = polygon.vertex(i);
			vec3 dv;
			vec3 pp = trimesh->project_to_mesh(to_geex_pnt(p), dv);
			all_points.push_back(MyPoint(p, to_cgal_pnt(pp), group_id));
			nb_pnts++;
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

	void RestrictedPolygonVoronoiDiagram::compute_clipped_VD(std::vector<Plane_3>& clipping_planes, std::vector<Point_3>& ref_pnts)
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
		smoothed_VD_regions.resize(nb_groups);
		typedef RDT_data_structure::Halfedge_around_vertex_circulator Edge_circulator;
		for (unsigned int i = 0; i < nb_groups; i++)
		{
			const VertGroup& vg = samp_pnts[i];
			const Plane_3& pln = clipping_planes[i];
			smoothed_VD_regions[i].clear();
			// for curvature correction optimization
			for (unsigned int j = 0; j < vg.size(); j++)
			{
				std::vector<bool> is_triple_pnt;
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
							nb_lost_facets++;
						else if (current_edge->facet()->is_delaunay)
						{
							Point_3 c = CGAL::circumcenter(v_pre->mp, v_nxt->mp, vg[j]->mp);
							if (v_nxt->group_id < 0 && v_pre->group_id < 0)
								vg[j]->vd_vertices.push_back( v_nxt->mp );
							else if (v_nxt->group_id < 0 && v_pre->group_id >= 0)
							{
								vg[j]->vd_vertices.push_back(v_nxt->mp);
							}
							else if (v_nxt->group_id >= 0 && v_pre->group_id < 0)
							{
								if (vg[j]->group_id != v_nxt->group_id)
									vg[j]->vd_vertices.push_back(c);
							}
							else
								vg[j]->vd_vertices.push_back( c );
							if (v_pre->group_id != vg[j]->group_id && v_nxt->group_id != vg[j]->group_id && v_pre->group_id != v_nxt->group_id
								|| v_pre->group_id < 0 && v_nxt->group_id < 0 )
								is_triple_pnt.push_back(true);
							else
								is_triple_pnt.push_back(false);
						}
						else
						{
							Point_3 c = CGAL::centroid(v_pre->mp, v_nxt->mp, vg[j]->mp);
							if (v_nxt->group_id < 0 && v_pre->group_id < 0)
								vg[j]->vd_vertices.push_back( v_nxt->mp );
							else if (v_nxt->group_id < 0 && v_pre->group_id >= 0)
							{
								vg[j]->vd_vertices.push_back(v_nxt->mp);
							}
							else if (v_nxt->group_id >= 0 && v_pre->group_id < 0)
							{
								if (vg[j]->group_id != v_nxt->group_id)
									vg[j]->vd_vertices.push_back(c);								
							}
							else
								vg[j]->vd_vertices.push_back( c );
							vg[j]->contain_non_delaunay_facet = true;
							if (v_pre->group_id != vg[j]->group_id && v_nxt->group_id != vg[j]->group_id && v_pre->group_id != v_nxt->group_id
								|| v_pre->group_id < 0 && v_nxt->group_id < 0)
								is_triple_pnt.push_back(true);
							else
								is_triple_pnt.push_back(false);
						}				
					}
					++current_edge;
				} while (current_edge != end);
				std::vector<Point_3>& vd_vertices = vg[j]->vd_vertices;
// #ifdef _DEMO_
 				vg[j]->no_prj_vd_vertices = vd_vertices;
// #endif
				//std::vector<bool>& is_triple_pnt = vg[j]->is_triple_pnt;
				for (unsigned int k = 0; k < vd_vertices.size(); k++)
				{
					vd_vertices[k] = pln.projection(vd_vertices[k]);
					if (is_triple_pnt[k])
						smoothed_VD_regions[i].push_back(vd_vertices[k]);
				}
			}
			std::vector<Point_2> conhull;
			std::vector<Point_2> pnt2d;
			pnt2d.reserve(smoothed_VD_regions[i].size());
			Vector_3 base1 = pln.base1(), base2 = pln.base2();
			cgal_vec_normalize(base1);
			cgal_vec_normalize(base2);
			//Point_3 o = pln.point();
			Point_3 o = ref_pnts[i];
			for (unsigned int k = 0; k < smoothed_VD_regions[i].size(); k++)
			{
				Vector_3 v(o, smoothed_VD_regions[i][k]);
				pnt2d.push_back(Point_2(v*base1, v*base2));
			}
			CGAL::ch_graham_andrew(pnt2d.begin(), pnt2d.end(), std::back_insert_iterator<std::vector<Point_2>>(conhull));
			smoothed_VD_regions[i].clear();
			for (unsigned int k = 0; k < conhull.size(); k++)
				smoothed_VD_regions[i].push_back(o + ( conhull[k].x()*base1 + conhull[k].y()*base2 ) );
		}
		std::cout<<"Number of lost facets (due to wrong RDT): "<<nb_lost_facets<<std::endl;
	}
	void RestrictedPolygonVoronoiDiagram::compute_mixed_VD(std::vector<Plane_3>& clipping_planes, std::vector<Point_3>& ref_pnts, std::vector<bool>& use_approx)
	{
		if (nb_groups == 0)
			return;
		assert(clipping_planes.size() == nb_groups);
		for (Vertex_iterator vit = rdt_ds.vertices_begin(); vit != rdt_ds.vertices_end(); ++vit)
		{
			vit->vd_vertices.clear();
			vit->contain_non_delaunay_facet = false;
			vit->penetration = true;
		}
		int nb_lost_facets = 0;
		smoothed_VD_regions.resize(nb_groups);
		typedef RDT_data_structure::Halfedge_around_vertex_circulator Edge_circulator;
		for (unsigned int i = 0; i < nb_groups; i++)
		{
			const VertGroup& vg = samp_pnts[i];
			const Plane_3& pln = clipping_planes[i];
			smoothed_VD_regions[i].clear();

			for (unsigned int j = 0; j < vg.size(); j++)
			{
				std::vector<bool> is_triple_pnt;
				Edge_circulator start_edge = vg[j]->vertex_begin();
				Edge_circulator current_edge = start_edge;
				do 
				{
					Vertex_handle v_pre = current_edge->prev()->vertex();
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
							nb_lost_facets++;
						else if (use_approx[i] || !current_edge->facet()->is_delaunay)
						{
							Point_3 c = CGAL::centroid(v_pre->mp, v_nxt->mp, vg[j]->mp);
							if (v_nxt->group_id < 0 && v_pre->group_id < 0)
								vg[j]->vd_vertices.push_back( v_nxt->mp );
							else if (v_nxt->group_id < 0 && v_pre->group_id >= 0)
							{
								vg[j]->vd_vertices.push_back(v_nxt->mp);
							}
							else if (v_nxt->group_id >= 0 && v_pre->group_id < 0)
							{
								if (vg[j]->group_id != v_nxt->group_id)
									vg[j]->vd_vertices.push_back(c);								
							}
							else
								vg[j]->vd_vertices.push_back( c );
							vg[j]->contain_non_delaunay_facet = true;
							if (v_pre->group_id != vg[j]->group_id && v_nxt->group_id != vg[j]->group_id && v_pre->group_id != v_nxt->group_id
								|| v_pre->group_id < 0 && v_nxt->group_id < 0)
								is_triple_pnt.push_back(true);
							else
								is_triple_pnt.push_back(false);
						}	
						else
						{
							Point_3 c = CGAL::circumcenter(v_pre->mp, v_nxt->mp, vg[j]->mp);
							if (v_nxt->group_id < 0 && v_pre->group_id < 0)
								vg[j]->vd_vertices.push_back( v_nxt->mp );
							else if (v_nxt->group_id < 0 && v_pre->group_id >= 0)
							{
								vg[j]->vd_vertices.push_back(v_nxt->mp);
							}
							else if (v_nxt->group_id >= 0 && v_pre->group_id < 0)
							{
								if (vg[j]->group_id != v_nxt->group_id)
									vg[j]->vd_vertices.push_back(c);
							}
							else
								vg[j]->vd_vertices.push_back( c );
							if (v_pre->group_id != vg[j]->group_id && v_nxt->group_id != vg[j]->group_id && v_pre->group_id != v_nxt->group_id
								|| v_pre->group_id < 0 && v_nxt->group_id < 0 )
								is_triple_pnt.push_back(true);
							else
								is_triple_pnt.push_back(false);
						}
					}
					++current_edge;
				} while (current_edge != end);
				std::vector<Point_3>& vd_vertices = vg[j]->vd_vertices;
				for (unsigned int k = 0; k < vd_vertices.size(); k++)
				{
					vd_vertices[k] = pln.projection(vd_vertices[k]);
					if (is_triple_pnt[k])
						smoothed_VD_regions[i].push_back(vd_vertices[k]);
				}
			}
			std::vector<Point_2> conhull;
			std::vector<Point_2> pnt2d;
			pnt2d.reserve(smoothed_VD_regions[i].size());
			Vector_3 base1 = pln.base1(), base2 = pln.base2();
			cgal_vec_normalize(base1);
			cgal_vec_normalize(base2);
			Point_3 o = ref_pnts[i];
			for (unsigned int k = 0; k < smoothed_VD_regions[i].size(); k++)
			{
				Vector_3 v(o, smoothed_VD_regions[i][k]);
				pnt2d.push_back(Point_2(v*base1, v*base2));
			}
			CGAL::ch_graham_andrew(pnt2d.begin(), pnt2d.end(), std::back_insert_iterator<std::vector<Point_2>>(conhull));
			smoothed_VD_regions[i].clear();
			smoothed_VD_regions[i].reserve(conhull.size());
			for (unsigned int k = 0; k < conhull.size(); k++)
				smoothed_VD_regions[i].push_back(o + ( conhull[k].x()*base1 + conhull[k].y()*base2 ) );
		}
		std::cout<<"Number of lost facets (due to wrong RDT): "<<nb_lost_facets<<std::endl;
	}

	void RestrictedPolygonVoronoiDiagram::compute_midpoint_VD(std::vector<Plane_3>& clipping_planes, std::vector<Point_3>& ref_pnts)
	{
		if (nb_groups == 0)
			return;
		assert(clipping_planes.size() == nb_groups);
		for (Vertex_iterator vit = rdt_ds.vertices_begin(); vit != rdt_ds.vertices_end(); ++vit)
		{
			vit->med_segs.clear();
			vit->contain_non_delaunay_facet = false;
			vit->penetration = true;
		}
		int nb_lost_facets = 0;
		smoothed_VD_regions.resize(nb_groups);
		typedef RDT_data_structure::Halfedge_around_vertex_circulator Edge_circulator;
		for (unsigned int i = 0; i < nb_groups; i++)
		{
			const VertGroup& vg = samp_pnts[i];
			const Plane_3& pln = clipping_planes[i];
			smoothed_VD_regions[i].clear();
			for (unsigned int j = 0; j < vg.size(); j++)
			{
				Edge_circulator start_edge = vg[j]->vertex_begin();
				Edge_circulator current_edge = start_edge;
				do 
				{						
					//Facet_handle f = current_edge->facet();
					//if (f == Facet_handle())
					//	nb_lost_facets++;
					//else 
					//{
						Vertex_handle v_nxt = current_edge->next()->vertex();
						Vertex_handle v_pre = current_edge->prev()->vertex();
						if (v_nxt->group_id != vg[j]->group_id || v_pre->group_id != vg[j]->group_id)
						{
							if (current_edge->facet() == Facet_handle())
								nb_lost_facets++;
							else if (v_pre->group_id < 0 && v_nxt->group_id < 0)
							{
								vg[j]->med_segs.push_back(Segment_3(v_pre->mp, v_nxt->mp));
								smoothed_VD_regions[i].push_back(v_pre->mp);
							}
							else if (v_pre->group_id < 0 && v_nxt->group_id != /*==*/ vg[j]->group_id )
							{
								vg[j]->med_segs.push_back(Segment_3(v_pre->mp, CGAL::midpoint(v_nxt->mp, vg[j]->mp)));
							}
							else if (v_nxt->group_id < 0 && v_pre->group_id != /*==*/ vg[j]->group_id)
							{
								vg[j]->med_segs.push_back(Segment_3(v_nxt->mp, CGAL::midpoint(vg[j]->mp, v_pre->mp)));
							}
							else if (v_pre->group_id < 0 && v_nxt->group_id == vg[j]->group_id )
							{
								//vg[j]->med_segs.push_back(Segment_3(v_pre->mp, CGAL::midpoint(v_nxt->mp, vg[j]->mp)));
							}
							else if (v_nxt->group_id < 0 && v_pre->group_id == vg[j]->group_id)
							{
								//vg[j]->med_segs.push_back(Segment_3(v_nxt->mp, CGAL::midpoint(vg[j]->mp, v_pre->mp)));
							}
							else if (v_pre->group_id == vg[j]->group_id && v_nxt->group_id != vg[j]->group_id)
							{
								Point_3 end0 = CGAL::midpoint(vg[j]->mp, v_nxt->mp), end1 = CGAL::midpoint(v_pre->mp, v_nxt->mp);
								vg[j]->med_segs.push_back(Segment_3(end0, end1));
							}
							else if (v_pre->group_id != vg[j]->group_id && v_nxt->group_id == vg[j]->group_id)
							{
								Point_3 end0 = CGAL::midpoint(v_nxt->mp, v_pre->mp), end1 = CGAL::midpoint(vg[j]->mp, v_pre->mp);
								vg[j]->med_segs.push_back(Segment_3(end0, end1));
							}
							else if ( v_pre->group_id != vg[j]->group_id && v_nxt->group_id != vg[j]->group_id && v_pre->group_id == v_nxt->group_id )
							{
								Point_3 end0 = CGAL::midpoint(vg[j]->mp, v_nxt->mp), end1 = CGAL::midpoint(v_pre->mp, vg[j]->mp);
								vg[j]->med_segs.push_back(Segment_3(end0, end1));
							}
							else if ( v_pre->group_id != vg[j]->group_id && v_nxt->group_id != vg[j]->group_id && v_pre->group_id != v_nxt->group_id )
							{
								Point_3 end0 = CGAL::midpoint(vg[j]->mp, v_nxt->mp), end1 = CGAL::midpoint(v_pre->mp, vg[j]->mp);
								Point_3 c = CGAL::centroid(vg[j]->mp, v_pre->mp, v_nxt->mp);
								vg[j]->med_segs.push_back(Segment_3(end0, c));
								vg[j]->med_segs.push_back(Segment_3(c, end1));
								smoothed_VD_regions[i].push_back(c);
							}
						}
						
					//}	
					--current_edge;
				} while (current_edge != start_edge);

				std::vector<Segment_3>& med_segs = vg[j]->med_segs;
// #ifdef _DEMO_
 				vg[j]->no_prj_med_segs = med_segs;
// #endif
				for (unsigned int k = 0; k < med_segs.size(); k++)
				{
					Point_3 src = pln.projection(med_segs[k].source());
					Point_3 tgt = pln.projection(med_segs[k].target());
					med_segs[k] = Segment_3(src, tgt);
				}
			}
			for (unsigned int k = 0; k < smoothed_VD_regions[i].size(); k++)
				smoothed_VD_regions[i][k] = pln.projection(smoothed_VD_regions[i][k]);
			std::vector<Point_2> conhull, pnt2d;
			pnt2d.reserve(smoothed_VD_regions[i].size());
			Vector_3 base1 = pln.base1(), base2 = pln.base2();
			cgal_vec_normalize(base1);
			cgal_vec_normalize(base2);
			Point_3 o = ref_pnts[i];
			for (unsigned int k = 0; k < smoothed_VD_regions[i].size(); k++)
			{
				Vector_3 v(o, smoothed_VD_regions[i][k]);
				pnt2d.push_back(Point_2(v*base1, v*base2));
			}
			CGAL::ch_graham_andrew(pnt2d.begin(), pnt2d.end(), std::back_insert_iterator<std::vector<Point_2>>(conhull));
			smoothed_VD_regions[i].clear();
			for (unsigned int k = 0; k < conhull.size(); k++)
				smoothed_VD_regions[i].push_back(o + ( conhull[k].x()*base1 + conhull[k].y()*base2 ) );
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

				if (e->vertex()->group_id == e->opposite()->vertex()->group_id)
					continue;
					//e->facet()->is_delaunay = e->opposite()->facet()->is_delaunay = false;
				if (vpedges.find(Edge(u, v)) != vpedges.end())
				{
					//std::cout<<"Unflippable edge!\n";
					nb_unflippable_edges++;
					//e->facet()->is_delaunay = oe->facet()->is_delaunay = false;
					continue;
				}

				//if (!e->facet()->is_delaunay && !oe->facet()->is_delaunay)
				//	continue;

				/** test validity of pointers after split and joint **/
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

				//fe->facet()->is_delaunay = fe->opposite()->facet()->is_delaunay = true;
			}
			else
			{
				//e->facet()->is_delaunay = true;
				//oe->facet()->is_delaunay = true;
			}

		}
		std::cout<<"Number of unflippable edges: "<<nb_unflippable_edges<<std::endl;
		// found out non-delaunay facets
		for (Facet_iterator fit = rdt_ds.facets_begin(); fit != rdt_ds.facets_end(); ++fit)
		{
			Halfedge_handle e = fit->halfedge();
			if (!is_delaunay_edge(e))
			{
				fit->is_delaunay = false;
				continue;
			}
			e = e->next();
			if (!is_delaunay_edge(e))
			{
				fit->is_delaunay = false;
				continue;
			}
			e = e->next();
			if (!is_delaunay_edge(e))
			{
				fit->is_delaunay = false;
				continue;
			}
			fit->is_delaunay = true;
		}
		// update facet normal
		for (Facet_iterator fit = rdt_ds.facets_begin(); fit != rdt_ds.facets_end(); ++fit)
		{
			RestrictedPolygonVoronoiDiagram::Halfedge_handle eh = fit->halfedge();
			Point_3 v0 = eh->vertex()->mp;
			eh = eh->next();
			Point_3 v1 = eh->vertex()->mp;
			eh = eh->next();
			Point_3 v2 = eh->vertex()->mp;
			Vector_3 v01(v0, v1), v12(v1, v2);
			cgal_vec_normalize(v01);
			cgal_vec_normalize(v12);
			fit->n = CGAL::cross_product(v01, v12);
			cgal_vec_normalize(fit->n);
		}
	}

	void RestrictedPolygonVoronoiDiagram::save_triangulation(const std::string fn) const
	{
		std::ofstream os(fn.c_str());
		std::vector<Point_3> pnts(rdt_ds.size_of_vertices());
		for (Vertex_const_iterator vit = rdt_ds.vertices_begin(); vit != vertices_end(); ++vit)
			pnts[vit->idx] = vit->mp;
		for (unsigned int i = 0; i < pnts.size() ; i++)
			os<<"v "<<pnts[i].x()<<' '<<pnts[i].y()<<' '<<pnts[i].z()<<'\n';
		for (Facet_const_iterator fit = rdt_ds.facets_begin(); fit != rdt_ds.facets_end(); ++fit)
		{
			Halfedge_const_handle eh = fit->halfedge();
			int i0 = eh->vertex()->idx;
			eh = eh->next();
			int i1 = eh->vertex()->idx;
			eh = eh->next();
			int i2 = eh->vertex()->idx;		
			os<<"f "<<i0+1<<' '<<i1+1<<' '<<i2+1<<'\n';			
		}
	}

}