#include "polygon_diagram.h"

namespace Geex
{
	
	RestrictedPolygonVoronoiDiagram::RestrictedPolygonVoronoiDiagram() 
	{
		mesh = 0;
		rvd = 0;
		pnt_coord = 0;
		open = false;
		del = Delaunay::create("CGAL");
	}

	RestrictedPolygonVoronoiDiagram::~RestrictedPolygonVoronoiDiagram()
	{
		delete rvd;
		delete del;
		delete pnt_coord;
		glDeleteLists(RDT_disp_list, 1);
	}

	void RestrictedPolygonVoronoiDiagram::set_mesh(const std::string& mesh_file_name)
	{
		assert(!mesh);
		mesh = new TopoPolyMesh;
		mesh->load(mesh_file_name);
		rvd = new RestrictedVoronoiDiagram(del, mesh);
	}
	
	void RestrictedPolygonVoronoiDiagram::insert_polygons(const Polygon_3& polygon, unsigned int group_id, unsigned int samp_nb)
	{
		//assert(open);
		samp_pnts.reserve(samp_pnts.size() + samp_nb + 1);
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
				RDT_data_structure::Vertex_handle vh = rdt_ds.create_vertex(RDT_data_structure::Vertex(p, -1));
				samp_pnts.push_back(vh);
			}
		}
		
	}

	void RestrictedPolygonVoronoiDiagram::begin_insert()
	{
		assert(!open);
		open = true;
		samp_pnts.clear();
		if (!pnt_coord) delete pnt_coord;
	}
	
	void RestrictedPolygonVoronoiDiagram::end_insert()
	{
		std::cout<<"Number of sample points: "<<samp_pnts.size()<<std::endl;
		pnt_coord = new double[samp_pnts.size()*3];
		double *ptr = pnt_coord;
		for ( unsigned int i = 0; i < samp_pnts.size(); i++ )
		{
			Point_3 p = samp_pnts[i]->point_3();
			*ptr++ = p.x();
			*ptr++ = p.y();
			*ptr++ = p.z();
		}
		del->set_vertices(samp_pnts.size(), pnt_coord);
		rdt_ds.set_dimension(2);
		rvd->for_each_primal_triangle(Construct_RDT_structure(rdt_ds, samp_pnts, edge_face_adjacency));
		open = false;
	}
	
	void RestrictedPolygonVoronoiDiagram::draw_DT()
	{
		if (mesh && rvd)
		{
			if (!glIsList(RDT_disp_list))
			{
				RDT_disp_list = glGenLists(1);
				glNewList(RDT_disp_list, GL_COMPILE);
				glLineWidth(1.5f);
				glColor3f(0.8f, 0.0f, 0.0f);
				glDisable(GL_LIGHTING);
				glBegin(GL_LINES);
				for (RDT_data_structure::Edge_iterator eit = rdt_ds.edges_begin(); eit != rdt_ds.edges_end(); eit++)
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
		find_adjacency(samp_pnts[i], samp_pnts[j], f->index(samp_pnts[k]), f);
		find_adjacency(samp_pnts[j], samp_pnts[k], f->index(samp_pnts[i]), f);
		find_adjacency(samp_pnts[k], samp_pnts[i], f->index(samp_pnts[j]), f);
	}

}