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
		which_group.reserve(which_group.size() + samp_nb + 1);
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
				//samp_pnts.push_back(p);
				samp_pnts.push_back(new Embedded_vertex_2<Embedded_point_2>(p));
				which_group.push_back(group_id);
			}
		}
		
	}

	void RestrictedPolygonVoronoiDiagram::begin_insert()
	{
		assert(!open);
		open = true;
		samp_pnts.clear();
		which_group.clear();
		if (!pnt_coord) delete pnt_coord;
	}
	
	void RestrictedPolygonVoronoiDiagram::end_insert()
	{
		std::cout<<"Number of sample points: "<<samp_pnts.size()<<std::endl;
		pnt_coord = new double[samp_pnts.size()*3];
		double *ptr = pnt_coord;
		for ( unsigned int i = 0; i < samp_pnts.size(); i++ )
		{
			*ptr++ = samp_pnts[i]->point().x();
			*ptr++ = samp_pnts[i]->point().y();
			*ptr++ = samp_pnts[i]->point().z();
		}
		del->set_vertices(samp_pnts.size(), pnt_coord);
		open = false;
	}
	//void RestrictedPolygonVoronoiDiagram::uniform_sample(const Polygon_3& polygon, unsigned int samp_nb)
	//{
	//	samp_pnts.reserve(samp_pnts.size() + samp_nb);
	//	double perimeter = 0.0;
	//	std::vector<double> edge_len(polygon.size());
	//	for ( unsigned int i = 0; i < polygon.size(); i++ )
	//	{
	//		Segment_3 e = polygon.edge(i);
	//		edge_len[i] = CGAL::sqrt(e.squared_length());
	//		perimeter += edge_len[i];
	//	}

	//	for ( unsigned int i = 0; i < polygon.size(); i++ )
	//	{
	//		Segment_3 e = polygon.edge(i);
	//		unsigned int n = edge_len[i] / perimeter * samp_nb;
	//		Point_3 src = e.source(), tgt = e.target();
	//		for ( unsigned int j = 0; j < n+1; j++ )
	//		{
	//			Point_3 p = CGAL::ORIGIN + ( (n+1-j)*(src - CGAL::ORIGIN) + j*(tgt - CGAL::ORIGIN) )/(n+1);
	//			samp_pnts.push_back(p);
	//		}
	//	}
	//}

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
				rvd->for_each_primal_triangle(draw_primal_triangles(pnt_coord));
				glEnable(GL_LIGHTING);
				glEndList();
			}
			glCallList(RDT_disp_list);
		}

	}
}