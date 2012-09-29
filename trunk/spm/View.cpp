#pragma  once

#include "View.h"

namespace Geex
{
	static const GLfloat surf_color_diff[] = {0.95f, 0.95f, 0.0f, 1.0f} ;
	static const GLfloat surf_color_spec[] = {0.8f, 0.8f, 0.8f, 1.0f} ;
	static const GLfloat shininess[] = {75.0f};
	static const GLfloat c1 = 0.35f ;
	static const GLfloat c2 = 0.5f ;
	static const GLfloat c3 = 1.0f ;
	static GLfloat color_table[12][3] = 
	{
		{c3, c2, c2},
		{c2, c3, c2},
		{c2, c2, c3},
		{c2, c3, c3},
		{c3, c2, c3},
		{c3, c3, c2},

		{c1, c2, c2},
		{c2, c1, c2},
		{c2, c2, c1},
		{c2, c1, c1},
		{c1, c2, c1},
		{c1, c1, c2}
	} ;

	static int random_color_index_ = 0 ;
	static void gl_random_color() {
		glColor3f(
			color_table[random_color_index_][0], 
			color_table[random_color_index_][1], 
			color_table[random_color_index_][2]
		) ;
		random_color_index_ = (random_color_index_ + 1) % 12 ;
	}
	static void gl_random_color(int index) {
		random_color_index_ = index % 12 ;
		gl_random_color() ;
	}
	static void gl_randomize_colors() {
		random_color_index_ = 0 ;
	}
	View::View(Doc* doc) : pDoc_(doc), show_domain_mesh_(false), show_polygons_(true), 
						show_bounding_spheres_(true), show_bisectors_(true), show_mesh_skeleton_(false),
						show_regular_triangulation_(false)
	{
	
	}
	void View::draw()
	{
		glUseProgramObjectARB(0);
		gl_randomize_colors();  
		if (show_domain_mesh_)
			draw_domain_mesh();
		if (show_mesh_skeleton_)
			draw_mesh_skeleton(pDoc_->getMeshDomain());
		if (show_polygons_)
			draw_polygons();
		if (show_bounding_spheres_)
			draw_bounding_spheres();
		if (show_bisectors_)
			draw_bisectors();
		if (show_regular_triangulation_)
			draw_triangulation();
		// draw the selected triangles on which we put polygons initially
		//draw_selected_triangles();

	}
	// set up for the mesh drawing
	void View::draw_domain_mesh()
	{
		glCullFace(GL_FRONT) ;
		glEnable(GL_LIGHTING);
		glPushMatrix();
		glMaterialfv(GL_FRONT, GL_DIFFUSE, surf_color_diff);
		glMaterialfv(GL_FRONT, GL_SPECULAR, surf_color_spec);
		glMaterialfv(GL_FRONT, GL_SHININESS, shininess);
		glPopMatrix();
		draw_trimesh(pDoc_->getMeshDomain()) ;
	}
	// draw purely the mesh
	void View::draw_trimesh(const TriMesh& m)
	{
		for(unsigned int i=0; i<m.size(); i++) 
		{
			vec3 n = -m[i].normal();
			glBegin(GL_TRIANGLES);
			glNormal(n);
			glVertex(m[i].vertex[0]) ;
			glVertex(m[i].vertex[1]) ;
			glVertex(m[i].vertex[2]) ;
			glEnd();
		}
	}
	bool View::check_edge_flag(const Facet& F, int e) {
		switch(F.edge_flag[e]) {
			case 1:
				glColor4f(0.2f, 0.2f, 0.2f, 0.5f) ;
				return true ;
			case 2:
				glColor3f(0.5f, 0.5f, 0.5f) ;
				return true ;
			case 3:
				glColor3f(0.1f, 0.1f, 0.9f) ;
				return true ;
		}
		return false ;
	}
	// draw purely the mesh edges
	void View::draw_mesh_skeleton(const TriMesh& m)
	{
		glDisable(GL_LIGHTING);
		glLineWidth(1.0f);
		glBegin(GL_LINES) ;
		for(unsigned int i=0; i<m.size(); i++) {
			const Facet& F = m[i] ;
			if(check_edge_flag(F,2)) { glVertex(F.vertex[0]) ; glVertex(F.vertex[1]) ; }
			if(check_edge_flag(F,0)) { glVertex(F.vertex[1]) ; glVertex(F.vertex[2]) ; }
			if(check_edge_flag(F,1)) { glVertex(F.vertex[2]) ; glVertex(F.vertex[0]) ; }
		}
		glEnd() ;
		glEnable(GL_LIGHTING);
	}
	void View::draw_polygons()
	{
		const vector<vector<Segment_3>>& polygons = pDoc_->getPolygons();
		const vector<vec3> normals = pDoc_->getPolygonNormals();
		glColor3f(0.0f, 0.0f, 0.0f);
		glEnable(GL_LIGHTING);
		for (unsigned int i = 0; i < polygons.size(); i++)
		{
			
			const vector<Segment_3>& ap = polygons[i];
			const vec3& n = normals[i];
			//vec3 center(0.0, 0.0, 0.0);
			glNormal(n);
			glBegin(GL_POLYGON);
			for (unsigned int j = 0; j < ap.size(); j++)
			{
				glVertex3d(ap[j].source().x(), ap[j].source().y(), ap[j].source().z());
				//center += to_geex(ap[j].source());
			}
			glEnd();
			// draw the normal of each polygons, debug usage
			//center /= ap.size();
			//vec3 end = center + n;
			//glDisable(GL_LIGHTING);
			//glColor3f(0.0f, 1.0f, 0.0f);
			//glLineWidth(2.0f);
			//glBegin(GL_LINES);
			//glVertex3d(center.x, center.y, center.z);
			//glVertex3d(end.x, end.y, end.z);
			//glEnd();
			//glLineWidth(1.0f);
			//glEnable(GL_LIGHTING);
		}
		// draw the normal of each facet on the mesh, debug use
		//const TriMesh& m = pDoc_->getMeshDomain();
		//for(unsigned int i=0; i<m.size(); i++) {
		//	vec3 n = -m[i].normal();
		//	glDisable(GL_LIGHTING);
		//	glBegin(GL_LINES);
		//	glLineWidth(4.0f);
		//	glColor3f(1.0f,0.0f,0.0f);
		//	vec3 c = (m[i].vertex[0] + m[i].vertex[1] + m[i].vertex[2])/3.0;
		//	vec3 e = c + n;
		//	glVertex3d(c.x, c.y, c.z);
		//	glVertex3d(e.x, e.y, e.z);
		//	glLineWidth(1.0f);
		//	glEnd();
		//	glEnable(GL_LIGHTING);
		//}
	}
	void View::draw_bounding_spheres()
	{
		const vector<SphericalBounder*>& bounders = pDoc_->getSphericalBounder();
		for (vector<SphericalBounder*>::const_iterator it = bounders.begin(); it!=bounders.end(); it++)
		{
			const vector<Weighted_point>& wps = (*it)->getMarginBounder();
			for (vector<Weighted_point>::const_iterator wit = wps.begin(); wit!=wps.end(); wit++)
			{
				Point_3 center = wit->point();
				double r = std::sqrt(wit->weight());
				glPushMatrix();
				glTranslated(center.x(), center.y(), center.z());
				glutSolidSphere(r, 10, 10);
				glPopMatrix();
			}

		}
	}
	inline void View::draw_segment(const Segment_3 *s)
	{
		glBegin(GL_LINES);
		glVertex3d(s->source().x(), s->source().y(), s->source().z());
		glVertex3d(s->target().x(), s->target().y(), s->target().z());
		glEnd();
	}
	inline void View::draw_ray(const Ray_3 *r)
	{
		Vector_3 emitDir(r->point(0), r->point(1));
		Point_3 src = r->point(0);
		Point_3 tgt = src + 1.0e2*emitDir;
		glBegin(GL_LINES);
		glVertex3d(src.x(), src.y(), src.z());
		glVertex3d(tgt.x(), tgt.y(), tgt.z());
		glEnd();
	}
	void View::draw_bisectors()
	{
		glColor3f(0.0f, 0.0f, 0.0f);
		glDisable(GL_LIGHTING);
		const vector<PowerDiagram_3::Vertex_handle>& verts = pDoc_->getPowerDiagram().getAllVertices();
		for (vector<PowerDiagram_3::Vertex_handle>::const_iterator it = verts.begin(); it!=verts.end(); it++)
		{
			vector<Segment_3>& bisectors = (*it)->tangentIntersection;
			vector<int>& vertGroupId = (*it)->tanIntersecDualVert;
			for (unsigned int i = 0; i < bisectors.size(); i++)
			{
				if (vertGroupId[i] == (*it)->group_id)
					continue;
				draw_segment(&bisectors[i]);
			}
		}
		glEnable(GL_LIGHTING);
	}
	//void View::draw_bisectors()
	//{
	//	glColor3f(0.0f, 0.0f, 0.0f);
	//	glDisable(GL_LIGHTING);
	//	const PowerDiagram_3& pd = pDoc_->getPowerDiagram();
	//	/************************** Original Version: Show all the bisectors ***********************************/
	//	
	//	//const vector<PowerDiagram_3::Vertex_handle>& seeds = pd.getAllVertices();
	//	//for (vector<PowerDiagram_3::Vertex_handle>::const_iterator it = seeds.begin(); it!=seeds.end(); it++)
	//	//{
	//	//	const vector<vector<Object>>& cell = (*it)->vCell;
	//	//	//gl_random_color((*it)->group_id);
	//	//	//draw each face of the cell
	//	//	//for (vector<vector<Object>>::const_iterator f = cell.begin(); f!=cell.end(); f++)
	//	//	for (unsigned int f = 0; f < cell.size(); f++)
	//	//	{
	//	//		if ((*it)->dualVert[f] == (*it)->group_id)
	//	//			continue;
	//	//		for (vector<Object>::const_iterator e = cell[f].begin(); e!=cell[f].end(); e++) 
	//	//		{
	//	//			if (const Segment_3 *s = object_cast<Segment_3>(&(*e)))
	//	//				draw_segment(s);
	//	//			else if (const Ray_3 *r = object_cast<Ray_3>(&(*e)))
	//	//				draw_ray(r);
	//	//		}
	//	//	}
	//	//}
	//	/************************************************************************/
	//	/****************** Only draw those that are bisectors between groups ************************/
	//	const vector<vector<Segment_3>>& bisectors = pd.getBisectors();
	//	unsigned int nbGroups = pd.getNbGroups();
	//	for (unsigned int i = 0; i < nbGroups; i++)
	//	{
	//		const vector<Segment_3>& gbisectors = bisectors[i];
	//		//glBegin(GL_LINES);
	//		gl_random_color(i);
	//		glDisable(GL_LIGHTING);
	//		glBegin(GL_POLYGON);
	//		for (unsigned int j = 0; j < gbisectors.size(); j++)
	//		{
	//			glPoint_3(gbisectors[j].source());
	//			//glPoint_3(gbisectors[j].target());
	//			//if (gbisectors[j].source().x()>=1.0 || gbisectors[j].source().y()>=1.0 || gbisectors[j].source().z()>=1.0)
	//			//{
	//			//	std::cout<<"("<<gbisectors[j].source()<<"), ";
	//			//	std::cout<<"("<<gbisectors[j].target()<<")"<<std::endl;
	//			//}
	//		}
	//		glEnd();
	//		glEnable(GL_LIGHTING);
	//	}
	//	//glEnable(GL_LIGHTING);
	//}

	void View::draw_selected_triangles()
	{
		glDisable(GL_LIGHTING);
		glColor3f(0.0f, 0.0f, 1.0f);
		const vector<int>& sf = pDoc_->getSelectedTriangles();
		const TriMesh& m = pDoc_->getMeshDomain();
		for (unsigned int i = 0; i < sf.size(); i++)
		{
			const Facet& f = m[sf[i]];
			glBegin(GL_TRIANGLES);
			glVertex(f.vertex[0]) ;
			glVertex(f.vertex[1]) ;
			glVertex(f.vertex[2]) ;
			glEnd();
			
		}
		glEnable(GL_LIGHTING);
	}

	void View::draw_triangulation()
	{
		const PowerDiagram_3& pd = pDoc_->getPowerDiagram();
		const vector<PowerDiagram_3::Cell_handle>& allCells = pd.getAllCells();
		glColor3f(1.0f, 0.0f, 0.0f);
		for (vector<PowerDiagram_3::Cell_handle>::const_iterator cit = allCells.begin(); cit != allCells.end(); cit++)
		{
			Point_3 verts[] = { (*cit)->vertex(0)->point(), (*cit)->vertex(1)->point(), (*cit)->vertex(2)->point(), (*cit)->vertex(3)->point()};
			if ((*cit)->vertex(0)->group_id !=0 ) continue;
			glBegin(GL_LINES);
			glPoint_3(verts[0]);	glPoint_3(verts[1]);
			glPoint_3(verts[0]);	glPoint_3(verts[2]);
			glPoint_3(verts[0]);	glPoint_3(verts[3]);
			glPoint_3(verts[1]);	glPoint_3(verts[2]);
			glPoint_3(verts[1]);	glPoint_3(verts[3]);
			glPoint_3(verts[2]);	glPoint_3(verts[3]);
			glEnd();
		}
	}
}