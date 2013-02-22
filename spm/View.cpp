#include "View.h"

namespace Geex
{
	static const GLfloat surf_color_diff[] = {0.95f, 0.95f, 0.0f, 1.0f} ;
	static const GLfloat surf_color_spec[] = {0.8f, 0.8f, 0.8f, 1.0f} ;
	static const GLfloat sphere_diff[] = {0.702f, 0.231f, 0.0f, 0.8f};
	static const GLfloat highlight_diff[] = {1.0f, 0.0f, 0.0f, 1.0f};
	static const GLfloat polygon_diff[] = {0.702f, 0.231f, 0.0f, 0.8f};
	static const GLfloat shininess[] = {75.0f};
	static const GLfloat c1 = 0.35f ;
	static const GLfloat c2 = 0.5f ;
	static const GLfloat c3 = 1.0f ;
	static GLfloat color_table[24][3] = 
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
		{c1, c1, c2},

		{c1, c3, c3},
		{c3, c1, c3},
		{c3, c3, c1},
		{c3, c1, c1},
		{c1, c3, c1},
		{c1, c1, c3},

		{c1, c2, c3}, 
		{c1, c3, c2}, 
		{c2, c1, c3}, 
		{c2, c3, c1}, 
		{c3, c1, c2}, 
		{c3, c2, c1}
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
	View::View(Doc* doc) : pDoc_(doc), show_domain_mesh_(true), show_polygons_(true), 
						show_bounding_spheres_(false), show_bisectors_(true), show_mesh_skeleton_(false),
						show_regular_triangulation_(false), show_holes_(false), show_boundary_(false), 
						show_features_(false), fill_polygon_(false), show_high_curvature_(true), 
						apply_texture_(false), highlighted(-1)
	{
		 Real_generator rg(boost::mt19937(time(0)), boost::random::uniform_real_distribution<>(-0.05, 0.05));
		 random_offset.resize(10);
		 for (unsigned int i = 0; i < random_offset.size(); i++)
			 random_offset[i] = rg();
	}
	View::~View()
	{
		if (texNames.size()>0)
			for (unsigned int i = 0; i < texNames.size(); i++)
				glDeleteTextures(1, &texNames[i]);
		if (bgTexNames.size()>0)
			for (unsigned int i = 0; i < bgTexNames.size(); i++)
				glDeleteTextures(1, &bgTexNames[i]);
		glDeleteTextures(1, &bgTexName);
	}
	void View::draw()
	{
		glUseProgramObjectARB(0);
		//glEnable(GL_DEPTH_TEST);
		//glClear(GL_DEPTH_BUFFER_BIT);
		gl_randomize_colors();
 		if (show_features_)
 			draw_features();
 		if (show_high_curvature_)
 			draw_high_curvature();

 		if (fill_polygon_)
 			fill_draw_polygons();
 		else if (apply_texture_)
 			polygon_draw_with_texture();
 		else
 			draw_polygons();

		if (show_domain_mesh_)
			draw_domain_mesh();		

		if (show_mesh_skeleton_)
			draw_mesh_skeleton(pDoc_->getMeshDomain());
		if (show_bounding_spheres_)
			draw_bounding_spheres();
		if (show_bisectors_)
		{
			draw_bisectors();
			//draw_power_cell();
		}
		if (show_holes_)
		{
			draw_holes();
			draw_undetected_holes();
		}
		if (show_boundary_)
			draw_domain_boundary();
		if (show_regular_triangulation_)
			draw_triangulation();
		//draw_undetected_holes();
		//draw_fitting_points();
		//draw_centroids();
	}
	//void View::draw_centroids()
	//{
	//	vector<vec3>& centroids = pDoc_->pgnCentroid;
	//	glPointSize(6.0f);
	//	glColor3f(1.0f, 0.0f, 0.0f);
	//	glBegin(GL_POINTS);
	//	for (unsigned int i = 0; i < centroids.size(); i++)
	//		glVertex3d(centroids[i].x, centroids[i].y, centroids[i].z);
	//	glEnd();
	//}
	void View::draw_undetected_holes()
	{
		
		//vector<PowerDiagram_3::SkeletonPoint>& hs = pDoc_->getPowerDiagram().holeSkeleton;
		//double scanRadius = pDoc_->getPowerDiagram().scanRadius;
		//glEnable(GL_LIGHTING);
		//for (unsigned int i = 0; i < hs.size(); i++)
		//{
		//	if (hs[i].cluster_id < 0)
		//		continue;
		//	glPushMatrix();
		//	glMaterialfv(GL_FRONT, GL_DIFFUSE, color_table[(hs[i].cluster_id+3)%24]);
		//	glTranslated(hs[i].p.x(), hs[i].p.y(), hs[i].p.z());
		//	glutSolidSphere(hs[i].v->radius/16.0, 10, 10);
		//	//glutSolidSphere(scanRadius, 10, 10);
		//	glPopMatrix();
		//}
		const vector<vector<PowerDiagram_3::Vertex_handle>>& udh = pDoc_->getPowerDiagram().holeBorders;
		for (unsigned int i = 0 ; i < udh.size(); i++)
		{
			const vector<PowerDiagram_3::Vertex_handle>& ahole = udh[i];
			for (unsigned int j = 0; j < ahole.size(); j++)
			{
				glPushMatrix();
				glMaterialfv(GL_FRONT, GL_DIFFUSE, color_table[i%24]);
				glTranslated(ahole[j]->point().x(), ahole[j]->point().y(), ahole[j]->point().z());
				glutSolidSphere(ahole[j]->radius, 10, 10);
				glPopMatrix();
			}
		}
	}
	void View::draw_fitting_points()
	{
		const vector<vector<Point_3>>& fpts = pDoc_->fitPoints();
		glDisable(GL_LIGHTING);
		if (highlighted >= 0)
		{
			glColor3f(1.0f, 0.0f, 0.0f);
			glPointSize(4.0f);
			glBegin(GL_POINTS);
			for (unsigned int i = 0;  i < fpts[highlighted].size(); i++)
				glPoint_3(fpts[highlighted][i]);
			glEnd();
		}
		glEnable(GL_LIGHTING);
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
		draw_trimesh(pDoc_->getMeshDomain()) ;
		glPopMatrix();
	}
	void View::draw_domain_boundary()
	{
		const TriMesh& m = pDoc_->getMeshDomain();
		//TriMesh::BoundaryEdgeSet boundary = m.getBoundary();
		//glLineWidth(4.0);
		//glColor3f(1.0, 0.0, 0.0);
		//for (TriMesh::BoundaryEdgeSet::const_iterator it = boundary.begin(); it != boundary.end(); it++)
		//{
		//	const MeshVertex& v0 = m.vertex(it->first);
		//	const MeshVertex& v1 = m.vertex(it->second);
		//	glBegin(GL_LINES);
		//	glVertex3d(v0.point().x, v0.point().y, v0.point().z);
		//	glVertex3d(v1.point().x, v1.point().y, v1.point().z);
		//	glEnd();
		//}
		//glLineWidth(1.0);
		glEnable(GL_LIGHTING);
		glPushMatrix();
		glMaterialfv(GL_FRONT, GL_DIFFUSE, sphere_diff);
		const std::vector<Weighted_point>& bp = pDoc_->getBoundaryPoints();
		glBegin(GL_POINTS);
		glPointSize(4.0f);
		glColor3f(1.0f, 0.0f, 0.0f);
		for ( std::vector<Weighted_point>::const_iterator it = bp.begin(); it != bp.end(); it++)
		{
			Point_3 center = it->point();
			//glPushMatrix();
			//glTranslated(center.x(), center.y(), center.z());
			//glutSolidSphere(0.001, 10, 10);
			//glPopMatrix();
			glPoint_3(center);
		}
		glEnd();
		glPopMatrix();
	}
	// draw purely the mesh
	void View::draw_trimesh(const TriMesh& m)
	{
		for(unsigned int i=0; i<m.size(); i++) 
		{
			vec3 n = m[i].normal();
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
		//glEnable(GL_LIGHTING);
		glDisable(GL_LIGHTING);
		glColor3f(0.0f, 0.0f, 0.0f);
		//glPushMatrix();
		//glMaterialfv(GL_FRONT, GL_DIFFUSE, sphere_diff);
		const vector<vector<Segment_3>>& polygons = pDoc_->getPolygons();
		const vector<vec3> normals = pDoc_->getPolygonNormals();
		//double offset = 1.0e-3;//offset the polygons out of mesh in order for better view
		glLineWidth(1.5f);
		for (unsigned int i = 0; i < polygons.size(); i++)
		{
			if (i == highlighted)
				continue;
			const vector<Segment_3>& ap = polygons[i];
			const vec3& n = normals[i];
			for (unsigned int j = 0; j < ap.size(); j++)
				draw_segment(&ap[j]);
		}
		//glPopMatrix();
		if ( highlighted > 0 && highlighted < polygons.size())
		{
			glPushMatrix();
			glMaterialfv(GL_FRONT, GL_DIFFUSE, highlight_diff);
			const vector<Segment_3>& ap = polygons[highlighted];
			const vec3& n = normals[highlighted];
			glNormal(n);
			glBegin(GL_POLYGON);
			for (unsigned int j = 0; j < ap.size(); j++)
				draw_segment(&ap[j]);
			glEnd();
			glPopMatrix();
		}
	}
	void View::fill_draw_polygons()
	{
		glCullFace(GL_FRONT) ;
		glEnable(GL_LIGHTING);
		glPushMatrix();
		glMaterialfv(GL_FRONT, GL_DIFFUSE, sphere_diff);
		//glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
		const vector<vector<Segment_3>>& polygons = pDoc_->getPolygons();
		const vector<vec3> normals = pDoc_->getPolygonNormals();
		double offset = 1.0e-3;//offset the polygons out of mesh in order for better view
		for (unsigned int i = 0; i < polygons.size(); i++)
		{
			if (i == highlighted)
				continue;
			const vec3& n = normals[i];
			//glNormal(n);
			glPushMatrix();
			glTranslated(offset*n.x, offset*n.y, offset*n.z);
			triangulate_draw(polygons[i], n);
			glPopMatrix();
		}
		glPopMatrix();
		if ( highlighted > 0 && highlighted < polygons.size())
		{
			glPushMatrix();
			glMaterialfv(GL_FRONT, GL_DIFFUSE, highlight_diff);
			const vector<Segment_3>& ap = polygons[highlighted];
			const vec3& n = normals[highlighted];
			glNormal(n);
			triangulate_draw(ap, n);
			glPopMatrix();
		}
		//glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
	}
	inline void View::triangulate_draw(const vector<Segment_3>& p, const vec3& n)
	{
		// compute center
		Vector_3 vc(0.0, 0.0, 0.0);
		for (unsigned int i = 0; i < p.size(); i++)
			vc = vc + (p[i].source() - CGAL::ORIGIN);
		Point_3 c = CGAL::ORIGIN + vc / p.size();		
		for (unsigned int i = 0;i < p.size(); i++)
		{
			glBegin(GL_TRIANGLES);	
			glNormal(n);
			glPoint_3(p[i].source());
			glNormal(n);
			glPoint_3(c);
			glNormal(n);
			glPoint_3(p[i].target());
			glEnd();
		}
		
	}
	void View::draw_bounding_spheres()
	{
		glEnable(GL_LIGHTING);
		glPushMatrix();
		glMaterialfv(GL_FRONT, GL_DIFFUSE, sphere_diff);
		const vector<PowerDiagram_3::VGroup>& vgs = pDoc_->getPowerDiagram().getGroupList();
		for (unsigned int i = 0; i < vgs.size(); i++)
		{
			if (i == highlighted)
				continue;
			for (vector<PowerDiagram_3::Vertex_handle>::const_iterator wit = vgs[i].begin(); wit!= vgs[i].end(); wit++)
			{
				Point_3 center = (*wit)->point();
				glPushMatrix();
				glTranslated(center.x(), center.y(), center.z());
				glutSolidSphere((*wit)->radius, 10, 10);
				glPopMatrix();
			}
		}
		glPopMatrix();
		if (highlighted >0 && highlighted < vgs.size())
		{
			glPushMatrix();
			glMaterialfv(GL_FRONT, GL_DIFFUSE, highlight_diff);
			for (vector<PowerDiagram_3::Vertex_handle>::const_iterator wit = vgs[highlighted].begin(); wit!= vgs[highlighted].end(); wit++)
			{
				Point_3 center = (*wit)->point();
				glPushMatrix();
				glTranslated(center.x(), center.y(), center.z());
				glutSolidSphere((*wit)->radius, 10, 10);
				glPopMatrix();
			}
			glPopMatrix();
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
		glLineWidth(1.0f);
		const vector<PowerDiagram_3::Vertex_handle>& verts = pDoc_->getPowerDiagram().getAllVertices();
		for (vector<PowerDiagram_3::Vertex_handle>::const_iterator it = verts.begin(); it!=verts.end(); it++)
		{
			vector<Segment_3>& bisectors = (*it)->tangentIntersection;
			//vector<int>& vertGroupId = (*it)->tanIntersecDualVert;
			for (unsigned int i = 0; i < bisectors.size(); i++)
			{
				//if (vertGroupId[i] == (*it)->group_id)
				//	continue;
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
	//	const vector<PowerDiagram_3::Vertex_handle>& seeds = pd.getAllVertices();
	//	int cnt = 0;
	//	for (vector<PowerDiagram_3::Vertex_handle>::const_iterator it = seeds.begin(); it!=seeds.end() && cnt < 2; it++, cnt++)
	//	{
	//		const vector<vector<Object>>& cell = (*it)->vCell;
	//		//gl_random_color((*it)->group_id);
	//		//draw each face of the cell
	//		//for (vector<vector<Object>>::const_iterator f = cell.begin(); f!=cell.end(); f++)
	//		for (unsigned int f = 0; f < cell.size(); f++)
	//		{
	//			if ((*it)->dualVert[f] == (*it)->group_id)
	//				continue;
	//			for (vector<Object>::const_iterator e = cell[f].begin(); e!=cell[f].end(); e++) 
	//			{
	//				if (const Segment_3 *s = object_cast<Segment_3>(&(*e)))
	//					draw_segment(s);
	//				else if (const Ray_3 *r = object_cast<Ray_3>(&(*e)))
	//					draw_ray(r);
	//			}
	//		}
	//	}
	//	/************************************************************************/
	//	/****************** Only draw those that are bisectors between groups ************************/
	//	//const vector<vector<Segment_3>>& bisectors = pd.getBisectors();
	//	//unsigned int nbGroups = pd.getNbGroups();
	//	//for (unsigned int i = 0; i < nbGroups; i++)
	//	//{
	//	//	const vector<Segment_3>& gbisectors = bisectors[i];
	//	//	//glBegin(GL_LINES);
	//	//	gl_random_color(i);
	//	//	glDisable(GL_LIGHTING);
	//	//	glBegin(GL_POLYGON);
	//	//	for (unsigned int j = 0; j < gbisectors.size(); j++)
	//	//	{
	//	//		glPoint_3(gbisectors[j].source());
	//	//		//glPoint_3(gbisectors[j].target());
	//	//		//if (gbisectors[j].source().x()>=1.0 || gbisectors[j].source().y()>=1.0 || gbisectors[j].source().z()>=1.0)
	//	//		//{
	//	//		//	std::cout<<"("<<gbisectors[j].source()<<"), ";
	//	//		//	std::cout<<"("<<gbisectors[j].target()<<")"<<std::endl;
	//	//		//}
	//	//	}
	//	//	glEnd();
	//	//	glEnable(GL_LIGHTING);
	//	//}
	//	glEnable(GL_LIGHTING);
	//}

	//void View::draw_selected_triangles()
	//{
	//	glDisable(GL_LIGHTING);
	//	glColor3f(0.0f, 0.0f, 1.0f);
	//	const vector<int>& sf = pDoc_->getSelectedTriangles();
	//	const TriMesh& m = pDoc_->getMeshDomain();
	//	for (unsigned int i = 0; i < sf.size(); i++)
	//	{
	//		const Facet& f = m[sf[i]];
	//		glBegin(GL_TRIANGLES);
	//		glVertex(f.vertex[0]) ;
	//		glVertex(f.vertex[1]) ;
	//		glVertex(f.vertex[2]) ;
	//		glEnd();
	//		
	//	}
	//	glEnable(GL_LIGHTING);
	//}

	void View::draw_triangulation()
	{
		const PowerDiagram_3& pd = pDoc_->getPowerDiagram();
		glColor3f(1.0f, 0.0f, 0.0f);
		glDisable(GL_LIGHTING);
		//for (PowerDiagram_3::Finite_cells_iterator cit = pd.finite_cells_begin(); cit != pd.finite_cells_end(); cit++)
		//{
		//	Point_3 verts[] = { cit->vertex(0)->point(), cit->vertex(1)->point(),
		//						cit->vertex(2)->point(), cit->vertex(3)->point()};
		//	if (cit->vertex(0)->group_id !=0 ) continue;
		//	glBegin(GL_LINES);
		//	glPoint_3(verts[0]);	glPoint_3(verts[1]);
		//	glPoint_3(verts[0]);	glPoint_3(verts[2]);
		//	glPoint_3(verts[0]);	glPoint_3(verts[3]);
		//	glPoint_3(verts[1]);	glPoint_3(verts[2]);
		//	glPoint_3(verts[1]);	glPoint_3(verts[3]);
		//	glPoint_3(verts[2]);	glPoint_3(verts[3]);
		//	glEnd();
		//}
		//glBegin(GL_LINES);
		//for (PowerDiagram_3::Finite_edges_iterator eit = pd.finite_edges_begin(); eit != pd.finite_edges_end(); eit++)
		//{
		//	PowerDiagram_3::Cell_handle ch = eit->first;
		//	PowerDiagram_3::Vertex_handle vh0 = ch->vertex(eit->second);
		//	PowerDiagram_3::Vertex_handle vh1 = ch->vertex(eit->third);
		//	glPoint_3(vh0->point());
		//	glPoint_3(vh1->point());
		//}
		//glEnd();
		glLineWidth(0.5f);
		glBegin(GL_LINES);
		if (highlighted >= 0)
		{
			if (highlighted >= pd.getNbGroups())
				highlighted = pd.getNbGroups();
			const PowerDiagram_3::VGroup& vg = pd.getGroup(highlighted);
			for (unsigned int i = 0; i < vg.size(); i++)
			{
				if (i != 0)
					continue;
				vector<PowerDiagram_3::Edge> adje;
				pd.finite_incident_edges(vg[i], back_inserter(adje));
				for (unsigned int j = 0; j < adje.size(); j++)
				{
					PowerDiagram_3::Cell_handle ch = adje[j].first;
					PowerDiagram_3::Vertex_handle vh0 = ch->vertex(adje[j].second);
					PowerDiagram_3::Vertex_handle vh1 = ch->vertex(adje[j].third);
					glPoint_3(vh0->point());
					glPoint_3(vh1->point());
				}
			}
		}
		glEnd();
		glEnable(GL_LIGHTING);
	}

	void View::draw_power_cell()
	{
		const PowerDiagram_3& pd = pDoc_->getPowerDiagram();
		//const vector<vector<Segment_3>>& clippedPolygons = pd.getClippedPolygons();
		const vector<vector<Bisector>>& clippedPolygons = pd.getClippedPolygons();
		glDisable(GL_LIGHTING);
		glLineWidth(2.5f);
		for (unsigned int i = 0; i < pd.getNbGroups(); i++)
		{
			gl_random_color(i);
			//glBegin(GL_POLYGON);
			glBegin(GL_LINES);
			for (unsigned int j = 0; j < clippedPolygons[i].size(); j++)
			{
				glPoint_3(clippedPolygons[i][j].source());
				glPoint_3(clippedPolygons[i][(j+1)%clippedPolygons[i].size()].source());
				//draw_segment(&clippedPolygons[i][j].segment());
			}
			glEnd();
		}
		glEnable(GL_LIGHTING);
		//for (unsigned int i = 0; i < pd.getNbGroups(); i++)
		//{
		//	const PowerDiagram_3::VGroup& vg = pd.getGroup(i);
		//	gl_random_color(i);
		//	for (unsigned int j = 0; j < vg.size(); j++)
		//	{
		//		//get the involved bisectors for each vertex
		//		vector<Segment_3>& bisec = vg[j]->tangentIntersection;
		//		vector<int>& dualGroupId = vg[j]->tanIntersecDualVert;			
		//		glBegin(GL_POLYGON);
		//		for (unsigned int k = 0; k < bisec.size(); k++)
		//		{
		//			if (dualGroupId[k] == vg[k]->group_id)
		//				continue;
		//			
		//		}
		//		glEnd();
		//	}
		//}
	}

	void View::draw_holes()
	{
		glDisable(GL_LIGHTING);	
		const vector<vector<Point_3>>& holes = pDoc_->getPowerDiagram().getHoles();
		for ( int i = 0; i < holes.size(); i++)
		{
			gl_random_color(i);
			glBegin(GL_LINE_LOOP);
			const vector<Point_3>& hole = holes[i];
			for (vector<Point_3>::const_iterator it = hole.begin(); it != hole.end(); it++)
			{
				glPoint_3(*it);
			}
			glEnd();
		}
		glEnable(GL_LIGHTING);
		
	}

	void View::draw_features()
	{
		static GLfloat feature_facet_color[] = {0.8, 0.0, 0.0};
		static GLfloat feature_edge_color[] = {0.0, 0.0, 0.0};
		//glEnable(GL_LIGHTING);
		//glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
		//glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, feature_diffuse);
		glDisable(GL_LIGHTING);
		glLineWidth(0.5f);
		const TriMesh& tm = pDoc_->getMeshDomain();
		const vector<pair<int, int>>& fearureEdges = tm.getFeatureFacets();
		for (unsigned int i = 0; i < fearureEdges.size(); i++)
		{
			const Facet& f0 = tm[fearureEdges[i].first];
			const Facet& f1 = tm[fearureEdges[i].second];
			glColor3fv(feature_facet_color);
			glBegin(GL_TRIANGLES);
			glVertex3d(f0.vertex[0].x, f0.vertex[0].y, f0.vertex[0].z);
			glVertex3d(f0.vertex[1].x, f0.vertex[1].y, f0.vertex[1].z);
			glVertex3d(f0.vertex[2].x, f0.vertex[2].y, f0.vertex[2].z);
			glVertex3d(f1.vertex[0].x, f1.vertex[0].y, f1.vertex[0].z);
			glVertex3d(f1.vertex[1].x, f1.vertex[1].y, f1.vertex[1].z);
			glVertex3d(f1.vertex[2].x, f1.vertex[2].y, f1.vertex[2].z);
			glEnd();
			glColor3fv(feature_edge_color);
			glBegin(GL_LINE_LOOP);
			glVertex3d(f0.vertex[0].x, f0.vertex[0].y, f0.vertex[0].z);
			glVertex3d(f0.vertex[1].x, f0.vertex[1].y, f0.vertex[1].z);
			glVertex3d(f0.vertex[2].x, f0.vertex[2].y, f0.vertex[2].z);
			glEnd();
			glBegin(GL_LINE_LOOP);
			glVertex3d(f1.vertex[0].x, f1.vertex[0].y, f1.vertex[0].z);
			glVertex3d(f1.vertex[1].x, f1.vertex[1].y, f1.vertex[1].z);
			glVertex3d(f1.vertex[2].x, f1.vertex[2].y, f1.vertex[2].z);
			glEnd();
		}
		glEnable(GL_LIGHTING);
		glLineWidth(1.0f);
		//glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
	}
	void View::draw_high_curvature()
	{
		static GLfloat curvature_facet_color[] = {0.0, 0.0, 0.8};
		static GLfloat curvature_edge_color[] = {0.0, 0.0, 0.0};
		glDisable(GL_LIGHTING);
		glLineWidth(0.5f);
		const TriMesh& tm = pDoc_->getMeshDomain();
		const set<int>& facetHighCurvature = tm.getHighCurvatureFacets();
		for (set<int>::const_iterator it = facetHighCurvature.begin(); it!=facetHighCurvature.end(); it++)
		{
			const Facet& f = tm[*it];
			glColor3fv(curvature_facet_color);
			glBegin(GL_TRIANGLES);
			glVertex3d(f.vertex[0].x, f.vertex[0].y, f.vertex[0].z);
			glVertex3d(f.vertex[1].x, f.vertex[1].y, f.vertex[1].z);
			glVertex3d(f.vertex[2].x, f.vertex[2].y, f.vertex[2].z);
			glEnd();
			glColor3fv(curvature_edge_color);
			glBegin(GL_LINE_LOOP);
			glVertex3d(f.vertex[0].x, f.vertex[0].y, f.vertex[0].z);
			glVertex3d(f.vertex[1].x, f.vertex[1].y, f.vertex[1].z);
			glVertex3d(f.vertex[2].x, f.vertex[2].y, f.vertex[2].z);
			glEnd();
		}
		glEnable(GL_LIGHTING);
		glLineWidth(1.0f);
	}
	void View::draw_project_facets()
	{
		static GLfloat facet_color[] = {0.0, 0.8, 0.0};
		const vector<int>& pf = pDoc_->getProjectFacets();
		const TriMesh& tm = pDoc_->getMeshDomain();
		glDisable(GL_LIGHTING);
		glColor3fv(facet_color);
		glBegin(GL_TRIANGLES);
		for (unsigned int i = 0; i < pf.size(); i++)
		{
			const Facet &f = tm[pf[i]];
			glVertex3d(f.vertex[0].x, f.vertex[0].y, f.vertex[0].z);
			glVertex3d(f.vertex[1].x, f.vertex[1].y, f.vertex[1].z);
			glVertex3d(f.vertex[2].x, f.vertex[2].y, f.vertex[2].z);
		}
		glEnd();
		glEnable(GL_LIGHTING);
	}

	void View::read_texture(const vector<string>& filenames)
	{
		//glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		texNames.resize(filenames.size());
		for (unsigned int i = 0; i < filenames.size(); i++)
		{
			Mat m = cv::imread(filenames[i]);
			if (!m.data)
			{
				std::cout<<"Fail to load texture image "<<filenames[i]<<std::endl;
				system("pause");
				exit(0);
			}
			glGenTextures(1, &texNames[i]);
			glBindTexture(GL_TEXTURE_2D, texNames[i]);
			int r, c;
			int e = 0;
			while ((1<<e) < m.rows)
				e++;
			r = (1<<e);
			e = 0;
			while ((1<<e) < m.cols)
				e++;
			c = (1<<e);
			//GLubyte *scaleTexels = new GLubyte[r*c*3*sizeof(GLubyte)];
			GLubyte *scaleTexels = (GLubyte*)malloc(r*c*3*sizeof(GLubyte));
			GLubyte *pixels = (GLubyte*)malloc(m.rows*m.cols*3*sizeof(GLubyte));
			GLubyte *ptr = pixels;
			for (int y = 0; y < m.rows; y++)
				for (int x = 0; x < m.cols; x++)
				{
					cv::Vec3b pixel = m.at<Vec3b>(y, x);
					// to rgb model
					*ptr++ = pixel.val[2];
					*ptr++ = pixel.val[1];
					*ptr++ = pixel.val[0];
				}
			if (gluScaleImage(GL_RGB, m.cols, m.rows, GL_UNSIGNED_BYTE, pixels, c, r, GL_UNSIGNED_BYTE, scaleTexels))
			{
				std::cout<<"Scale image failed\n";
				system("pause");
				exit(0);
			}
			//gluBuild2DMipmaps(GL_TEXTURE_2D,GL_RGB,c,r,GL_RGB,GL_UNSIGNED_BYTE,scaleTexels);
			glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
			glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
			glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, c, r, 0, GL_RGB, GL_UNSIGNED_BYTE, scaleTexels);
			free(scaleTexels);
			free(pixels);
		}
	}
	void View::polygon_draw_with_texture()
	{
		//glEnable(GL_DEPTH_TEST);
		glCullFace(GL_FRONT);
		glEnable(GL_CULL_FACE);
		//glEnable(GL_AUTO_NORMAL)
		//glClear(GL_DEPTH_BUFFER_BIT);
		glDisable(GL_LIGHTING);
		glEnable(GL_TEXTURE_2D);
		const vector<int>& pgnLibIdx = pDoc_->getPolygonLibIndex();
		const vector<vector<Segment_3>>& pgns = pDoc_->getPolygons();
		const vector<Polygon_2>& texCoords = pDoc_->getTextCoords();
		const vector<vec3>& pgnNormals = pDoc_->getPolygonNormals();
		const vector<bool>& fbmasks = pDoc_->getFBMasks();
		bool hasBackground = pDoc_->hasBackground();
		//std::cout<<"number of polygons: "<<pgns.size()<<std::endl;
		for (unsigned int i = 0; i < pgns.size(); i++)
		{
			const vector<Segment_3>& p = pgns[i];
			Vector_3 vc(0.0, 0.0, 0.0);
			for (unsigned int j = 0; j < p.size(); j++)
				vc = vc + (p[j].source() - CGAL::ORIGIN);
			Point_3 c = CGAL::ORIGIN + vc / p.size();
			//const Polygon_2& tp = texCoords[pgnLibIdx[i]];
			const Polygon_2& tp = texCoords[i];
			Vector_2 tvc(0.0, 0.0);
			for (unsigned int j = 0; j < tp.size(); j++)
				tvc = tvc + (tp.vertex(j) - CGAL::ORIGIN);
			Point_2 tc = CGAL::ORIGIN + tvc/tp.size();;

			if (hasBackground && fbmasks[i] || !hasBackground)
			{
				glBindTexture(GL_TEXTURE_2D, texNames[pgnLibIdx[i]]);
			}
			else if (hasBackground && !fbmasks[i] && bgTexNames.size() > 0)
				glBindTexture(GL_TEXTURE_2D, bgTexNames[pgnLibIdx[i]]);
			else
				glBindTexture(GL_TEXTURE_2D, bgTexName);

			glBegin(GL_TRIANGLES);
			for (unsigned int j = 0; j < p.size(); j++)
			{			
				glNormal(pgnNormals[i]);
				glTexCoord2f(tp[j].x(), tp[j].y());  
				glPoint_3(p[j].source());
				glNormal(pgnNormals[i]);
				glTexCoord2f(tc.x(), tc.y()); 
				glPoint_3(c);
				glNormal(pgnNormals[i]);
				glTexCoord2f(tp[(j+1)%p.size()].x(), tp[(j+1)%p.size()].y()); 
				glPoint_3(p[j].target());
			}
			glEnd();

		}
		glDisable(GL_TEXTURE_2D);
		glEnable(GL_LIGHTING);
		glDisable(GL_CULL_FACE);	
		//glDisable(GL_DEPTH_TEST);
	}
	void View::read_backgroud_texture(const string& bg_tex_filename)
	{
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		Mat m = cv::imread(bg_tex_filename);
		if (!m.data)
		{
			std::cout<<"Fail to load backgroud texture image "<<bg_tex_filename<<std::endl;
			system("pause");
			exit(0);
		}
		glGenTextures(1, &bgTexName);
		glBindTexture(GL_TEXTURE_2D, bgTexName);
		int r, c;
		int e = 0;
		while ((1<<e) < m.rows)
			e++;
		r = (1<<e);
		e = 0;
		while ((1<<e) < m.cols)
			e++;
		c = (1<<e);
		//GLubyte *scaleTexels = new GLubyte[r*c*3*sizeof(GLubyte)];
		GLubyte *scaleTexels = (GLubyte*)malloc(r*c*3*sizeof(GLubyte));
		GLubyte *pixels = (GLubyte*)malloc(m.rows*m.cols*3*sizeof(GLubyte));
		GLubyte *ptr = pixels;
		for (int y = 0; y < m.rows; y++)
			for (int x = 0; x < m.cols; x++)
			{
				cv::Vec3b pixel = m.at<Vec3b>(y, x);
				// to rgb model
				*ptr++ = pixel.val[2];
				*ptr++ = pixel.val[1];
				*ptr++ = pixel.val[0];
			}
		if (gluScaleImage(GL_RGB, m.cols, m.rows, GL_UNSIGNED_BYTE, pixels, c, r, GL_UNSIGNED_BYTE, scaleTexels))
		{
			std::cout<<"Scale image failed\n";
			system("pause");
			exit(0);
		}
		gluBuild2DMipmaps(GL_TEXTURE_2D,GL_RGB,c,r,GL_RGB,GL_UNSIGNED_BYTE,scaleTexels);
		free(scaleTexels);
		free(pixels);
	}
	void View::read_backgroud_texture(const vector<string>& bg_tex_filenames)
	{
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		bgTexNames.resize(bg_tex_filenames.size());
		for (unsigned int i = 0; i < bg_tex_filenames.size(); i++)
		{
			Mat m = cv::imread(bg_tex_filenames[i]);
			if (!m.data)
			{
				std::cout<<"Fail to load texture image "<<bg_tex_filenames[i]<<std::endl;
				system("pause");
				exit(0);
			}
			glGenTextures(1, &bgTexNames[i]);
			glBindTexture(GL_TEXTURE_2D, bgTexNames[i]);
			int r, c;
			int e = 0;
			while ((1<<e) < m.rows)
				e++;
			r = (1<<e);
			e = 0;
			while ((1<<e) < m.cols)
				e++;
			c = (1<<e);
			//GLubyte *scaleTexels = new GLubyte[r*c*3*sizeof(GLubyte)];
			GLubyte *scaleTexels = (GLubyte*)malloc(r*c*3*sizeof(GLubyte));
			GLubyte *pixels = (GLubyte*)malloc(m.rows*m.cols*3*sizeof(GLubyte));
			GLubyte *ptr = pixels;
			for (int y = 0; y < m.rows; y++)
				for (int x = 0; x < m.cols; x++)
				{
					cv::Vec3b pixel = m.at<Vec3b>(y, x);
					// to rgb model
					*ptr++ = pixel.val[2];
					*ptr++ = pixel.val[1];
					*ptr++ = pixel.val[0];
				}
				if (gluScaleImage(GL_RGB, m.cols, m.rows, GL_UNSIGNED_BYTE, pixels, c, r, GL_UNSIGNED_BYTE, scaleTexels))
				{
					std::cout<<"Scale image failed\n";
					system("pause");
					exit(0);
				}
				//std::cout<<"rows = "<<r<<", columns = "<<c<<", channels = "<<channels<<", depth = "<<depth<<std::endl;
				gluBuild2DMipmaps(GL_TEXTURE_2D,GL_RGB,c,r,GL_RGB,GL_UNSIGNED_BYTE,scaleTexels);
				free(scaleTexels);
				free(pixels);
		}
	}
	void View::draw_fractal_polygons()
	{
		using CGAL::ORIGIN;
		glCullFace(GL_FRONT);
		glEnable(GL_CULL_FACE);
		glDisable(GL_LIGHTING);
		glEnable(GL_TEXTURE_2D);
		const vector<int>& pgnLibIdx = pDoc_->getPolygonLibIndex();
		glPushMatrix();
		glMaterialfv(GL_FRONT, GL_DIFFUSE, sphere_diff);
		//glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
		const vector<vector<Segment_3>>& polygons = pDoc_->getPolygons();
		const vector<vec3> normals = pDoc_->getPolygonNormals();
		const unsigned int noiseNb = random_offset.size();
		for (unsigned int i = 0; i < polygons.size(); i++)
		{
			const vector<Segment_3>& p = polygons[i];
			glNormal(normals[i]);
			glBegin(GL_TRIANGLES);
			Vector_3 c(0.0, 0.0, 0.0);
			for (unsigned int j = 0; j < p.size(); j++)
				c = c + (p[j].source() - CGAL::ORIGIN);
			c = c/p.size();
			Point_3 cent = CGAL::ORIGIN + c;
			for (unsigned int j = 0; j < p.size(); j++)
			{
				Point_3 pre = p[j].source();
				for (unsigned int k = 1; k < noiseNb; k++)
				{
					Point_3 sp = ORIGIN + ( (noiseNb - k)*(p[j].source() - ORIGIN) + k*(p[j].target()- ORIGIN) ) / noiseNb;
					double t = random_offset[(j+k)%noiseNb];
					Point_3 cur = ORIGIN + (t*(cent-ORIGIN) + (1-t)*(sp-ORIGIN));
					// draw this thin triangle
					glPoint_3(cent);
					glPoint_3(pre);
					glPoint_3(cur);
					pre = cur;
				}
				glPoint_3(cent);
				glPoint_3(pre);
				glPoint_3(p[j].target());
			}
			glEnd();
		}
		glPopMatrix();
	}
	void View::thicken_polygons()
	{

	}
}
