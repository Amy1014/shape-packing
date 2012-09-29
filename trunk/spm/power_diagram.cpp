#include "power_diagram.h"

namespace Geex
{
	PowerDiagram_3::PowerDiagram_3(TriMesh* m): bdMesh(m), open(false), threshold(1.0e-10){}

	void PowerDiagram_3::clear()
	{
		baseclass::clear();
		all_vertices_.clear();
		clipped_tangent_plane_.clear();
	}
	void PowerDiagram_3::begin_insert()	
	{
		open = true;
		groupCount = 0;
		vertIdx = 0;
		clear();
	}

	PowerDiagram_3::Vertex_handle PowerDiagram_3::insert(const Weighted_point& wp)
	{
		Point_3 c = wp.point();
		double sw = wp.weight();
		if (is_nan(c.x()) || is_nan(c.y()) || is_nan(c.z()) || is_nan(sw))
		{
			std::cerr << "Nan !" << std::endl ;
			return 0 ;
		}
		return baseclass::insert(wp);
	}
	void PowerDiagram_3::insert_grouped_points( const vector<Weighted_point>& vlist)
	{
		gx_assert(open);
		//VGroup vg;
		for (unsigned int i = 0; i < vlist.size(); i++)
		{
			Vertex_handle v = insert(vlist[i]);
			v->index = vertIdx;
			v->group_id = groupCount;
			v->radius = CGAL::sqrt(vlist[i].weight());
			//vg.push_back(v);
			vertIdx++;
		}
		//groupList.push_back(vg);
		groupCount++;
	}

	void PowerDiagram_3::end_insert(bool redraw /* = true */)
	{
		assert(is_valid() && open && dimension() == 3);
		all_vertices_.clear();
		all_cells_.clear();
		groupList.clear();
		groupList.resize(groupCount);
		// find a point inside the triangulation
		vec3 midpoint(0.0, 0.0, 0.0);
		for (Finite_vertices_iterator fvit = finite_vertices_begin(); fvit != finite_vertices_end(); fvit++)
			midpoint += to_geex(fvit->point());
		midpoint /= number_of_vertices();
		Point_3 interiorPoint = to_cgal(midpoint);
		for (All_vertices_iterator it = all_vertices_begin(); it!=all_vertices_end(); it++)
		{
			it->dual_intersects_boundary = false;
			if (!is_infinite(it))
			{
				all_vertices_.push_back(it);
				groupList[it->group_id].push_back(it);
				vector<Edge> incidentEdges;
				back_insert_iterator<vector<Edge>> bit(incidentEdges);
				finite_incident_edges(it, bit);
				for (vector<Edge>::const_iterator eit = incidentEdges.begin(); eit!=incidentEdges.end(); eit++)
				{
					vector<Object> face;//the face of the polyhedron encompassing this edge
					Vertex_handle v0 = eit->first->vertex(eit->second), v1 = eit->first->vertex(eit->third);
					Cell_circulator cellStart = incident_cells(*eit),
									cellEnd = cellStart;
					Cell_circulator cit = cellStart;
					do {
						cit->visited = false;
						cit++;
					}while (cit!=cellEnd);
					Cell_handle cur = cellStart, pre;
					bool loopEnd = false;
					unsigned int coplane_case = 0;
					while (!loopEnd)
					{
						int start_infidx = -1, end_infidx = -1;// index of infinite points, 0, 1, 2 or 3
						Triangle start_tri, end_tri;// only for infinite vertex case
						Point_3 startPoint, endPoint;
						if (is_infinite(cur))
						{
							start_infidx = cur->index(infinite_vertex());
							start_tri = triangle(cur, start_infidx);
							coplane_case = 0;
						}
						else if (coplane(cur))
						{
							start_infidx = 0;
							//vec3 p = to_geex(cur->vertex(0)->point()), q = to_geex(cur->vertex(1)->point()),
							//	 r = to_geex(cur->vertex(2)->point()), s = to_geex(cur->vertex(3)->point());
							//vec3 cp0 = cross(q-p, r-p), cp1 = cross(q-p, s-p);
							//if (std::fabs(cp0.length()) < 0.00000001 || std::fabs(cp1.length()) < 0.00000001)
								//std::cout<<cp0.length()<<" , "<<cp1.length()<<std::endl;
							if (!collinear(cur->vertex(0)->point(), cur->vertex(1)->point(), cur->vertex(2)->point()))
							{
								start_tri = Triangle(cur->vertex(0)->point(), cur->vertex(1)->point(), cur->vertex(2)->point());
								//if (std::fabs(cp0.length()) < 0.00000001)
								//	std::cout<<cp0.length()<<std::endl;
							}
							else 
							{
								start_tri = Triangle(cur->vertex(0)->point(), cur->vertex(1)->point(), cur->vertex(3)->point()); 
								//if (std::fabs(cp1.length()) < 0.00000001)
								//	std::cout<<cp1.length()<<std::endl;
							}
							coplane_case++;
						}
						else
						{
							startPoint = baseclass::dual(cur);
							coplane_case = 0;
						}
						cur->visited = true;
						bool isEdgeVertex[] = {false, false, false, false};
						int idx = cur->index(v0);
						isEdgeVertex[idx] = true;
						idx = cur->index(v1);
						isEdgeVertex[idx] = true;
						vector<Cell_handle> candidate;
						if (!isEdgeVertex[0])	candidate.push_back(cur->neighbor(0));
						if (!isEdgeVertex[1])	candidate.push_back(cur->neighbor(1));
						if (!isEdgeVertex[2])	candidate.push_back(cur->neighbor(2));
						if (!isEdgeVertex[3])	candidate.push_back(cur->neighbor(3));
						//pre = cur;
						if ( (candidate[0]->visited ^ candidate[1]->visited) || (!candidate[0]->visited && !candidate[1]->visited))
							cur = candidate[0]->visited ? candidate[1]:candidate[0];
						else
						{
							loopEnd = true;
							cur = cellEnd;// == start
						}
						if (is_infinite(cur))
						{
							end_infidx = cur->index(infinite_vertex());
							end_tri = triangle(cur, end_infidx);
							coplane_case=0;
						}
						else if (coplane(cur))
						{
							end_infidx = 0;
							//vec3 p = to_geex(cur->vertex(0)->point()), q = to_geex(cur->vertex(1)->point()),
							//	r = to_geex(cur->vertex(2)->point()), s = to_geex(cur->vertex(3)->point());
							//vec3 cp0 = cross(q-p, r-p), cp1 = cross(q-p, s-p);
							//if (std::fabs(cp0.length()) < 0.00000001 || std::fabs(cp1.length()) < 0.00000001)
								//std::cout<<cp0.length()<<" , "<<cp1.length()<<std::endl;
							if (!collinear(cur->vertex(0)->point(), cur->vertex(1)->point(), cur->vertex(2)->point()))
							{
								end_tri = Triangle(cur->vertex(0)->point(), cur->vertex(1)->point(), cur->vertex(2)->point());
								//if (std::fabs(cp0.length()) < 0.00000001)
								//	std::cout<<cp0.length()<<std::endl;
							}
							else
							{
								end_tri = Triangle(cur->vertex(0)->point(), cur->vertex(1)->point(), cur->vertex(3)->point());
								//if (std::fabs(cp1.length()) < 0.00000001)
								//	std::cout<<cp1.length()<<std::endl;
							}
							coplane_case++;
						}
						else 
						{
							endPoint = baseclass::dual(cur);
							coplane_case=0;
						}
						//if (coplane_case == 3)
						//	system("pause");
						if (start_infidx >= 0 && end_infidx >= 0) // the two cells are both infinite
							continue;
						else if ( start_infidx >= 0 && end_infidx < 0) 
						{
							Plane_3 sp = start_tri.supporting_plane();	
							Ray_3 r;
							if (sp.has_on_negative_side(interiorPoint))
								r = Ray_3(endPoint, sp.orthogonal_vector());
							else 
								r = Ray_3(endPoint, -sp.orthogonal_vector());
							face.push_back(make_object(r));
						}
						else if (start_infidx < 0 && end_infidx >= 0) 
						{
							Ray_3 r;
							Plane_3 sp = end_tri.supporting_plane();
							if (sp.has_on_negative_side(interiorPoint))
								r = Ray_3(startPoint, sp.orthogonal_vector());
							else
								r = Ray_3(startPoint, -sp.orthogonal_vector());
							face.push_back(make_object(r));
						}		
						else // normal case
						{
							Segment_3 s(startPoint, endPoint);
							face.push_back(make_object(s));
						}
					}// while (!loopEnd)
					it->vCell.push_back(face);	
					if (it == v0)	it->dualVert.push_back(v1->group_id);
					else			it->dualVert.push_back(v0->group_id);
				}// for every edge
			}//if (not infinite vertex)
		}// for every vertex
		if (tangentPlanes.size() > 0) // if we specify the tangent planes
			clipped_tangent_plane();
		open = false;
		//if (redraw)
		//	glut_viewer_redraw();
	}

	void PowerDiagram_3::single_clip( const Plane_3& p, const VGroup& vg/*, vector<Segment_3> *clippedRegion*/ )
	{
		//clippedRegion->clear();
		for (unsigned int i = 0; i < vg.size(); i++)
		{
			const vector<vector<Object>>& vCell = vg[i]->vCell;
			for (unsigned int j = 0; j < vCell.size(); j++)
			{
				//if (vg[i]->dualVert[j] == vg[i]->group_id)//bisector between the same group
					//continue;
				// compute the intersection with all bisectors
				const vector<Object>& face = vCell[j];
				vector<Point_3> intersectPoints; //we assume the face of the cell is convex, so number of intersection points must be 0, 1, 2
				vector<Point_3> db_src;
				vector<char> db_type;
				for (unsigned int k = 0; k < face.size(); k++)
				{
					if (const Segment_3* s = object_cast<Segment_3>(&face[k]))
					{
						db_src.push_back(s->source());
						db_src.push_back(s->target());
						db_type.push_back('s');
						if (!do_intersect(p, *s))
							continue;
						Object o = intersection(p, *s);
						if (const Point_3 *intp = object_cast<Point_3>(&o))
							intersectPoints.push_back(*intp);
						else if (const Segment_3 *ints = object_cast<Segment_3>(&o))
						{
							intersectPoints.push_back(ints->source());
							intersectPoints.push_back(ints->target());
						}
					}
					else if (const Ray_3 *r = object_cast<Ray_3>(&face[k]))
					{
						db_src.push_back(r->source());
						db_type.push_back('r');
						if (!do_intersect(p, *r))
							continue;
						Object o = intersection(p, *r);
						if (const Point_3 *intp = object_cast<Point_3>(&o))
							intersectPoints.push_back(*intp);
						else if (const Ray_3 *intr = object_cast<Ray_3>(&o))
							std::cerr<<"Impossible case at line "<<__LINE__<<" in file "<<__FILE__<<std::endl;
					}
				}
				//eliminate duplicate intersection points
				if (intersectPoints.size() >= 2)
				{
					vector<Point_3> duplicate;
					vector<Point_3>::iterator it = intersectPoints.begin();
					while (it != intersectPoints.end())
					{
						//if (std::find(duplicate.begin(), duplicate.end(), *it) == duplicate.end())
						if (approxFind(duplicate.begin(), duplicate.end(), *it) == duplicate.end())
						{
							duplicate.push_back(*it);
							it++;
						}
						else
							it = intersectPoints.erase(it);
					}
				}
				vector<double> onPlane;
				for (unsigned int db_i = 0; db_i < db_src.size(); db_i++)
				{
					Point_3 pp = db_src[db_i];
					onPlane.push_back(p.a()*pp.x() + p.b()*pp.y() + p.c()*pp.z() + p.d());
				}
				assert(intersectPoints.size() <= 2);
				if (intersectPoints.size()==2)	
				{
					//clippedRegion->push_back(Segment_3(intersectPoints[0], intersectPoints[1]));
					vg[i]->tangentIntersection.push_back(Segment_3(intersectPoints[0], intersectPoints[1]));
					vg[i]->tanIntersecDualVert.push_back(vg[i]->dualVert[j]);
				}
			}			
		}
#if 0
		if (clippedRegion->size() < 2)
			return;
		//std::cout<<"========================\n";
		//for (unsigned int n = 0; n < clippedRegion->size(); n++)
			//std::cout<<"("<<clippedRegion->at(n).source()<<"),"<<"("<<clippedRegion->at(n).target()<<"),"<<std::endl;
		//std::cout<<std::endl;
		// connect the segments in the clipped region computed above, if possible
		//double threshold = 1.0e-12;
		list<Segment_3> loop;
		vector<bool> determined(clippedRegion->size(), false);
		bool allDetermined = false;
		list<Segment_3>::iterator lit = loop.insert(loop.end(),clippedRegion->at(0));
		//for (list<Segment_3>::iterator it = loop.begin(); it!=loop.end(); it++)
		//	std::cout<<"("<<it->source()<<"),("<<it->target()<<")\n";
		determined[0] = true;
		std::cout << "Entering the merging process\n";
		while (!allDetermined)
		{
			allDetermined = true;
			for (lit = loop.begin(); lit!=loop.end(); lit++)
				for (unsigned int n = 0; n < clippedRegion->size(); n++)
				{
					allDetermined &= determined[n];
					if (determined[n])
						continue;
					const Segment_3& current = clippedRegion->at(n);
					//if (current.source() == lit->source() && current.target() != lit->target())
					//	lit = loop.insert(lit, current.opposite());
					//else if (current.source() == lit->target() && current.target() != lit->source())
					//	lit = loop.insert(lit, current);
					//else if (current.target() == lit->source() && current.source() != lit->target())
					//{
					//	lit++;
					//	lit = loop.insert(lit, current);
					//}
					//else if (current.target() == lit->target() && current.source() != lit->source())
					//{
					//	lit++;
					//	lit = loop.insert(lit, current.opposite());
					//}
					if (approxEqual(current.source(), lit->source(), threshold) && !approxEqual(current.target(), lit->target(), threshold))
						lit = loop.insert(lit, current.opposite());
					else if ( approxEqual(current.source(), lit->target(), threshold) && !approxEqual(current.target(), lit->source(), threshold))
					{	
						lit++;
						lit = loop.insert(lit, current);
					}
					else if (approxEqual(current.target(), lit->source(), threshold) && !approxEqual(current.source(), lit->target(), threshold))
					{
						lit = loop.insert(lit, current);
						//lit--;
					}
					else if (approxEqual(current.target(), lit->target(), threshold ) && !approxEqual(current.source(), lit->source(), threshold) )
					{
						lit++;
						lit = loop.insert(lit, current.opposite());
						//lit--;
					}
					else //this edge does not match
						continue;
					determined[n] = true;
 				//	std::cout<<"*************************\n";
 				//	for (list<Segment_3>::iterator it = loop.begin(); it!=loop.end(); it++)
 				//		std::cout<<"("<<it->source()<<"),("<<it->target()<<")\n";
					break;
				}
		}
		std::cout<<"Exiting the merging process\n";
		//if ( loop.front().source() != loop.back().target() )
		if (!approxEqual(loop.front().source() , loop.back().target(), threshold))
			loop.push_back(Segment_3(loop.back().target(), loop.front().source()));
		clippedRegion->resize(loop.size());
		unsigned int idx = 0;
		//std::cout<<"=================================\n";
		Segment_3 pre = loop.back();
		for (std::list<Segment_3>::iterator it = loop.begin(); it!=loop.end(); it++)
		{	
			//std::cout<<"("<<it->source()<<"),("<<it->target()<<")\n";
 			if ( !approxEqual(pre.target(), it->source(), threshold) )
			{
				std::cout<<"Not consistent here!\n";
				system("pause");
			}
			clippedRegion->at(idx) = *it;
			idx++;
			pre = *it;
		}
#endif
	}
	void PowerDiagram_3::accurate_single_clip(const Plane_3& p, const VGroup& vg)
	{
		Ext_plane_3 ep = to_ext(p);
		for (unsigned int i = 0; i < vg.size(); i++)
		{
			const vector<vector<Object>>& vCell = vg[i]->vCell;
			for (unsigned int j = 0; j < vCell.size(); j++)
			{
				// compute the intersection with all bisectors
				const vector<Object>& face = vCell[j];
				vector<Ext_point_3> intersectPoints; //we assume the face of the cell is convex, so number of intersection points must be 0, 1, 2
				for (unsigned int k = 0; k < face.size(); k++)
				{
					if (const Ext_segment_3* s = object_cast<Ext_segment_3>(&face[k]))
					{
						if (!do_intersect(ep, *s))
							continue;
						Object o = intersection(ep, *s);
						if (const Ext_point_3 *intp = object_cast<Ext_point_3>(&o))
							intersectPoints.push_back(*intp);
						else if (const Ext_segment_3 *ints = object_cast<Ext_segment_3>(&o))
						{
							intersectPoints.push_back(ints->source());
							intersectPoints.push_back(ints->target());
						}
					}
					else if (const Ray_3 *r = object_cast<Ray_3>(&face[k]))
					{
						if (!do_intersect(p, *r))
							continue;
						Object o = intersection(p, *r);
						if (const Ext_point_3 *intp = object_cast<Ext_point_3>(&o))
							intersectPoints.push_back(*intp);
						else if (const Ext_ray_3 *intr = object_cast<Ext_ray_3>(&o))
							std::cerr<<"Impossible case at line "<<__LINE__<<" in file "<<__FILE__<<std::endl;
					}
				}
				//eliminate duplicate intersection points
				if (intersectPoints.size() >= 2)
				{
					vector<Ext_point_3> duplicate;
					vector<Ext_point_3>::iterator it = intersectPoints.begin();
					while (it != intersectPoints.end())
					{
						if (std::find(duplicate.begin(), duplicate.end(), *it) == duplicate.end())
						//if (approxFind(duplicate.begin(), duplicate.end(), *it) == duplicate.end())
						{
							duplicate.push_back(*it);
							it++;
						}
						else
							it = intersectPoints.erase(it);
					}
				}
				assert(intersectPoints.size() <= 2);
				if (intersectPoints.size()==2)	
				{
					//clippedRegion->push_back(Segment_3(intersectPoints[0], intersectPoints[1]));
					vg[i]->tangentIntersection.push_back(Segment_3(inv_to_ext(intersectPoints[0]), inv_to_ext(intersectPoints[1])));
					vg[i]->tanIntersecDualVert.push_back(vg[i]->dualVert[j]);
				}
			}			
		}
	}
	//void PowerDiagram_3::accurate_single_clip( const Plane_3& p, const VGroup& vg/*, vector<Segment_3> *clippedRegion*/ )
	//{
	//	//clippedRegion->clear();
	//	vector<Ext_segment_3> extClippedRegion;
	//	Ext_plane_3 ep = to_ext(p);
	//	for (unsigned int i = 0; i < vg.size(); i++)
	//	{
	//		const vector<vector<Object>>& vCell = vg[i]->vCell;
	//		for (unsigned int j = 0; j < vCell.size(); j++)
	//		{
	//			if (vg[i]->dualVert[j] == vg[i]->group_id)
	//				continue;
	//			const vector<Object>& face = vCell[j];
	//			vector<Ext_point_3> intersectionPoints;
	//			//convert to MP_Float representation, then do the computation
	//			for (unsigned int k = 0; k < face.size(); k++)
	//			{
	//				if (const Segment_3 *s = object_cast<Segment_3>(&face[k]))
	//				{
	//					Ext_segment_3 es = to_ext(*s);
	//					if (!do_intersect(ep, es))
	//						continue;
	//					Object o = intersection(ep, es);
	//					if (const Ext_point_3 *intp = object_cast<Ext_point_3>(&o))
	//						intersectionPoints.push_back(*intp);
	//					else if (const Ext_segment_3 *ints = object_cast<Ext_segment_3>(&o))
	//					{
	//						intersectionPoints.push_back(ints->source());
	//						intersectionPoints.push_back(ints->target());
	//					}
	//				}
	//				else if (const Ray_3 *r = object_cast<Ray_3>(&face[k]))
	//				{
	//					Ext_ray_3 er = to_ext(*r);
	//					if (!do_intersect(ep, er))
	//						continue;
	//					Object o = intersection(ep, er);
	//					if (const Ext_point_3 *intp = object_cast<Ext_point_3>(&o))
	//						intersectionPoints.push_back(*intp);
	//					else if (const Ext_ray_3 *intr = object_cast<Ext_ray_3>(&o))
	//						std::cerr<<"Impossible case at line "<<__LINE__<<" in file "<<__FILE__<<std::endl;
	//				}
	//			}
	//			// eliminate duplicate intersection points
	//			if (intersectionPoints.size() >= 2)
	//			{
	//				vector<Ext_point_3> duplicate;
	//				vector<Ext_point_3>::iterator it = intersectionPoints.begin();
	//				while (it != intersectionPoints.end())
	//				{
	//					if (std::find(duplicate.begin(), duplicate.end(), *it) == duplicate.end())
	//					{
	//						duplicate.push_back(*it);
	//						it++;
	//					}
	//					else
	//						it = intersectionPoints.erase(it);
	//				}
	//			}
	//			assert(intersectionPoints.size()<=2);
	//			if (intersectionPoints.size()==2)
	//				extClippedRegion.push_back(Ext_segment_3(intersectionPoints[0], intersectionPoints[1]));
	//		}// end of a face
	//	}// end of a vertex
	//	if (extClippedRegion.size() < 2)
	//		return;
	//	list<Ext_segment_3> loop;
	//	vector<bool> determined(extClippedRegion.size(), false);
	//	bool allDetermined = false;
	//	list<Ext_segment_3>::iterator lit = loop.insert(loop.end(), extClippedRegion[0]);
	//	determined[0] = true;
	//	std::cout << "Enter while loop\n";
	//	while (!allDetermined)
	//	{
	//		allDetermined = true;
	//		for (lit = loop.begin(); lit != loop.end(); lit++)
	//			for (unsigned int n = 0; n < extClippedRegion.size(); n++)
	//			{
	//				allDetermined &= determined[n];
	//				if (determined[n])
	//					continue;
	//				const Ext_segment_3& current = extClippedRegion[n];
	//				if (current.source() == lit->source() && current.target() != lit->target())
	//					lit = loop.insert(lit, current.opposite());
	//				else if (current.source() == lit->target() && current.target() != lit->source())
	//				{
	//					lit++;
	//					lit = loop.insert(lit, current);
	//				}
	//				else if (current.target() == lit->source() && current.source() != lit->target())
	//					lit = loop.insert(lit, current);
	//				else if (current.target() == lit->target() && current.source() != lit->source())
	//				{
	//					lit++;
	//					lit = loop.insert(lit, current.opposite());
	//				}
	//				else
	//					continue;
	//				determined[n] = true;
	//				break;
	//			}
	//	}
	//	std::cout<<"Exit while loop\n";
	//	if (loop.front().source() != loop.back().target())
	//		loop.push_back(Ext_segment_3(loop.back().target(), loop.front().source()));
	//	extClippedRegion.resize(loop.size());
	//	unsigned int idx = 0;
	//	for (std::list<Ext_segment_3>::iterator it = loop.begin(); it != loop.end(); it++, idx++)
	//		extClippedRegion[idx] = *it;
	//	//convert back to double precision 
	//	for (vector<Ext_segment_3>::const_iterator it = extClippedRegion.begin(); it != extClippedRegion.end(); it++)
	//	{
	//		Point_3 s(to_double(it->source().x()), to_double(it->source().y()), to_double(it->source().z()));
	//		Point_3 t(to_double(it->target().x()), to_double(it->target().y()), to_double(it->target().z()));
	//		//clippedRegion->push_back(Segment_3(s, t));
	//	}
	//}
	void PowerDiagram_3::clipped_tangent_plane()
	{
		clipped_tangent_plane_.resize(groupCount);
		for (unsigned int i = 0; i < groupList.size(); i++)
			//accurate_single_clip(tangentPlanes[i], groupList[i]/*, &clipped_tangent_plane_[i]*/);
			single_clip(tangentPlanes[i], groupList[i]/*, &clipped_tangent_plane_[i]*/);
	}
}