#include "power_diagram.h"

namespace Geex
{
	BisectorEdge::SPLIT_STATE BisectorEdge::split(const vector<Bisector>& bisecList, list<Vertex_handle> *left, list<Vertex_handle>* right)
	{
		left->clear();
		right->clear();
		list<int>::const_iterator cpyit, minit, it;
		double minDist = DBL_MAX, maxDist = DBL_MIN;
		// find the end of the left list
		for (it = segIdx.begin(); it != segIdx.end(); it++)
		{
			double dist = CGAL::squared_distance(bisecList[*it].sphere->point(), bisecList[*it].seg);
			if ( dist < minDist)
			{
				minDist = dist;
				minit = it;
			}
			if ( dist > maxDist)
				maxDist = dist;
		}
		if (maxDist < CGAL::squared_distance(bisecList[*minit].sphere->point(), bisecList[*minit].seg)*4.0)
			return EMPTY;
		if (minit != segIdx.end() && minDist > CGAL::squared_distance(bisecList[*minit].sphere->point(), bisecList[*minit].seg)*16)
			minit = segIdx.end();
		if (minit != segIdx.end())
			minit++;//append one more, make left and right connected
		for ( cpyit = segIdx.begin(); cpyit != minit; cpyit++ )
			left->push_back(bisecList[*cpyit].sphere);
		
		if ( minit == segIdx.end() ) // all left list
		{
			*right = *left;
			return FULL;
		}	
		minit--;
		for (cpyit = minit; cpyit != segIdx.end(); cpyit++)
			right->push_back(bisecList[*cpyit].sphere);
		if ( right->size() == 0 && left->size() == 0)
			return EMPTY;
		return HALF;
	}
	//BisectorEdge::SPLIT_STATE BisectorEdge::split(const vector<Bisector>& bisecList, list<Vertex_handle> *left, list<Vertex_handle>* right)
	//{
	//	left->clear();
	//	right->clear();
	//	list<int>::const_iterator cpyit, minit, it;
	//	double minDist = DBL_MAX, maxDist = DBL_MIN;
	//	double w = 0.25;
	//	// find the end of the left list
	//	for (it = segIdx.begin(); it != segIdx.end(); it++)
	//	{
	//		double sqdist = bisecList[*it].sqDist;
	//		if ( sqdist < minDist)
	//		{
	//			minDist = sqdist;
	//			minit = it;
	//		}
	//		if ( sqdist > maxDist)
	//			maxDist = sqdist;
	//	}
	//	if (maxDist < (1+w)*(1+w)*bisecList[*minit].sphere->point().weight())
	//		return EMPTY;
	//	if (minit != segIdx.end() && minDist > bisecList[*minit].sphere->point().weight()*4.0)
	//		minit = segIdx.end();
	//	for ( cpyit = segIdx.begin(); cpyit != minit; cpyit++ )
	//		left->push_back(bisecList[*cpyit].sphere);

	//	if ( minit == segIdx.end() ) // all left list
	//	{
	//		*right = *left;
	//		return FULL;
	//	}	
	//	for (cpyit = minit; cpyit != segIdx.end(); cpyit++)
	//		right->push_back(bisecList[*cpyit].sphere);
	//	if ( right->size() == 0 && left->size() == 0)
	//		return EMPTY;
	//	return HALF;
	//}
	void BisectorEdge::rightSplitAndVisit(list<Vertex_handle>* l)
	{
		if (visited!=FULLY_VISITED && visited!=RIGHT_VISITED)
		{
			if (ss == FULL || ss == EMPTY)
				visited = FULLY_VISITED;
			else //if (rlist.size() > 0)
				setRightVisited();			
			l->insert(l->end(), rlist.begin(), rlist.end());
		}
	}
	void BisectorEdge::leftSplitAndVisit(list<Vertex_handle>* l)
	{
		if (visited!=FULLY_VISITED && visited!=LEFT_VISITED)
		{
			if (ss == FULL || ss == EMPTY)
				visited = FULLY_VISITED;
			else //if (llist.size() > 0)
				setLeftVisited();				
			l->insert(l->end(), llist.begin(), llist.end());
		}
	}
	PowerDiagram_3::PowerDiagram_3(TriMesh* m): bdMesh(m), open(false), threshold(1.0e-12), use_omp(0){}

	void PowerDiagram_3::clear()
	{
		baseclass::clear();
		all_vertices_.clear();
		//clipped_tangent_plane_.clear();
		groupList.clear();
		tangentPlanes.clear();
		all_vertices_.clear();
		holes.clear();
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
			if (v == 0) //hidden point
				continue;
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
		groupList.clear();
		groupList.resize(groupCount);
		for (unsigned int i = 0; i < clippedPolygons.size(); i++)
			clippedPolygons[i].clear();
		clippedPolygons.resize(groupCount);
		for (unsigned int i = 0; i < holes.size(); i++)
			holes[i].clear();
		holes.clear();
		holeArea.clear();
		rawHoles.clear();
		bisecEdgeSet.clear();
		// find a point inside the triangulation
		vec3 midpoint(0.0, 0.0, 0.0);
		int nbValid = 0;
		for (Finite_vertices_iterator fvit = finite_vertices_begin(); fvit != finite_vertices_end(); fvit++)
		{
			if (fvit->group_id >= 0)
			{
				midpoint += to_geex(fvit->point());
				nbValid++;
			}
		}
		//midpoint /= number_of_vertices();
		midpoint /= nbValid;
		Point_3 interiorPoint = to_cgal(midpoint);
		//baseclass::insert(Weighted_point(interiorPoint, 0.0));
		///////////////	added for use of OpenMP	///////////////
		int nbThread = 1;
		All_vertices_iterator *starts, *ends;
		//if (use_omp)
		//{
		//	nbThread = std::min<int>(int(number_of_vertices()), omp_get_num_procs());
		//	omp_set_num_threads(nbThread);
		//}
#ifdef CILK
		nbThread = __cilkrts_get_nworkers();
#endif
		starts = new All_vertices_iterator[nbThread];
		ends = new All_vertices_iterator[nbThread];
		int step = number_of_vertices()/nbThread;
		starts[0] = all_vertices_begin();
		for (int i = 1; i < nbThread; i++)
		{
			starts[i] = starts[i-1];
			std::advance(starts[i], step);
		}
		for (int i = 0; i < nbThread-1; i++)
			ends[i] = starts[i+1];
		ends[nbThread-1] = all_vertices_end();
		// construct the visit mask of each cell
		for (All_cells_iterator cit = all_cells_begin(); cit!=all_cells_end(); cit++)
			cit->visited.resize(nbThread);//each thread use an independent mask on the others
		////////////////////////////////////////////////////////
//#pragma omp parallel for schedule(static) if(use_omp)
#ifdef CILK
		cilk_for (int i = 0; i < nbThread; i++)
#else
		for (int i = 0; i < nbThread; i++)
#endif
			for (All_vertices_iterator it = starts[i]; it!=ends[i]; it++)
			//for (All_vertices_iterator it = all_vertices_begin(); it!=all_vertices_end(); it++)
			{
				if (!is_infinite(it) && it->group_id >= 0)
				{
					//all_vertices_.push_back(it);
					//groupList[it->group_id].push_back(it);
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
							cit->visited[i] = false;
							//cit->visited = false;
							cit++;
						}while (cit!=cellEnd);
						/** traversal for an edge loop start here**/
						// while loop variable
						Cell_handle cur = cellStart, pre;
						int start_infidx = -1, end_infidx = -1;// index of infinite points, 0, 1, 2 or 3
						Triangle start_tri, end_tri;// only for infinite vertex case
						Point_3 startPoint, endPoint;
						bool loopEnd = false;
						bool startExt = false, endExt = false;
						Ext_triangle_3 start_tri_ext, end_tri_ext;
						// check the property of the first cell
						if (is_infinite(cur))
						{
							start_infidx = cur->index(infinite_vertex());
							start_tri = triangle(cur, start_infidx);
							//startInf = true;
						}
						else if (coplane(cur))
						{
							start_infidx = 0;
							if (!collinear(cur->vertex(0)->point(), cur->vertex(1)->point(), cur->vertex(2)->point()))
								start_tri = Triangle(cur->vertex(0)->point(), cur->vertex(1)->point(), cur->vertex(2)->point());
							else 
								start_tri = Triangle(cur->vertex(0)->point(), cur->vertex(1)->point(), cur->vertex(3)->point()); 
							startExt = false;
							//assert(!collinear(start_tri.vertex(0), start_tri.vertex(1), start_tri.vertex(2)));
							if ( collinear(start_tri.vertex(0), start_tri.vertex(1), start_tri.vertex(2)) )
							{
								Vector_3 det0 = CGAL::cross_product(Vector_3(cur->vertex(0)->point(), cur->vertex(1)->point()),
																	Vector_3(cur->vertex(0)->point(), cur->vertex(2)->point()));
								Vector_3 det1 = CGAL::cross_product(Vector_3(cur->vertex(0)->point(), cur->vertex(1)->point()),
																	Vector_3(cur->vertex(0)->point(), cur->vertex(3)->point()));
								if ( det0.squared_length() >= det1.squared_length() )
									start_tri_ext = to_ext(Triangle(cur->vertex(0)->point(), cur->vertex(1)->point(), cur->vertex(2)->point()));
								else
									start_tri_ext = to_ext(Triangle(cur->vertex(0)->point(), cur->vertex(1)->point(), cur->vertex(3)->point()));
								startExt = true;
							}
						}
						else
						{
							start_infidx = -1;
							startPoint = baseclass::dual(cur);
						}
						while (!loopEnd)
						{
							// get the next neighboring cell
							cur->visited[i] = true;
							//cur->visited = true;
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
							pre = cur;
							if ( (candidate[0]->visited[i] ^ candidate[1]->visited[i]) || (!candidate[0]->visited[i] && !candidate[1]->visited[i]))
								cur = candidate[0]->visited[i] ? candidate[1]:candidate[0];
							else
							{
								loopEnd = true;
								cur = cellEnd;// == start
							}
							if (is_infinite(cur))
							{
								end_infidx = cur->index(infinite_vertex());
								end_tri = triangle(cur, end_infidx);
								//endExt = false;
							}
							else if (coplane(cur))
							{
								end_infidx = 0;
								if (!collinear(cur->vertex(0)->point(), cur->vertex(1)->point(), cur->vertex(2)->point()))
									end_tri = Triangle(cur->vertex(0)->point(), cur->vertex(1)->point(), cur->vertex(2)->point());
								else
									end_tri = Triangle(cur->vertex(0)->point(), cur->vertex(1)->point(), cur->vertex(3)->point());
								endExt = false;
								//assert(!collinear(end_tri.vertex(0), end_tri.vertex(1), end_tri.vertex(2)));
								//endInf = false;
								if ( collinear(end_tri.vertex(0), end_tri.vertex(1), end_tri.vertex(2)) )
								{
									Vector_3 det0 = CGAL::cross_product(Vector_3(cur->vertex(0)->point(), cur->vertex(1)->point()),
																		Vector_3(cur->vertex(0)->point(), cur->vertex(2)->point()));
									Vector_3 det1 = CGAL::cross_product(Vector_3(cur->vertex(0)->point(), cur->vertex(1)->point()),
																		Vector_3(cur->vertex(0)->point(), cur->vertex(3)->point()));
									if ( det0.squared_length() >= det1.squared_length() )
										end_tri_ext = to_ext(Triangle(cur->vertex(0)->point(), cur->vertex(1)->point(), cur->vertex(2)->point()));
									else
										end_tri_ext = to_ext(Triangle(cur->vertex(0)->point(), cur->vertex(1)->point(), cur->vertex(3)->point()));
									endExt = true;
								}
							}
							else 
							{
								end_infidx = -1;
								endPoint = baseclass::dual(cur);
								//endInf = false;
							}
							if (start_infidx >= 0 && end_infidx >= 0) // the two cells are both infinite
							{
								//go to the end of the loop
							}
							else if ( start_infidx >= 0 && end_infidx < 0) 
							{
								Ray_3 r;
								if (!startExt)
								{
									Plane_3 sp = start_tri.supporting_plane();	
									if (sp.has_on_negative_side(interiorPoint))
										r = Ray_3(endPoint, sp.orthogonal_vector());
									else 
										r = Ray_3(endPoint, -sp.orthogonal_vector());
								}
								else
								{
									Ext_plane_3 sp = start_tri_ext.supporting_plane();
									if (sp.has_on_negative_side(to_ext(interiorPoint)))
										r = Ray_3(endPoint, inv_to_ext(sp.orthogonal_vector()));
									else
										r = Ray_3(endPoint, -inv_to_ext(sp.orthogonal_vector()));
								}
								face.push_back(make_object(r));
							}
							else if (start_infidx < 0 && end_infidx >= 0) 
							{
								Ray_3 r;
								if (!endExt)
								{
									Plane_3 sp = end_tri.supporting_plane();
									if (sp.has_on_negative_side(interiorPoint))
										r = Ray_3(startPoint, sp.orthogonal_vector());
									else
										r = Ray_3(startPoint, -sp.orthogonal_vector());
								}
								else
								{
									Ext_plane_3 sp = end_tri_ext.supporting_plane();
									if (sp.has_on_negative_side(to_ext(interiorPoint)))
										r = Ray_3(startPoint, inv_to_ext(sp.orthogonal_vector()));
									else
										r = Ray_3(startPoint, -inv_to_ext(sp.orthogonal_vector()));
								}
								face.push_back(make_object(r));
							}		
							else // normal case
							{
								Segment_3 s(startPoint, endPoint);
								face.push_back(make_object(s));
							}
							// the end of the previous becomes the start of the next
							start_infidx = end_infidx;
							start_tri = end_tri;
							startPoint = endPoint;
							if (startExt = endExt)
								start_tri_ext = end_tri_ext;
							
							//startInf = endInf;
						}// while (!loopEnd)
						it->vCell.push_back(face);	
						if (it == v0)	{ it->dualVert.push_back(v1->group_id); it->dualVertHandle.push_back(v1); }
						else			{ it->dualVert.push_back(v0->group_id); it->dualVertHandle.push_back(v0); }
					}// for every edge
				}//if (not infinite vertex)
			}// for every vertex
		delete []starts;
		delete []ends;
		for (All_vertices_iterator it = all_vertices_begin(); it!=all_vertices_end(); it++)
		//for ( unsigned int i = 0; i < groupList.size(); i++)
		{
			if (!is_infinite(it) && it->group_id >= 0)
			{
				all_vertices_.push_back(it);
				groupList[it->group_id].push_back(it);
			}
		}
 		if (tangentPlanes.size() > 0) // if we specify the tangent planes
 			clipped_tangent_plane();
		open = false;
	}
	void PowerDiagram_3::preDetectHoles()
	{
		for (int group_id = 0; group_id < getNbGroups(); group_id++)
		{
			vector<Bisector>& clippedPolygon = clippedPolygons[group_id];
			int idx = 0;
			int preNeighbor = -1, curNeighbor;
			list<BisectorEdge> tmp;
			curNeighbor = clippedPolygon[idx].dualVertGroupId;
			tmp.push_back(BisectorEdge(group_id, curNeighbor, preNeighbor, -1, &clippedPolygon));
			preNeighbor = curNeighbor;
			tmp.back().segIdx.push_back(idx);
			idx++;
			for (; idx < clippedPolygon.size(); idx++)
			{	
				curNeighbor = clippedPolygon[idx].dualVertGroupId;
				if (curNeighbor != preNeighbor)
				{
					tmp.back().end.second = curNeighbor;// complete the previous edge
					if (tmp.size() > 0 && curNeighbor == tmp.front().base.second) // we have been back
					{
						tmp.front().end.first = preNeighbor;
						break;
					}
					tmp.push_back(BisectorEdge(group_id, curNeighbor, preNeighbor, -1, &clippedPolygon));
					preNeighbor = curNeighbor;
				}
				tmp.back().segIdx.push_back(idx);
			}
			list<int>::iterator insertPos = tmp.front().segIdx.begin();
			// consider the remaining part of the first edge
			if (curNeighbor != tmp.front().base.second)
				tmp.front().end.first = curNeighbor;
			for (; idx < clippedPolygon.size(); idx++)
				tmp.front().segIdx.insert(insertPos, idx);
			if (tmp.back().end.second < 0)
				tmp.back().end.second = curNeighbor;
			for (list<BisectorEdge>::const_iterator cpyit = tmp.begin(); cpyit != tmp.end(); cpyit++)
			{
				pair<int, int> edgeBase = make_pair(cpyit->base.first, cpyit->base.second);
				pair<map<pair<int, int>, BisectorEdge, BisectorEdgeCmp>::iterator, bool> e = bisecEdgeSet.insert(make_pair(edgeBase, *cpyit));
				e.first->second.split();
			}
		}
	}
	void PowerDiagram_3::preFinalDetectHoles(const vector<vec3>& normals)
	{
		for_each(sholes.begin(), sholes.end(), mem_fun_ref(&vector<Weighted_point>::clear));
		sholes.clear();
		// compute the 1-ring neighbor for each polygon
		vector<set<int>> ring_1(groupList.size());
		for (unsigned int i = 0; i < groupList.size(); i++)
		{
			VGroup& vg = groupList[i];
			for (unsigned int j = 0; j < vg.size(); j++)
			{
				Vertex_handle& v = vg[j];
				vector<int>& dualVertGroupId = v->dualVert;
				for (unsigned int k = 0; k < dualVertGroupId.size(); k++)
				{
					if (dualVertGroupId[k] >= 0)
						ring_1[i].insert(dualVertGroupId[k]);
				}
			}
		}
		// compute 2-ring
		vector<set<int>> ring_2(groupList.size());
		for (unsigned int i = 0; i < groupList.size(); i++)
		{
			ring_2[i].insert(ring_1[i].begin(), ring_1[i].end());
			for (set<int>::const_iterator it = ring_1[i].begin(); it != ring_1[i].end(); it++)
				ring_2[i].insert(ring_1[*it].begin(), ring_1[*it].end());//!!! note: this may include the group i itself
		}
		// compute distance from spheres to bisector edges
		vector<DistBisectorTriplet> dbt;
		static const double PI = 3.14159265358979323846;
		double normalAngleThreshold = std::cos(PI/6.0);
		for (unsigned int i = 0; i < groupList.size(); i++)
		{
			VGroup& vg = groupList[i];
			//const vec3& nthis = normals[i];
			for (unsigned int j = 0; j < vg.size(); j++)
			{
				vector<Segment_3>& bisec = vg[j]->tangentIntersection;
				vector<int>& dualVertGroupId = vg[j]->dualVert;
				for (int k = 0; k < bisec.size(); k++)
				{				
					//const vec3& ndual = normals[dualVertGroupId[k]];
					//if (dot(nthis, ndual)/(nthis.length()*ndual.length()) <= normalAngleThreshold)
					//	continue;
					double d = CGAL::squared_distance(vg[j]->point(), bisec[k]);
					//if (d >= 8*8*vg[j]->radius*vg[j]->radius)
					//	continue;
					dbt.push_back(DistBisectorTriplet(d, vg[j], k));
				}
				vg[j]->included = false;
			}
		}
		std::sort(dbt.begin(), dbt.end(), DistBisectorCmp());
		// from the longest, construct spheres at holes
		//vector<vector<Weighted_point>> holes;
		for (unsigned int i = 0; i < dbt.size(); i++)
		{
			if (dbt[i].v->included)
				continue;
			if (dbt[i].distance <= 3*dbt[i].v->radius*dbt[i].v->radius)//too small, not considered any longer
				break;
			Segment_3& b = dbt[i].getBisector();
			//double distance as radius
			Weighted_point s(CGAL::midpoint(b.source(), b.target()), dbt[i].distance*2.5);
			sholes.push_back(vector<Weighted_point>());
			//search 2-ring neighbor
			set<int>& ngh = ring_2[dbt[i].v->group_id];
			int nbGroupIncluded = 0;
			for (set<int>::const_iterator it = ngh.begin(); it != ngh.end(); it++)
			{
				bool groupIncluded = false;
				VGroup& neighGroup = groupList[*it];
				int nbIncluded = 0;
				for (unsigned int j = 0; j < neighGroup.size(); j++)
				{
					Point_3 c = neighGroup[j]->point();
					double sr = neighGroup[j]->radius*neighGroup[j]->radius;
					if (sphereContain(s.point(), s.weight(), c, sr))
					{
						groupIncluded = true;
						sholes.back().push_back(Weighted_point(c, sr));
						neighGroup[j]->included = true;
						nbIncluded++;
					}
				}
				if (nbIncluded >= neighGroup.size()/2.0)
				{
					sholes.pop_back();
					break;
				}
				if (groupIncluded)
					nbGroupIncluded++;
			}
			if (nbGroupIncluded == 1)
				sholes.pop_back();
			//holespheres.push_back(s);
			dbt[i].v->included = true;
		}
	}
	/*  traverse along the edge list until either of the following situations
	1. the end of the list is reached; 2. the distance to the corresponding sphere less than threshold
	if (1), turn to its neighbors at the end;(at most one neighbor will be considered)
	if (2), turn to its opposite neighbor
	*/
	inline void PowerDiagram_3::recoverBisectorEdges(const list<BisectorEdge*>& rcvList)
	{
		for (list<BisectorEdge*>::const_iterator it = rcvList.begin(); it != rcvList.end(); it++)
			if ((*it)->reversed)
				(*it)->reverse();
	}
	void PowerDiagram_3::simplifyHole(list<Vertex_handle>* hole)
	{
		list<Vertex_handle> simHole;
		list<Vertex_handle>::const_iterator pre = hole->end(), cur, nxt;
		pre--;//point to the last element
		nxt = cur = hole->begin();
		nxt++;
		while ( cur != hole->end() )
		{
			Point_3 p0 = (*pre)->point();
			Point_3 p1 = (*cur)->point();
			Point_3 p2 = (*nxt)->point();
			if (!collinear(p0, p1, p2))
				simHole.push_back(*cur);
			pre = cur;
			cur++;
			nxt++;
			if (nxt == hole->end())
				nxt = hole->begin();
		}
		if (simHole.size() == hole->size())
			return;
		hole->assign(simHole.begin(), simHole.end());
	}
	void PowerDiagram_3::simplifyPolygon(vector<Point_3> *pgn)
	{
		vector<Point_3> simPgn;
		for ( int i = 0; i < pgn->size(); i++)
		{
			const Point_3& p0 = pgn->at((i+pgn->size()-1)%pgn->size());
			const Point_3& p1 = pgn->at(i);
			const Point_3& p2 = pgn->at((i+1)%pgn->size());
			if (!collinear(p0, p1, p2))
				simPgn.push_back(p1);
		}
		if (simPgn.size() == pgn->size())
			return;
		pgn->clear();
		for (vector<Point_3>::const_iterator it = simPgn.begin(); it!=simPgn.end(); it++)
			pgn->push_back(*it);
	}
	void PowerDiagram_3::detectHoles()
	{
		/* for each bisector edge */
		std::cout<<"Start detecting holes"<<std::endl;
		preDetectHoles();
		list<BisectorEdge*> rcvList;
		for (BisectorEdgeIterator it = bisecEdgeSet.begin(); it != bisecEdgeSet.end(); it++, recoverBisectorEdges(rcvList))
		{
			if (it->second.visited == BisectorEdge::FULLY_VISITED )
				continue;
			rcvList.clear();//bisector edges that need recovering from reverse
			BisectorEdge& sideAB = it->second;
			BisectorEdgeIterator findIter;
			/***  The first end neighbor ***/
			findIter = bisecEdgeSet.find(sideAB.firstEndNeighbor());
			if (findIter == bisecEdgeSet.end())	continue;
			BisectorEdge& sideAC = findIter->second;
			/*** bisector edge CA on C's side ***/
			findIter = bisecEdgeSet.find(sideAC.oppositeSideIndex());
			if (findIter == bisecEdgeSet.end())
			{
				findIter = bisecEdgeSet.find(sideAC.firstEndNeighbor());
				if (findIter == bisecEdgeSet.end())
				{ 
					//std::cout<<"Cannot find at "<<__LINE__<<std::endl; 
					continue;
				}
				BisectorEdge& tmp = findIter->second;
				findIter = bisecEdgeSet.find(tmp.oppositeSideIndex());
				if (findIter == bisecEdgeSet.end()){ /*std::cout<<"Cannot find at "<<__LINE__<<std::endl;*/ continue;}
				if (findIter->second.sameDirection(tmp))
					findIter->second.reverse();
			}
			else
			{
				if (findIter->second.sameDirection(sideAC))
					findIter->second.reverse();
			}
			BisectorEdge& sideCA = findIter->second;
			rcvList.push_back(&sideCA);
			/***  The second end neighbor ***/
			findIter = bisecEdgeSet.find(sideAB.secondEndNeighbor());
			if (findIter == bisecEdgeSet.end())
			{
				//std::cout<<"Cannot find at "<<__LINE__<<std::endl; 
				continue;
			}
			BisectorEdge& sideAD = findIter->second;
			/*** bisector edge AD on D's side ***/
			findIter = bisecEdgeSet.find(sideAD.oppositeSideIndex());
			if (findIter == bisecEdgeSet.end())//avoid the case where the opposite does not exist
			{
				findIter = bisecEdgeSet.find(sideAD.secondEndNeighbor());
				if (findIter == bisecEdgeSet.end())
				{
					//std::cout<<"Cannot find at "<<__LINE__<<std::endl;
					continue;
				}
				BisectorEdge& tmp = findIter->second;
				findIter = bisecEdgeSet.find(tmp.oppositeSideIndex());
				if (findIter == bisecEdgeSet.end())
				{
					//std::cout<<"Cannot find at "<<__LINE__<<std::endl;
					continue;
				}
				if (findIter->second.sameDirection(tmp))
					findIter->second.reverse();
			}
			else
			{
				if (findIter->second.sameDirection(sideAD))
					findIter->second.reverse();
			}
			BisectorEdge& sideDA = findIter->second;
			rcvList.push_back(&sideDA);
			/*** bisector edge BD on D's side ***/
			findIter = bisecEdgeSet.find(sideDA.secondEndNeighbor());
			if (findIter == bisecEdgeSet.end())
			{
				//std::cout<<"Cannot find at "<<__LINE__<<std::endl;
				continue;
			}
			BisectorEdge& sideDB = findIter->second;
			if (sideDA.reversed)
				sideDB.reverse();
			rcvList.push_back(&sideDB);
			/*** opposite side of sideAB ***/
			findIter = bisecEdgeSet.find(sideAB.oppositeSideIndex());
			BisectorEdgeIterator sideBA = findIter;
			if (sideBA != bisecEdgeSet.end() && sideBA->second.sameDirection(sideAB))
			{
				sideBA->second.reverse();
				rcvList.push_back(&sideBA->second);
			}
			/*** bisector edge BD on B's side ***/
			if (sideBA != bisecEdgeSet.end())
			{
				if (sideBA->second.ss == BisectorEdge::EMPTY || sideBA->second.visited == BisectorEdge::FULLY_VISITED)
				{
					sideAB.visited = BisectorEdge::FULLY_VISITED;
					sideBA->second.visited = BisectorEdge::FULLY_VISITED;
					continue;
				}
				findIter = bisecEdgeSet.find(sideBA->second.firstEndNeighbor());
				if (findIter == bisecEdgeSet.end())
				{
					//std::cout<<"Cannot find at "<<__LINE__<<std::endl;
					continue;
				}	
				if (sideBA->second.reversed)
					findIter->second.reverse();
			}
			else
			{
				findIter = bisecEdgeSet.find(sideDB.oppositeSideIndex());
				if (findIter == bisecEdgeSet.end())
				{
					//std::cout<<"Cannot find at "<<__LINE__<<std::endl;
					continue;
				}
				if (findIter->second.sameDirection(sideDB))
					findIter->second.reverse();
			}
			BisectorEdge& sideBD = findIter->second;
			rcvList.push_back(&sideBD);
			/*** bisector edge BC on B's side ***/
			if (sideBA != bisecEdgeSet.end())
			{
				findIter = bisecEdgeSet.find(sideBA->second.secondEndNeighbor());
				if (findIter == bisecEdgeSet.end())
				{
					//std::cout<<"Cannot find at "<<__LINE__<<std::endl;
					continue;
				}
				if (sideBA->second.reversed)
					findIter->second.reverse();
			}
			else
			{
				findIter = bisecEdgeSet.find(sideCA.firstEndNeighbor());
				if (findIter == bisecEdgeSet.end())
				{
					//std::cout<<"Cannot find at "<<__LINE__<<std::endl;
					continue;
				}
				BisectorEdge& tmp = findIter->second;
				findIter = bisecEdgeSet.find(tmp.oppositeSideIndex());
				if (findIter == bisecEdgeSet.end())
				{
					//std::cout<<"Cannot find at "<<__LINE__<<std::endl;
					continue;
				}
				if (findIter->second.sameDirection(tmp))
					findIter->second.reverse();
			}
			BisectorEdge& sideBC = findIter->second;
			rcvList.push_back(&sideBC);
			/*** bisector edge BC on C's side ***/
			findIter = bisecEdgeSet.find(sideCA.firstEndNeighbor());
			if (findIter == bisecEdgeSet.end())
			{
				//std::cout<<"Cannot find at "<<__LINE__<<std::endl;
				continue;
			}
			BisectorEdge& sideCB = findIter->second;
			if (sideCA.reversed)
				sideCB.reverse();
			rcvList.push_back(&sideCB);
			
			// validate
			bool leftBranch = sideAC.rightVacant() && sideAB.leftVacant() && 
								sideBC.leftVacant() && sideCB.rightVacant() && sideCA.leftVacant();
			bool rightBranch = sideAB.rightVacant() && sideAD.leftVacant() && 
								sideDA.rightVacant() && sideDB.leftVacant() && sideBD.rightVacant();
			bool middleBranch = sideAB.rightVacant() || sideAB.leftVacant();
			if (sideBA != bisecEdgeSet.end())
			{
				leftBranch = leftBranch && sideBA->second.rightVacant();
				rightBranch = rightBranch && sideBA->second.leftVacant();
				middleBranch = middleBranch && (sideBA->second.rightVacant() || sideBA->second.leftVacant());
			}
			if ( !leftBranch && !rightBranch && !middleBranch)
				continue;
			// split starts
			rawHoles.push_back(list<Vertex_handle>());
			list<Vertex_handle>& hl = rawHoles.back();
			if (sideAB.visited == BisectorEdge::NO_VISITED)
			{
				BisectorEdge::SPLIT_STATE ss = sideAB.ss;
				// if the principle axis is not broken
				if (ss == BisectorEdge::FULL)
				{
					sideAC.rightSplitAndVisit(&hl);
					sideAB.rightSplitAndVisit(&hl);
					sideAD.leftSplitAndVisit(&hl);
					sideDA.rightSplitAndVisit(&hl);
					sideDB.leftSplitAndVisit(&hl);
					sideBD.rightSplitAndVisit(&hl);
					if (sideBA != bisecEdgeSet.end())
					{
						if (sideBA->second.ss != BisectorEdge::FULL)
						{
							//std::cout<<"Assertion fail!!!"<<std::endl;
							hl.insert(hl.end(), sideBA->second.llist.begin(), sideBA->second.llist.end());
							hl.insert(hl.end(), sideBA->second.rlist.begin(), sideBA->second.rlist.end());
							sideBA->second.visited = BisectorEdge::FULLY_VISITED;
						}
						else
							sideBA->second.rightSplitAndVisit(&hl);
					}
					sideBC.leftSplitAndVisit(&hl);
					sideCB.rightSplitAndVisit(&hl);
					sideCA.leftSplitAndVisit(&hl);
				}
				else if (ss == BisectorEdge::HALF)
				{
					sideAC.rightSplitAndVisit(&hl);
					sideAB.leftSplitAndVisit(&hl);
					if (sideBA != bisecEdgeSet.end())
						sideBA->second.rightSplitAndVisit(&hl);
					sideBC.leftSplitAndVisit(&hl);
					sideCB.rightSplitAndVisit(&hl);
					sideCA.leftSplitAndVisit(&hl);
				}
				else// ss == EMPTY
				{
					sideAB.visited = BisectorEdge::FULLY_VISITED;
					if (sideBA != bisecEdgeSet.end())
						sideBA->second.visited = BisectorEdge::FULLY_VISITED;
				}
			}
			else if (sideAB.visited == BisectorEdge::LEFT_VISITED)	
			{
				sideAB.rightSplitAndVisit(&hl);
				sideAD.leftSplitAndVisit(&hl);
				sideDA.rightSplitAndVisit(&hl);
				sideDB.leftSplitAndVisit(&hl);
				sideBD.rightSplitAndVisit(&hl);
				if (sideBA != bisecEdgeSet.end())
					sideBA->second.leftSplitAndVisit(&hl);
			}
			else if (sideAB.visited == BisectorEdge::RIGHT_VISITED)
			{
				sideAC.rightSplitAndVisit(&hl);
				sideAB.leftSplitAndVisit(&hl);
				if (sideBA != bisecEdgeSet.end())
					sideBA->second.rightSplitAndVisit(&hl);
				sideBC.leftSplitAndVisit(&hl);
				sideCB.rightSplitAndVisit(&hl);
				sideCA.leftSplitAndVisit(&hl);
			}
			// find the centroid of all the involved spheres
			//simplifyHole(&hl);
			if (hl.size() > 2)
			{
				hl.unique();
				vec3 c(0.0, 0.0, 0.0), n, prjc;
				for (list<Vertex_handle>::const_iterator vit = hl.begin(); vit != hl.end(); vit++)
					c += to_geex((*vit)->point());
				c /= hl.size();
				prjc = bdMesh->project_to_mesh(c, n);
				Point_3 o(to_cgal(prjc));
				Plane_3 medialPlane(o, to_cgal(n)-CGAL::ORIGIN);
				vector<Point_3> ahole;
				for (list<Vertex_handle>::const_iterator vit = hl.begin(); vit != hl.end(); vit++)
					ahole.push_back(medialPlane.projection((*vit)->point()));
				//check simple
				//Vector_3 uvz = medialPlane.orthogonal_vector();
				Vector_3 uvz(n.x, n.y, n.z);//point to the outside of mesh
				uvz = uvz / sqrt(uvz.squared_length());
				Vector_3 uvx(ahole[0], ahole[1]);
				uvx = uvx / sqrt(uvx.squared_length());
				Vector_3 uvy = CGAL::cross_product(uvz, uvx); 
				Polygon_2 pgn_2;
				for (unsigned int i = 0; i < ahole.size(); i++)
				{
					pgn_2.push_back(to_uv(uvx, uvy, o, ahole[i]));
				}
				//simplifyPolygon(&ahole);
				if (!pgn_2.is_simple() || ahole.size() < 3)
				{
					rawHoles.pop_back();
					continue;	
				}		
				if (pgn_2.is_clockwise_oriented())
				{
					//pgn_2.reverse_orientation();
					rawHoles.back().reverse();
					vector<Point_3> temp = ahole;
					ahole.clear();
					for (vector<Point_3>::reverse_iterator rit = temp.rbegin(); rit != temp.rend(); rit++)
						ahole.push_back(*rit);
				}
				holes.push_back(ahole);
				simplifyPolygon(&ahole);
				holeArea.push_back(make_pair(std::fabs(pgn_2.area()), holes.size()-1));
			}
			else
				rawHoles.pop_back();
		}// for each bisector edge
		std::cout<<"End detecting holes"<<std::endl;
	}
	void PowerDiagram_3::findTwoLongestEdges(list<Vertex_handle>& rh, 
							list<Vertex_handle>::iterator& maxStart, list<Vertex_handle>::iterator& maxEnd, 
							list<Vertex_handle>::iterator& subMaxStart, list<Vertex_handle>::iterator& subMaxEnd)
	{
		double maxLen = CGAL::squared_distance(rh.back()->point(), rh.front()->point()), 
				subMaxLen = DBL_MIN;
		maxStart = rh.end();
		maxStart--;
		maxEnd = rh.begin();
		list<Vertex_handle>::iterator cur, nxt;
		nxt = cur = rh.begin();
		nxt++;
		for ( ; nxt!=rh.end(); cur++, nxt++)
		{
			double len = CGAL::squared_distance((*cur)->point(), (*nxt)->point());
			if (len > maxLen)
			{
				subMaxStart = maxStart;
				subMaxEnd = maxEnd;
				subMaxLen = maxLen;
				maxStart = cur;
				maxEnd = nxt;
				maxLen = len;
			}
			else if (len > subMaxLen) // subMaxLen < len < maxLen
			{
				subMaxStart = cur;
				subMaxEnd = nxt;
				subMaxLen = len;
			}
		}	
	}
	bool PowerDiagram_3::merge(list<Vertex_handle>& rh0, list<Vertex_handle>::iterator& start0, list<Vertex_handle>::iterator& end0, 
							   list<Vertex_handle>& rh1, list<Vertex_handle>::iterator& start1, list<Vertex_handle>::iterator& end1)
	{
		if ( *start0 == *start1 && *end0 == *end1) //same direction
		{
			rh1.reverse();
			if (*start1 == rh1.front())
				rh0.insert(end0, start1, end1);
			else
			{
				rh0.insert(end0, start1, rh1.end());
				rh0.insert(end0, rh1.begin(), end1);
			}
			rh0.erase(start0);
			return true;
		}
		else if (*start0 == *end1 && *end0 == *start1)//opposite direction
		{
			if (*end1 == rh1.front())
				rh0.insert(end0, end1, start1);
			else
			{
				rh0.insert(end0, end1, rh1.end());
				rh0.insert(end0, rh1.begin(), start1);
			}
			rh0.erase(start0);
			return true;
		}
		return false;
		
	}
	void PowerDiagram_3::mergeHoles()
	{
		std::sort(holeArea.begin(), holeArea.end(), holeAreaCmp);
		//take the first 10% largest 
		double f = 0.2;
		int nbLargeHoles(std::floor(holeArea.size()*f + 0.5));
		int nbMerge = nbLargeHoles*2;
		vector<bool> visited(nbMerge, false);
		int i = 0;
		while (i < nbLargeHoles)
		{
			if (visited[i])
			{
				i++;
				continue;
			}
			int j;
			for (j = i+1; j < nbMerge; j++)
			{
				if (visited[j])
					continue;
				list<Vertex_handle>& rh0 = rawHoles[holeArea[i].second];
				list<Vertex_handle>& rh1 = rawHoles[holeArea[j].second];
				// find the two longest edge in the first polygon
				list<Vertex_handle>::iterator maxStart[2], maxEnd[2], subMaxStart[2], subMaxEnd[2];
				findTwoLongestEdges(rh0, maxStart[0], maxEnd[0], subMaxStart[0], subMaxEnd[0]);
				findTwoLongestEdges(rh1, maxStart[1], maxEnd[1], subMaxStart[1], subMaxEnd[1]);
				bool merged =  merge(rh0, maxStart[0], maxEnd[0], rh1, maxStart[1], maxEnd[1])
							|| merge(rh0, maxStart[0], maxEnd[0], rh1, subMaxStart[1], subMaxEnd[1])
							|| merge(rh0, subMaxStart[0], subMaxEnd[0], rh1, maxStart[1], maxEnd[1])
							|| merge(rh0, subMaxStart[0], subMaxEnd[0], rh1, subMaxStart[1], subMaxEnd[1]);
				if (merged)
				{
					visited[j] = true;
					break;
				}
			}
			//if ( j == nbLargeHoles ) // no match. no need to consider it any more.
			if ( j == nbMerge )
				i++;
		}
		for (i = 0; i < holes.size(); i++)
			holes[i].clear();
		holes.clear();
		for (i = 0; i < holeArea.size(); i++)
		{
			//if ( i < nbLargeHoles && !visited[i] || i >= nbLargeHoles )
			if ( i < nbMerge && !visited[i] || i >= nbMerge )
			{
				vec3 c(0.0, 0.0, 0.0), n, prjc;
				int vi = holeArea[i].second;
				for (list<Vertex_handle>::const_iterator vit = rawHoles[vi].begin(); vit!=rawHoles[vi].end(); vit++)
					c+=to_geex((*vit)->point());
				c /= rawHoles[vi].size();
				prjc = bdMesh->project_to_mesh(c, n);
				Point_3 o = to_cgal(prjc);
				Plane_3 medialPlane(o, to_cgal(n)-CGAL::ORIGIN);
				holes.push_back(vector<Point_3>());
				for (list<Vertex_handle>::const_iterator vit = rawHoles[vi].begin(); vit!=rawHoles[vi].end(); vit++)
					holes.back().push_back(medialPlane.projection((*vit)->point()));
				//check simple
				//Vector_3 nm = medialPlane.orthogonal_vector();
				Vector_3 uvz(n.x, n.y, n.z);//point to outside mesh
				uvz = uvz / sqrt(uvz.squared_length());
				vector<Point_3>& ahole = holes.back();
				Vector_3 uvx(ahole[0], ahole[1]);
				uvx = uvx / sqrt(uvx.squared_length());
				Vector_3 uvy = cross_product(uvz, uvx);
				Polygon_2 pgn_2;
				for (unsigned int j = 0; j < ahole.size(); j++)
					pgn_2.push_back(to_uv(uvx, uvy, o, ahole[j]));
				if (!pgn_2.is_simple() || holes.back().size() < 3)
				{
					holes.pop_back();	
					continue;
				}
				if (pgn_2.is_clockwise_oriented())
				{
					//pgn_2.reverse_orientation();
					vector<Point_3> temp = ahole;
					ahole.clear();
					for (vector<Point_3>::reverse_iterator rit = temp.rbegin(); rit != temp.rend(); rit++)
						ahole.push_back(*rit);
				}
				simplifyPolygon(&holes.back());
			}
		}
	}
	void PowerDiagram_3::single_clip(int group_id)
	{
		Plane_3& p = tangentPlanes[group_id];
		const VGroup& vg = groupList[group_id];
		//vector<Segment_3>& clippedPolygon = clippedPolygons[group_id];
		vector<Bisector>& clippedPolygon = clippedPolygons[group_id];
		clippedPolygon.clear();
		vector<int> dualVertGroupId;
		vector<Vertex_handle> dualVertHandle;
		for (unsigned int i = 0; i < vg.size(); i++)
		{
			const vector<vector<Object>>& vCell = vg[i]->vCell;
			for (unsigned int j = 0; j < vCell.size(); j++)
			{
				// compute the intersection with all bisectors
				if (vg[i]->group_id == vg[i]->dualVert[j])
					continue;
				const vector<Object>& face = vCell[j];
				vector<Point_3> intersectPoints; //we assume the face of the cell is convex, so number of intersection points must be 0, 1, 2
				for (unsigned int k = 0; k < face.size(); k++)
				{
					if (const Segment_3* s = object_cast<Segment_3>(&face[k]))
					{
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
				if (intersectPoints.size() > 2)
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
				assert(intersectPoints.size() <= 2);
				if (intersectPoints.size()==2)	
				{
					//clippedRegion->push_back(Segment_3(intersectPoints[0], intersectPoints[1]));
					vg[i]->tangentIntersection.push_back(Segment_3(intersectPoints[0], intersectPoints[1]));
					clippedPolygon.push_back(Bisector(vg[i]->dualVert[j], Segment_3(intersectPoints[0], intersectPoints[1]), vg[i]/*, sqDist*/));
					dualVertGroupId.push_back(vg[i]->dualVert[j]);
					dualVertHandle.push_back(vg[i]->dualVertHandle[j]);
				}
			}	
			vg[i]->dualVert.assign(dualVertGroupId.begin(), dualVertGroupId.end());
			vg[i]->dualVertHandle.assign(dualVertHandle.begin(), dualVertHandle.end());
			vg[i]->vCell.clear();
		}

		/* Compute the close loop(clipped polygon) */
		list<Bisector> loop;
		vector<bool> determined(clippedPolygon.size(), false);
		bool allDetermined = false;
		list<Bisector>::iterator lit = loop.insert(loop.end(),clippedPolygon[0]);
		determined[0] = true;
		double newThreshold = 1e-14;
		//std::cout << "Entering the merging process\n";
		int iterTimes = 0, iterLim = clippedPolygon.size()*clippedPolygon.size();
		while (!allDetermined)
		{
			allDetermined = true;
			if (iterTimes > iterLim)
			{
				//std::cout<<"Stuck in single clip loop!\n";
				newThreshold = 2*newThreshold;
				iterTimes = 0;
			}
			for (lit = loop.begin(); lit!=loop.end(); lit++)
			{
				unsigned int n;
				for ( n = 0; n < clippedPolygon.size(); n++)
				{
					allDetermined &= determined[n];
					if (determined[n])
						continue;
					const Bisector& current = clippedPolygon[n];
					if (approxEqual(current.source(), lit->source(), newThreshold) && !approxEqual(current.target(), lit->target(), newThreshold))
						lit = loop.insert(lit, current.opposite());
					else if ( approxEqual(current.source(), lit->target(), newThreshold) && !approxEqual(current.target(), lit->source(), newThreshold))
					{	
						lit++;
						lit = loop.insert(lit, current);
					}
					else if (approxEqual(current.target(), lit->source(), newThreshold) && !approxEqual(current.source(), lit->target(), newThreshold))
						lit = loop.insert(lit, current);
					else if (approxEqual(current.target(), lit->target(), newThreshold ) && !approxEqual(current.source(), lit->source(), newThreshold) )
					{
						lit++;
						lit = loop.insert(lit, current.opposite());
					}
					else //this edge does not match
					{
						iterTimes++;
						continue;
					}
					determined[n] = true;
					iterTimes = 0;
					newThreshold = 1.0e-14;
					break;
				}
			}

		}
		//std::cout<<"Exiting the merging process\n";
		clippedPolygon.resize(loop.size());
		int idx = 0;
		for (list<Bisector>::const_iterator cpyit = loop.begin(); cpyit != loop.end(); cpyit++, idx++)
			clippedPolygon[idx] = *cpyit;
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
					//vg[i]->tanIntersecDualVert.push_back(vg[i]->dualVert[j]);
				}
			}			
		}
	}
	void PowerDiagram_3::clipped_tangent_plane()
	{
		//clipped_tangent_plane_.resize(groupCount);
		//if (use_omp)
		//{
		//	int nbThreads = std::min<int>(int(groupList.size()), omp_get_num_procs())*2;
		//	omp_set_num_threads(nbThreads);
		//}
//#pragma omp parallel for if(use_omp)	
#ifdef CILK
		cilk_for (int i = 0 ; i < (int)groupList.size(); i++)
#else
		for (int i = 0; i < int(groupList.size()); i++)
#endif
			single_clip(/*tangentPlanes[i], */ i /*, &clipped_tangent_plane_[i]*/);
	}
	Transformation_3 PowerDiagram_3::setConstraints(unsigned int group_id, vector<Constraint> *constraintList, const vec3& nm)
	{
		constraintList->clear();
		const VGroup& vg = groupList[group_id];
		// compute the matrix that transforms this group of vertices on 2D plane
		vec3 c(0.0, 0.0, 0.0);
		for (unsigned int i = 0; i < vg.size(); i++)
			c += to_geex(vg[i]->point());
		c /= vg.size();
		Transformation_3 move_to_origin(CGAL::TRANSLATION, Vector_3(-c.x, -c.y, -c.z));
		//compute the rotation
		double cos_theta = nm.z;
		double sin_theta = sqrt(nm.x*nm.x + nm.y*nm.y);
		Transformation_3 rot(CGAL::IDENTITY);
		if (nm.x != 0.0 || nm.y != 0.0)
		{
			double cos_phi = nm.x/sin_theta;
			double sin_phi = nm.y/sin_theta;
			rot = Transformation_3(cos_theta*cos_phi, cos_theta*sin_phi, -sin_theta,
								   -sin_phi,		   cos_phi,			  0.0,
								   sin_theta*cos_phi, sin_theta*sin_phi, cos_theta);
		}
		Transformation_3 t = rot*move_to_origin;
		for (unsigned int i = 0; i < vg.size(); i++)
		{
			//get the involved bisectors for each vertex
			vector<Segment_3>& bisec = vg[i]->tangentIntersection;
			//vector<int>& dualGroupId = vg[i]->tanIntersecDualVert;
			Point_3 p3d = t(vg[i]->point());
			Point_2 p(p3d.x(), p3d.y());
			for (unsigned int j = 0; j < bisec.size(); j++)
			{
				//if (dualGroupId[j] == vg[i]->group_id)
				//	continue;
				Segment_3 s3d = bisec[j].transform(t);
				Point_2 src(s3d.source().x(), s3d.source().y()), tgt(s3d.target().x(), s3d.target().y());//we simply omit the third coordinate of the transformed segments
				constraintList->push_back(Constraint(p, Segment_2(src, tgt), vg[i]->radius));
			}
		}
		return t;
	}
	void PowerDiagram_3::setConstraints(unsigned int group_id, vector<Constraint> *constraintList, 
													const Vector_3& uvx, const Vector_3& uvy, Point_3 *rotCent)
	{
		constraintList->clear();
		const VGroup& vg = groupList[group_id];
		// compute the matrix that transforms this group of vertices on 2D plane
		vec3 c(0.0, 0.0, 0.0);
		for (unsigned int i = 0; i < vg.size(); i++)
			c += to_geex(vg[i]->point());
		c /= vg.size();
		*rotCent = to_cgal(c);
		for (unsigned int i = 0; i < vg.size(); i++)
		{
			//get the involved bisectors for each vertex
			vector<Segment_3>& bisec = vg[i]->tangentIntersection;
			Point_2 p = to_uv(uvx, uvy, *rotCent, vg[i]->point());
			for (unsigned int j = 0; j < bisec.size(); j++)
			{
				Point_2 src = to_uv(uvx, uvy, *rotCent, bisec[j].source());
				Point_2 tgt = to_uv(uvx, uvy, *rotCent, bisec[j].target());//we simply omit the third coordinate of the transformed segments
				constraintList->push_back(Constraint(p, Segment_2(src, tgt), vg[i]->radius));
			}
		}
	}
	void PowerDiagram_3::setConstraints(unsigned int group_id, vector<Constraint> *constraintList, const Vector_3& uvx, 
										const Vector_3& uvy, Point_3 *rotCent, const vector<Segment_3>& pgn)
	{
		constraintList->clear();
		const VGroup& vg = groupList[group_id];
		// compute the matrix that transforms this group of vertices on 2D plane
		vec3 c(0.0, 0.0, 0.0);
		for (unsigned int i = 0; i < vg.size(); i++)
			c += to_geex(vg[i]->point());
		c /= vg.size();
		*rotCent = to_cgal(c);
		vector<Segment_2> outer;
		for (unsigned int i = 0; i < vg.size(); i++)
		{
			vector<Segment_3>& bisec = vg[i]->tangentIntersection;
			for (unsigned int j = 0; j < bisec.size(); j++)
			{
				Point_2 src = to_uv(uvx, uvy, *rotCent, bisec[j].source());
				Point_2 tgt = to_uv(uvx, uvy, *rotCent, bisec[j].target());
				outer.push_back(Segment_2(src, tgt));
			}
		}
		vector<Point_2> inner;
		for (unsigned int i = 0; i < pgn.size(); i++)
			inner.push_back(to_uv(uvx, uvy, *rotCent, pgn[i].source()));
		for (unsigned int i = 0; i < inner.size(); i++)
		{
			// find the two most nearest edges
			//double minDist = DBL_MAX, subMinDist;
			//unsigned int minDistIdx, subMinDistIdx;
			for (unsigned int j = 0; j < outer.size(); j++)
			{
			//	double d = CGAL::squared_distance(inner[i], outer[j]);
			//	if ( d < minDist )
			//	{
			//		subMinDist = minDist;
			//		subMinDistIdx = minDistIdx;
			//		minDist = d;
			//		minDistIdx = j;
			//	}
			//	else if ( d < subMinDist )
			//	{
			//		subMinDist = d;
			//		subMinDistIdx = j;
			//	}
				constraintList->push_back(Constraint(inner[i], outer[j]));
			}
			//constraintList->push_back(Constraint(inner[i], outer[minDistIdx]));
			//constraintList->push_back(Constraint(inner[i], outer[subMinDistIdx]));
		}
	}
	void PowerDiagram_3::clusterHoleDetect()
	{
		holeSkeleton.clear();
		for_each(holeBorders.begin(), holeBorders.end(), mem_fun_ref(&(vector<Vertex_handle>::clear)));
		holeBorders.clear();
		//const int samplenb = 4;
		const int dim = 3;
		//int nbedges = 0;
		//double epsilon = DBL_MIN;
		double avgr = 0.0;//to determine the sample rate and neigbhor threshold used in DBScan
		map<Vertex_handle, vector<bool>> dualVertHandleMask;
		for (Finite_vertices_iterator vit = finite_vertices_begin(); vit != finite_vertices_end(); vit++)
		{
			if (vit->group_id < 0)
				continue;
			dualVertHandleMask[vit].resize(vit->dualVertHandle.size(), false);
			vit->included = false;
			avgr += vit->radius;
		}
		avgr /= number_of_vertices();
		int nbPnts = 0;//number of points sampled and used in DBScan
		double sample_rate = 0.2*avgr;
		//double maxs = 9.0*9.0, mins = 3.0*3.0;
		double srmin = 3.5*3.5*avgr*avgr, srmax = 16.0*16.0*avgr*avgr;
		// count expected bisect segments and compute number of sampling points
		for (Finite_vertices_iterator vit = finite_vertices_begin(); vit != finite_vertices_end(); vit++)
		{
			if (vit->group_id < 0)
				continue;
			//double srmin = mins*vit->radius*vit->radius, srmax = maxs*vit->radius*vit->radius;
			const vector<Segment_3>& bisec = vit->tangentIntersection;
			for (unsigned int i = 0; i < bisec.size(); i++)
			{
				Point_3 mp = CGAL::midpoint(bisec[i].source(), bisec[i].target());
				double d = CGAL::squared_distance(mp, vit->point());
				if ( d >= srmin && d <= srmax )
				{
					double segLen = sqrt(bisec[i].squared_length());
					int sn(segLen/sample_rate);
					//if (sn == 0)
					//	sn = 1;
					sn += 2;//the end points of the segments must be sampled
					nbPnts += sn;
					dualVertHandleMask[vit][i] = true;
					//nbedges++;
				}
			} 
		} 
		//int nbPnts = nbedges*samplenb;
		double *ptsCoord = new double[nbPnts*dim];
		double *ptr = ptsCoord;
		for (Finite_vertices_iterator vit = finite_vertices_begin(); vit != finite_vertices_end(); vit++)
		{
			if (vit->group_id < 0/* || vit->included*/)
				continue;
			//double srmin = mins*vit->radius*vit->radius, srmax = maxs*vit->radius*vit->radius;
			const vector<Segment_3>& bisec = vit->tangentIntersection;
			for (vector<Segment_3>::const_iterator sit = bisec.begin(); sit != bisec.end(); sit++)
			{
				Point_3 mp = CGAL::midpoint(sit->source(), sit->target());
				double d = CGAL::squared_distance(mp, vit->point());
				if ( d >= srmin && d <= srmax )
				{
					double segLen = sqrt(sit->squared_length());
					int sn(segLen/sample_rate);
					//if (sn == 0)
					//	sn = 1;//the middle point is sampled
					sn += 2; // end points must be sampled
					for (int i = 0; i < sn; i++)
					{
						Point_3 sp = CGAL::ORIGIN + ((sn-1-i)*(sit->source()-CGAL::ORIGIN) + i*(sit->target()-CGAL::ORIGIN))/(sn-1);
						vec3 dummyn;
						vec3 prjsp = bdMesh->project_to_mesh(to_geex(sp), dummyn);
						holeSkeleton.push_back(SkeletonPoint(to_cgal(prjsp), vit));
						//vit->included = true;
						*ptr++ = prjsp.x; *ptr++ = prjsp.y; *ptr++ = prjsp.z;
					}
				}
			} 
		}
		// do clustering
		int min_neigbhor_n = 0, max_neighbor_n = 0;
		double epsilon = 0.732*avgr;
		scanRadius = epsilon;
		std::cout<<"Scan radius "<<scanRadius<<std::endl;
		DBSCAN dbscan(ptsCoord, nbPnts, dim, 0.8*epsilon, 0.60*epsilon, 6, min_neigbhor_n, max_neighbor_n);
		std::cout<<"min neighbor number = "<<min_neigbhor_n<<", max neighbor number "<<max_neighbor_n<<std::endl;
		//assign cluster id
		//count number of groups in each cluster
		vector<set<int>> clusterGroupCount(dbscan.label_n_);
		vector<bool> clusterValid(dbscan.label_n_, true);
		for (int i = 0; i < nbPnts; i++)
		{
			if (dbscan.internal_label_[i] >= 0)
				clusterGroupCount[dbscan.internal_label_[i]].insert(holeSkeleton[i].v->group_id);
		}
		for (unsigned int i = 0; i < clusterGroupCount.size(); i++)
			if (clusterGroupCount[i].size() <= 1)
				clusterValid[i] = false;
		//assert(holeSkeleton.size() == nbPnts);
		for (unsigned int i = 0; i < holeSkeleton.size(); i++)
		{
			if (dbscan.internal_label_[i] < 0 || !clusterValid[dbscan.internal_label_[i]])
			{
				holeSkeleton[i].cluster_id = -1;
				//std::cout<<"noise point!\n";
			}
			else
				holeSkeleton[i].cluster_id = dbscan.internal_label_[i];
		}
		delete []ptsCoord;
		// outline bounding spheres 
		vector<set<Vertex_handle>> hb(dbscan.label_n_);//hole borders
		for (unsigned int i = 0; i < holeSkeleton.size(); i++)
		{
			int cid = dbscan.internal_label_[i];
			if (clusterValid[cid] && cid >= 0)
				hb[cid].insert(holeSkeleton[i].v);
		}
		// further filter, eliminate possible isolated points
		for (unsigned int i = 0; i < hb.size(); i++)
		{
			set<Vertex_handle>::iterator vit0 = hb[i].begin(); 
			while ( vit0 != hb[i].end() )
			{
				set<Vertex_handle>::iterator vit1;
				for ( vit1 = hb[i].begin(); vit1 != hb[i].end(); vit1++)
					if (*vit0 != *vit1)
					{
						Cell_handle dummycell;
						int dummyi, dummyj;
						if (is_edge(*vit0, *vit1, dummycell, dummyi, dummyj) && (*vit0)->group_id != (*vit1)->group_id)
							break;
					}
				if (vit1 == hb[i].end())
				{
					set<Vertex_handle>::iterator tmp = vit0;
					vit0++;
					hb[i].erase(tmp);
				}
				else
					vit0++;
			}
		}
		// add the dual vertices if they are not included in any hole
		// mark vertices in holes
		//for (unsigned int i = 0; i < hb.size(); i++)
		//{
		//	for (set<Vertex_handle>::iterator it = hb[i].begin(); it != hb[i].end(); it++)
		//		(*it)->included = true;
		//}
		// check dual vertices
		//for (unsigned int i = 0; i < hb.size(); i++)
		//{
		//	if (hb[i].size() < 3)
		//		continue;
		//	for (set<Vertex_handle>::iterator it = hb[i].begin(); it != hb[i].end(); it++)
		//	{
		//		vector<Vertex_handle>& dualVertHandle = (*it)->dualVertHandle;
		//		for (unsigned int j = 0; j < dualVertHandle.size(); j++)
		//			if (!dualVertHandle[j]->included && dualVertHandleMask[*it][j])
		//			{
		//				hb[i].insert(dualVertHandle[j]);
		//				dualVertHandle[j]->included = true;
		//			}
		//	}
		//}
		for (unsigned int i = 0; i < hb.size(); i++)
			if (hb[i].size() >= 8)
			{
				holeBorders.push_back(vector<Vertex_handle>());
				for (set<Vertex_handle>::const_iterator it = hb[i].begin(); it != hb[i].end(); it++)
					holeBorders.back().push_back(*it);
			}
	}
}