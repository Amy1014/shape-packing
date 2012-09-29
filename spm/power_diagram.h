/*
 *  _____ _____ ________  _
 * /  __//  __//  __/\  \//
 * | |  _|  \  |  \   \  / 
 * | |_//|  /_ |  /_  /  \ 
 * \____\\____\\____\/__/\\
 *
 * Graphics Environment for EXperimentations.
 *  Copyright (C) 2006 INRIA - Project ALICE
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: 
 *
 *     ALICE Project - INRIA
 *     INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 *  Note that the GNU General Public License does not permit incorporating
 *  the Software into proprietary programs. 
 */

#ifndef __POWERDIAGRAM__
#define __POWERDIAGRAM__

#include <cfloat>
#include <vector>
#include <list>
#include <iterator>
#include "spm_cgal.h"
#include "meshes.h"
#include "Containment.h"
#include <GL/gl.h>

using std::vector;
using std::list;
using std::back_insert_iterator;
using Geex::Numeric::is_nan;
using CGAL::object_cast;
using CGAL::intersection;
using CGAL::do_intersect;
using CGAL::squared_distance;
using CGAL::to_double;
//using CGAL::collinear;

namespace Geex {

	struct GroupVertex;

//______________________________________________________________________________________

    class PowerDiagram_3 : public RegularBase_3 
	{

		typedef RegularBase_3 baseclass ;
		typedef vector<PowerDiagram_3::Vertex_handle> VGroup;
		typedef PowerDiagram_3::Edge	Edge;
		typedef PowerDiagram_3::Cell_circulator Cell_circulator;
		typedef PowerDiagram_3::Triangle	Triangle;

	public:

		typedef PowerDiagram_3::Vertex_handle Vertex_handle;
		typedef PowerDiagram_3::Cell_handle		Cell_handle;

		PowerDiagram_3(TriMesh* m);
		/* insert methods */
		// insert a single weighed point
		void clear();
		void begin_insert() ;
		Vertex_handle insert(const vec3& p, const double w) ;
		Vertex_handle insert(const Weighted_point& wp);
		void end_insert(bool redraw = true) ;
		// insert a group of list
		void insert_weighted_points( int nb, double* x, const double * w );
		void insert_weighted_points( const vector<GroupVertex>& vlist );
		void insert_grouped_points( const vector<Weighted_point>& vlist);
		void insert_grouped_points( const vector<Weighted_point>& vlist, const Plane_3& tangentPlane)
		{
			insert_grouped_points(vlist);
			tangentPlanes.push_back(tangentPlane);
		}
		// access methods
        const vector<Vertex_handle>& getAllVertices() const { return all_vertices_ ; }
        size_t nb_vertices() const { return all_vertices_.size() ; }
		unsigned int getNbGroups() const { return groupCount;}
		const TriMesh& get_boundary_mesh()	{ return *bdMesh; }
		const vector<Cell_handle>& getAllCells() const { return all_cells_; }
		const vector<vector<Segment_3>>& getBisectors()	const {return clipped_tangent_plane_;}
		/* computation */
		// projection of one group of weighed points to the mesh boundary
		//vec3 project_to_mesh(const vec3& pt, vec3& norm) { return boundary_.project_to_mesh(pt, norm); } 
		void clipped_tangent_plane();//clip the tangent plane with bisectors
	private:
		void single_clip(const Plane_3& p, const VGroup& vg/*, vector<Segment_3> *clippedRegion*/);
		void accurate_single_clip(const Plane_3& p, const VGroup& vg/*, vector<Segment> *clippedRegion*/);
		inline bool approxEqual(const Point_3& p, const Point_3& q, double td)
		{
			return squared_distance(p, q) <= td;
		}
		inline vector<Point_3>::iterator approxFind(vector<Point_3>::iterator first, vector<Point_3>::iterator last, const Point_3& p )
		{
			for (vector<Point_3>::iterator it = first; it != last; it++)
				if (approxEqual(*it, p, threshold))
					return it;
			return last;
		}
		inline Ext_point_3 to_ext(const Point_3& p)
		{
			return Ext_point_3(NT(p.x()), NT(p.y()), NT(p.z()));
		}
		inline Point_3 inv_to_ext(const Ext_point_3& p)
		{
			return Point_3(to_double(p.x()), to_double(p.y()), to_double(p.z()));
		}
		inline Ext_segment_3 to_ext(const Segment_3& s)
		{
			return Ext_segment_3(to_ext(s.source()), to_ext(s.target()));
		}
		inline Segment_3 inv_to_ext(const Ext_segment_3& s)
		{
			return Segment_3(inv_to_ext(s.source()), inv_to_ext(s.target()));
		}
		inline Ext_vector_3 to_ext(const Vector_3& v)
		{
			return Ext_vector_3(NT(v.x()), NT(v.y()), NT(v.z()));
		}
		inline Ext_ray_3 to_ext(const Ray_3& r)
		{
			return Ext_ray_3(to_ext(r.source()), to_ext(r.to_vector()));
		}
		inline Ext_plane_3 to_ext(const Plane_3& p)
		{
			return Ext_plane_3(NT(p.a()), NT(p.b()), NT(p.c()), NT(p.d()));
		}
		inline Ext_triangle to_ext(const Triangle& t)
		{
			return Ext_triangle(to_ext(t.vertex(0)), to_ext(t.vertex(1)), to_ext(t.vertex(2)));
		}
		bool collinear(const Point_3& P, const Point_3& Q, const Point_3& R )
		{
			vec3 p = to_geex(P), q = to_geex(Q), r = to_geex(R);
			vec3 cp = cross(q - p, r - p);
			return (cp.length2() <= 1.0e-12);
		}
		bool coplane(Cell_handle c)
		{
			// for numerical stability, check co-plane in advance
			double num_x, num_y, num_z, den;
			Vertex_handle p = c->vertex(0);
			Vertex_handle q = c->vertex(1);
			Vertex_handle r = c->vertex(2);
			Vertex_handle s = c->vertex(3);
			if ( p->group_id == q->group_id && p->group_id == r->group_id && p->group_id == s->group_id 
			   && q->group_id == r->group_id && q->group_id == s->group_id 
			   && r->group_id == s->group_id)
			   return true;
			CGAL::determinants_for_weighted_circumcenterC3(p->point().x(), p->point().y(), p->point().z(), p->point().weight(),
				q->point().x(), q->point().y(), q->point().z(), q->point().weight(),
				r->point().x(), r->point().y(), r->point().z(), r->point().weight(),
				s->point().x(), s->point().y(), s->point().z(), s->point().weight(),
				num_x,  num_y, num_z,den);
			//return (CGAL::is_zero(den));
			return (std::fabs(den) <= threshold);

		}
	private:
		TriMesh* bdMesh;
		vector<Vertex_handle> all_vertices_;
		vector<Cell_handle> all_cells_;
		vector<vector<Segment_3>> clipped_tangent_plane_;
		vector<VGroup> groupList;
		vector<Plane_3> tangentPlanes;
		bool open;
		unsigned int groupCount;
		int vertIdx;//the group ID and vertex index currently being inserted
		double threshold;
    } ;
	struct GroupVertex
	{
		int group_id;
		RegularBase_3::Vertex_handle vert;
		double radius;
	};
}

#endif
