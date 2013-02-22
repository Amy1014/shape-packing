#ifndef __POWERDIAGRAM__
#define __POWERDIAGRAM__


#include <cfloat>
#include <cmath>
#include <vector>
#include <list>
#include <map>
#include <utility>
#include <iterator>
#include "spm_cgal.h"
#include "meshes.h"
#include "Containment.h"
#include <exception>
#include <GL/gl.h>
//#include <boost/shared_ptr.hpp>
//#include "dbscan/db_scan.h"
//#include "dbscan/distance.h"
#include "DBSCAN.h"
#include <CGAL/create_offset_polygons_2.h>

#ifdef CILK
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/reducer_opadd.h>
#include <cilk/reducer_min.h>
#endif

using std::vector;
using std::list;
using std::map;
using std::pair;
using std::make_pair;
using std::back_insert_iterator;
using Geex::Numeric::is_nan;
using CGAL::object_cast;
using CGAL::intersection;
using CGAL::do_intersect;
using CGAL::squared_distance;
using CGAL::to_double;
using CGAL::sqrt;
//using CGAL::collinear;

namespace Geex {

	struct Bisector;
	class BisectorEdge;
	//class BisectorEdgeCmp;
	struct GroupVertex;
	class BisectorEdgeCmp
	{
	public:
		inline bool operator()(const pair<int, int>& p0, const pair<int, int>& p1) const
		{
			if (p0.first != p1.first)
				return (p0.first < p1.first);
			else
				return (p0.second < p1.second);
		}
	};
	inline bool holeAreaCmp(const pair<double, unsigned int>& ha0, const pair<double, unsigned int>& ha1)
	{
		return (ha0.first > ha1.first);
	}
//______________________________________________________________________________________

    class PowerDiagram_3 : public RegularBase_3 
	{

		typedef RegularBase_3 baseclass ;	
		
		typedef PowerDiagram_3::Cell_circulator Cell_circulator;
		typedef PowerDiagram_3::Triangle	Triangle;
		typedef map<pair<int, int>, BisectorEdge, BisectorEdgeCmp>::iterator BisectorEdgeIterator;
		typedef boost::shared_ptr<Geex::Polygon_2> PolygonPrt;
		typedef std::vector<PolygonPrt> PolygonPtrVector;	
	public:
		typedef PowerDiagram_3::Edge	Edge;
		typedef vector<PowerDiagram_3::Vertex_handle> VGroup;
		typedef PowerDiagram_3::Vertex_handle Vertex_handle;
		typedef PowerDiagram_3::Cell_handle		Cell_handle;
		typedef pair<list<int>, list<Vertex_handle>> Hole;

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
		/* access methods */
        const vector<Vertex_handle>& getAllVertices() const { return all_vertices_ ; }
		const vector<VGroup>& getGroupList() const { return groupList; }
		const VGroup& getGroup(unsigned int i) const { return groupList[i]; }
        size_t nb_vertices() const { return all_vertices_.size() ; }
		unsigned int getNbGroups() const { return groupCount;}
		const TriMesh& get_boundary_mesh()	{ return *bdMesh; }
		const vector<vector<Bisector>>& getClippedPolygons() const { return clippedPolygons; }

		/* computation */
		// projection of one group of weighed points to the mesh boundary
		void clipped_tangent_plane();//clip the tangent plane with bisectors
		Transformation_3 setConstraints(unsigned int group_id, vector<Constraint> *constraintList, const vec3& nm);
		void setConstraints(unsigned int group_id, vector<Constraint> *constraintList,  
										const Vector_3& uvx, const Vector_3& uvy, Point_3 *rotCent);
		void setConstraints(unsigned int group_id, vector<Constraint> *constraintList, const Vector_3& uvx, 
										const Vector_3& uvy, Point_3 *rotCent, const vector<Segment_3>& pgn);
		void detectHoles();//find all the holes of the power diagram formed by polygons
		void mergeHoles(); // merge the small holes detected by detect holes
		const vector<vector<Point_3>>& getHoles() const { return holes; }
		void preFinalDetectHoles(const vector<vec3>&);//hole detection with sphere bounding
		void clusterHoleDetect();
		//debug
		const vector<vector<Weighted_point>>& getUndetectedHoles() const { return sholes; }
	private:
		void preDetectHoles();//build the data structure to detect holes
		
		void simplifyHole(list<Vertex_handle>* hole);//remove the collinear points in a raw hole
		void simplifyPolygon(vector<Point_3> *pgn);//eliminate the collinear vertices in a polygon
		//find the two longest edges in a hole, only used in mergeHoles routine
		void findTwoLongestEdges(list<Vertex_handle>& rh, 
								list<Vertex_handle>::iterator& maxStart, list<Vertex_handle>::iterator& maxEnd,
								list<Vertex_handle>::iterator& subMaxStart, list<Vertex_handle>::iterator& subMaxEnd);
		bool merge(list<Vertex_handle>& rh0, list<Vertex_handle>::iterator& start0, list<Vertex_handle>::iterator& end0,
				   list<Vertex_handle>& rh1, list<Vertex_handle>::iterator& start1, list<Vertex_handle>::iterator& end1);
		//void preDetectHoles();//build the data structure to detect holes
		inline void recoverBisectorEdges(const list<BisectorEdge*>& rcvList);
		void single_clip(int group_id/* const VGroup& vg*/);
		void clip_polygons(int group_id);
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
		inline Vector_3 inv_to_ext(const Ext_vector_3& v)
		{
			 return Vector_3(to_double(v.x()), to_double(v.y()), to_double(v.z()));
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
		inline Ext_triangle_3 to_ext(const Triangle& t)
		{
			return Ext_triangle_3(to_ext(t.vertex(0)), to_ext(t.vertex(1)), to_ext(t.vertex(2)));
		}
		bool collinear(const Point_3& P, const Point_3& Q, const Point_3& R )
		{
			vec3 p = to_geex(P), q = to_geex(Q), r = to_geex(R);
			vec3 cp = cross(q - p, r - p);
			return (cp.length2() <= 1e-13);// for egg, 1.0e-16
			//return CGAL::is_zero(cp.length2());
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
			return (std::fabs(den) <= 1e-13/*threshold*/);// for egg, 1.0e-16

		}
		inline bool sphereContain(const Point_3& bigCent, double bigSqRadius, const Point_3& smallCent, double smallSqRadius)
		{
			// square of sphere center distance
			double cendist = CGAL::squared_distance(bigCent, smallCent);
			// square of difference of radius, || r1 - r2 ||^2
			double sqdiff = bigSqRadius - 2*sqrt(bigSqRadius*smallSqRadius)+smallSqRadius;
			if (cendist <= sqdiff)
				return true;
			else
				return false;
		}
	private:
		TriMesh* bdMesh;
		vector<Vertex_handle> all_vertices_;
		//vector<vector<Segment_3>> clipped_tangent_plane_;
		vector<VGroup> groupList;
		vector<Plane_3> tangentPlanes;
		// vector<vector<Segment_3>> clippedPolygons;// for each group, the polygon generated by clipping power cells with tangent planes
		vector<vector<Bisector>> clippedPolygons;
		map<pair<int, int>, BisectorEdge, BisectorEdgeCmp> bisecEdgeSet;//contain all the bisector segments between polygons in this power diagram
		vector<vector<Point_3>> holes;
		vector<vector<Weighted_point>> sholes;
		vector<list<Vertex_handle>> rawHoles;
		vector<pair<double, unsigned int>> holeArea;
		bool open;
		unsigned int groupCount;
		int vertIdx;//the group ID and vertex index currently being inserted
		double threshold;
		int use_omp;
		//////////////////////////////////////////////////////////////////////////
		struct DistBisectorTriplet
		{
			double distance;
			Vertex_handle v;
			int whichBisec;
			DistBisectorTriplet(double d, Vertex_handle vh, int bisecIdx) : distance(d), v(vh), whichBisec(bisecIdx) {}
			Segment_3& getBisector() { return v->tangentIntersection[whichBisec];}
		};
		struct DistBisectorCmp
		{
		public:
			inline bool operator()(const DistBisectorTriplet& p0, const DistBisectorTriplet& p1)
			{
				return p0.distance > p1.distance;
			}
		};
		//debug
		vector<Segment_3> undetectedHoles;
		public:
			struct SkeletonPoint
			{
				Point_3 p;
				Vertex_handle v;
				int cluster_id;
				SkeletonPoint(const Point_3& pos, Vertex_handle vh) : p(pos), v(vh), cluster_id(-1) {}
			};
		vector<SkeletonPoint> holeSkeleton;
		vector<vector<Vertex_handle>> holeBorders;
		double scanRadius;
    } ;
	struct GroupVertex
	{
		int group_id;
		PowerDiagram_3::Vertex_handle vert;
		double radius;
	};

	/////////////////////Bisectors////////////////////////////
	struct Bisector
	{
		typedef PowerDiagram_3::Vertex_handle Vertex_handle;

	public:
		Bisector() : dualVertGroupId(-1){}
		Bisector(int groupId, const Segment_3& s, const Vertex_handle& vh) : 
									dualVertGroupId(groupId), seg(s), sphere(vh) {}
		//Bisector(int groupId, const Segment_3& s, const Vertex_handle& vh, double sd)
		//				:dualVertGroupId(groupId), seg(s), sphere(vh), sqDist(sd) {}
		Point_3 source() const { return seg.source(); }
		Point_3 target() const { return seg.target(); }
		Bisector opposite() const { return Bisector(dualVertGroupId, seg.opposite(), sphere); }
		const Segment_3& segment() const { return seg; }
		int dualVertGroupId;
		Segment_3 seg;
		Vertex_handle sphere;//the sphere that form this bisector
		//double sqDist;//distance from the center of sphere to bisector plane
	};



	/*
					\					 /
					 \   A(base.frist)	/
		C(end.frist)  \________________/  D(end.second)
					  /				   \
					 /	 B(base.second)	\
					/					 \
	*/
	class BisectorEdge // clipped segments between two polygons formed by intersecting tangent planes with power cells
	{
		friend class PowerDiagram_3; 
		friend class BisectorEdgeCmp;
		typedef PowerDiagram_3::Vertex_handle Vertex_handle;

	private:// only for internal use in power diagram class
		typedef enum{EMPTY, HALF, FULL} SPLIT_STATE;//return state defined for split function
		typedef enum{NO_VISITED, FULLY_VISITED, LEFT_VISITED, RIGHT_VISITED} VISIT_STATE;

		BisectorEdge(int b0, int b1, int e0, int e1): 
					base(b0, b1), end(e0, e1), visited(NO_VISITED), reversed(false){}
		BisectorEdge(int b0, int b1, int e0, int e1, vector<Bisector>* bl): 
					base(b0, b1), end(e0, e1), visited(NO_VISITED), reversed(false), bisecList(bl) {}

		inline bool sameDirection(const BisectorEdge& be) const
		{
			return (end.first == be.end.first || end.second == be.end.second);
		}
		void reverse()
		{
			std::swap(end.first, end.second);
			//segIdx.reverse();
			reversed = !reversed;
			if (visited == LEFT_VISITED)
				visited = RIGHT_VISITED;
			else if (visited == RIGHT_VISITED)
				visited = LEFT_VISITED;
			llist.reverse();
			rlist.reverse();
			std::swap(llist, rlist);
		}
		inline pair<int, int> oppositeSideIndex() const
		{
			return pair<int, int>(base.second, base.first);
		}
		inline pair<int, int> firstEndNeighbor() const 
		{
			return pair<int, int>(base.first, end.first);
		}
		inline pair<int, int> secondEndNeighbor() const
		{
			return pair<int, int>(base.first, end.second);
		}
		inline bool rightVacant() const
		{
			return (visited == NO_VISITED || visited == LEFT_VISITED)/* && ss != EMPTY*/;
		}
		inline bool leftVacant() const
		{
			return (visited == NO_VISITED || visited == RIGHT_VISITED)/* && ss != EMPTY*/;
		}
		inline void setLeftVisited()
		{
			if (visited == RIGHT_VISITED)
				visited = FULLY_VISITED;
			else
				visited = LEFT_VISITED;
		}
		inline void setRightVisited()
		{
			if (visited == LEFT_VISITED)
				visited = FULLY_VISITED;
			else
				visited = RIGHT_VISITED;
		}
		// split from the end[0] to end[1], according to the distance to corresponding spheres
		SPLIT_STATE split(const vector<Bisector>& bisecList, list<Vertex_handle> *left, list<Vertex_handle>* right);
		void leftSplitAndVisit(list<Vertex_handle>* l);
		void rightSplitAndVisit(list<Vertex_handle>* l);
		void split()
		{
			if (segIdx.size() == 0)
				ss = EMPTY;
			else
				ss = split(*bisecList, &llist, &rlist);
		}
		/* data field */
		pair<int, int> base;
		pair<int, int> end;
		list<int> segIdx; // indices in polygons[base0]
		list<Vertex_handle> llist;
		list<Vertex_handle> rlist;
		VISIT_STATE visited;
		bool reversed;
		vector<Bisector>* bisecList;
		SPLIT_STATE ss;
	};	


}

#endif
