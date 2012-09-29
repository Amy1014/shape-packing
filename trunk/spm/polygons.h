#ifndef __POLYGONS__
#define __POLYGONS__


#include <Geex/mathematics/glsl_linear.h>
#include "geometry.h"
#include <set>
#include "spm_cgal.h"
#include <algorithm>

namespace Geex {

    class Convex : public std::vector< Line<real> > 
	{
    public:
        inline bool contains(const vec2& p) 
		{
            for(unsigned int i=0; i<size(); i++) 
			{
                if((*this)[i].side(p) < 0) 
                    return false ;
            }
            return true ;
        }
    } ;


    template <class T> inline void sets_intersect(
        const std::set<T>& S1, const std::set<T>& S2, std::set<T>& I
    ) 
	{
        I.clear() ;
        std::set_intersection(S1.begin(), S1.end(), S2.begin(), S2.end(), std::inserter(I, I.begin())) ;
    }

    class PolygonVertex : public vec2 {
    public:
        PolygonVertex() : on_boundary(false), id(-1) { }
        PolygonVertex(real x, real y) : vec2(x,y), on_boundary(false),id(-1) { }
        PolygonVertex(const vec2& V) : vec2(V), on_boundary(false),id(-1) { }
        
        std::set<int> bisectors ;
        std::set<int> boundary_edges ;
        bool on_boundary ;
		int id;
    } ;


    class Polygon : public std::vector<PolygonVertex> {
    public:
        void load(const std::string& filename) ;
        void normalize() ;
        bool contains(const vec2& p) ;
        Line<real> edge_line(unsigned int i) const ;
        void convex_clip(Polygon& to, const Line<real>& L, int E) ;
    } ;

    class PolygonEdge 
	{
    public:
        PolygonEdge() {vh[0]=NULL; vh[1]=NULL;  }
        PolygonEdge(const PolygonVertex& v1, const PolygonVertex& v2) 
		{
            vertex[0] = v1 ;
            vertex[1] = v2 ;
			vh[0]=NULL; vh[1]=NULL;
        }
        Line<real> line() const {
            vec2 V = vertex[1] - vertex[0] ;
            return Line<real>(vertex[0], vec2(-V.y, V.x)) ;
        }
        PolygonVertex vertex[2] ;
		RegularBase::Vertex_handle vh[2]; // dual vertices

    } ;

    class Polygon2 : public std::vector<PolygonEdge> 
	{
    public:
        void load(const std::string& filename) ;
        void normalize() ;
        bool contains(const vec2& p) ;
        void convex_clip(Polygon2& to, const Line<real>& L, int E, bool close = true) ;
		void convex_clip(Polygon2& to, const Line<real>& L, RegularBase::Edge E, bool close);
		size_t  nb_vertices(){return vertices.size();}
		vec2 mid_point();
	public:
		std::vector<vec2> vertices;
    } ;
}

#endif
