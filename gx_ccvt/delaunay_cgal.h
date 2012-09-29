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

#ifndef __DELAUNAY_CGAL__
#define __DELAUNAY_CGAL__

#include <Geex/CVT/geometry.h>

//______________________________________________________________________________________
// CGAL stuff

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
//#include <CGAL/Triangulation_hierarchy_3.h>
#include <CGAL/Timer.h>
#include <CGAL/Tetrahedron_3.h>

namespace Geex {


// ------------- My stuff to extend CGAL --------------------------

	// three types of vertex. 
	// V_CORNER: vertex's VD contains a feature point, this vertex is fixed during the optimization
	// V_CREASE: vertex's VD contains a connected feature line
	// V_SMOOTH: vertex's VD contains a smooth surface patch
	enum VertexType {V_CORNER, V_CREASE, V_SMOOTH} ;

    class VertexDecoration {
	public:
        VertexDecoration() {
            id = -1 ;
            weight_ = 1. ;
            vtype = V_SMOOTH ;
        }

    public:
        int    id ;
        real   weight_ ;
        VertexType vtype ; 
    public:
        int&   idx()       { return id ; }
        VertexType& type() { return vtype ; }
		real& weight()     { return weight_ ; }
    } ;

    class CellDecoration {
    public:

        CellDecoration() {
            clear_marks() ;
        }

        void clear_marks() {
            for(int i=0; i<6; i++) {
                edge_is_marked[i] = false ;
            }
            infinite = false ;
            aux      = false ;
			weight_  = 1. ;
			cell_index = -1 ;
        }

        int& idx()     { return cell_index ; }
		real& weight() { return weight_ ; }

        bool infinite ;
        vec3 dual ;
        bool dual_outside ;
        bool edge_is_marked[6] ;
        int  cell_index ; // index of cell
		bool aux;    // 
		real weight_ ;
    } ;

// ------------- CGAL stuff ---------------------------------

    // ----------------------- A CGAL::Vertex with decoration ------------------
    template < class Gt, class Vb = CGAL::Triangulation_vertex_base_3<Gt> >
    class Vertex : public  Vb, public VertexDecoration  {
        typedef Vb superclass;
    public:
        typedef typename Vb::Vertex_handle      Vertex_handle;
        typedef typename Vb::Cell_handle        Cell_handle;
        typedef typename Vb::Point              Point;
        
        template < typename TDS2 >
        struct Rebind_TDS {
            typedef typename Vb::template Rebind_TDS<TDS2>::Other Vb2;
            typedef Vertex<Gt,Vb2> Other;
        } ;
        
    public:
        Vertex() : superclass() {}
        Vertex(const Point & p) : superclass(p) {}
        Vertex(const Point & p, Cell_handle f) : superclass(f,p) {}
        Vertex(Cell_handle f) : superclass(f) {}
    } ;


    // ----------------------- A CGAL::Cell with decoration ------------------
    template < class Gt, class Cb = CGAL::Triangulation_cell_base_3<Gt> >
    class Cell : public Cb, public CellDecoration {
        typedef Cb superclass;
    public:
        typedef typename Cb::Vertex_handle      Vertex_handle;
        typedef typename Cb::Cell_handle        Cell_handle;
        template < typename TDS2 >
        struct Rebind_TDS {
            typedef typename Cb::template Rebind_TDS<TDS2>::Other Cb2;
            typedef Cell<Gt,Cb2> Other;
        };


        Cell() : superclass() {  }
        Cell(
            Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3
        ) : superclass(v0,v1,v2,v3) { }
            
        Cell(
            Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3,
            Cell_handle n0, Cell_handle n1, Cell_handle n2, Cell_handle n3
        ) : superclass(v0,v1,v2,v3,n0,n1,n2,n3) { }

    } ;
    
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;
    typedef Vertex<K> Vb ;
    //typedef CGAL::Triangulation_hierarchy_vertex_base_3<Vb> Vbh ;
    typedef Vb Vbh ;
    typedef Cell<K>   Cb ;
    typedef CGAL::Triangulation_data_structure_3<Vbh,Cb> TDS;
    typedef CGAL::Delaunay_triangulation_3<K, TDS> TRI ;
	//typedef CGAL::Triangulation_hierarchy_3<TRI>   Dh ;

	typedef TDS::Vertex_handle Vertex_handle ;
    typedef K::Point_3    Point ;
    typedef K::Vector_3   CGAL_Vector;
    typedef K::Segment_3  Segment ;
    typedef K::Ray_3      Ray ;
	typedef K::Tetrahedron_3 Tetrahedron ;

    inline vec3 to_geex(const Point& p) { return vec3(p.x(), p.y(), p.z()) ; }
    inline Point to_cgal(const vec3& p) { return Point(p.x, p.y, p.z) ; }


    //class DelaunayBase : public CGAL::Triangulation_hierarchy_3<TRI> {
	  class DelaunayBase : public TRI {
    public:
		virtual ~DelaunayBase(){}
    } ;

}


#endif
