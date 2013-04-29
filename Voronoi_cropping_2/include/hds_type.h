// ------------------------------------------------------------------------------
// Copyright (c) 2013 GeometryFactory (FRANCE)
// Distributed under the Boost Software License, Version 1.0. (See accompany-
// ing file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
// ------------------------------------------------------------------------------

#include <CGAL/HalfedgeDS_default.h>

// A face type with a color member variable.
template <class Refs>
class Face_with_int : public CGAL::HalfedgeDS_face_base<Refs> {
  typedef typename Refs::Halfedge_handle Halfedge_handle;
  int m_index;
public:
  Face_with_int(int c=-1):m_index(c) {}
  int& color() {return m_index; }
  std::vector<Halfedge_handle> holes;
};

// A face type with a ptr as a color member variable.
template <class Refs>
class Face_with_int_ptr : public CGAL::HalfedgeDS_face_base<Refs> {
  typedef typename Refs::Halfedge_handle Halfedge_handle;
  std::pair<int,int>* m_index_and_color_ptr;
public:
  Face_with_int_ptr(std::pair<int,int>* ptr=NULL):m_index_and_color_ptr(ptr) {}
  int& color() {return m_index_and_color_ptr->second; }
  int  index() {return m_index_and_color_ptr->first; }
  std::vector<Halfedge_handle> holes;
};

// A vertex type with a boolean
template <class Refs, class Point>
class Vertex_with_bool : public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point> {
  bool m_bool;
  typedef CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point> Vertex_base;
public:
  Vertex_with_bool():m_bool(false) {}
  Vertex_with_bool(const Point& p):Vertex_base(p), m_bool(false) {}
  bool& is_on_border() {return m_bool; }
};


struct HDS_Item_extra : public CGAL::HalfedgeDS_items_2 {
    template <class Refs, class Traits>
    struct Face_wrapper {
        typedef Face_with_int<Refs> Face;
    };

    template <class Refs, class Traits>
    struct Vertex_wrapper {
        typedef Vertex_with_bool<Refs, typename Traits::Point_2 > Vertex;
    };
};

struct HDS_Item_extra_ptr : public CGAL::HalfedgeDS_items_2 {
    template <class Refs, class Traits>
    struct Face_wrapper {
        typedef Face_with_int_ptr<Refs> Face;
    };

    template <class Refs, class Traits>
    struct Vertex_wrapper {
        typedef Vertex_with_bool<Refs, typename Traits::Point_2 > Vertex;
    };
};
