// ------------------------------------------------------------------------------
// Copyright (c) 2013 GeometryFactory (FRANCE)
// Distributed under the Boost Software License, Version 1.0. (See accompany-
// ing file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
// ------------------------------------------------------------------------------


#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/HalfedgeDS_decorator.h>
#include <iterator>
#include <map>
#include <cmath>

#ifdef DEBUG_JOIN_FACES
#include <sstream>
#include <fstream>
#endif

/*! returns k if `pt=bbox[k]` or `(k + k+1)/2` if `pt` lies inbetween bbox[k] and bbox[k+1] */
template <class Point_2, class Iso_rectangle_2>
double index_on_bbox(const Point_2& pt, const Iso_rectangle_2& bbox)
{
  CGAL_precondition( bbox[0].x()==bbox.xmin() && bbox[0].y()==bbox.ymin() );

  if ( pt.y() == bbox.ymin() ){
    if ( pt.x() == bbox.xmin() ) return 0;
    if ( pt.x() == bbox.xmax() ) return 1;
    return 0.5;
  }

  if ( pt.y() == bbox.ymax() ){
    if ( pt.x() == bbox.xmin() ) return 3;
    if ( pt.x() == bbox.xmax() ) return 2;
    return 2.5;
  }

  if ( pt.x() == bbox.xmin() ) return 3.5;
  CGAL_assertion ( pt.x() == bbox.xmax() );
  return 1.5;
}

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <bitset>

/// a class to store the dual of a Delaunay face and the orientation
///  of the dual wrt to each boundary segment of the iso-rectangle
template <class Input_kernel,class Exact_kernel>
class Face_info_for_DT2
{
  std::bitset<4> m_orientation_vs_bbox;
  std::bitset<4> m_on_boundary;
  std::size_t m_dual_index;

  typedef std::vector<typename Exact_kernel::Point_2> Voronoi_vertices;

public:
  void set_orientation(int i, CGAL::Orientation o)
  {
    switch (o)
    {
      case CGAL::CLOCKWISE:
        m_orientation_vs_bbox.set(i);
        break;
      case CGAL::COLLINEAR:
        m_on_boundary.set(i);
      default:
        //do nothing
        break;
    }
  }

  std::size_t dual_index() const
  {
    return m_dual_index;
  }

  void set_dual_index(std::size_t i)
  {
    m_dual_index=i;
  }

  const typename Exact_kernel::Point_2&
  dual(const Voronoi_vertices& voronoi_vertices) const
  {
    return voronoi_vertices[m_dual_index];
  }

  const typename Input_kernel::Point_2
  inexact_dual(const Voronoi_vertices& voronoi_vertices) const
  {
    CGAL::Cartesian_converter<Exact_kernel, Input_kernel> to_input;
    return to_input( voronoi_vertices[m_dual_index] );
  }

  bool is_dual_on_bbox_boundary() const {
    return m_on_boundary.any();
  }

  bool is_dual_in_bbox() const{
    return m_orientation_vs_bbox.none();
  }
};

template <class T>
std::pair<T,T> make_sorted_pair(T t1, T t2)
{
  return t1 < t2 ? std::make_pair(t1, t2):std::make_pair(t2, t1);
}




template < class Input_kernel,
           class Exact_kernel,
           class HDS,
           class DT2,
           class Create_hds_face >
class Cropped_Voronoi_hds_creator{
/// typedefs
/// ========
  typedef typename HDS::Halfedge_handle Halfedge_handle;
  typedef typename HDS::Face_handle HDS_Face_handle;
  typedef typename HDS::Vertex_handle HDS_Vertex_handle;
  typedef typename DT2::Face_handle DT2_Face_handle;
  typedef typename DT2::Vertex_handle DT2_Vertex_handle;

  typedef std::map< std::size_t, HDS_Vertex_handle > Vertex_map;
  typedef std::map< std::pair< DT2_Face_handle, DT2_Face_handle >, Halfedge_handle > Dual_hedge_map;
/// data members
/// ============
  /// this map is used to keep track of the halfedges incident to
  /// a vertex on the iso-rectangle boundary, in order to link the
  /// border halfedges together. Filling the map is done using
  /// set_border_hedge_source and set_border_hedge_target
  std::map<HDS_Vertex_handle, std::pair<Halfedge_handle, Halfedge_handle> > border_hedges;
  /// map used to match a Delaunay edge to a dual cropped edge
  Dual_hedge_map dual_hedge_map;
  /// map used to match a Delaunay face to a dual vertex
  Vertex_map vertex_map;
  /// a vector containing all dual vertices
  std::vector< typename Exact_kernel::Point_2 > exact_voronoi_vertices;
  /// decorator
  CGAL::HalfedgeDS_items_decorator<HDS> item_decorator;

  /// convert from input to exact
  CGAL::Cartesian_converter<Input_kernel, Exact_kernel> to_exact;
  /// reference to the Delaunay triangulation with the input points
  const DT2& dt2;
  /// iso rectangle to which we crop the Voronoi diagram
  typename Input_kernel::Iso_rectangle_2 input_box;
  /// the exact version
  typename Exact_kernel::Iso_rectangle_2 exact_box;
  /// reference to the output hds
  HDS& hds;
  /// reference to the functor for creating faces
  const Create_hds_face& create_face;

/// private functions
/// =================
  /// see border_hedges
  void set_border_hedge_source( Halfedge_handle h, HDS_Vertex_handle v )
  {
    border_hedges[v].first=h;
  }

  /// same border_hedges
  void set_border_hedge_target( Halfedge_handle h, HDS_Vertex_handle v )
  {
    border_hedges[v].second=h;
  }

  /// This function adds the missing halfedges and vertices to close a Voronoi
  /// face intersected by the iso-rectangle boundary
  void
  add_halfedge_on_iso_rectangle(
    HDS_Face_handle new_face,
    Halfedge_handle previous,
    Halfedge_handle next )
  {
    //ensure target of previous is on the rectangle boundary
    CGAL_assertion(
      input_box.xmin()==previous->vertex()->point().x() ||  input_box.xmax()==previous->vertex()->point().x() ||
      input_box.ymin()==previous->vertex()->point().y() ||  input_box.ymax()==previous->vertex()->point().y()
    );

    //ensure source of next is on the rectangle boundary
    CGAL_assertion(
      input_box.xmin()==next->opposite()->vertex()->point().x() ||  input_box.xmax()==next->opposite()->vertex()->point().x() ||
      input_box.ymin()==next->opposite()->vertex()->point().y() ||  input_box.ymax()==next->opposite()->vertex()->point().y()
    );

    //ensure the points are different
    CGAL_assertion( next->opposite()->vertex()->point()!=previous->vertex()->point() );

    //input_box[0], input_box[1], input_box[2], input_box[3] are counterclockwise oriented
    double index_target=index_on_bbox(previous->vertex()->point(),input_box);
    double index_next=index_on_bbox(next->opposite()->vertex()->point(),input_box);

    if (index_next < index_target) index_next+=4;
    int k = std::floor( index_target);

  //add the missing vertices
    while (++k < index_next)
    {
      HDS_Vertex_handle new_vertex =
        hds.vertices_push_back( typename HDS::Vertex( input_box[k] ) );

      Halfedge_handle new_hh=hds.edges_push_back(typename HDS::Halfedge(), typename HDS::Halfedge());
      item_decorator.set_face(new_hh,new_face);
      item_decorator.set_vertex(new_hh,new_vertex);
      set_border_hedge_target(new_hh->opposite(), previous->vertex());
      set_border_hedge_source(new_hh->opposite(), new_vertex);
      item_decorator.set_vertex(new_hh->opposite(),previous->vertex());

      previous->set_next(new_hh);
      item_decorator.set_prev(new_hh,previous);

      previous=new_hh;
    }

    //add the final halfedge closing the face
    Halfedge_handle new_hh=hds.edges_push_back(typename HDS::Halfedge(), typename HDS::Halfedge());
    item_decorator.set_face(new_hh,new_face);
    item_decorator.set_vertex(new_hh,next->opposite()->vertex());
    item_decorator.set_vertex(new_hh->opposite(),previous->vertex());

    set_border_hedge_target(new_hh->opposite(), previous->vertex());
    set_border_hedge_source(new_hh->opposite(), new_hh->vertex());

    previous->set_next(new_hh);
    item_decorator.set_prev(new_hh,previous);

    new_hh->set_next(next);
    item_decorator.set_prev(next,new_hh);
  }

  /// utility function to ensure that each dual of Delaunay vertex is unique
  HDS_Vertex_handle
  get_dual_vertex( DT2_Face_handle face )
  {
    CGAL_precondition( face->info().is_dual_in_bbox() );
    std::pair<typename Vertex_map::iterator, bool> insert_res =
      vertex_map.insert( std::make_pair(face->info().dual_index(), HDS_Vertex_handle()) );
    if (insert_res.second)
    {
      insert_res.first->second =
        hds.vertices_push_back(
          typename HDS::Vertex( face->info().inexact_dual(exact_voronoi_vertices) ) );
      insert_res.first->second->is_on_border()=face->info().is_dual_on_bbox_boundary();
    }
    return insert_res.first->second;
  }

  /// return true is `primitive` intersects exact_box as a segment. The result is put
  /// into `res`. The orientation of the segment is the same as that of `bisector`
  template <class Primitive, class Bisector>
  bool crop_to_iso_rectangle( const Primitive& primitive,
                              const Bisector& bisector,
                              typename Exact_kernel::Segment_2& res)
  {
    CGAL::Object obj =
      typename Exact_kernel::Intersect_3()( primitive, exact_box );
    const typename Exact_kernel::Segment_2* s=
      CGAL::object_cast<typename Exact_kernel::Segment_2>(&obj);

    if(!s) return false;

    //the bisector is correctly oriented
    if ( s->to_vector() * bisector.to_vector() > 0)
      res=*s;
    else
      res=s->opposite();
    return true;
  }

  /// function to convert exact vertices into the input kernel while ensuring that points on the iso_rectangle
  /// boundary are correctly rounded into the input type.
  typename Input_kernel::Point_2
  convert_to_input(const typename Exact_kernel::Point_2& p)
  {
    CGAL::NT_converter<typename Exact_kernel::FT, typename Input_kernel::FT> to_input;
    typename Input_kernel::FT x,y;

    //handle x coordinate
    if ( p.x() == exact_box.xmin() )
      x=to_input(exact_box.xmin());
    else{
      if ( p.x() == exact_box.xmax() ) x=to_input(exact_box.xmax());
      else x=to_input(p.x());
    }

    //handle y coordinate
    if ( p.y() == exact_box.ymin() )
      y=to_input(exact_box.ymin());
    else{
      if ( p.y() == exact_box.ymax() ) y=to_input(exact_box.ymax());
      else y=to_input(p.y());
    }

    return typename Input_kernel::Point_2(x,y);
  }

  /// given an oriented pair of incident Delaunay faces, returns the dual halfedge cropped to
  /// `exact_box` if the intersection is not restricted to a point and Halfedge_handle() otherwise.
  /// Edge are computed thanks to a map
  Halfedge_handle
  get_dual_halfedge_cropped(  DT2_Face_handle source_face,
                              DT2_Face_handle target_face,
                              int edge_index_in_source_face)
  {
    typedef Halfedge_handle Halfedge_handle;
    typedef typename HDS::Vertex_handle Vertex_handle;

    std::pair<DT2_Face_handle, DT2_Face_handle> sorted_pair=
      make_sorted_pair(source_face, target_face);

    std::pair< typename Dual_hedge_map::iterator, bool> insert_res =
      dual_hedge_map.insert( std::make_pair( sorted_pair, Halfedge_handle() ) );

    HDS_Vertex_handle source_vertex, target_vertex;

    if ( insert_res.second )
    {
      if ( dt2.is_infinite(source_face) )
      {
        if (   dt2.is_infinite( source_face->vertex(DT2::cw(edge_index_in_source_face)) )
            || dt2.is_infinite( source_face->vertex(DT2::ccw(edge_index_in_source_face)) ) )
        {
          return Halfedge_handle();
        }

        typename Exact_kernel::Line_2 bisector =
          typename Exact_kernel::Construct_bisector_2()(
              to_exact( source_face->vertex( DT2::cw(edge_index_in_source_face) )->point() ),
              to_exact( source_face->vertex( DT2::ccw(edge_index_in_source_face) )->point() )
        );
        CGAL_assertion( bisector.has_on_positive_side(
          to_exact( source_face->vertex( DT2::cw(edge_index_in_source_face) )->point() )  ) );

        if ( dt2.is_infinite(target_face) )
        {
          //the dual is a line
          typename Exact_kernel::Segment_2 inter;
          if ( !crop_to_iso_rectangle(bisector, bisector, inter) )
            return Halfedge_handle();
          source_vertex=hds.vertices_push_back( typename HDS::Vertex( convert_to_input( inter.source() )) );
          target_vertex=hds.vertices_push_back( typename HDS::Vertex( convert_to_input( inter.target() )) );
          source_vertex->is_on_border()=true;
          target_vertex->is_on_border()=true;
        }
        else
        {
          //the dual is a ray
          typename Exact_kernel::Ray_2 ray( target_face->info().dual(exact_voronoi_vertices), bisector.opposite() );
          typename Exact_kernel::Segment_2 inter;
          if ( !crop_to_iso_rectangle(ray, bisector, inter) )
            return Halfedge_handle();
          if ( target_face->info().is_dual_in_bbox() )
          {
            target_vertex=get_dual_vertex(target_face);
            source_vertex=hds.vertices_push_back( typename HDS::Vertex( convert_to_input( inter.source() )) );
            source_vertex->is_on_border()=true;
          }
          else{
            source_vertex=hds.vertices_push_back( typename HDS::Vertex( convert_to_input( inter.source() )) );
            target_vertex=hds.vertices_push_back( typename HDS::Vertex( convert_to_input( inter.target() )) );
            source_vertex->is_on_border()=true;
            target_vertex->is_on_border()=true;
          }
        }
      }
      else{
        if ( dt2.is_infinite(target_face) )
        {
          if (   dt2.is_infinite( source_face->vertex(DT2::cw(edge_index_in_source_face)) )
              || dt2.is_infinite( source_face->vertex(DT2::ccw(edge_index_in_source_face)) ) )
          {
            return Halfedge_handle();
          }

          typename Exact_kernel::Line_2 bisector =
            typename Exact_kernel::Construct_bisector_2()(
                to_exact( source_face->vertex( DT2::cw(edge_index_in_source_face) )->point() ),
                to_exact( source_face->vertex( DT2::ccw(edge_index_in_source_face) )->point() ) );
          CGAL_assertion( bisector.has_on_positive_side(
            to_exact( source_face->vertex( DT2::cw(edge_index_in_source_face) )->point() )  ) );

          //the dual is a ray
          typename Exact_kernel::Ray_2 ray( source_face->info().dual(exact_voronoi_vertices), bisector );
          typename Exact_kernel::Segment_2 inter;
          if ( !crop_to_iso_rectangle(ray, bisector, inter) )
            return Halfedge_handle();
          if ( source_face->info().is_dual_in_bbox() )
          {
            source_vertex=get_dual_vertex(source_face);
            target_vertex=hds.vertices_push_back( typename HDS::Vertex( convert_to_input( inter.target() )) );
            target_vertex->is_on_border()=true;
          }
          else{
            source_vertex=hds.vertices_push_back( typename HDS::Vertex( convert_to_input( inter.source() )) );
            target_vertex=hds.vertices_push_back( typename HDS::Vertex( convert_to_input( inter.target() )) );
            source_vertex->is_on_border()=true;
            target_vertex->is_on_border()=true;
          }
        }
        else
        {
          //in case of co-circular input, a halfedge might be a point
          if ( source_face->info().dual_index() == target_face->info().dual_index() ) return Halfedge_handle();

          //due to the approximation into double, the embedded into double edge can be null or swapped
          //thus we check that the interval box around the points do not overlap

          //the dual is a segment
          if ( source_face->info().is_dual_in_bbox() &&
               target_face->info().is_dual_in_bbox() )
          {
            source_vertex = get_dual_vertex(source_face);
            target_vertex = get_dual_vertex(target_face);
          }
          else{
            typename Exact_kernel::Segment_2 segment( source_face->info().dual(exact_voronoi_vertices),
                                                      target_face->info().dual(exact_voronoi_vertices) );
            typename Exact_kernel::Segment_2 inter;
            if ( !crop_to_iso_rectangle(segment, segment, inter) )
              return Halfedge_handle();


            if ( source_face->info().is_dual_in_bbox() )
              source_vertex = get_dual_vertex(source_face);
            else
            {
              source_vertex = hds.vertices_push_back( typename HDS::Vertex( convert_to_input( inter.source() )) );
              source_vertex->is_on_border()=true;
            }
            if ( target_face->info().is_dual_in_bbox() )
              target_vertex = get_dual_vertex(target_face);
            else
            {
              target_vertex = hds.vertices_push_back( typename HDS::Vertex( convert_to_input( inter.target() )) );
              target_vertex->is_on_border()=true;
            }
          }
        }
      }

      CGAL_assertion( source_vertex != HDS_Vertex_handle() );
      CGAL_assertion( target_vertex != HDS_Vertex_handle() );

      Halfedge_handle hh=hds.edges_push_back(typename HDS::Halfedge(), typename HDS::Halfedge());
      if (source_face == sorted_pair.first)
      {
        item_decorator.set_vertex(hh,target_vertex);
        item_decorator.set_vertex(hh->opposite(),source_vertex);
        item_decorator.set_vertex_halfedge(target_vertex, hh);
        item_decorator.set_vertex_halfedge(source_vertex, hh->opposite());
      }
      else
      {
        item_decorator.set_vertex(hh,source_vertex);
        item_decorator.set_vertex(hh->opposite(),target_vertex);
        item_decorator.set_vertex_halfedge(source_vertex, hh);
        item_decorator.set_vertex_halfedge(target_vertex, hh->opposite());
      }
      insert_res.first->second=hh;
    }

    return ( insert_res.first->second == Halfedge_handle() ||
           source_face == sorted_pair.first ) ?
      insert_res.first->second : insert_res.first->second->opposite();
  }

  void init_dt_face_info()
  {
    typename Exact_kernel::Construct_circumcenter_2 circumcenter;

    exact_voronoi_vertices.reserve( dt2.number_of_faces() );

    typedef std::map<typename Exact_kernel::Point_2, std::size_t> Point_to_index_map;
    Point_to_index_map point_to_index;

    std::size_t index=0;
    for (typename DT2::Finite_faces_iterator  fit=dt2.finite_faces_begin(),
                                              fit_end=dt2.finite_faces_end();
                                              fit!=fit_end; ++fit )
    {
      std::pair<typename Point_to_index_map::iterator, bool> insert_res=
        point_to_index.insert(
          std::make_pair(
                  circumcenter( to_exact( fit->vertex(0)->point() ),
                                to_exact( fit->vertex(1)->point() ),
                                to_exact( fit->vertex(2)->point() ) ),
                  index )
        );

      if ( insert_res.second ){
        exact_voronoi_vertices.push_back( insert_res.first->first );
        ++index;
      }

      fit->info().set_dual_index(insert_res.first->second);
    }
  }

  /// init input_box and exact_box so that all Voronoi vertices are included
  void init_boxes()
  {
    std::size_t nb_element = exact_voronoi_vertices.size();
    CGAL::Bbox_2 bbox;
    if ( nb_element <4 )
    {
      typename DT2::Finite_vertices_iterator vit=dt2.finite_vertices_begin();
      typename DT2::Finite_vertices_iterator vit_end=dt2.finite_vertices_end();
      bbox=vit->point().bbox();
      for(;++vit!=vit_end;)
        bbox=bbox+vit->point().bbox();
    }
    else
      bbox=exact_voronoi_vertices[0].bbox();
    for(std::size_t i=0; i< nb_element; ++i)
      bbox=bbox+exact_voronoi_vertices[i].bbox();

    double dx = 0.0025 * ( bbox.xmax() - bbox.xmin() );
    double dy = 0.0025 * ( bbox.ymax() - bbox.ymin() );

    double xmin=bbox.xmin()-dx;
    double xmax=bbox.xmax()+dx;
    double ymin=bbox.ymin()-dy;
    double ymax=bbox.ymax()+dy;

    if (dx==0) { xmin=ymin; xmax=ymax; }
    if (dy==0) { ymin=xmin; ymax=xmax; }

    input_box=typename Input_kernel::Iso_rectangle_2(xmin, ymin, xmax, ymax);
    exact_box=to_exact(input_box);
  }

  void copy_the_box_in_the_hds(HDS_Face_handle new_face)
  {
    HDS_Vertex_handle vertices[4] = {
      hds.vertices_push_back( typename HDS::Vertex( input_box[0] ) ),
      hds.vertices_push_back( typename HDS::Vertex( input_box[1] ) ),
      hds.vertices_push_back( typename HDS::Vertex( input_box[2] ) ),
      hds.vertices_push_back( typename HDS::Vertex( input_box[3] ) )
    };

    Halfedge_handle new_hedges[4] = {
      hds.edges_push_back(typename HDS::Halfedge(), typename HDS::Halfedge()),
      hds.edges_push_back(typename HDS::Halfedge(), typename HDS::Halfedge()),
      hds.edges_push_back(typename HDS::Halfedge(), typename HDS::Halfedge()),
      hds.edges_push_back(typename HDS::Halfedge(), typename HDS::Halfedge())
    };

    for (int k=0;k<4; ++k)
    {
      item_decorator.set_face(new_hedges[k], new_face);
      //set hedge target vertex
      item_decorator.set_vertex(new_hedges[k], vertices[k]);
      item_decorator.set_vertex_halfedge(vertices[k], new_hedges[k]);
      item_decorator.set_vertex(new_hedges[(k+1)%4]->opposite(), vertices[k]);
      //set next-prev for inner hedges
      new_hedges[k]->set_next( new_hedges[ (k+1)%4 ] );
      item_decorator.set_prev( new_hedges[ (k+1)%4 ], new_hedges[k] );
      //set next-prev for border hedges
      new_hedges[k]->opposite()->set_next( new_hedges[ (k+3)%4 ]->opposite() );
      item_decorator.set_prev( new_hedges[ (k+3)%4 ]->opposite(), new_hedges[k]->opposite() );
    }

    item_decorator.set_face_halfedge(new_face, new_hedges[0]);
  }

  void handle_new_face(std::vector<Halfedge_handle>& new_hedges, DT2_Vertex_handle vh )
  {
    if ( !new_hedges.empty() ) // some edges have been created
    {
      std::size_t nb_hedges=new_hedges.size();
      if (nb_hedges == 1) //if there is one edge
      {
        HDS_Vertex_handle start_vertex = new_hedges[0]->vertex();
        HDS_Vertex_handle goal_vertex = new_hedges[0]->opposite()->vertex();

        CGAL_assertion( start_vertex->is_on_border() && goal_vertex->is_on_border() );

        double start_index=index_on_bbox( start_vertex->point(), input_box );
        double goal_index =index_on_bbox( goal_vertex->point(), input_box );

        if (goal_index<start_index) goal_index+=4;

        if ( std::ceil(goal_index)-std::floor(start_index) == 1 )
        {
          set_border_hedge_source( new_hedges[0], new_hedges[0]->opposite()->vertex() );
          set_border_hedge_target( new_hedges[0], new_hedges[0]->vertex() );
          return; // they are on the same input_box edge, oriented outside of the box
        }
      }

      new_hedges.push_back(new_hedges[0]);
      HDS_Face_handle new_face=hds.faces_push_back( create_face(vh) );
      CGAL_assertion( new_hedges[0]!=Halfedge_handle() );
      item_decorator.set_face_halfedge(new_face,new_hedges[0]); //set one halfedge for the face

      for (std::size_t i=0;i<nb_hedges;++i)
      {
        item_decorator.set_face(new_hedges[i],new_face);

        if ( new_hedges[i]->vertex()->is_on_border() &&
             new_hedges[i]->vertex()!=new_hedges[i+1]->opposite()->vertex() )
        {
          CGAL_assertion( new_hedges[i+1]->opposite()->vertex()->is_on_border() );
          //we need one to several hedges on the iso-rectangle boundary to link the two vertices
          add_halfedge_on_iso_rectangle(new_face, new_hedges[i], new_hedges[i+1]);
        }
        else
        {
          new_hedges[i]->set_next(new_hedges[i+1]);
          item_decorator.set_prev(new_hedges[i+1], new_hedges[i]);
        }
      }
    }
  }

  struct Sort_along_line
  {
    typename Exact_kernel::Vector_2 base;
    CGAL::Cartesian_converter<Input_kernel, Exact_kernel> to_exact;
    Sort_along_line(DT2_Vertex_handle v1, DT2_Vertex_handle v2)
      :base( to_exact( v1->point() ),
             to_exact( v2->point() ) ) {}

    bool operator() (DT2_Vertex_handle v1, DT2_Vertex_handle v2) const
    {
      return ( base * (to_exact( v2->point() ) - to_exact( v1->point() )) ) > 0;
    }
  };

  void build_hds_dimension_1()
  {

    std::vector<DT2_Vertex_handle> vertices;
    for (typename DT2::Finite_vertices_iterator vit=dt2.finite_vertices_begin();
                                                vit!=dt2.finite_vertices_end(); ++vit)
    {
      vertices.push_back(vit);
    }

    /// sort the points along the line
    std::sort(vertices.begin(), vertices.end(), Sort_along_line(vertices[0], vertices[1]) );

    std::size_t nb_pts=vertices.size();

    Halfedge_handle previous_halfedge=Halfedge_handle();
    std::size_t k=0;
    //find the first bisector intersecting the rectangle
    do{
      typename Exact_kernel::Line_2 bisector =
        typename Exact_kernel::Construct_bisector_2()(
            to_exact( vertices[k]->point() ),
            to_exact( vertices[k+1]->point() )
      );
      CGAL_assertion( bisector.has_on_positive_side( to_exact( vertices[k]->point() )  ) );

      typename Exact_kernel::Segment_2 inter;
      if ( crop_to_iso_rectangle(bisector, bisector, inter) )
      {
        //compute the halfedge cropped and set it
        Halfedge_handle new_halfedge = hds.edges_push_back(typename HDS::Halfedge(), typename HDS::Halfedge());
        HDS_Vertex_handle source_vertex=hds.vertices_push_back( typename HDS::Vertex( convert_to_input( inter.source() )) );
        HDS_Vertex_handle target_vertex=hds.vertices_push_back( typename HDS::Vertex( convert_to_input( inter.target() )) );
        source_vertex->is_on_border()=true;
        target_vertex->is_on_border()=true;
        item_decorator.set_vertex(new_halfedge,target_vertex);
        item_decorator.set_vertex(new_halfedge->opposite(),source_vertex);
        item_decorator.set_vertex_halfedge(target_vertex, new_halfedge);
        item_decorator.set_vertex_halfedge(source_vertex, new_halfedge->opposite());

        //set the new face
        std::vector<Halfedge_handle> new_hedges;
        if( previous_halfedge != Halfedge_handle() )
          new_hedges.push_back( previous_halfedge->opposite() );

        new_hedges.push_back( new_halfedge);
        handle_new_face(new_hedges, vertices[k] );

        previous_halfedge=new_halfedge;
      }
    }while (++k!=nb_pts-1);

    if( previous_halfedge != Halfedge_handle() )
    {
      //set the last face
      std::vector<Halfedge_handle> new_hedges;
      new_hedges.push_back( previous_halfedge->opposite() );
      handle_new_face(new_hedges, vertices[k] );
    }
  }

public:

/// constructors
/// ============
  Cropped_Voronoi_hds_creator(
    const DT2& dt2_,
    const typename Input_kernel::Iso_rectangle_2& iso_rect_,
    HDS& hds_,
    const Create_hds_face& create_face_)
  :dt2(dt2_), input_box(iso_rect_), exact_box(to_exact(input_box)), hds(hds_), create_face(create_face_)
  {
    init_dt_face_info();
  }

  Cropped_Voronoi_hds_creator(
    const DT2& dt2_,
    HDS& hds_,
    const Create_hds_face& create_face_)
  :dt2(dt2_), hds(hds_), create_face(create_face_)
  {
    init_dt_face_info();
    init_boxes();
  }


/// main function
/// =============
  void run()
  {
    //set orientation of dual vertex wrt rectangle
    typename Exact_kernel::Orientation_2 orientation;
    for (typename DT2::Finite_faces_iterator  fit=dt2.finite_faces_begin(),
                                              fit_end=dt2.finite_faces_end();
                                              fit!=fit_end; ++fit )
    {
      for (int i=0;i<4;++i)
        fit->info().set_orientation(
                i,
                orientation(
                  fit->info().dual(exact_voronoi_vertices),
                  exact_box[i],
                  exact_box[i+1] ) );
    }

    Halfedge_handle null_halfedge=Halfedge_handle();

    switch( dt2.dimension() )
    {
      case -1:
      {
        HDS_Face_handle new_face=hds.faces_push_back( create_face() );
        copy_the_box_in_the_hds(new_face);
        return;
      }

      case 0:
      {
        HDS_Face_handle new_face=hds.faces_push_back(
          create_face(dt2.finite_vertices_begin()) );
        copy_the_box_in_the_hds(new_face);
        return;
      }

      case 1:
        build_hds_dimension_1();
        break;

      default:
        for (typename DT2::Finite_vertices_iterator  vit=dt2.finite_vertices_begin(),
                                                     vit_end=dt2.finite_vertices_end();
                                                     vit!=vit_end; ++vit )
        {
          DT2_Face_handle current_face=vit->face(), first_face=current_face;
          std::vector<Halfedge_handle> new_hedges;

          do{
            int other_vertex_index=DT2::ccw( current_face->index(vit) );
            DT2_Face_handle neighbor_face=current_face->neighbor(other_vertex_index);

            Halfedge_handle hedge =
              get_dual_halfedge_cropped(current_face, neighbor_face,
                                        other_vertex_index );

            if (hedge!=null_halfedge) new_hedges.push_back(hedge);

            current_face=neighbor_face;
          }
          while(current_face!=first_face);

          handle_new_face( new_hedges, vit );
        }
    }

    //link border halfedges
    for (
      typename std::map<HDS_Vertex_handle, std::pair<Halfedge_handle, Halfedge_handle> >::iterator
      link_it=border_hedges.begin(), link_it_end=border_hedges.end();link_it!=link_it_end; ++link_it)
    {
      CGAL_assertion(
        link_it->second.first->opposite()->vertex()==link_it->second.second->vertex() );
      link_it->second.second->set_next(link_it->second.first);
      item_decorator.set_prev(link_it->second.first, link_it->second.second);
      item_decorator.set_vertex_halfedge(link_it->first, link_it->second.second);
    }

    if (border_hedges.empty())
    {
      HDS_Face_handle new_face=hds.faces_push_back(
        create_face(
          typename Input_kernel::Point_2(
            (input_box.xmax()+input_box.xmin())/2,
            (input_box.ymax()+input_box.ymin())/2 ),
          dt2) );
      copy_the_box_in_the_hds(new_face);
    }
  }

  const typename Input_kernel::Iso_rectangle_2& get_input_box() const
  { return input_box; }
};

/// function calling the functor
template <  class Input_kernel,
            class Exact_kernel,
            class HDS,
            class DT2,
            class Create_hds_face >
void create_hds_for_cropped_voronoi_diagram(
  const DT2& dt2,
  const typename Input_kernel::Iso_rectangle_2& iso_rect,
  HDS& hds,
  const Create_hds_face& create_face)
{
  Cropped_Voronoi_hds_creator<Input_kernel, Exact_kernel, HDS, DT2, Create_hds_face> cropping(dt2, iso_rect, hds, create_face);
  cropping.run();
}

/// function calling the functor
template <  class Input_kernel,
            class Exact_kernel,
            class HDS,
            class DT2,
            class Create_hds_face >
const typename Input_kernel::Iso_rectangle_2&
create_hds_for_cropped_voronoi_diagram(
  const DT2& dt2,
  HDS& hds,
  const Create_hds_face& create_face)
{
  Cropped_Voronoi_hds_creator<Input_kernel, Exact_kernel, HDS, DT2, Create_hds_face> cropping(dt2, hds, create_face);
  cropping.run();
  return cropping.get_input_box();
}

template< class DT2, class HDS>
struct Default_create_face{
  typename HDS::Face
  operator()( typename DT2::Vertex_handle ) const
  {
    return typename HDS::Face();
  }

  typename HDS::Face
  operator()( const typename DT2::Point& ,const DT2&) const
  {
    return typename HDS::Face();
  }

  typename HDS::Face
  operator()() const
  {
    return typename HDS::Face();
  }
};

template< class DT2, class HDS>
struct Create_face_from_info{
  typename HDS::Face
  operator()( typename DT2::Vertex_handle vh) const
  {
    return typename HDS::Face(vh->info());
  }

  typename HDS::Face
  operator()( const typename DT2::Point& p,const DT2& dt2) const
  {
    return typename HDS::Face(	dt2.nearest_vertex(p)->info() );
  }
  typename HDS::Face
  operator()() const
  {
    return typename HDS::Face();
  }
};

template < class DT2, class PointIterator, class InfoIterator>
void dt2_insert_with_info(
  DT2& dt2,
  PointIterator point_begin,
  PointIterator point_end,
  InfoIterator info_begin,
  InfoIterator /*info_end*/)
{
  std::vector<std::ptrdiff_t> indices;
  std::vector<typename DT2::Point> points;
  std::vector<typename DT2::Triangulation_data_structure::Vertex::Info> infos;
  std::ptrdiff_t index=0;
  InfoIterator iit=info_begin;
  for ( PointIterator pit=point_begin;
        pit!=point_end;++pit, ++iit)
  {
    points.push_back( *pit );
    infos.push_back ( *iit );
    indices.push_back(index++);
  }

  typedef CGAL::Spatial_sort_traits_adapter_2<typename DT2::Geom_traits,typename DT2::Point*> Search_traits;

  CGAL::spatial_sort(indices.begin(),indices.end(),Search_traits(&(points[0]),dt2.geom_traits()));

  typename DT2::Vertex_handle v_hint;
  typename DT2::Face_handle hint;
  for (typename std::vector<std::ptrdiff_t>::const_iterator
    it = indices.begin(), end = indices.end();
    it != end; ++it){
    v_hint = dt2.insert(points[*it], hint);
    if (v_hint!=typename DT2::Vertex_handle()){
      v_hint->info()=infos[*it];
      hint=v_hint->face();
    }
  }

}


/*!
Version taking a range of points and a range of infos and creating the HDS cropped to the iso-rectangle.
The size of both ranges are the same and the info of each point is associated to the corresponding face
\todo Create_face_from_info should be a template parameter
*/
template <  class Input_kernel,
            class Exact_kernel,
            class PointIterator,
            class InfoIterator,
            class HDS >
void create_hds_for_cropped_voronoi_diagram(
  PointIterator point_begin,
  PointIterator point_end,
  InfoIterator info_begin,
  InfoIterator info_end,
  const typename Input_kernel::Iso_rectangle_2& iso_rect,
  HDS& hds )
{
  typedef typename std::iterator_traits<InfoIterator>::value_type Vertex_info;
  typedef CGAL::Triangulation_vertex_base_with_info_2<Vertex_info, Input_kernel>     Vb;
  typedef Face_info_for_DT2<Input_kernel, Exact_kernel> Face_info;
  typedef CGAL::Triangulation_face_base_with_info_2<Face_info, Input_kernel> Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                       Tds;
  typedef CGAL::Delaunay_triangulation_2<Input_kernel, Tds>                  DT2;

  DT2 dt2;
  //workaround to fix a bug in boost 1.46 with g++ 4.12
  dt2_insert_with_info(dt2, point_begin, point_end, info_begin, info_end);
  //dt2.insert(
  //  boost::make_zip_iterator( boost::make_tuple(point_begin, info_begin) ),
  //  boost::make_zip_iterator( boost::make_tuple(point_end, info_end) ) );

    create_hds_for_cropped_voronoi_diagram<Input_kernel, Exact_kernel>(dt2,
                                                                     iso_rect,
                                                                     hds,
                                                                     Create_face_from_info<DT2, HDS>() );
}

/*!
Version taking a range of points and creating the HDS cropped to the default iso-rectangle
*/
template <  class Input_kernel,
            class Exact_kernel,
            class PointIterator,
            class HDS >
void create_hds_for_cropped_voronoi_diagram(
  PointIterator point_begin,
  PointIterator point_end,
  const typename Input_kernel::Iso_rectangle_2& iso_rect,
  HDS& hds)
{
  typedef CGAL::Triangulation_vertex_base_2<Input_kernel>     Vb;
  typedef Face_info_for_DT2<Input_kernel, Exact_kernel> Face_info;
  typedef CGAL::Triangulation_face_base_with_info_2<Face_info, Input_kernel> Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                       Tds;
  typedef CGAL::Delaunay_triangulation_2<Input_kernel, Tds>                  DT2;

  DT2 dt2( point_begin, point_end );

  create_hds_for_cropped_voronoi_diagram<Input_kernel, Exact_kernel>( dt2,
                                                                      iso_rect,
                                                                      hds,
                                                                      Default_create_face<DT2,HDS>());
}

/*!
Version taking a range of points and a range of infos and creating the HDS cropped to the default iso-rectangle.
The size of both ranges are the same and the info of each point is associated to the corresponding face
\todo Create_face_from_info should be a template parameter
*/
template <  class Input_kernel,
            class Exact_kernel,
            class PointIterator,
            class InfoIterator,
            class HDS >
const typename Input_kernel::Iso_rectangle_2&
create_hds_for_cropped_voronoi_diagram(
  PointIterator point_begin,
  PointIterator point_end,
  InfoIterator info_begin,
  InfoIterator info_end,
  HDS& hds )
{
  typedef typename std::iterator_traits<InfoIterator>::value_type Vertex_info;
  typedef CGAL::Triangulation_vertex_base_with_info_2<Vertex_info, Input_kernel>     Vb;
  typedef Face_info_for_DT2<Input_kernel, Exact_kernel> Face_info;
  typedef CGAL::Triangulation_face_base_with_info_2<Face_info, Input_kernel> Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                       Tds;
  typedef CGAL::Delaunay_triangulation_2<Input_kernel, Tds>                  DT2;

  DT2 dt2;
  //workaround to fix a bug in boost 1.46 with g++ 4.12
  dt2_insert_with_info(dt2, point_begin, point_end, info_begin, info_end);
  //dt2.insert(
  //  boost::make_zip_iterator( boost::make_tuple(point_begin, info_begin) ),
  //  boost::make_zip_iterator( boost::make_tuple(point_end, info_end) ) );

  return
    create_hds_for_cropped_voronoi_diagram<Input_kernel, Exact_kernel>( dt2,
                                                                          hds,
                                                                          Create_face_from_info<DT2, HDS>() );
}

/*!
Version taking a range of points and creating the HDS cropped to the default iso-rectangle
*/
template <  class Input_kernel,
            class Exact_kernel,
            class PointIterator,
            class HDS >
const typename Input_kernel::Iso_rectangle_2&
create_hds_for_cropped_voronoi_diagram(
  PointIterator point_begin,
  PointIterator point_end,
  HDS& hds)
{
  typedef CGAL::Triangulation_vertex_base_2<Input_kernel>     Vb;
  typedef Face_info_for_DT2<Input_kernel, Exact_kernel> Face_info;
  typedef CGAL::Triangulation_face_base_with_info_2<Face_info, Input_kernel> Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                       Tds;
  typedef CGAL::Delaunay_triangulation_2<Input_kernel, Tds>                  DT2;

  DT2 dt2( point_begin, point_end );

  return
    create_hds_for_cropped_voronoi_diagram<Input_kernel, Exact_kernel>( dt2,
                                                                        hds,
                                                                        Default_create_face<DT2,HDS>());
}

#ifdef DEBUG_JOIN_FACES
template <class HDS, class MapHandle, class UF>
void write_hds_debug(HDS& hds, MapHandle& faces_to_handles, UF& faces_uf, const char* fname)
{
  std::ofstream output(fname);
  for( typename HDS::Face_iterator fit=hds.faces_begin();fit!=hds.faces_end();++fit)
  {
    std::stringstream sstream;

    typename HDS::Face_handle face_master=*faces_uf.find( faces_to_handles.find(fit)->second );
    if ( fit!=face_master ) continue; //do it only once

    //the halfedge pointer of the face might be incorrect, we look for a valid one.
    typename HDS::Halfedge_iterator hit=hds.halfedges_begin();
    for( ;hit!=hds.halfedges_end();++hit)
    {
      if ( hit->face()==typename HDS::Face_handle() ) continue; //hedge to be deleted
      typename HDS::Face_handle hface =
        *faces_uf.find( faces_to_handles.find(hit->face())->second ) ;
      if (hface==fit) break;
    }
    CGAL_assertion( hit!=hds.halfedges_end() );
    typename HDS::Halfedge_handle curr=hit, end=curr;

    int i=0;
    do{
      CGAL_assertion( curr!= typename HDS::Halfedge_handle() );
      //CGAL_assertion( curr->face()!=typename HDS::Face_handle() );
      //CGAL_assertion( curr->face()==fit );
      ++i;
      sstream << curr->vertex()->point() << " 0\n";
      curr=curr->next();
    }while(curr!=end);
    sstream << curr->vertex()->point() << " 0\n";
    output << i+1 << " " << sstream.str() << std::endl;
  }
  output.close();
}
#endif

/// this function remove halfedge of length 0
template <class HDS>
void clean_up_hds(HDS&hds)
{
  std::vector<typename HDS::Halfedge_handle> hedges_to_remove;
  for (typename HDS::Halfedge_iterator hit=hds.halfedges_begin(),
                                       hit_end=hds.halfedges_end(); hit!=hit_end;
                                       std::advance(hit,2))
  {
    if ( hit->vertex()->point() == hit->opposite()->vertex()->point() )
      hedges_to_remove.push_back(hit);
  }

  CGAL::HalfedgeDS_decorator<HDS> D(hds);

  std::size_t nb_hedges=hedges_to_remove.size();
  for (std::size_t i=0; i<nb_hedges; ++i)
  {
    D.join_vertex(hedges_to_remove[i]);
  }
}


template<class Halfedge_handle>
CGAL::Bbox_2
bbox_of_halfedges_sequences(Halfedge_handle h)
{
  CGAL::Bbox_2 bbox=h->vertex()->point().bbox();
  Halfedge_handle hit=h->next();
  while(h!=hit)
  {
    bbox=bbox+hit->vertex()->point().bbox();
    hit=hit->next();
  }
  return bbox;
}

#include <CGAL/Union_find.h>

template <class HDS>
void join_faces_with_same_color(HDS& hds, bool do_clean_hds=false)
{
  // since the voronoi vertices are created from exact versions,
  // it might happen during the rounding that two vertices that are
  // close become identical
  if (do_clean_hds) clean_up_hds(hds);

  // collect the set of halfedges that need to be removed: incident faces have the same color
  std::vector<typename HDS::Halfedge_handle> hedges_to_remove;
  for (typename HDS::Halfedge_iterator hit=hds.halfedges_begin(),
                                       hit_end=hds.halfedges_end(); hit!=hit_end;
                                       std::advance(hit,2))
  {
    if ( hit->is_border() || hit->opposite()->is_border() ) continue;
    if ( hit->face()->color()==-1 || hit->opposite()->face()->color()==-1 ) continue;
    if (hit->face()->color() == hit->opposite()->face()->color() )
      hedges_to_remove.push_back(hit);
  }

  // we use a union-find structure to keep track of the merge of faces
  typedef CGAL::Union_find< typename HDS::Face_handle > Faces_UF;
  Faces_UF faces_uf;
  std::map<typename HDS::Face_handle, typename Faces_UF::handle> faces_to_handles;

  //each face is placed into one component
  for (typename HDS::Face_iterator fit=hds.faces_begin(), fit_end=hds.faces_end();fit!=fit_end;++fit)
    faces_to_handles[fit]=faces_uf.make_set(fit);

  //save the creation of the null halfedge
  const typename HDS::Halfedge_handle null_halfedge = typename HDS::Halfedge_handle();

  //set decorators for hds editing
  CGAL::HalfedgeDS_decorator<HDS> D(hds);
  CGAL::HalfedgeDS_items_decorator<HDS> item_decorator;

  //main loop removing halfedges one by one
  std::size_t nb_hedges=hedges_to_remove.size();
  for (std::size_t i=0; i<nb_hedges; ++i)
  {
    #ifdef DEBUG_JOIN_FACES
    std::stringstream ss;
    if (i<10)
      ss << "debug/hds-0"<<i<<".cgal";
    else
      ss << "debug/hds-"<<i<<".cgal";
    write_hds_debug(hds, faces_to_handles, faces_uf , ss.str().c_str());

    std::cout << "current edge ("<<i<<") is " << hedges_to_remove[i]->vertex()->point() << " 0 " << hedges_to_remove[i]->opposite()->vertex()->point() << " 0\n";
    #endif

    //check if the halfedge removal has already been taken care of:
    //we set incident faces to null in that case
    if ( hedges_to_remove[i]->is_border() &&
         hedges_to_remove[i]->opposite()->is_border() )
    {
      hds.edges_erase(hedges_to_remove[i]);
      continue;
    }

    //what we call a dangling edge is a sequence of edges inside a face. It can be attached to one
    //or two connected component of the boundary of the face
    bool is_dangling_edge = faces_uf.same_set(
      faces_to_handles[ hedges_to_remove[i]->face()],
      faces_to_handles[ hedges_to_remove[i]->opposite()->face() ]);

    //what we call an isolated cycle is a hole in a given face
    bool is_isolated_cycle=false;


    //these two halfedges are used to bound the sequence of halfedges incident to the same
    //faces and that need to be removed together.
    //After the two while-loops, the target vertex of target is either connected
    //to 1 or more that 3 edges. The same holds for the source vertex of source;
    typename HDS::Halfedge_handle target = hedges_to_remove[i], source=target;

    //check if the valence of the target vertex is 2
    while ( target->opposite() == target->next()->opposite()->next() &&
            ( !is_dangling_edge || target->next()!=target->opposite() ) )//if the dangling edge is attached at one vertex only
    {
      target=target->next();
      if (source==target)
      {
        is_isolated_cycle=true;
        break;
      }
    }

    //check if the valence of the source vertex is 2
    if (!is_isolated_cycle)
      while ( source == source->opposite()->next()->opposite()->next() &&
             ( !is_dangling_edge || source!=source->opposite()->next() ) )//if the dangling edge is attached at one vertex only
        source=source->opposite()->next()->opposite();

    //in this case it is safe to remove the edge
    if (source==target && !is_dangling_edge && !is_isolated_cycle)
    {
      #ifdef DEBUG_JOIN_FACES
      std::cout << "remove face\n";
      #endif
      //copy-paste of join_face function from the HalfedgeDS_decorator
      //we just skip the face removal as we need it for the union-find
      typename HDS::Halfedge_handle h=hedges_to_remove[i];
      typename HDS::Halfedge_handle hprev = D.find_prev( h);
      typename HDS::Halfedge_handle gprev = D.find_prev( h->opposite());
      D.remove_tip( hprev);
      D.remove_tip( gprev);

      //put the two adjacent faces in the same component
      faces_uf.unify_sets(
        faces_to_handles[h->face()],
        faces_to_handles[h->opposite()->face()] );
      //hds.faces_erase( D.get_face( gprev));

      hds.edges_erase( h);
      h = hprev;
      // 'half' of the halfedges have their correct faces.
      // Here we do the remaining halfedges.
      while ( h != gprev) h = h->next();

      D.set_vertex_halfedge( hprev);
      D.set_vertex_halfedge( gprev);


      //D.join_face( hedges_to_remove[i] );
    }
    else
    {
      #ifdef DEBUG_JOIN_FACES
      std::cout << "vertex " << target->opposite()->vertex()->point() << " 0\n";
      std::cout << "2 " << target->vertex()->point() << " 0 " << source->opposite()->vertex()->point() << " 0\n";
      std::cout << "target " << target->opposite()->vertex()->point() << " 0 " << target->vertex()->point() << " 0\n";
      std::cout << "source " << source->opposite()->vertex()->point() << " 0 " << source->vertex()->point() << " 0\n";
      #endif

      if (!is_isolated_cycle)
      {
        if( target->vertex() == source->opposite()->vertex() )
        {
          //in this case we have an attached circle
          #ifdef DEBUG_JOIN_FACES
          std::cout << "FOUND ONE CYCLE ATTACHED\n";
          #endif
          if (target->next() != source)
          {
            typename HDS::Halfedge_handle tmp=source->opposite();
            source=target->opposite();
            target=tmp;
          }

          CGAL_assertion( target->next() == source );
          //hds.faces_erase( target->face() );
          faces_uf.unify_sets(
            faces_to_handles[target->face()],
            faces_to_handles[target->opposite()->face()] );

          item_decorator.get_prev(target->opposite())->set_next( source->opposite()->next() );
          item_decorator.set_prev( source->opposite()->next(),
                                   item_decorator.get_prev(target->opposite()) );
          //update vertex hedge pointer
          item_decorator.set_vertex_halfedge( item_decorator.get_prev(target->opposite()) );
          //update face pointer
          //item_decorator.set_face_halfedge( source->opposite()->next() );
        }
        else{
          typename HDS::Halfedge_handle target_next = target->next();
          const typename HDS::Halfedge_handle target_opp_prev = item_decorator.get_prev(target->opposite());
          const typename HDS::Halfedge_handle source_prev = item_decorator.get_prev(source);
          const typename HDS::Halfedge_handle source_opp_next = source->opposite()->next();

          //update next/prev around target
          if( target_next != target->opposite() ) //not a pending edge
          {
            item_decorator.set_prev(target_next, target_opp_prev);
            target_opp_prev->set_next( target_next );
            item_decorator.set_vertex_halfedge( target_opp_prev );
            #ifdef DEBUG_JOIN_FACES
            std::cout << "Linking " << target_opp_prev->opposite()->vertex()->point() << " 0 " << target_opp_prev->vertex()->point() << " 0"
                      << " to "     << target_next->opposite()->vertex()->point() << " 0 " << target_next->vertex()->point() << " 0" << std::endl;
            #endif
            //update face halfedge pointer in case it was removed
            //item_decorator.set_face_halfedge(face_kept, target_next);
          }
          else
          {
            #ifdef DEBUG_JOIN_FACES
            std::cout << "delete vertex a\n";
            #endif
            hds.vertices_erase(target->vertex());
          }


          //update next/prev around source
          if(source != source_opp_next)
          {
            source_prev->set_next( source_opp_next );
            item_decorator.set_vertex_halfedge( source_prev );
            item_decorator.set_prev(source_opp_next, source_prev);
            #ifdef DEBUG_JOIN_FACES
            std::cout << "Linking " << source_prev->opposite()->vertex()->point() << " 0 " << source_prev->vertex()->point() << " 0"
                      << " to "     << source_opp_next->opposite()->vertex()->point() << " 0 " << source_opp_next->vertex()->point() << " 0" << std::endl;
            #endif
            //update face halfedge pointer in case it was removed
            //item_decorator.set_face_halfedge(face_kept, source_prev);
          }
          else
          {
            #ifdef DEBUG_JOIN_FACES
            std::cout << "delete vertex b\n";
            #endif
            hds.vertices_erase(source->opposite()->vertex());
          }

          if ( !is_dangling_edge )
          {
            //hds.faces_erase( target->face() );
            faces_uf.unify_sets(
              faces_to_handles[target->face()],
              faces_to_handles[target->opposite()->face()] );

            //update face pointer of edges of the face removed
            while(target_next!=source_opp_next)
            {
              //item_decorator.set_face(target_next, face_kept);
              target_next=target_next->next();
            }
          }
        }
      }
      else
      {
        #ifdef DEBUG_JOIN_FACES
        std::cout << "FOUND A CYCLE\n";
        std::cout << "delete vertex d\n";
        #endif
        CGAL_assertion(source==target);
        //hds.faces_erase( target->face() );
        faces_uf.unify_sets(
          faces_to_handles[target->face()],
          faces_to_handles[target->opposite()->face()] );
        source=source->next();

        hds.vertices_erase(target->vertex()); //remove target vertex not removed in the loop below
      }

      //set to null the pointer of edges to remove
      //and remove vertices in the middle
      while(source!=target)
      {
        item_decorator.set_face(source, typename HDS::Face_handle());
        item_decorator.set_face(source->opposite(), typename HDS::Face_handle());
        hds.vertices_erase(source->vertex());
        #ifdef DEBUG_JOIN_FACES
        std::cout << "delete vertex c\n";
        #endif
        source=source->next();
      }
      item_decorator.set_face(target, typename HDS::Face_handle());
      item_decorator.set_face(target->opposite(), typename HDS::Face_handle());

      //do remove current edge
      hds.edges_erase(hedges_to_remove[i]);
    }
  }

  //once all edges have been removed, we need to update the face pointers
  for (typename HDS::Face_iterator fit=hds.faces_begin(), fit_end=hds.faces_end();
                                   fit!=fit_end; ++fit)
  {
    //reset face halfedge to indicate that they have not been processed
    item_decorator.set_face_halfedge(fit, null_halfedge );
  }

  typename HDS::Face_iterator next_available_face=hds.faces_begin();

  std::set<typename HDS::Halfedge_handle> hedge_visited;
  std::vector< std::pair<int, typename HDS::Halfedge_handle> > face_info;

  for (typename HDS::Halfedge_iterator hit=hds.halfedges_begin(), hit_end=hds.halfedges_end();
                                       hit!=hit_end; ++hit)
  {
    if (hit->is_border()) continue;
    if ( !hedge_visited.insert(hit).second ) continue; //already visited

    typename HDS::Face_handle face=*faces_uf.find( faces_to_handles[hit->face()] );
    item_decorator.set_face(hit, face); //set the face of the halfedge
    if (face->halfedge() == null_halfedge)
      item_decorator.set_face_halfedge(face, hit); //set the first as the outer ccb
    else
      face->holes.push_back(hit); //this is a hole, the face "master" halfedge is already set

    //set the face of the halfedges
    typename HDS::Halfedge_handle h=hit->next();
    while(h!=hit){
      hedge_visited.insert(h);
      item_decorator.set_face(h, face);
      h=h->next();
    }
  }

  for (typename HDS::Face_iterator fit=hds.faces_begin(); fit!=hds.faces_end(); )
  {
    if (fit->halfedge() == null_halfedge)
      hds.faces_erase(fit++);
    else
    {
      if ( !fit->holes.empty() )
      {
        /// Ensure that the connected component of edges (CCE) of the halfedge
        /// associated to the face is the one "containing" all the holes.
        /// To check this, we set the CCE to be the one associated to the face
        /// if its bounding box contains that of all other CCE.
        typename HDS::Halfedge_handle halfedge = fit->halfedge();
        CGAL::Bbox_2 max_bbox=bbox_of_halfedges_sequences(halfedge);
        std::size_t nb_holes = fit->holes.size();
        for( std::size_t i=0; i<nb_holes; ++i)
        {
          CGAL::Bbox_2 bbox =
            bbox_of_halfedges_sequences( fit->holes[i] );
          //there is only one bbox that contains them all
          if (bbox+max_bbox == bbox)
          {
            std::swap(halfedge, fit->holes[i]);
            max_bbox=bbox;
          }
        }
        item_decorator.set_face_halfedge(fit, halfedge);
      }
      ++fit;
    }
  }
}
