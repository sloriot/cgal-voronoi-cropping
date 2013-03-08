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

template<class HDS>
void set_border_hedge_source(
  typename HDS::Halfedge_handle h,
  typename HDS::Vertex_handle v,
  std::map<typename HDS::Vertex_handle, std::pair<typename HDS::Halfedge_handle, typename HDS::Halfedge_handle> >& border_hedges
)
{
  border_hedges[v].first=h;
}

template<class HDS>
void set_border_hedge_target(
  typename HDS::Halfedge_handle h,
  typename HDS::Vertex_handle v,
  std::map<typename HDS::Vertex_handle, std::pair<typename HDS::Halfedge_handle, typename HDS::Halfedge_handle> >& border_hedges
)
{
  border_hedges[v].second=h;
}

template <class HDS, class Iso_rectangle_2>
void
add_halfedge_on_iso_rectangle(
  HDS& hds,
  typename HDS::Face_handle new_face,
  typename HDS::Halfedge_handle previous,
  typename HDS::Halfedge_handle next,
  const Iso_rectangle_2& bbox,
  std::map<typename HDS::Vertex_handle, std::pair<typename HDS::Halfedge_handle, typename HDS::Halfedge_handle> >& border_hedges
)
{
  //ensure target of previous is on the bbox boundary
  CGAL_assertion(
    bbox.xmin()==previous->vertex()->point().x() ||  bbox.xmax()==previous->vertex()->point().x() ||
    bbox.ymin()==previous->vertex()->point().y() ||  bbox.ymax()==previous->vertex()->point().y()
  );

  //ensure source of next is on the bbox boundary
  CGAL_assertion(
    bbox.xmin()==next->opposite()->vertex()->point().x() ||  bbox.xmax()==next->opposite()->vertex()->point().x() ||
    bbox.ymin()==next->opposite()->vertex()->point().y() ||  bbox.ymax()==next->opposite()->vertex()->point().y()
  );

  //ensure the points are different
  /// \todo remove this as it can be true for a voronoi vertex on the boundary
  /// \todo add something to make the vertices identical if the points are exactly equal: this can be done by checking the orientation of the dual voronoi vertex and storing the hds vertex in the hds vertex
  /// \todo and also do a link of hedges and return in case of identical vertex
  CGAL_assertion( next->opposite()->vertex()->point()!=previous->vertex()->point() );

  //bbox[0], bbox[1], bbox[2], bbox[3] are counterclockwise oriented
  double index_target=index_on_bbox(previous->vertex()->point(),bbox);
  double index_next=index_on_bbox(next->opposite()->vertex()->point(),bbox);

  if (index_next < index_target) index_next+=4;
  int k = std::floor( index_target);

  CGAL::HalfedgeDS_items_decorator<HDS> item_decorator;
  while (++k <= index_next)
  {
    typename HDS::Vertex_handle new_vertex =
      hds.vertices_push_back( typename HDS::Vertex( bbox[k] ) );

    typename HDS::Halfedge_handle new_hh=hds.edges_push_back(typename HDS::Halfedge(), typename HDS::Halfedge());
    item_decorator.set_face(new_hh,new_face);
    item_decorator.set_vertex(new_hh,new_vertex);
    set_border_hedge_target<HDS>(new_hh->opposite(), previous->vertex(), border_hedges);
    set_border_hedge_source<HDS>(new_hh->opposite(), new_vertex, border_hedges);
    item_decorator.set_vertex(new_hh->opposite(),previous->vertex());

    previous->set_next(new_hh);
    item_decorator.set_prev(new_hh,previous);

    previous=new_hh;
  }

  typename HDS::Halfedge_handle new_hh=hds.edges_push_back(typename HDS::Halfedge(), typename HDS::Halfedge());
  item_decorator.set_face(new_hh,new_face);
  item_decorator.set_vertex(new_hh,next->opposite()->vertex());
  item_decorator.set_vertex(new_hh->opposite(),previous->vertex());

  set_border_hedge_target<HDS>(new_hh->opposite(), previous->vertex(), border_hedges);
  set_border_hedge_source<HDS>(new_hh->opposite(), new_hh->vertex(), border_hedges);

  previous->set_next(new_hh);
  item_decorator.set_prev(new_hh,previous);

  new_hh->set_next(next);
  item_decorator.set_prev(next,new_hh);
}



#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <bitset>

template <class Input_kernel,class Exact_kernel>
class Face_info_for_DT2
{
  std::bitset<4> m_orientation_vs_bbox;
  std::bitset<4> m_on_boundary;
  typename Exact_kernel::Point_2 m_dual;
public:
  void set_orientation(int i, CGAL::Orientation o)
  {
    switch (o)
    {
      case CGAL::CLOCKWISE:
        m_orientation_vs_bbox.set(i);
        break;
      case CGAL::COLLINEAR:
        //~ m_orientation_vs_bbox.set(i+4);
        std::cout << "inexact_dual " << inexact_dual() << std::endl;
        m_on_boundary.set(i);
      default:
        //do nothing
        break;
    }
  }

  template <class Face_handle>
  void compute_dual(Face_handle fh)
  {
    CGAL::Cartesian_converter<Input_kernel,Exact_kernel> to_exact;
    m_dual = typename Exact_kernel::Construct_circumcenter_2()(
                        to_exact( fh->vertex(0)->point() ),
                        to_exact( fh->vertex(1)->point() ),
                        to_exact( fh->vertex(2)->point() )
    );
  }

  const typename Exact_kernel::Point_2&
  dual() const
  {
    return m_dual;
  }

  const typename Input_kernel::Point_2
  inexact_dual() const
  {
    CGAL::Cartesian_converter<Exact_kernel, Input_kernel> to_input;
    return to_input( m_dual );
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

template <class HDS, class Face_handle, class Vertex_map>
typename HDS::Vertex_handle
get_dual_vertex(  Face_handle face,
                  HDS& hds,
                  Vertex_map& vertex_map )
{
  CGAL_precondition( face->info().is_dual_in_bbox() );
  std::pair<typename Vertex_map::iterator, bool> insert_res =
    vertex_map.insert( std::make_pair(face, typename HDS::Vertex_handle()) );
  if (insert_res.second)
  {
    insert_res.first->second =
      hds.vertices_push_back(
        typename HDS::Vertex( face->info().inexact_dual() ) );
  }
  return insert_res.first->second;
}

template <class Input_kernel, class Exact_kernel, class Primitive, class Bisector>
bool crop_to_iso_rectangle( const Primitive& primitive,
                            const Bisector& bisector,
                            const typename Exact_kernel::Iso_rectangle_2& exact_iso_rect,
                            typename Exact_kernel::Segment_2& res)
{
  CGAL::Object obj =
    typename Exact_kernel::Intersect_3()( primitive, exact_iso_rect );
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

template <class Input_kernel, class Exact_kernel>
typename Input_kernel::Point_2
convert_to_input(const typename Exact_kernel::Point_2& p, const typename Exact_kernel::Iso_rectangle_2& bbox)
{
  CGAL::NT_converter<typename Exact_kernel::FT, typename Input_kernel::FT> to_input;
  typename Input_kernel::FT x,y;

  //handle x coordinate
  if ( p.x() == bbox.xmin() )
    x=to_input(bbox.xmin());
  else{
    if ( p.x() == bbox.xmax() ) x=to_input(bbox.xmax());
    else x=to_input(p.x());
  }

  //handle y coordinate
  if ( p.y() == bbox.ymin() )
    y=to_input(bbox.ymin());
  else{
    if ( p.y() == bbox.ymax() ) y=to_input(bbox.ymax());
    else y=to_input(p.y());
  }

  return typename Input_kernel::Point_2(x,y);
}

template<class Input_kernel, class Exact_kernel, class HDS, class DT2, class Hedge_map, class Vertex_map>
typename HDS::Halfedge_handle
get_dual_halfedge_cropped(  typename DT2::Face_handle source_face,
                            typename DT2::Face_handle target_face,
                            int edge_index_in_source_face,
                            const typename CGAL::Iso_rectangle_2<Exact_kernel>& exact_iso_rect,
                            HDS& hds,
                            Hedge_map& hedge_map,
                            Vertex_map& vertex_map,
                            const DT2& dt2)
{
  typedef typename HDS::Halfedge_handle Halfedge_handle;
  typedef typename HDS::Vertex_handle Vertex_handle;
  typedef typename DT2::Face_handle Face_handle;
  CGAL::HalfedgeDS_items_decorator<HDS> item_decorator;

  std::pair<Face_handle, Face_handle> sorted_pair=
    make_sorted_pair(source_face, target_face);

  std::pair< typename Hedge_map::iterator, bool> insert_res =
    hedge_map.insert( std::make_pair( sorted_pair, Halfedge_handle() ) );

  Vertex_handle source_vertex, target_vertex;
  CGAL::Cartesian_converter<Input_kernel, Exact_kernel> to_exact;

  if ( insert_res.second )
  {
    if ( dt2.is_infinite(source_face) )
    {
      if (   dt2.is_infinite( source_face->vertex(DT2::cw(edge_index_in_source_face)) )
          || dt2.is_infinite( source_face->vertex(DT2::ccw(edge_index_in_source_face)) ) )
      {
        return typename HDS::Halfedge_handle();
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
        if ( !crop_to_iso_rectangle<Input_kernel, Exact_kernel>(bisector, bisector,exact_iso_rect, inter) )
          return typename HDS::Halfedge_handle();
        source_vertex=hds.vertices_push_back( typename HDS::Vertex( convert_to_input<Input_kernel, Exact_kernel>( inter.source(), exact_iso_rect )) );
        target_vertex=hds.vertices_push_back( typename HDS::Vertex( convert_to_input<Input_kernel, Exact_kernel>( inter.target(), exact_iso_rect )) );
        source_vertex->is_on_border()=true;
        target_vertex->is_on_border()=true;
      }
      else
      {
        //the dual is a ray
        typename Exact_kernel::Ray_2 ray( target_face->info().dual(), bisector.opposite() );
        typename Exact_kernel::Segment_2 inter;
        if ( !crop_to_iso_rectangle<Input_kernel, Exact_kernel>(ray, bisector,exact_iso_rect, inter) )
          return typename HDS::Halfedge_handle();
        if ( target_face->info().is_dual_in_bbox() )
        {
          target_vertex=get_dual_vertex(target_face, hds, vertex_map);
          source_vertex=hds.vertices_push_back( typename HDS::Vertex( convert_to_input<Input_kernel, Exact_kernel>( inter.source(), exact_iso_rect )) );
          source_vertex->is_on_border()=true;
        }
        else{
          source_vertex=hds.vertices_push_back( typename HDS::Vertex( convert_to_input<Input_kernel, Exact_kernel>( inter.source(), exact_iso_rect )) );
          target_vertex=hds.vertices_push_back( typename HDS::Vertex( convert_to_input<Input_kernel, Exact_kernel>( inter.target(), exact_iso_rect )) );
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
          return typename HDS::Halfedge_handle();
        }

        typename Exact_kernel::Line_2 bisector =
          typename Exact_kernel::Construct_bisector_2()(
              to_exact( source_face->vertex( DT2::cw(edge_index_in_source_face) )->point() ),
              to_exact( source_face->vertex( DT2::ccw(edge_index_in_source_face) )->point() ) );
        CGAL_assertion( bisector.has_on_positive_side(
          to_exact( source_face->vertex( DT2::cw(edge_index_in_source_face) )->point() )  ) );

        //the dual is a ray
        typename Exact_kernel::Ray_2 ray( source_face->info().dual(), bisector );
        typename Exact_kernel::Segment_2 inter;
        if ( !crop_to_iso_rectangle<Input_kernel, Exact_kernel>(ray, bisector,exact_iso_rect, inter) )
          return typename HDS::Halfedge_handle();
        if ( source_face->info().is_dual_in_bbox() )
        {
          source_vertex=get_dual_vertex(source_face, hds, vertex_map);
          target_vertex=hds.vertices_push_back( typename HDS::Vertex( convert_to_input<Input_kernel, Exact_kernel>( inter.target(), exact_iso_rect )) );
          target_vertex->is_on_border()=true;
        }
        else{
          source_vertex=hds.vertices_push_back( typename HDS::Vertex( convert_to_input<Input_kernel, Exact_kernel>( inter.source(), exact_iso_rect )) );
          target_vertex=hds.vertices_push_back( typename HDS::Vertex( convert_to_input<Input_kernel, Exact_kernel>( inter.target(), exact_iso_rect )) );
          source_vertex->is_on_border()=true;
          target_vertex->is_on_border()=true;
        }
      }
      else
      {
        //the dual is a segment
        if ( source_face->info().is_dual_in_bbox() &&
             target_face->info().is_dual_in_bbox() )
        {
          source_vertex = get_dual_vertex(source_face, hds, vertex_map);
          target_vertex = get_dual_vertex(target_face, hds, vertex_map);
        }
        else{
          typename Exact_kernel::Segment_2 segment( source_face->info().dual(), target_face->info().dual() );
          typename Exact_kernel::Segment_2 inter(typename Exact_kernel::Point_2(66,66), typename Exact_kernel::Point_2(66,66));
          if ( !crop_to_iso_rectangle<Input_kernel, Exact_kernel>(segment, segment,exact_iso_rect, inter) )
            return typename HDS::Halfedge_handle();


          if ( source_face->info().is_dual_in_bbox() )
            source_vertex = get_dual_vertex(source_face, hds, vertex_map);
          else
          {
            source_vertex = hds.vertices_push_back( typename HDS::Vertex( convert_to_input<Input_kernel, Exact_kernel>( inter.source(), exact_iso_rect )) );
            source_vertex->is_on_border()=true;
          }
          if ( target_face->info().is_dual_in_bbox() )
            target_vertex = get_dual_vertex(target_face, hds, vertex_map);
          else
          {
            target_vertex = hds.vertices_push_back( typename HDS::Vertex( convert_to_input<Input_kernel, Exact_kernel>( inter.target(), exact_iso_rect )) );
            target_vertex->is_on_border()=true;
          }
        }
      }
    }

    CGAL_assertion( source_vertex != typename HDS::Vertex_handle() );
    CGAL_assertion( target_vertex != typename HDS::Vertex_handle() );

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
  /// \todo handle the case of a degenerate input set
  CGAL::Cartesian_converter<Input_kernel, Exact_kernel> to_exact;
  typename Exact_kernel::Orientation_2 orientation;
  typename Exact_kernel::Iso_rectangle_2 exact_iso_rect = to_exact( iso_rect );

  for (typename DT2::Finite_faces_iterator  fit=dt2.finite_faces_begin(),
                                            fit_end=dt2.finite_faces_end();
                                            fit!=fit_end; ++fit )
  {
    fit->info().compute_dual(fit);
    for (int i=0;i<4;++i)
      fit->info().set_orientation(
              i,
              orientation(
                fit->info().dual(),
                exact_iso_rect[i],
                exact_iso_rect[i+1] ) );
  }

  CGAL::HalfedgeDS_items_decorator<HDS> item_decorator;

  std::map< std::pair< typename DT2::Face_handle,
                       typename DT2::Face_handle >,
            typename HDS::Halfedge_handle > dual_hedge_map;
  std::map< typename DT2::Face_handle, typename HDS::Vertex_handle > vertex_map;

  typename HDS::Halfedge_handle null_halfedge=typename HDS::Halfedge_handle();

  std::map<typename HDS::Vertex_handle, std::pair<typename HDS::Halfedge_handle, typename HDS::Halfedge_handle> > border_hedges;

  for (typename DT2::Finite_vertices_iterator  vit=dt2.finite_vertices_begin(),
                                               vit_end=dt2.finite_vertices_end();
                                               vit!=vit_end; ++vit )
  {
    typename DT2::Face_handle current_face=vit->face(), first_face=current_face;
    std::vector<typename HDS::Halfedge_handle> new_hedges;

    do{
      int other_vertex_index=DT2::ccw( current_face->index(vit) );
      typename DT2::Face_handle neighbor_face=current_face->neighbor(other_vertex_index);

      typename HDS::Halfedge_handle hedge =
        get_dual_halfedge_cropped<Input_kernel>(current_face, neighbor_face,
                                                other_vertex_index, exact_iso_rect,
                                                hds, dual_hedge_map, vertex_map, dt2);

      if (hedge!=null_halfedge) new_hedges.push_back(hedge);

      current_face=neighbor_face;
    }
    while(current_face!=first_face);

    if ( !new_hedges.empty() )
    {
      std::size_t nb_hedges=new_hedges.size();
      new_hedges.push_back(new_hedges[0]);
      typename HDS::Face_handle new_face=hds.faces_push_back( create_face(vit) );
      CGAL_assertion( new_hedges[0]!=typename HDS::Halfedge_handle() );
      item_decorator.set_face_halfedge(new_face,new_hedges[0]); //set one halfedge for the face

      for (std::size_t i=0;i<nb_hedges;++i)
      {
        item_decorator.set_face(new_hedges[i],new_face);

        if ( new_hedges[i]->vertex()->is_on_border() )
        {
          CGAL_assertion( new_hedges[i+1]->opposite()->vertex()->is_on_border() );
          //we need one to several hedges on the iso-rectangle boundary to link the two vertices
          add_halfedge_on_iso_rectangle(hds, new_face, new_hedges[i], new_hedges[i+1], iso_rect, border_hedges);
        }
        else
        {
          new_hedges[i]->set_next(new_hedges[i+1]);
          item_decorator.set_prev(new_hedges[i+1], new_hedges[i]);
        }
      }
    }
  }

  //link border halfedges
  for (
    typename std::map<typename HDS::Vertex_handle, std::pair<typename HDS::Halfedge_handle, typename HDS::Halfedge_handle> >::iterator
    link_it=border_hedges.begin(), link_it_end=border_hedges.end();link_it!=link_it_end; ++link_it)
  {
    CGAL_assertion(
      link_it->second.first->opposite()->vertex()==link_it->second.second->vertex() );
    link_it->second.second->set_next(link_it->second.first);
    item_decorator.set_prev(link_it->second.first, link_it->second.second);
    item_decorator.set_vertex_halfedge(link_it->first, link_it->second.second);
  }

  /// \todo test me + add the vertex info!!!!
  if (border_hedges.empty())
  {
    typename DT2::Vertex_handle vh;
    typename HDS::Face_handle new_face=hds.faces_push_back( create_face(vh) );
    typename HDS::Vertex_handle vertices[4] = {
      hds.vertices_push_back( typename HDS::Vertex( iso_rect[0] ) ),
      hds.vertices_push_back( typename HDS::Vertex( iso_rect[1] ) ),
      hds.vertices_push_back( typename HDS::Vertex( iso_rect[2] ) ),
      hds.vertices_push_back( typename HDS::Vertex( iso_rect[3] ) )
    };

    typename HDS::Halfedge_handle new_hedges[4] = {
      hds.edges_push_back(typename HDS::Halfedge(), typename HDS::Halfedge()),
      hds.edges_push_back(typename HDS::Halfedge(), typename HDS::Halfedge()),
      hds.edges_push_back(typename HDS::Halfedge(), typename HDS::Halfedge()),
      hds.edges_push_back(typename HDS::Halfedge(), typename HDS::Halfedge())
    };

    for (int k=0;k<4; ++k)
    {
      item_decorator.set_face(new_hedges[k], new_face);
      item_decorator.set_vertex(new_hedges[k], vertices[ (k+1)%4 ]);
      item_decorator.set_vertex_halfedge(vertices[ (k+1)%4 ], new_hedges[k]);
      item_decorator.set_vertex(new_hedges[k]->opposite(), vertices[ k ]);
      new_hedges[k]->set_next( new_hedges[ (k+1)%4 ] );
      item_decorator.set_prev( new_hedges[ (k+1)%4 ], new_hedges[k] );
    }

    item_decorator.set_face_halfedge(new_face, new_hedges[0]);
  }

}

/// \todo in case when no iso_rectangle is provided,  consider a sufficently large iso_rectangle to contain the result and mark
///  halfedge as fake and vertices as infinite.

template< class DT2, class HDS>
struct Default_create_face{
  typename HDS::Face
  operator()( typename DT2::Vertex_handle ) const
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
};


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
  dt2.insert(
    boost::make_zip_iterator( boost::make_tuple(point_begin, info_begin) ),
    boost::make_zip_iterator( boost::make_tuple(point_end, info_end) ) );

    create_hds_for_cropped_voronoi_diagram<Input_kernel, Exact_kernel>(dt2,
                                                                     iso_rect,
                                                                     hds,
                                                                     Create_face_from_info<DT2, HDS>() );
}

/*!
Version taking a range of points and creating the HDS cropped to the iso-rectangle
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

#include <CGAL/Union_find.h>

template <class HDS>
void join_faces_with_same_color(HDS& hds)
{
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

        hds.vertices_erase(target->vertex()); //remove target vertex not remove in the loop below
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
    {
      const typename HDS::Vertex::Point& p=hit->opposite()->vertex()->point();
      const typename HDS::Vertex::Point& q=hit->vertex()->point();
      typename HDS::Halfedge_handle tmp=hit;
      CGAL::Orientation orient =
        orientation( p, q, tmp->next()->vertex()->point() );
      while (orient==CGAL::COLLINEAR)
      {
        tmp=tmp->next();
        orient = orientation( p, q, tmp->next()->vertex()->point() );
      }

      if( orient==CGAL::COUNTERCLOCKWISE )
        item_decorator.set_face_halfedge(face, hit); //this is the outer ccb
      else
        face->holes.push_back(hit); //this is a hole
    }
    else
      face->holes.push_back(hit); //this is a hole

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
      ++fit;
  }
}
