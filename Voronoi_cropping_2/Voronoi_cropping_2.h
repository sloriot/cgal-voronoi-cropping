// ------------------------------------------------------------------------------
// Copyright (c) 2013 GeometryFactory (FRANCE)
// Distributed under the Boost Software License, Version 1.0. (See accompany-
// ing file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
// ------------------------------------------------------------------------------

#ifndef SWIG_CGAL_VORONOI_CROPPING_2_VORONOI_CROPPING_2_H
#define SWIG_CGAL_VORONOI_CROPPING_2_VORONOI_CROPPING_2_H

#include <SWIG_CGAL/Common/Macros.h>
#include <SWIG_CGAL/Common/Iterator.h>
#include <SWIG_CGAL/Common/Wrapper_iterator_helper.h>
#include <SWIG_CGAL/Kernel/Point_2.h>
#include <SWIG_CGAL/Kernel/Iso_rectangle_2.h>
#include <voronoi_cropping.h>

#include <boost/shared_ptr.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>


template <class HDS_wrapper>
class Voronoi_cropping_wrapper_2{
  boost::shared_ptr< std::list< std::pair<int, int> > > m_indices_and_colors_ptr;
  std::vector<typename Point_2::cpp_base> m_points;

  struct Get_pointer: std::unary_function<std::pair<int, int>, std::pair<int, int>*>
  {
      std::pair<int, int>* operator()(std::pair<int, int>& p) const { return &p; }
  };

public:

  Voronoi_cropping_wrapper_2(): m_indices_and_colors_ptr( new std::list< std::pair<int, int> >() ){}

  void insert( const Point_2& pt, int color=-1) {
    m_indices_and_colors_ptr->push_back( std::make_pair(m_points.size(), color) );
    m_points.push_back( pt.get_data() );
  }
  void insert(Wrapper_iterator_helper<Point_2>::input points, boost::shared_ptr<std::vector<int> > colors)
  {
    std::size_t nb_elements=colors->size();
    std::size_t capacity=m_points.size()+nb_elements;
    m_points.reserve(capacity);

    int index = (int) m_points.size();
    for( std::vector<int>::iterator it=colors->begin(), it_end=colors->end(); it!=it_end; ++it)
      m_indices_and_colors_ptr->push_back( std::make_pair(index++, *it) );
    std::copy( SWIG_CGAL::get_begin(points), SWIG_CGAL::get_end(points), std::back_inserter(m_points) );
  }

  void voronoi_diagram(HDS_wrapper& hds, const Iso_rectangle_2& iso_rectangle)
  {
    hds.clear();
    hds.set_info_shared_ptr(m_indices_and_colors_ptr);

    typedef CGAL::Simple_cartesian< CGAL::Lazy_exact_nt<CGAL::Gmpq> > Exact_kernel;
    create_hds_for_cropped_voronoi_diagram<EPIC_Kernel, Exact_kernel>(
      m_points.begin(),
      m_points.end(),
      boost::make_transform_iterator(m_indices_and_colors_ptr->begin(), Get_pointer()),
      boost::make_transform_iterator(m_indices_and_colors_ptr->end(), Get_pointer()),
      iso_rectangle.get_data(),
      hds.get_data()
    );
  }

  /// \todo shall we keep an internal dt3 to avoid rebuilding it in case of a few additionnal points inserted?

  Iso_rectangle_2 voronoi_diagram(HDS_wrapper& hds)
  {
    hds.clear();
    hds.set_info_shared_ptr(m_indices_and_colors_ptr);

    typedef CGAL::Simple_cartesian< CGAL::Lazy_exact_nt<CGAL::Gmpq> > Exact_kernel;
    return Iso_rectangle_2(
      create_hds_for_cropped_voronoi_diagram<EPIC_Kernel, Exact_kernel>(
        m_points.begin(),
        m_points.end(),
        boost::make_transform_iterator(m_indices_and_colors_ptr->begin(), Get_pointer()),
        boost::make_transform_iterator(m_indices_and_colors_ptr->end(), Get_pointer()),
        hds.get_data() )
    );
  }

  void set_point(int i, const Point_2& p)
  {
    m_points[i]=p.get_data();
  }

  void set_point(int i, double x, double y)
  {
    m_points[i]=EPIC_Kernel::Point_2(x,y);
  }

  Point_2 get_point(int i)
  {
    return Point_2(m_points[i]);
  }

  void get_point(int i, Point_2& p)
  {
    p.get_data()=m_points[i];
  }

  void clear()
  {
    m_indices_and_colors_ptr=boost::shared_ptr<std::list< std::pair<int,int> > >( new std::list< std::pair<int,int> >() );
    m_points.clear();
  }

};



#endif //SWIG_CGAL_VORONOI_CROPPING_2_VORONOI_CROPPING_2_H
