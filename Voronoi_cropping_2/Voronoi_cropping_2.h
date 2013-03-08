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
  boost::shared_ptr< std::list<int> > m_colors_ptr;
  std::vector<typename Point_2::cpp_base> m_points;

  struct Get_pointer: std::unary_function<int, int*>
  {
      int* operator()(int& i) const { return &i; }
  };

public:

  Voronoi_cropping_wrapper_2(): m_colors_ptr( new std::list<int>() ){}

  void insert( const Point_2& pt, int color=-1) {
    m_points.push_back( pt.get_data() );
    m_colors_ptr->push_back(color);
  }
  void insert(Wrapper_iterator_helper<Point_2>::input points, boost::shared_ptr<std::vector<int> > colors)
  {
    std::size_t nb_elements=colors->size();
    std::size_t capacity=m_points.size()+nb_elements;
    m_points.reserve(capacity);

    std::copy( colors->begin(), colors->end(), std::back_inserter(*m_colors_ptr) );
    std::copy( SWIG_CGAL::get_begin(points), SWIG_CGAL::get_end(points), std::back_inserter(m_points) );
  }

  void voronoi_diagram(HDS_wrapper& hds, const Iso_rectangle_2& iso_rectangle)
  {
    hds.clear();
    hds.set_info_shared_ptr(m_colors_ptr);

    typedef CGAL::Simple_cartesian< CGAL::Lazy_exact_nt<CGAL::Gmpq> > Exact_kernel;
    create_hds_for_cropped_voronoi_diagram<EPIC_Kernel, Exact_kernel>(
      m_points.begin(),
      m_points.end(),
      boost::make_transform_iterator(m_colors_ptr->begin(), Get_pointer()),
      boost::make_transform_iterator(m_colors_ptr->end(), Get_pointer()),
      iso_rectangle.get_data(),
      hds.get_data()
    );
  }

  //~ Iso_rectangle_2 get_bonding_iso_rectangle()
  //~ {return Iso_rectangle_2();}

  //~ Iso_rectangle_2 voronoi_diagram(HDS_wrapper& hds)
  //~ {return Iso_rectangle_2();}

  void clear()
  {
    m_colors_ptr=boost::shared_ptr<std::list<int> >( new std::list<int>() );
    m_points.clear();
  }

};



#endif //SWIG_CGAL_VORONOI_CROPPING_2_VORONOI_CROPPING_2_H
