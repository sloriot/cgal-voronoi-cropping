// ------------------------------------------------------------------------------
// Copyright (c) 2013 GeometryFactory (FRANCE)
// Distributed under the Boost Software License, Version 1.0. (See accompany-
// ing file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
// ------------------------------------------------------------------------------

#ifndef SWIG_CGAL_VORONOI_CROPPING_2_HALFEDGEDS_H
#define SWIG_CGAL_VORONOI_CROPPING_2_HALFEDGEDS_H

#include <SWIG_CGAL/Common/Macros.h>
#include <SWIG_CGAL/Common/Iterator.h>
#include <SWIG_CGAL/HalfedgeDS/HalfedgeDS.h>


template < class HDS_cpp>
class HalfedgeDS_vc_wrapper: public HalfedgeDS_wrapper<HDS_cpp>
{
  typedef HalfedgeDS_wrapper<HDS_cpp> Base;
  boost::shared_ptr<std::list<std::pair<int,int> > > m_data_ptr;
public:
  HalfedgeDS_vc_wrapper():Base(){};
  HalfedgeDS_vc_wrapper(int v, int h, int f):Base(v,h,f){}
  HalfedgeDS_vc_wrapper( const HalfedgeDS_vc_wrapper& hds2):Base(hds2){}
  #ifndef SWIG
  void set_info_shared_ptr(boost::shared_ptr<std::list< std::pair<int,int> > > ptr)
  {m_data_ptr=ptr;}
  #endif
};

#endif //SWIG_CGAL_VORONOI_CROPPING_2_HALFEDGEDS_H