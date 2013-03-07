// ------------------------------------------------------------------------------
// Copyright (c) 2013 GeometryFactory (FRANCE)
// Distributed under the Boost Software License, Version 1.0. (See accompany-
// ing file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
// ------------------------------------------------------------------------------

//#define DEBUG_JOIN_FACES
#include "../Voronoi_cropping_2/include/voronoi_cropping.h"
#include "../Voronoi_cropping_2/include/hds_type.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>
#include <fstream>
#include <sstream>

typedef CGAL::Simple_cartesian< CGAL::Lazy_exact_nt<CGAL::Gmpq> > Exact_kernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Creator_uniform_2<double,K::Point_2>  Creator;

typedef CGAL::HalfedgeDS_default<K,HDS_Item_extra> HDS;

void write_hds(const HDS& hds, const char* fname)
{
  std::ofstream output(fname);
  for( HDS::Face_const_iterator fit=hds.faces_begin();fit!=hds.faces_end();++fit)
  {
    std::stringstream sstream;
    HDS::Halfedge_const_handle curr=fit->halfedge(), end=curr;
    int i=0;
    CGAL_assertion( curr!= HDS::Halfedge_const_handle() );
    do{
      CGAL_assertion( curr!= HDS::Halfedge_const_handle() );
      CGAL_assertion( curr->face()!=HDS::Face_handle() );
      CGAL_assertion( curr->face()==fit );
      ++i;
      sstream << curr->vertex()->point() << " 0\n";
      curr=curr->next();
    }while(curr!=end);
    sstream << curr->vertex()->point() << " 0\n";
    output << i+1 << " " << sstream.str() << std::endl;
  }
  output.close();
}

int main(int argc, char** argv)
{
  int n=100;
  if (argc==2)
    n=atoi(argv[1]);

  HDS hds;

  std::vector< K::Point_2 > points;
  points.reserve(n);
  std::vector< int > colors(n,0);

  CGAL::Random rng(0);
  CGAL::Random_points_in_disc_2<K::Point_2,Creator> g( 1,rng);
  CGAL::cpp11::copy_n( g, n, std::back_inserter(points));

  std::ofstream output("points.xyz");
  for (int i=0; i <n; ++i)
    output << points[i] << " 0\n";
  output.close();

  K::Iso_rectangle_2 bbox(-1, -1, 1.5, 1.5);

  output.open("bbox.cgal");
  output << "5 ";
  for (int i=0;i<5;++i)
    output << bbox[i]<< " 0 ";
  output << "\n";
  output.close();

  //create_hds_for_cropped_voronoi_diagram<K, Exact_kernel>(points.begin(), points.end(), bbox, hds);
  create_hds_for_cropped_voronoi_diagram<K, Exact_kernel>(points.begin(), points.end(), colors.begin(), colors.end(), bbox, hds);
  write_hds(hds, "hds.cgal");


  CGAL::HalfedgeDS_decorator<HDS> d(hds);
  CGAL_assertion( d.is_valid(false,3) );

  std::cout << hds.size_of_vertices() << " " << hds.size_of_halfedges() << " " << hds.size_of_faces() << std::endl;
  join_faces_with_same_color(hds);
  std::cout << hds.size_of_vertices() << " " << hds.size_of_halfedges() << " " << hds.size_of_faces() << std::endl;

  CGAL_assertion( d.is_valid(false,3) );

  write_hds(hds, "joint.cgal");
}

