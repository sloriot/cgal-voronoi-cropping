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

bool check(const HDS& hds, std::size_t v, std::size_t e, std::size_t f)
{
  bool res= hds.size_of_vertices()==v &&
            hds.size_of_halfedges()==e &&
            hds.size_of_faces()==f;
  if (!res)
    std::cerr << "Error " << hds.size_of_vertices() << " " << hds.size_of_halfedges() << " " << hds.size_of_faces() << std::endl;

  return res;
}

void test1()
{
  std::cout << "test1" <<std::endl;
  std::vector<K::Point_2> points;
  points.push_back( K::Point_2(0,0) );
  points.push_back( K::Point_2(0,1) );
  points.push_back( K::Point_2(1,1) );
  points.push_back( K::Point_2(1,0) );
  points.push_back( K::Point_2(0.5,0.5) );

  {
    std::cout << "  test1-1\n";
    HDS hds;
    K::Iso_rectangle_2 bbox(0, 0, 1, 1);
    create_hds_for_cropped_voronoi_diagram<K, Exact_kernel>(points.begin(), points.end(), bbox, hds);
    CGAL_assertion ( check(hds,8, 24, 5) );
    CGAL::HalfedgeDS_decorator<HDS> d(hds);
    CGAL_assertion( d.is_valid(false,3) );
  }

  {
    std::cout << "  test1-2\n";
    HDS hds;
    K::Iso_rectangle_2 bbox(0.5, 0, 1, 1);
    create_hds_for_cropped_voronoi_diagram<K, Exact_kernel>(points.begin(), points.end(), bbox, hds);
    CGAL_assertion ( check(hds, 5, 14, 3) );
    CGAL::HalfedgeDS_decorator<HDS> d(hds);
    CGAL_assertion( d.is_valid(false,3) );
  }

  {
    std::cout << "  test1-3\n";
    HDS hds;
    K::Iso_rectangle_2 bbox(0, 0, 0.5, 1);
    create_hds_for_cropped_voronoi_diagram<K, Exact_kernel>(points.begin(), points.end(), bbox, hds);
    CGAL_assertion ( check(hds, 5, 14, 3) );
    CGAL::HalfedgeDS_decorator<HDS> d(hds);
    CGAL_assertion( d.is_valid(false,3) );
  }
}

void test2()
{
  std::cout << "test2" <<std::endl;
  std::vector<K::Point_2> points;
  points.push_back( K::Point_2(1,0) );
  points.push_back( K::Point_2(1,1) );
  points.push_back( K::Point_2(1,2) );
  points.push_back( K::Point_2(0,1) );
  points.push_back( K::Point_2(2,1) );

  HDS hds;
  K::Iso_rectangle_2 bbox(0.5, 0.5, 1.5, 1.5);
  create_hds_for_cropped_voronoi_diagram<K, Exact_kernel>(points.begin(), points.end(), bbox, hds);
  CGAL_assertion ( check(hds, 4, 8, 1) );

  CGAL::HalfedgeDS_decorator<HDS> d(hds);
  CGAL_assertion( d.is_valid(false,3) );
}

void test3()
{
  std::cout << "test3" <<std::endl;
  std::vector<K::Point_2> points;
  points.push_back( K::Point_2(0,0) );
  points.push_back( K::Point_2(2,0) );

  {
    std::cout << "  test3-1\n";
    HDS hds;
    K::Iso_rectangle_2 bbox(0.5, 0.5, 1.5, 1.5);
    create_hds_for_cropped_voronoi_diagram<K, Exact_kernel>(points.begin(), points.end(), bbox, hds);
    CGAL_assertion ( check(hds, 6, 14, 2) );

    CGAL::HalfedgeDS_decorator<HDS> d(hds);
    CGAL_assertion( d.is_valid(false,3) );
  }

  {
    std::cout << "  test3-2\n";
    HDS hds;
    K::Iso_rectangle_2 bbox(0, 0, 0.5, 0.5);
    create_hds_for_cropped_voronoi_diagram<K, Exact_kernel>(points.begin(), points.end(), bbox, hds);
    CGAL_assertion ( check(hds, 4, 8, 1) );

    CGAL::HalfedgeDS_decorator<HDS> d(hds);
    CGAL_assertion( d.is_valid(false,3) );
  }

  {
    std::cout << "  test3-3\n";
    points.push_back( K::Point_2(1,0) );
    HDS hds;
    K::Iso_rectangle_2 bbox(0, 0, 2, 2);
    create_hds_for_cropped_voronoi_diagram<K, Exact_kernel>(points.begin(), points.end(), bbox, hds);
    CGAL_assertion ( check(hds, 8, 20, 3) );

    CGAL::HalfedgeDS_decorator<HDS> d(hds);
    CGAL_assertion( d.is_valid(false,3) );
  }

  {
    std::cout << "  test3-4\n";
    HDS hds;
    K::Iso_rectangle_2 bbox(-2, -2, -1, -1);
    create_hds_for_cropped_voronoi_diagram<K, Exact_kernel>(points.begin(), points.end(), bbox, hds);
    CGAL_assertion ( check(hds, 4, 8, 1) );

    CGAL::HalfedgeDS_decorator<HDS> d(hds);
    CGAL_assertion( d.is_valid(false,3) );
  }
}

void test4()
{
  std::cout << "test4" <<std::endl;
  std::vector<K::Point_2> points;
  points.push_back( K::Point_2(0,0) );

  {
    std::cout << "  test4-1\n";
    HDS hds;
    K::Iso_rectangle_2 bbox(0.5, 0.5, 1.5, 1.5);
    create_hds_for_cropped_voronoi_diagram<K, Exact_kernel>(points.begin(), points.end(), bbox, hds);
    CGAL_assertion ( check(hds, 4, 8, 1) );

    CGAL::HalfedgeDS_decorator<HDS> d(hds);
    CGAL_assertion( d.is_valid(false,3) );
  }

  {
    std::cout << "  test4-2\n";
    HDS hds;
    K::Iso_rectangle_2 bbox(0.5, 0.5, 1.5, 1.5);
    create_hds_for_cropped_voronoi_diagram<K, Exact_kernel>(points.begin(), points.begin(), bbox, hds);
    CGAL_assertion ( check(hds, 4, 8, 1) );

    CGAL::HalfedgeDS_decorator<HDS> d(hds);
    CGAL_assertion( d.is_valid(false,3) );
  }
}

void test5()
{
  std::cout << "test5" <<std::endl;
  std::vector<K::Point_2> points;
  points.push_back( K::Point_2(10,0) );
  points.push_back( K::Point_2(20,10) );
  points.push_back( K::Point_2(10,20) );
  points.push_back( K::Point_2(0,10) );
  points.push_back( K::Point_2(10,10) );

  int colors[]={1,2,3,4,5};

  HDS hds;
  K::Iso_rectangle_2 bbox(5, 1, 15, 2);
  create_hds_for_cropped_voronoi_diagram<K, Exact_kernel>(points.begin(), points.end(), colors, colors+5, bbox, hds);
  CGAL_assertion ( check(hds, 4, 8, 1) );

  CGAL_assertion( hds.faces_begin()->color()==1 );

  CGAL::HalfedgeDS_decorator<HDS> d(hds);
  CGAL_assertion( d.is_valid(false,3) );
}

void test6()
{
  std::cout << "test6" <<std::endl;
  std::vector<K::Point_2> points;
  points.push_back( K::Point_2(1,0) );
  points.push_back( K::Point_2(1,1) );
  points.push_back( K::Point_2(0,1) );
  points.push_back( K::Point_2(0,0) );

  int colors[]={1,2,3,4};

  HDS hds;
  K::Iso_rectangle_2 bbox(-1, -1, 2, 2);
  create_hds_for_cropped_voronoi_diagram<K, Exact_kernel>(points.begin(), points.end(), colors, colors+4, bbox, hds);
  CGAL_assertion ( check(hds, 9, 24, 4) );

  CGAL::HalfedgeDS_decorator<HDS> d(hds);
  CGAL_assertion( d.is_valid(false,3) );
}

void test6bis()
{
  std::cout << "test6bis" <<std::endl;
  std::vector<K::Point_2> points;
 
  points.push_back( K::Point_2(0.3780546928833294,10.077223947566674)); 
  points.push_back( K::Point_2(0.28861197378333725,10.256109385766658));
  points.push_back( K::Point_2(0.7113880262166623,10.243890614233342));
  points.push_back( K::Point_2(0.621945307116671,10.422776052433326));
  points.push_back( K::Point_2(1.0346946247734936,10.406211498509538));
  points.push_back( K::Point_2(0.9457859736717265,10.584028800713073));
  points.push_back( K::Point_2(0.9864771388046294,10.599081442384993));
  points.push_back( K::Point_2(1.3575868958369663,10.4863190833188));
  points.push_back( K::Point_2(1.3090797708297002,10.680347583347867));
  points.push_back( K::Point_2(1.6909202291703,10.569652416652133));
  points.push_back( K::Point_2(1.642413104163033,10.7636809166812));
  // comment the next line makes it work...
  points.push_back( K::Point_2(2.0242535625036333,10.652985749985467));

  int colors[]={2,1,2,1,2,1,1,2,1,2,1,2};

  std::ofstream output("input.xyz");
  for (std::size_t i=0; i <points.size(); ++i)
    output << points[i] << " 0\n";
  output.close();
  
  HDS hds;
  K::Iso_rectangle_2 bbox( K::Point_2(-17.49, -13.6), K::Point_2(30.39, 25.6) );
  create_hds_for_cropped_voronoi_diagram<K, Exact_kernel>(points.begin(), points.end(), colors, colors+12, bbox, hds);
  write_hds(hds, "before.cgal");
  join_faces_with_same_color(hds);
  write_hds(hds, "after.cgal");
}

int main(int argc, char** argv)
{
  test1();
  test2();
  test3();
  test4();
  test5();
  test6();
  test6bis();

  int n=100;
  if (argc==2)
    n=atoi(argv[1]);

  HDS hds;

  std::vector< K::Point_2 > points;
  points.reserve(n);
  std::vector< int > colors(n,0);

  CGAL::Random rng(0);
  CGAL::Random_points_in_disc_2<K::Point_2,Creator> g( 1,rng);
  CGAL::cpp11::copy_n( g, n, std::back_inserter(points) );

  std::ofstream output("points.xyz");
  for (int i=0; i <n; ++i)
    output << points[i] << " 0\n";
  output.close();

  K::Iso_rectangle_2 bbox(0.5, 0.5, 1.5, 1.5);


  output.open("bbox.cgal");
  output << "5 ";
  for (int i=0;i<5;++i)
    output << bbox[i]<< " 0 ";
  output << "\n";
  output.close();

  //create_hds_for_cropped_voronoi_diagram<K, Exact_kernel>(points.begin(), points.end(), bbox, hds);
  create_hds_for_cropped_voronoi_diagram<K, Exact_kernel>(points.begin(), points.end(), colors.begin(), colors.end(), bbox, hds);
  //create_hds_for_cropped_voronoi_diagram<K, Exact_kernel>(points.begin(), points.end(), colors.begin(), colors.end(), hds);
  write_hds(hds, "hds.cgal");


  CGAL::HalfedgeDS_decorator<HDS> d(hds);
  CGAL_assertion( d.is_valid(false,3) );

  std::cout << hds.size_of_vertices() << " " << hds.size_of_halfedges() << " " << hds.size_of_faces() << std::endl;
  join_faces_with_same_color(hds);
  std::cout << hds.size_of_vertices() << " " << hds.size_of_halfedges() << " " << hds.size_of_faces() << std::endl;

  CGAL_assertion( d.is_valid(false,3) );

  write_hds(hds, "joint.cgal");
}

