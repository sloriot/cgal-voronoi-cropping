// ------------------------------------------------------------------------------
// Copyright (c) 2013 GeometryFactory (FRANCE)
// Distributed under the Boost Software License, Version 1.0. (See accompany-
// ing file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
// ------------------------------------------------------------------------------

%module CGAL_Voronoi_cropping_2

%include "SWIG_CGAL/common.i"
Decl_void_type()

SWIG_CGAL_add_java_loadLibrary(CGAL_Voronoi_cropping_2)

%import  "SWIG_CGAL/Common/Macros.h"
%import  "SWIG_CGAL/Kernel/CGAL_Kernel.i"
%include  "SWIG_CGAL/Common/Iterator.h"

%{
#include <SWIG_CGAL/User_packages/Voronoi_cropping_2/typedefs.h>
#include <SWIG_CGAL/HalfedgeDS/all_includes.h>
#include <SWIG_CGAL/User_packages/Voronoi_cropping_2/HalfedgeDS.h>
#include <SWIG_CGAL/User_packages/Voronoi_cropping_2/Voronoi_cropping_2.h>
%}

//include files

%extend HDSFace_wrapper{
  int color(){
    return $self->get_data()->color();
  }

  void set_color(int c){
    $self->get_data()->color()=c;
  }

  int point_index(){
    return $self->get_data()->index();
  }

  bool has_holes(){
    return !$self->get_data()->holes.empty();
  }

  SWIG_CGAL_Iterator<std::vector<HDS_d::Halfedge_handle>::iterator,HDSHalfedge_wrapper<HDS_d> >
  holes()
  {
    return
      SWIG_CGAL_Iterator<std::vector<HDS_d::Halfedge_handle>::iterator,HDSHalfedge_wrapper<HDS_d> >
      ( $self->get_data()->holes.begin(), $self->get_data()->holes.end() );
  }

}

%pragma(java) jniclassimports=%{
  import java.util.Iterator;
  import java.util.Collection;
  import CGAL.Kernel.Point_2;
  import CGAL.Kernel.Iso_rectangle_2;
%}

//definitions
%include "SWIG_CGAL/HalfedgeDS/HalfedgeDS.h"
%include "SWIG_CGAL/User_packages/Voronoi_cropping_2/HalfedgeDS.h"
%include "SWIG_CGAL/HalfedgeDS/HalfedgeDS_handles.h"
%include "SWIG_CGAL/HalfedgeDS/HalfedgeDS_decorator.h"
%include "SWIG_CGAL/HalfedgeDS/General_modifier.h"

//template instantiations of handles
SWIG_CGAL_declare_identifier_of_template_class(HDS_Halfedge_handle,HDSHalfedge_wrapper<HDS_d>)
SWIG_CGAL_declare_identifier_of_template_class(HDS_Face_handle,HDSFace_wrapper<HDS_d>)
%typemap(javaimports)                       HDSVertex_wrapper %{import CGAL.Kernel.Point_2;%}
SWIG_CGAL_declare_identifier_of_template_class(HDS_Vertex_handle,HDSVertex_wrapper<HDS_d>)

//Iterators
SWIG_CGAL_set_as_java_iterator(SWIG_CGAL_Iterator,HDS_Vertex_handle,)
SWIG_CGAL_declare_identifier_of_template_class(HDS_Vertex_iterator,SWIG_CGAL_Iterator< HDS_d::Vertex_iterator, HDSVertex_wrapper<HDS_d> >)
SWIG_CGAL_set_as_java_iterator(SWIG_CGAL_Iterator,HDS_Halfedge_handle,)
SWIG_CGAL_declare_identifier_of_template_class(HDS_Halfedge_iterator,SWIG_CGAL_Iterator< HDS_d::Halfedge_iterator, HDSHalfedge_wrapper<HDS_d> >)
SWIG_CGAL_set_as_java_iterator(SWIG_CGAL_Iterator,HDS_Face_handle,)
SWIG_CGAL_declare_identifier_of_template_class(HDS_Face_iterator,SWIG_CGAL_Iterator< HDS_d::Face_iterator, HDSFace_wrapper<HDS_d> >)
SWIG_CGAL_set_as_java_iterator(SWIG_CGAL_Iterator,HDS_Halfedge_handle,)
SWIG_CGAL_declare_identifier_of_template_class(HDS_Face_hole_iterator,SWIG_CGAL_Iterator<std::vector<HDS_d::Halfedge_handle>::iterator,HDSHalfedge_wrapper<HDS_d> >)
//general modifier for convenience
%typemap(javaimports)       General_modifier<HDS_d> %{import CGAL.Kernel.Point_2;%}
SWIG_CGAL_declare_identifier_of_template_class(HalfedgeDS_modifier,General_modifier<HDS_d>)
// template instantiation of HDS class
%typemap(javaimports)                       HalfedgeDS_wrapper %{import CGAL.Kernel.Point_2;%}
SWIG_CGAL_declare_identifier_of_template_class(HalfedgeDS_vc_base,HalfedgeDS_wrapper<HDS_d>)
SWIG_CGAL_declare_identifier_of_template_class(HalfedgeDS,HalfedgeDS_vc_wrapper<HDS_d>)
// template instantiation of HDS decorator class
%typemap(javaimports)                       HalfedgeDS_decorator_wrapper %{import CGAL.Kernel.Point_2;%}
SWIG_CGAL_declare_identifier_of_template_class(HalfedgeDS_decorator,HalfedgeDS_decorator_wrapper<HDS_d>)

//typemap for point input iterator
SWIG_CGAL_set_wrapper_iterator_helper_input(Point_2)
%include "SWIG_CGAL/Java/typemaps.i"
SWIG_CGAL_array_of_int_to_vector_of_int_typemap_in

%pragma(java) moduleimports=%{import CGAL.Kernel.Iso_rectangle_2; import CGAL.Kernel.Point_2; import java.util.Collection; import java.util.Iterator;%}

void join_faces(HalfedgeDS_vc_wrapper<HDS_d>&);

%{
  void join_faces(HalfedgeDS_vc_wrapper<HDS_d>& hds)
  {
    join_faces_with_same_color(hds.get_data(), true);
  }
%}

%typemap(javaimports)  Voronoi_cropping_wrapper_2 %{ import CGAL.Kernel.Point_2; import CGAL.Kernel.Iso_rectangle_2; import java.util.Iterator; %}
%include "SWIG_CGAL/User_packages/Voronoi_cropping_2/Voronoi_cropping_2.h"
%template (Voronoi_cropping_2) Voronoi_cropping_wrapper_2<HalfedgeDS_vc_wrapper<HDS_d> >;
