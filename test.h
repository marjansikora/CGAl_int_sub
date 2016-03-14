/////////////////////////////////////////////////////////////////////////////
//
// test.h - Declarations for test.cpp
//
// Author: Marjan Sikora (sikora@fesb.hr)
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __KERNEL_CGAL_H_INCLUDED__
#define __KERNEL_CGAL_H_INCLUDED__

#define NOITER 10000

//#include <CGAL/Simple_cartesian.h>
//typedef CGAL::Simple_cartesian<double> Kernel;

//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

#include <CGAL/Polygon_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Boolean_set_operations_2.h>

typedef Kernel::Point_2 Point_2_CGAL;
typedef Kernel::Line_2 Line_2_CGAL;
typedef Kernel::Ray_2 Ray_2_CGAL;
typedef Kernel::Segment_2 Segment_2_CGAL;
typedef Kernel::Vector_2 Vector_2_CGAL;
typedef Kernel::Triangle_2 Triangle_2_CGAL;
typedef CGAL::Polygon_2<Kernel> Polygon_2_CGAL;
typedef CGAL::Polygon_with_holes_2<Kernel>                Polygon_with_holes_2_CGAL;
typedef std::list<Polygon_with_holes_2_CGAL>                   Pwh_list_2_CGAL;

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel_EE;
typedef Kernel_EE::Point_2 Point_2_CGAL_EE;
typedef Kernel_EE::Line_2 Line_2_CGAL_EE;
typedef Kernel_EE::Ray_2 Ray_2_CGAL_EE;
typedef Kernel_EE::Segment_2 Segment_2_CGAL_EE;
typedef Kernel_EE::Vector_2 Vector_2_CGAL_EE;
typedef Kernel_EE::Triangle_2 Triangle_2_CGAL_EE;
typedef CGAL::Polygon_2<Kernel_EE> Polygon_2_CGAL_EE;
typedef CGAL::Delaunay_triangulation_2<Kernel_EE> Delaunay_triangulation_2_CGAL_EE;

#endif

