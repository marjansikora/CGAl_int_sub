/////////////////////////////////////////////////////////////////////////////
//
// test.cpp - Program for testing of triangle2dCGAL class
//
// Author: Marjan Sikora (sikora@fesb.hr)
//
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <sstream>
#include <chrono>
#include <random>
#include <fstream>
#include <cmath>

#include "test.h"
#include "triangle2dCGAL.h"

using namespace std;
using namespace std::chrono;

// global variables with generated triangles and intersection and subtraction results
Point_2_CGAL a[NOITER][3], b[NOITER][3];
Point_2_CGAL a0, a1, a2, b0, b1, b2;
vector<CGALTriangle2d*> colIntCGTri;
vector<CGALTriangle2d*> colSubCGTri;

// used for random number generation
unsigned seed = (unsigned)chrono::system_clock::now().time_since_epoch().count();
minstd_rand0 generator(seed);
unsigned mx = generator.max();
unsigned mn = generator.min();

// prototype
void GenerateTriangles();

int main() {
// for testing CGAL triangle class intersection operations
// vector with pairs of triangles A,B
	vector<Triangle_2_CGAL_EE> t2ee;
// for testing CGAL polygon class subtraction operations 
	vector<Polygon_2_CGAL> p2ee;

// for measuring area checking for errors
	vector<double> t2ee_area_A;
	double t2ee_area_A_sum;
	vector<double> t2ee_area_int;
	double t2ee_area_int_sum = 0.0;
	vector<double> p2ee_area_subt; 
	double p2ee_area_subt_sum = 0.0;
	vector<double> t2ee_area_delta;
	double t2ee_area_delta_sum;

// for testing my intersection and subtraction algorithms
	vector<CGALTriangle2d> t2my;

	vector<double> t2my_A_area;
	double t2my_A_area_sum;
	vector<double> t2my_int_area;
	double t2my_int_area_sum;
	vector<double> t2my_sub_area;
	double t2my_sub_area_sum;
	vector<double> t2my_delta_area;
	double t2my_delta_area_sum;

// for triangulating intersection results of CGAL function
	Delaunay_triangulation_2_CGAL_EE dtee;

// for measuring time
	high_resolution_clock::time_point t1;
	duration<double> time_span;

// generating triangles and filling CGAL 2d triangle array (for intersection base results)
	for (int i = 0; i < NOITER; i++) {
		GenerateTriangles();
		a[i][0] = a0; a[i][1] = a1; a[i][2] = a2;
		b[i][0] = b0; b[i][1] = b1; b[i][2] = b2;
	}

// filling CGAL 2d triangle array (for intersection base results)
	t1 = high_resolution_clock::now();
	for (int i = 0; i < NOITER; i++) {
		t2ee.push_back(Triangle_2_CGAL_EE(a[i][0], a[i][1], a[i][2]));
		t2ee_area_A.push_back( CGAL::to_double( t2ee.back().area() ) );
		t2ee.push_back(Triangle_2_CGAL_EE(b[i][0], b[i][1], b[i][2]));
	}
	time_span = duration_cast<duration<double>>(high_resolution_clock::now() - t1);
	cout << "CGAL triangles generation: " << time_span.count() << " s" << endl;

// filling CGAL 2d polygon array (for subtraction base results)
	t1 = high_resolution_clock::now();
	for (int i = 0; i < NOITER; i++) {
		Polygon_2_CGAL p1, p2;
		p1.push_back(a[i][0]);
		p1.push_back(a[i][1]);
		p1.push_back(a[i][2]);
		p2.push_back(b[i][0]);
		p2.push_back(b[i][1]);
		p2.push_back(b[i][2]);
		p2ee.push_back(p1);
		p2ee.push_back(p2);
	}
	time_span = duration_cast<duration<double>>(high_resolution_clock::now() - t1);
	cout << "CGAL polygones generation: " << time_span.count() << " s" << endl;

// filling my triangle array
	t1 = high_resolution_clock::now();
	for (int i = 0; i < NOITER; i++) {
		CGALTriangle2d a(a[i][0], a[i][1], a[i][2]);
		t2my.push_back(a);
		t2my_A_area.push_back( CGAL::to_double( a.GetTriangle().area() ));
		CGALTriangle2d b(b[i][0], b[i][1], b[i][2]);
		t2my.push_back(b);
	}
	time_span = duration_cast<duration<double>>(high_resolution_clock::now() - t1);
	cout << "triangle2dCGAL generation: " << time_span.count() << " s" << endl << endl;

// intersecting with CGAL triangle algorithm
	t1 = high_resolution_clock::now();
	for (int i = 0; i < NOITER*2 - 1; i=i+2) {
		auto result = CGAL::intersection(t2ee[i], t2ee[i + 1]);
		double temp = 0.0;
		if (result) {
			if (const Triangle_2_CGAL_EE *t = boost::get<Triangle_2_CGAL_EE>(&*result)) {
				temp = -(CGAL::to_double( t->area() ) );
			} else if (const Point_2_CGAL_EE *t = boost::get<Point_2_CGAL_EE>(&*result)) {
			} else if (const Segment_2_CGAL_EE *t = boost::get<Segment_2_CGAL_EE>(&*result)) {
			} else if (const vector<Point_2_CGAL_EE> *t = boost::get<vector<Point_2_CGAL_EE>>(&*result)) {
				for(auto i : *t) {
					dtee.insert(i);
				}
				for (Delaunay_triangulation_2_CGAL_EE::Finite_faces_iterator fit = dtee.finite_faces_begin();
					fit != dtee.finite_faces_end();
					++fit) {
					Delaunay_triangulation_2_CGAL_EE::Face_handle face = fit;
					temp = temp + CGAL::to_double((dtee.triangle(face)).area());
				}
				dtee.clear();
			}
		} else {
			/* empty intersection */
		}
		t2ee_area_int.push_back( temp );
	}
	time_span = duration_cast<duration<double>>(high_resolution_clock::now() - t1);
	cout << "CGAL triangle intersection: " << time_span.count() << " s" << endl << endl;

// intersecting with my intersection algorithm
	vector<CGALTriangle2d*> colCGALTri;
	t1 = high_resolution_clock::now();
	for (int i = 0; i < NOITER * 2 - 1; i = i + 2) {
		double temp = 0.0;
		t2my[i].IntersectAndTriangulateSH(t2my[i + 1], colCGALTri);
		for (auto t : colCGALTri) {
			temp = temp + CGAL::to_double(t->GetTriangle().area());
		}
		for (auto it : colCGALTri) {
			delete it;
		}
		colCGALTri.clear();
		t2my_int_area.push_back(temp);
	}
	time_span = duration_cast<duration<double>>(high_resolution_clock::now() - t1);
	cout << "Triangle2dCGAL intersection: " << time_span.count() << " s" << endl;
	cout << "Errors:                      " << CGALTriangle2d::lNoExceptionsIntersect << endl << endl;

// subracting with CGAL polygon algorithm
	t1 = high_resolution_clock::now();
	for (int i = 0; i < NOITER * 2 - 1; i = i + 2) {
		Pwh_list_2_CGAL diffL;
		Pwh_list_2_CGAL::const_iterator  it;
		double temp = 0.0;
		CGAL::difference(p2ee[i], p2ee[i + 1], std::back_inserter(diffL));
		for (it = diffL.begin(); it != diffL.end(); ++it) {
			auto po = it->outer_boundary();
			temp = temp + CGAL::to_double(po.area());
		}
		p2ee_area_subt.push_back(temp);
	}
	time_span = duration_cast<duration<double>>(high_resolution_clock::now() - t1);
	cout << "CGAL polygones subtraction: " << time_span.count() << " s" << endl << endl;

// subtracting with my subtraction algorithm
	colCGALTri.clear();
	t1 = high_resolution_clock::now();
	for (int i = 0; i < NOITER*2 - 1; i=i+2) {
		double temp = 0.0;
		t2my[i].SubtractAndTriangulateWA(t2my[i + 1], colCGALTri);
		for (auto t : colCGALTri) {
			temp = temp + CGAL::to_double(t->GetTriangle().area());
		}
		for (auto it : colCGALTri) {
			delete it;
		}
		colCGALTri.clear();
		t2my_sub_area.push_back(temp);
	}
	time_span = duration_cast<duration<double>>(high_resolution_clock::now() - t1);
	cout << "Triangle2dCGAL subtraction: " << time_span.count() << " s" << endl;
	cout << "Errors:                     " << CGALTriangle2d::lNoExceptionsSubtract << endl << endl;

// checking my algorithms for error
// area of triangle A has to be equal to sum of areas of intersection and subtraction
	t2my_delta_area_sum = t2my_A_area_sum = t2my_int_area_sum = t2my_sub_area_sum = 0.0;
	t2ee_area_delta_sum = t2ee_area_A_sum = t2ee_area_int_sum = p2ee_area_subt_sum = 0.0;
	int errorMy = 0;
	for (unsigned int i = 0; i < t2my_A_area.size(); ++i){
		double delta = abs( t2my_A_area[i] - t2my_int_area[i] - t2my_sub_area[i] );
		if (delta > 1e-10) {
			cout << "triangle2d error #" << ++errorMy << ":" << t2my_A_area[i] << " " << t2my_int_area[i]
				<< " " << t2my_sub_area[i] << " " << delta << endl;
		}
		t2my_delta_area_sum = t2my_delta_area_sum + delta;
		t2my_A_area_sum = t2my_A_area_sum + t2my_A_area[i];
		t2my_int_area_sum = t2my_int_area_sum + t2my_int_area[i];
		t2my_sub_area_sum = t2my_sub_area_sum + t2my_sub_area[i];
	}

// checking CGAL algorithms for error
// area of triangle A has to be equal to sum of areas of intersection and subtraction
	int errorCGAL = 0;
	for (unsigned int i = 0; i < t2ee_area_A.size(); ++i) {
		double delta = abs( t2ee_area_A[i] - t2ee_area_int[i] - p2ee_area_subt[i] );
		if (delta > 1e-10) {
			cout << "cgal error #" << ++errorCGAL << ":" << t2ee_area_A[i] << " " << t2ee_area_int[i]
				<< " " << p2ee_area_subt[i] << " " << delta << endl;
		}
		t2ee_area_delta_sum += delta;
		t2ee_area_A_sum += t2ee_area_A[i];
		t2ee_area_int_sum += t2ee_area_int[i];
		p2ee_area_subt_sum += p2ee_area_subt[i];
	}

	cout << endl << "Number of pairs of triangles:       " << NOITER << endl << endl;
	cout << "triangle2dCGAL A triangle area sum: " << t2my_A_area_sum << endl;
	cout << "triangle2dCGAL subt + int area sum: " << t2my_int_area_sum + t2my_sub_area_sum << endl;
	cout << "triangle2dCGAL delta area sum     : " << t2my_delta_area_sum << endl;
	cout << "triangle2dCGAL #errors            : " << errorMy << endl << endl;

	cout << "CGAL A triangle area sum: " << t2ee_area_A_sum << endl;
	cout << "CGAL subt + int area sum: " << t2ee_area_int_sum + p2ee_area_subt_sum << endl;
	cout << "CGAL delta area sum     : " << t2ee_area_delta_sum << endl << endl;
	cout << "CGAL #errors            : " << errorCGAL << endl << endl;


	return 0;
}

// function generates two triangles that have the same vertex 
// pre:		-
// post:	- triangle vertices in a0, a1, a2, b0, b1, b2
void GenerateTrianglesColVerts() {
	high_resolution_clock::time_point time1;
	duration<double> time_span;
	double x, y;
	unsigned val;

// ColVerts
	Line_2_CGAL l;
	Triangle_2_CGAL t;
	do {
// generating triangle A
		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		a0 = Point_2_CGAL(x, y);
		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		a1 = Point_2_CGAL(x, y);
		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		a2 = Point_2_CGAL(x, y);
		t = Triangle_2_CGAL(a0, a1, a2);
		if (t.orientation() == CGAL::CLOCKWISE) {
			Point_2_CGAL temp;
			temp = a2;
			a2 = a1;
			a1 = temp;
		}

// generating triangle B - first common vertex
		b0 = a0;

// second vertex
		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		b1 = Point_2_CGAL(x, y);
		if ((generator() % 3 == 0)) {
			if (CGAL::squared_distance(b1, Line_2_CGAL(a1, a2)) < 0.01) {
				b1 = (Line_2_CGAL(a1, a2)).projection(b1);
			}
			else if (CGAL::squared_distance(b1, Line_2_CGAL(a2, a0)) < 0.01) {
				b1 = (Line_2_CGAL(a2, a0)).projection(b1);
			}
			if (CGAL::squared_distance(b1, a1) < 0.01) {
				b1 = a1;
			}
			if (CGAL::squared_distance(b1, a2) < 0.01) {
				b1 = a2;
			}
		}
// third vertex
		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		b2 = Point_2_CGAL(x, y);
		if ((generator() % 3 == 0)) {
			if (CGAL::squared_distance(b2, Line_2_CGAL(a1, a2)) < 0.01) {
				b2 = (Line_2_CGAL(a1, a2)).projection(b2);
			}
			else if (CGAL::squared_distance(b2, Line_2_CGAL(a2, a0)) < 0.01) {
				b2 = (Line_2_CGAL(a2, a0)).projection(b2);
			}
			if (CGAL::squared_distance(b2, a1) < 0.01) {
				b2 = a1;
			}
			if (CGAL::squared_distance(b2, a2) < 0.01) {
				b2 = a2;
			}
		}
		t = Triangle_2_CGAL(b0, b1, b2);
	} while (t.is_degenerate());

	if (t.orientation() == CGAL::CLOCKWISE) {
		Point_2_CGAL temp;
		temp = b2;
		b2 = b1;
		b1 = temp;
	}

// randomizing the order of vertices...
	val = generator() % 3;
	for (unsigned int i = 0; i < val; ++i) {
		Point_2_CGAL temp;
		temp = a0;
		a0 = a1;
		a1 = a2;
		a2 = temp;
	}
	val = generator() % 3;
	for (unsigned int i = 0; i < val; ++i) {
		Point_2_CGAL temp;
		temp = b0;
		b0 = b1;
		b1 = b2;
		b2 = temp;
	}

// make sure the order is positive
	t = Triangle_2_CGAL(a0, a1, a2);
	if (t.orientation() == CGAL::CLOCKWISE) {
		Point_2_CGAL temp;
		temp = a2;
		a2 = a1;
		a1 = temp;
	}
	t = Triangle_2_CGAL(b0, b1, b2);
	if (t.orientation() == CGAL::CLOCKWISE) {
		Point_2_CGAL temp;
		temp = b2;
		b2 = b1;
		b1 = temp;
	}
}

// function generates two triangles A and B, 
//		such that one vertex of B lies on the edge of A
// pre:		-
// post:	- triangle vertices in a0, a1, a2, b0, b1, b2
void GenerateTrianglesColVertBEdgeA() {
	high_resolution_clock::time_point time1;
	duration<double> time_span;
	double x, y;
	double coeff;
	unsigned val;

	Line_2_CGAL l;
	Triangle_2_CGAL t;
	do {
		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		a0 = Point_2_CGAL(x, y);
		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		a1 = Point_2_CGAL(x, y);
		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		a2 = Point_2_CGAL(x, y);
		t = Triangle_2_CGAL(a0, a1, a2);
		if (t.orientation() == CGAL::CLOCKWISE) {
			Point_2_CGAL temp;
			temp = a2;
			a2 = a1;
			a1 = temp;
		}

// triangle B - the vertex that lies on edge of A
		val = generator();
		Vector_2_CGAL vec;
		vec = a1 - a0;
		coeff = double(val) / (mx - mn);
		vec = vec * coeff;
		b0 = a0 + vec;
		b0 = (Line_2_CGAL(a0, a1)).projection(b0);

		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		b1 = Point_2_CGAL(x, y);
		if ((generator() % 3 == 0)) {
			if (CGAL::squared_distance(b1, Line_2_CGAL(a1, a2)) < 0.01) {
				b1 = (Line_2_CGAL(a1, a2)).projection(b1);
			}
			else if (CGAL::squared_distance(b1, Line_2_CGAL(a2, a0)) < 0.01) {
				b1 = (Line_2_CGAL(a2, a0)).projection(b1);
			}
			if (CGAL::squared_distance(b1, a1) < 0.01) {
				b1 = a1;
			}
			if (CGAL::squared_distance(b1, a2) < 0.01) {
				b1 = a2;
			}
		}

		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		b2 = Point_2_CGAL(x, y);
		if ((generator() % 3 == 0)) {
			if (CGAL::squared_distance(b2, Line_2_CGAL(a1, a2)) < 0.01) {
				b2 = (Line_2_CGAL(a1, a2)).projection(b2);
			}
			else if (CGAL::squared_distance(b2, Line_2_CGAL(a2, a0)) < 0.01) {
				b2 = (Line_2_CGAL(a2, a0)).projection(b2);
			}
			if (CGAL::squared_distance(b2, a1) < 0.01) {
				b2 = a1;
			}
			if (CGAL::squared_distance(b2, a2) < 0.01) {
				b2 = a2;
			}
		}
		t = Triangle_2_CGAL(b0, b1, b2);
	} while (t.is_degenerate());

	if (t.orientation() == CGAL::CLOCKWISE) {
		Point_2_CGAL temp;
		temp = b2;
		b2 = b1;
		b1 = temp;
	}

	val = generator() % 3;
	for (unsigned int i = 0; i < val; ++i) {
		Point_2_CGAL temp;
		temp = a0;
		a0 = a1;
		a1 = a2;
		a2 = temp;
	}
	val = generator() % 3;
	for (unsigned int i = 0; i < val; ++i) {
		Point_2_CGAL temp;
		temp = b0;
		b0 = b1;
		b1 = b2;
		b2 = temp;
	}

	t = Triangle_2_CGAL(a0, a1, a2);
	if (t.orientation() == CGAL::CLOCKWISE) {
		Point_2_CGAL temp;
		temp = a2;
		a2 = a1;
		a1 = temp;
	}
	t = Triangle_2_CGAL(b0, b1, b2);
	if (t.orientation() == CGAL::CLOCKWISE) {
		Point_2_CGAL temp;
		temp = b2;
		b2 = b1;
		b1 = temp;
	}

}

// function generates two triangles A and B, 
//		such that one vertex of A lies on the edge of B
// pre:		-
// post:	- triangle vertices in a0, a1, a2, b0, b1, b2
void GenerateTrianglesColVertAEdgeB() {
	high_resolution_clock::time_point time1;
	duration<double> time_span;
	double x, y;
	double coeff;
	unsigned val;

	Line_2_CGAL l;
	Triangle_2_CGAL t;
	do {
		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		a0 = Point_2_CGAL(x, y);
		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		a1 = Point_2_CGAL(x, y);
		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		a2 = Point_2_CGAL(x, y);
		t = Triangle_2_CGAL(a0, a1, a2);
		if (t.orientation() == CGAL::CLOCKWISE) {
			Point_2_CGAL temp;
			temp = a2;
			a2 = a1;
			a1 = temp;
		}

		val = generator();
		Vector_2_CGAL vec;
		vec = a1 - a0;
		coeff = double(val) / (mx - mn);
		vec = vec * coeff;
		b0 = a0 + vec;
		b0 = (Line_2_CGAL(a0, a1)).projection(b0);

		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		b1 = Point_2_CGAL(x, y);
		if ((generator() % 3 == 0)) {
			if (CGAL::squared_distance(b1, Line_2_CGAL(a1, a2)) < 0.01) {
				b1 = (Line_2_CGAL(a1, a2)).projection(b1);
			}
			else if (CGAL::squared_distance(b1, Line_2_CGAL(a2, a0)) < 0.01) {
				b1 = (Line_2_CGAL(a2, a0)).projection(b1);
			}
			if (CGAL::squared_distance(b1, a1) < 0.01) {
				b1 = a1;
			}
			if (CGAL::squared_distance(b1, a2) < 0.01) {
				b1 = a2;
			}
		}
		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		b2 = Point_2_CGAL(x, y);
		if ((generator() % 3 == 0)) {
			if (CGAL::squared_distance(b2, Line_2_CGAL(a1, a2)) < 0.01) {
				b2 = (Line_2_CGAL(a1, a2)).projection(b2);
			}
			else if (CGAL::squared_distance(b2, Line_2_CGAL(a2, a0)) < 0.01) {
				b2 = (Line_2_CGAL(a2, a0)).projection(b2);
			}
			if (CGAL::squared_distance(b2, a1) < 0.01) {
				b2 = a1;
			}
			if (CGAL::squared_distance(b2, a2) < 0.01) {
				b2 = a2;
			}
		}
		t = Triangle_2_CGAL(b0, b1, b2);
	} while (t.is_degenerate());

	if (t.orientation() == CGAL::CLOCKWISE) {
		Point_2_CGAL temp;
		temp = b2;
		b2 = b1;
		b1 = temp;
	}

	val = generator() % 3;
	for (unsigned int i = 0; i < val; ++i) {
		Point_2_CGAL temp;
		temp = a0;
		a0 = a1;
		a1 = a2;
		a2 = temp;
	}
	val = generator() % 3;
	for (unsigned int i = 0; i < val; ++i) {
		Point_2_CGAL temp;
		temp = b0;
		b0 = b1;
		b1 = b2;
		b2 = temp;
	}
	t = Triangle_2_CGAL(a0, a1, a2);
	if (t.orientation() == CGAL::CLOCKWISE) {
		Point_2_CGAL temp;
		temp = a2;
		a2 = a1;
		a1 = temp;
	}
	t = Triangle_2_CGAL(b0, b1, b2);
	if (t.orientation() == CGAL::CLOCKWISE) {
		Point_2_CGAL temp;
		temp = b2;
		b2 = b1;
		b1 = temp;
	}
	Point_2_CGAL temp;
	temp = a0;
	a0 = b0;
	b0 = temp;
	temp = a1;
	a1 = b1;
	b1 = temp;
	temp = a2;
	a2 = b2;
	b2 = temp;
}

// function generates two triangles A and B, 
//		that have one colinear edge
// pre:		-
// post:	- triangle vertices in a0, a1, a2, b0, b1, b2
void GenerateTrianglesColEdges() {
	high_resolution_clock::time_point time1;
	duration<double> time_span;
	double x, y;
	unsigned val;
	
	Line_2_CGAL l;
	Triangle_2_CGAL t;
	do {
		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		a0 = Point_2_CGAL(x, y);
		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		a1 = Point_2_CGAL(x, y);
		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		a2 = Point_2_CGAL(x, y);
		t = Triangle_2_CGAL(a0, a1, a2);
		if (t.orientation() == CGAL::CLOCKWISE) {
			Point_2_CGAL temp;
			temp = a2;
			a2 = a1;
			a1 = temp;
		}

// triangle B, colinear edge
		l = Line_2_CGAL(a1, a2);
		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		b1 = Point_2_CGAL(x, y);
		b1 = l.projection(b1);

// snap
		if ((generator() % 3 == 0)) {
			if (CGAL::squared_distance(b1, a1) < 0.01) {
				b1 = a1;
			}
		}
		if ((generator() % 3 == 0)) {
			if (CGAL::squared_distance(b1, a2) < 0.01) {
				b1 = a2;
			}
		}

		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		b2 = Point_2_CGAL(x, y);
		b2 = l.projection(b2);

// snap
		if ((generator() % 3 == 0)) {
			if (CGAL::squared_distance(b2, a1) < 0.01) {
				b2 = a1;
			}
		}
		if ((generator() % 3 == 0)) {
			if (CGAL::squared_distance(b2, a2) < 0.01) {
				b2 = a2;
			}
		}

// triangle B, third vertex
		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		b0 = Point_2_CGAL(x, y);
		if ((generator() % 3 == 0)) {
			if (CGAL::squared_distance(b0, a0) < 0.01) {
				b0 = a0;
			}
			if (CGAL::squared_distance(b0, Line_2_CGAL(a0, a1))< 0.01) {
				b0 = (Line_2_CGAL(a0, a1)).projection(b0);
			}
			else if (CGAL::squared_distance(b0, Line_2_CGAL(a0, a2))<0.01) {
				b0 = (Line_2_CGAL(a0, a2)).projection(b0);
			}
		}
		t = Triangle_2_CGAL(b0, b1, b2);


	} while (t.is_degenerate());

	if (t.orientation() == CGAL::CLOCKWISE) {
		Point_2_CGAL temp;
		temp = b2;
		b2 = b1;
		b1 = temp;
	}

	val = generator() % 3;
	for (unsigned int i = 0; i < val; ++i) {
		Point_2_CGAL temp;
		temp = a0;
		a0 = a1;
		a1 = a2;
		a2 = temp;
	}
	val = generator() % 3;
	for (unsigned int i = 0; i < val; ++i) {
		Point_2_CGAL temp;
		temp = b0;
		b0 = b1;
		b1 = b2;
		b2 = temp;
	}
	t = Triangle_2_CGAL(a0, a1, a2);
	if (t.orientation() == CGAL::CLOCKWISE) {
		Point_2_CGAL temp;
		temp = a2;
		a2 = a1;
		a1 = temp;
	}	
	t = Triangle_2_CGAL(b0, b1, b2);
	if (t.orientation() == CGAL::CLOCKWISE) {
		Point_2_CGAL temp;
		temp = b2;
		b2 = b1;
		b1 = temp;
	}
}

// function generates two triangles A and B, trivial case
// pre:		-
// post:	- triangle vertices in a0, a1, a2, b0, b1, b2
void GenerateTrianglesTrivial() {
	high_resolution_clock::time_point time1;
	duration<double> time_span;
	double x, y;
	unsigned val;

	Triangle_2_CGAL ta;
	do {
		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		a0 = Point_2_CGAL(x, y);
		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		a1 = Point_2_CGAL(x, y);
		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		a2 = Point_2_CGAL(x, y);
		ta = Triangle_2_CGAL(a0, a1, a2);
		if (ta.orientation() == CGAL::CLOCKWISE) {
			Point_2_CGAL temp;
			temp = a2;
			a2 = a1;
			a1 = temp;
		}
		ta = Triangle_2_CGAL(a0, a1, a2);
	} while (ta.is_degenerate());

	Triangle_2_CGAL tb;
	do {
		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		b0 = Point_2_CGAL(x, y);
		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		b1 = Point_2_CGAL(x, y);
		val = generator();
		x = val * 20.0 / (mx - mn) - 10.0;
		val = generator();
		y = val * 20.0 / (mx - mn) - 10.0;
		b2 = Point_2_CGAL(x, y);
		tb = Triangle_2_CGAL(b0, b1, b2);
		if (tb.orientation() == CGAL::CLOCKWISE) {
			Point_2_CGAL temp;
			temp = b2;
			b2 = b1;
			b1 = temp;
		}
		tb = Triangle_2_CGAL(b0, b1, b2);
	} while (tb.is_degenerate());

	val = generator() % 3;
	for (unsigned int i = 0; i < val; ++i) {
		Point_2_CGAL temp;
		temp = a0;
		a0 = a1;
		a1 = a2;
		a2 = temp;
	}
	val = generator() % 3;
	for (unsigned int i = 0; i < val; ++i) {
		Point_2_CGAL temp;
		temp = b0;
		b0 = b1;
		b1 = b2;
		b2 = temp;
	}
	ta = Triangle_2_CGAL(a0, a1, a2);
	if (ta.orientation() == CGAL::CLOCKWISE) {
		Point_2_CGAL temp;
		temp = a2;
		a2 = a1;
		a1 = temp;
	}
	tb = Triangle_2_CGAL(b0, b1, b2);
	if (tb.orientation() == CGAL::CLOCKWISE) {
		Point_2_CGAL temp;
		temp = b2;
		b2 = b1;
		b1 = temp;
	}
}

// function generates two triangles A and B, 
//		randomly using trivial or special cases
// pre:		-
// post:	- triangle vertices in a0, a1, a2, b0, b1, b2
void GenerateTriangles() {
	high_resolution_clock::time_point time1;
	duration<double> time_span;
	unsigned val;

	val = generator() % 5;
	switch(val) {
		case 0:
			GenerateTrianglesTrivial();			
			break;
		case 1:
			GenerateTrianglesColEdges();		
			break;
		case 2:
			GenerateTrianglesColVerts();		
			break;
		case 3:
			GenerateTrianglesColVertBEdgeA();	
			break;
		case 4:
			GenerateTrianglesColVertAEdgeB();	
			break;
	}
}

