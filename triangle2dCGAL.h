/////////////////////////////////////////////////////////////////////////////
//
// Triangle2dCGAL.h - Declarations for CGALTriangle2d class with CGAL kernel
//
// Author: Marjan Sikora (sikora@fesb.hr)
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __TRIANGLE2CGAL_H_INCLUDED__
#define __TRIANGLE2CGAL_H_INCLUDED__

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <algorithm>
#include <deque>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::Point_2 Point_2_CGAL;
typedef Kernel::Line_2 Line_2_CGAL;
typedef Kernel::Ray_2 Ray_2_CGAL;
typedef Kernel::Segment_2 Segment_2_CGAL;
typedef Kernel::Vector_2 Vector_2_CGAL;
typedef Kernel::Triangle_2 Triangle_2_CGAL;

using namespace std;

/////////////////////////////////////////////////////////////////////////////
//
// class CTriangle2d	- a simple 2D triangle class 
//						  made for testing of intersection and 
//						  subtraction of two triangles
//                 
/////////////////////////////////////////////////////////////////////////////

class CGALTriangle2d {
// traingle geometry
	Point_2_CGAL _point[3];
	Segment_2_CGAL _segment[3];
	Triangle_2_CGAL _mCGALTr;

public:
// constructors used for initializing triangle geometry
	CGALTriangle2d() {};
	CGALTriangle2d(Point_2_CGAL p0, Point_2_CGAL p1, Point_2_CGAL p2) {
		_mCGALTr = Triangle_2_CGAL(p0, p1, p2);
		_point[0] = Point_2_CGAL(p0.x(), p0.y());
		_point[1] = Point_2_CGAL(p1.x(), p1.y());
		_point[2] = Point_2_CGAL(p2.x(), p2.y());
		_segment[0] = Segment_2_CGAL(Point_2_CGAL(p0.x(), p0.y()), Point_2_CGAL(p1.x(), p1.y()));
		_segment[1] = Segment_2_CGAL(Point_2_CGAL(p1.x(), p1.y()), Point_2_CGAL(p2.x(), p2.y()));
		_segment[2] = Segment_2_CGAL(Point_2_CGAL(p2.x(), p2.y()), Point_2_CGAL(p0.x(), p0.y()));
	}
	CGALTriangle2d(const CGALTriangle2d& temp) {
		_point[0] = Point_2_CGAL(temp._mCGALTr[0].x(), temp._mCGALTr[0].y());
		_point[1] = Point_2_CGAL(temp._mCGALTr[1].x(), temp._mCGALTr[1].y());
		_point[2] = Point_2_CGAL(temp._mCGALTr[2].x(), temp._mCGALTr[2].y());
		_mCGALTr = Triangle_2_CGAL(_point[0], _point[1], _point[2]);
		_segment[0] = Segment_2_CGAL(_point[0], _point[1]);
		_segment[1] = Segment_2_CGAL(_point[1], _point[2]);
		_segment[2] = Segment_2_CGAL(_point[2], _point[0]);
	}
	CGALTriangle2d(const CGALTriangle2d& temp, int i) {
		_point[0] = Point_2_CGAL(temp._mCGALTr[i].x(), temp._mCGALTr[i].y());
		_point[1] = Point_2_CGAL(temp._mCGALTr[(i + 1) % 3].x(), temp._mCGALTr[(i + 1) % 3].y());
		_point[2] = Point_2_CGAL(temp._mCGALTr[(i + 2) % 3].x(), temp._mCGALTr[(i + 2) % 3].y());
		_mCGALTr = Triangle_2_CGAL(_point[0], _point[1], _point[2]);
		_segment[0] = Segment_2_CGAL(_point[0], _point[1]);
		_segment[1] = Segment_2_CGAL(_point[1], _point[2]);
		_segment[2] = Segment_2_CGAL(_point[2], _point[0]);
	}
	~CGALTriangle2d() {};

// operator =
	CGALTriangle2d& operator=(const CGALTriangle2d& temp) {
		_point[0] = Point_2_CGAL(temp._mCGALTr[0].x(), temp._mCGALTr[0].y());
		_point[1] = Point_2_CGAL(temp._mCGALTr[1].x(), temp._mCGALTr[1].y());
		_point[2] = Point_2_CGAL(temp._mCGALTr[2].x(), temp._mCGALTr[2].y());
		_mCGALTr = Triangle_2_CGAL(_point[0], _point[1], _point[2]);
		_segment[0] = Segment_2_CGAL(_point[0], _point[1]);
		_segment[1] = Segment_2_CGAL(_point[1], _point[2]);
		_segment[2] = Segment_2_CGAL(_point[2], _point[0]);
		return *this;
	}

// get function for private variable with triangle geometry
	const Triangle_2_CGAL& GetTriangle() {
		return _mCGALTr;
	}


//// MAIN ALGORTIHM FUNCTIONS ////
// intersect triangles that have positive orientation using Southerland-Hodgson algorithm
//    and triangulate the result
	void IntersectAndTriangulateSH(const CGALTriangle2d& other, vector<CGALTriangle2d*>& rColTri);
// subtract triangle other from this (both with positive orientation) with Wilson algorithm
//    and triangulate the result
	void SubtractAndTriangulateWA(const CGALTriangle2d& other, vector<CGALTriangle2d*>& rColTri);


// exception counters
	static long lNoExceptionsIntersect;
	static long lNoExceptionsSubtract;

// operation counters
	static long lNoOperationsSubtraction;
	static long lNoOperationsIntersersection;

private:
// private methods used to triangulate results of intersection
//    and subtraction operations
	void TriangulateFourConvexPoints(const Point_2_CGAL& rA, const Point_2_CGAL& rB,
		const Point_2_CGAL& rC, const Point_2_CGAL& rD, vector<CGALTriangle2d*>& rColTri);
	void TriangulateMoreThanFourConvexPoints(const vector<Point_2_CGAL>& arrPoints, vector<CGALTriangle2d*>& rColTri);
	void TriangulateFiveConcavePoints(const vector<Point_2_CGAL>& arrPoints, vector<CGALTriangle2d*>& rColTri);
	void TriangulateSixConcavePoints(const vector<Point_2_CGAL>& arrPoints, vector<CGALTriangle2d*>& rColTri);
	void TriangulateSevenConcavePoints(const vector<Point_2_CGAL>& arrPoints, vector<CGALTriangle2d*>& rColTri);
	void TriangulateHole(Point_2_CGAL arrThis[], Point_2_CGAL arrRTri[], vector<CGALTriangle2d*>& rColTri);

// helper function used to check if triangulation is a Delaunay one
	bool IsDealunay(const Point_2_CGAL& rA, const Point_2_CGAL& rB, const Point_2_CGAL& rC, const Point_2_CGAL& rD);
};


// helper classes used in subtraction algorithm
class thisTriangle {
public:
	Point_2_CGAL point;		// point position
	bool processed;			// has this point been processed yet
	bool outside;			// is point outside the other triangle
	bool boundary;			// is point on the boundary of the other triangle
	char cross;				// is this crossing point to the other triangle
	char noTriVert;			// the index if it is the vertex of the original this triangle

	thisTriangle(const Point_2_CGAL& tempPoint, bool tempProc, char tempNoTriVert) {
		point = tempPoint;
		processed = tempProc;
		noTriVert = tempNoTriVert;
	}
	thisTriangle(const Point_2_CGAL& tempPoint, bool tempProc, bool tempOut,
		bool tempBo, char tempCr, char tempNoTriVert) {
		point = tempPoint;
		processed = tempProc;
		outside = tempOut;
		boundary = tempBo;
		cross = tempCr;
		noTriVert = tempNoTriVert;
	}
};

class otherTriangle {
public:
	Point_2_CGAL point;		// point position
	bool processed;			// has this point been processed yet
	bool outside;			// is point outside the other triangle
	bool boundary;			// is point on the boundary of the other triangle
	char cross;				// is this crossing point to the other triangle

	otherTriangle(const Point_2_CGAL& tempPoint, bool tempProc) {
		point = tempPoint;
		processed = tempProc;
	}
	otherTriangle(const Point_2_CGAL& tempPoint, bool tempProc, char tempCr,
		bool tempOut, bool tempBo) {
		point = tempPoint;
		processed = tempProc;
		outside = tempOut;
		boundary = tempBo;
		cross = tempCr;
	}
};

class otherTemporary {
public:
	Point_2_CGAL point;		// point position
	char cross;				// is this crossing point to the other triangle
	char segNo;				// on which segment of other does intersection lie

	otherTemporary(const Point_2_CGAL& tempPoint, char tempCr, char tempSegNo) {
		point = tempPoint;
		cross = tempCr;
		segNo = tempSegNo;
	}
};

#endif