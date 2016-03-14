/////////////////////////////////////////////////////////////////////////////
//
// Triangle2dCGAL.cpp - CGALTriangle2d class with CGAL kernel
//
// Author: Marjan Sikora (sikora@fesb.hr)
//
/////////////////////////////////////////////////////////////////////////////


#include "triangle2dCGAL.h"

using namespace std;

long CGALTriangle2d::lNoExceptionsIntersect = 0;
long CGALTriangle2d::lNoExceptionsSubtract = 0;

// function for intersection of two triangles, intersects this with other
//   based on Southerland-Hodgson algorithm for polygon intersection
//
// pre:		- this is a valid CGAL 2D triangle
//			- other is a valid CGAL 2D triangle
// post:	- function adds resulting triangles to rColTri
void CGALTriangle2d::IntersectAndTriangulateSH(const CGALTriangle2d& other, vector<CGALTriangle2d*>& rColTri) {
	vector<Point_2_CGAL> oL;
	vector<Point_2_CGAL> iL;
	Line_2_CGAL edge[3];

	for (int i = 0; i < 3; ++i) {
		edge[i] = Line_2_CGAL(other._segment[i].vertex(0), other._segment[i].vertex(1));
	}
	oL.push_back(this->_point[0]);
	oL.push_back(this->_point[1]);
	oL.push_back(this->_point[2]);
	for (int i = 0; i < 3; ++i) {
		iL.clear();
		for (auto it : oL) {
			iL.push_back(it);
		}
		oL.clear();
		if (!iL.empty()) {
			Point_2_CGAL S = iL.back();
			for (auto E : iL) {
				if (CGAL::orientation(E, other._segment[i].vertex(0), other._segment[i].vertex(1)) == CGAL::LEFT_TURN){
					if (CGAL::orientation(S, other._segment[i].vertex(0), other._segment[i].vertex(1)) != CGAL::LEFT_TURN){
						auto v = CGAL::intersection(Segment_2_CGAL(S, E), edge[i]);
						if (v) {
							if (const Point_2_CGAL *p = boost::get<Point_2_CGAL>(&*v)) {
								oL.push_back(*p);
							} else {
								const Segment_2_CGAL *s = boost::get<Segment_2_CGAL>(&*v);
							}
						} else {
							/* empty intersection */
						}
					}
					oL.push_back(E);
				} else if (CGAL::orientation(S, other._segment[i].vertex(0), other._segment[i].vertex(1)) == CGAL::LEFT_TURN) {
					auto v = CGAL::intersection(Segment_2_CGAL(S, E), edge[i]);
					if (v) {
						if (const Point_2_CGAL *p = boost::get<Point_2_CGAL>(&*v)) {
							oL.push_back(*p);
						} else {
							const Segment_2_CGAL *s = boost::get<Segment_2_CGAL>(&*v);
						}
					} else {
						/* empty intersection */
					}
				}
				S = E;
			}
		}
	}
	if (oL.size() == 3) {
		rColTri.push_back(new CGALTriangle2d(oL[0], oL[1], oL[2]));
	} else if (oL.size() == 4){
		TriangulateFourConvexPoints(oL[0], oL[1], oL[2], oL[3], rColTri);
	} else if (oL.size() > 4){
		TriangulateMoreThanFourConvexPoints(oL, rColTri);
	} else if (oL.size() == 2 || oL.size() == 1 ) {
		++lNoExceptionsIntersect;
	}
}


// function for subtraction of two triangles, subtracts other from this
//   algorithm is based on technical report: Polygon Subtraction in Two or Three Dimensions, October 2013, JE Wilson
//   http://www.pnnl.gov/main/publications/external/technical_reports/PNNL-SA-97135.pdf
//
// pre:		- this is a valid CGAL 2D triangle
//			- other is a valid CGAL 2D triangle
// post:	- function adds resulting triangles to rColTri
void CGALTriangle2d::SubtractAndTriangulateWA(const CGALTriangle2d& other, vector<CGALTriangle2d*>& rColTri) {
// first table - stores this triangle data
	vector<thisTriangle> thisTr;
// which other vertex is colocated with vertex of this
	int thisColVertOther[3] = { -1, -1, -1 };
// where would I return other (using cross)
	int thisColVertOtherCross[3] = { -1, -1, -1 };
// which other edge has vertex of this
	int thisIsOnOtherSeg[3] = {-1, -1, -1};
// whivh vertex of this is colocated with other vertex
	int otherColVertThis[3] = { -1, -1, -1 };

// temporary intersections collection (used for transfer from first table (thisTr) 
//    to second one (otherTr)
	deque<otherTemporary> otherTemp;

// indexes: i is used for the first table (thisTr); o is for the second (otherTr) one
	int i = 0, o = 0;


/////////////////////////////////////////////////////
// 1. phase - procesing of vertices of both triangles
//            in order to create two triangle tables 
//			  containig vertices and intersections

///////// check for colocated vertices
	for (i = 0; i < 3; ++i) {
		for (o = 0; o < 3; ++o) {
			if (this->_point[i] == other._point[o]) {
				thisColVertOther[i] = o;
				otherColVertThis[o] = i;
			}
		}
	}

///////// loop that insert this vertices and checks for boundary situations
	for (i = 0; i < 3; ++i) {
// insert this vertex in thisTr
		thisTr.push_back(thisTriangle(this->_mCGALTr[i], false, i));
// is this vertex #i inside other?
		if (other._mCGALTr.has_on_bounded_side(this->_point[i])) {
			thisTr.back().outside = false;
			thisTr.back().boundary = false;
		} else if (other._mCGALTr.has_on_boundary(this->_point[i])) {
			thisTr.back().outside = false;
			thisTr.back().boundary = true;
		} else {
			thisTr.back().outside = true;
			thisTr.back().boundary = false;
		}
// is vertex #i colocated with any of vertices of edge #o
		if (thisColVertOther[i] != -1) {
			thisTr.back().cross = -2;
// to save the return (cross) position
			thisColVertOtherCross[i] = thisTr.size() - 1;
		} else {
// does any of other edges contain this vertex 
			if (thisTr.back().boundary) {
// check which other edge contains this vertex
				for (o = 0; 0 < 3; ++o) {
					if (other._segment[o].has_on(thisTr.back().point)) {
// record which edge contains this vertex, and cross data
						thisIsOnOtherSeg[i] = o;
						thisTr.back().cross = -3;
						break;
					}
				}
			} else {
				thisTr.back().cross = -1;
			}
		}
	}

//////////// loop that calculates intersections

// position on this where intersections will be inserted
	int thisInsertPos = 0;
// temporary variables for determinint the order of two intersections
	const Point_2_CGAL* tempInt1 = nullptr;
	const Point_2_CGAL* tempInt2 = nullptr;
// temporary indexes with edges on which intersections lie
	char tempIntRTriSegNo1, tempIntRTriSegNo2;

	for (i = 0; i < 3; ++i) {
		++thisInsertPos;
		for (o = 0; o < 3; o++) {

// special cases - vertex #i+1 is colocated with vertices #o or #o+1
			if( ( thisColVertOther[(i + 1) % 3] != o ) &&
				(thisColVertOther[(i + 1) % 3] != ( (o+1)%3) ) &&
// special cases - is vertex #i colocated with any of vertices of edge #o of other
				(thisColVertOther[i] != o) &&
				(thisColVertOther[i] != ((o + 1) % 3) ) ) {
// special case - check if vertex #i is on the edge #o of other
				if (thisIsOnOtherSeg[i] == o) {
					otherTemp.push_back(otherTemporary(this->_point[i], thisInsertPos - 1, o ));
// special case - skip if vertex #i+1 is on the edge #o+1 of other
				} else if ( thisIsOnOtherSeg[(i + 1) % 3] != o) {

//// trivial case - first part
// check if edge #i of this intersect with edge #o of other
					auto v = CGAL::intersection(this->_segment[i], other._segment[o]);
					if (v) {
						if (const Point_2_CGAL* p = boost::get<Point_2_CGAL>(&*v)) {
// if it is the first intersection
							if (tempInt1 == nullptr) {
								tempInt1 = new Point_2_CGAL(*p);
								tempIntRTriSegNo1 = o;
// if it is the second intersection
							} else {
// be sure it is not another special case:
//	 - if the second intersection is the same as the first one, then the vertex of other lies on the edge of this
								if (*tempInt1 != *p) {
// check which intersection is the nearest to the beggining of the edge and switch intersections if necessary
									if (CGAL::squared_distance(this->_point[i], *p) < CGAL::squared_distance(this->_point[i], *tempInt1)) {
										tempInt2 = tempInt1;
										tempIntRTriSegNo2 = tempIntRTriSegNo1;
										tempInt1 = new Point_2_CGAL(*p);
										tempIntRTriSegNo1 = o;
									} else {
										tempInt2 = new Point_2_CGAL(*p);
										tempIntRTriSegNo2 = o;
									}
									break;
								}
							}
						} else {
// if the intersection is segment, then edges are collinear
// and this is special case, taken care elswhere
						}
					}
				} 
// special cases - vertex #i or #i+1 lie on #o edge of other
			} else {
// special case - vertex #i lie on #o edge of other
				if (thisIsOnOtherSeg[i] == o) {
					otherTemp.push_back(otherTemporary(this->_point[i], thisInsertPos - 1, o ));
// special case - skip if vertex #i+1 lie on #o edge of other
				}
			}
		}

//// trivial case - second part
// if intersection exists add it to thisTr table
		if (tempInt1 != nullptr) {
			thisTr.insert(thisTr.begin() + thisInsertPos, thisTriangle(*tempInt1, false, false, true, 0, -1));
			++thisInsertPos;
// record the intersection in otherTemp, in order to later add it to otherTr table
			otherTemp.push_back(otherTemporary( *tempInt1, thisInsertPos - 1, tempIntRTriSegNo1));
// if there is second intersection, add it to thisTr table
			if (tempInt2 != nullptr) {
// be sure that this is not special case - check if second intersection is different from the first one
				if ( *tempInt1 != *tempInt2 ) {
					thisTr.insert(thisTr.begin() + thisInsertPos, thisTriangle(*tempInt2, false, false, true, 0, -1));
					++thisInsertPos;
// record the second intersection in otherTemp, in order to later add it to otherTr table
					otherTemp.push_back(otherTemporary(*tempInt2, thisInsertPos - 1, tempIntRTriSegNo2));
				}
// dealocation of temporary intersection
				if (tempInt2 != nullptr) {
					delete tempInt2;
					tempInt2 = nullptr;
				}
			}
// dealocation of temporary intersection
			if (tempInt1 != nullptr) {
				delete tempInt1;
				tempInt1 = nullptr;
			}
		}
	}

/////////// sort otherTemp table by the edge on which intersection lies

// deque index
	unsigned int tempI = 0;
// edge number of previous intersection
	int indLast = -1;
// loop over the otherTemp table
	while (tempI < otherTemp.size()) {
// switch intersections if their edge numbers are not in order
		if (otherTemp[tempI].segNo < indLast) {
// go to the end of container - pop from back and push to the beggining
			unsigned int i = tempI;
			while (i < otherTemp.size() ) {
				otherTemp.push_front(otherTemp.back());
				otherTemp.pop_back();
				++i;
			}
			break;
		} else {
			indLast = otherTemp[tempI].segNo;
		}
		++tempI;
	}

///////////// sort by distance if there are two intersections on the same edge of this

	tempI = 0;
	indLast = -1;
	while (tempI < otherTemp.size()) {
// if next intersection lies on the same edge, check which is closer to the beggining of the edge
		if (otherTemp[tempI].segNo == indLast) {
// if the second intersection is closer, switch them
			if (CGAL::squared_distance(other._point[indLast], otherTemp[tempI].point) < CGAL::squared_distance(other._point[indLast], otherTemp[tempI - 1].point)) {
				swap(otherTemp[tempI - 1], otherTemp[tempI]);
			}
// skip two
			++tempI;
			if (tempI >= otherTemp.size()){
				break;
			}
			indLast = otherTemp[tempI].segNo;
		} else {
			indLast = otherTemp[tempI].segNo;
		}
		++tempI;
	}

///////////// forming otherTr table from other vertices and intersections stored in otherTemp

// table with other data 
	vector<otherTriangle> otherTr;
// auxilary variable for special case
	int intNextCross = -1;

// creating otherTr table
	for (o = 0; o < 3; ++o) {
// insert the #o other vertex
		otherTr.push_back(otherTriangle(other._mCGALTr[o], false));
// be sure that it is not special case
		if (intNextCross == -1) {
// if it is not colocated vertex
			if (otherColVertThis[o] == -1) {
				otherTr.back().cross = -1;
			} else {
// if it is colocated vertex, I have to set cross indexes
// find thisPoint that vas a vertex with index equal to thisColVertOtherCross[otherColVertThis[o]]
				for (unsigned int i = 0; i < thisTr.size(); ++i) {
					if (thisTr[i].noTriVert == thisColVertOtherCross[otherColVertThis[o]]) {
						otherTr.back().cross = i;
						thisTr[i].cross = otherTr.size() - 1;
						break;
					}
				}
			}
// special case - intersection from last iteration
		} else {
			otherTr.back().cross = intNextCross;
			thisTr[intNextCross].cross = otherTr.size() - 1;
			intNextCross = -1;
		}
// is the vertex inside this or on its boundary?
		if (this->_mCGALTr.has_on_bounded_side(other._point[o]) ) {
			otherTr.back().outside = false;
			otherTr.back().boundary = false;
		} else if(	this->_mCGALTr.has_on_boundary(other._point[o]) ) {
			otherTr.back().outside = false;
			otherTr.back().boundary = true;
		} else {
			otherTr.back().outside = true;
			otherTr.back().boundary = false;
		}
// search otherTemp for index #o
		tempI = 0;
		while (tempI != otherTemp.size() && otherTemp[tempI].segNo != o) {
			++tempI;
		}
// check if index #o is found
		if (tempI != otherTemp.size()) {
// check if start or end vertex of edge is the same point as the intersection 
			bool s, e;
			s = (otherTemp[tempI].point == other._mCGALTr[o]);
			e = (otherTemp[tempI].point == other._mCGALTr[(o + 1) % 3]);
// not the same - it is real intersection
			if ( !s && !e ) {
				otherTr.push_back(otherTriangle(otherTemp[tempI].point, false, otherTemp[tempI].cross, false, true));
				thisTr[ otherTemp[tempI].cross ].cross = otherTr.size() - 1;
// same as start point - it does not intersect, but touches the edge
			} else if( s ) {
				otherTr.back().cross = otherTemp[tempI].cross;
				thisTr[ otherTemp[tempI].cross ].cross = otherTr.size() - 1;
// same as end point - it does not intersect, but touches the edge
// remember the index of intersection, and process it as a special case in the next iteration
			} else if( e ) {
				intNextCross = otherTemp[tempI].cross;
			}
			++tempI;
		}
// check for second intersection with this edge
		if (tempI != otherTemp.size() && otherTemp[tempI].segNo == o) {
// check if start or end vertex of edge is the same point as the intersection 
			bool s, e;
			s = (otherTemp[tempI].point == other._mCGALTr[o]);
			e = (otherTemp[tempI].point == other._mCGALTr[(o + 1) % 3]);
// not the same - it is real intersection
			if (!s && !e) {
				otherTr.push_back(otherTriangle(otherTemp[tempI].point, false, otherTemp[tempI].cross, false, true));
				thisTr[otherTemp[tempI].cross].cross = otherTr.size() - 1;
// same as start point - it does not intersect, but touches the edge
			} else if (s) {
				otherTr.back().cross = otherTemp[tempI].cross;
				thisTr[otherTemp[tempI].cross].cross = otherTr.size() - 1;
// same as end point - it does not intersect, but touches the edge
// remember the index of intersection, and process it as a special case in the next iteration
			} else if (e) {
				intNextCross = otherTemp[tempI].cross;
			}
			++tempI;
		}
	}
// special case - the intersection left for next iteration is in the last iteration
// it should be the first point
	if (intNextCross != -1) {
		otherTr[0].cross = intNextCross;
		thisTr[intNextCross].cross = 0;
		intNextCross = -1;
	}


/////////////////////////////////////////////////////////////
// 2. phase - determining the subtraction polygones 

///////////////// special boundary cases - the triangles don't intersect

// if otherTr has only three points, then this does not intersect 
	if (otherTr.size() == 3) {
		if (otherTr[0].outside == true || otherTr[1].outside == true || otherTr[2].outside == true) {
			bool thisTotOutside = true;
			unsigned int i = 0;
			while (i < thisTr.size()) {
				if (thisTr[i].outside == false && thisTr[i].boundary == false) {
					thisTotOutside = false;
					break;
				}
				i++;
			}
// if all vertices of otherTr are outside this 
			if (thisTotOutside) {
// if none of this vertices touch otherTr, then return this
				if ((thisTr[0].outside == true || thisTr[1].outside == true || thisTr[2].outside == true)) {
					rColTri.push_back(new CGALTriangle2d(*this));
					return;
				} else {
// at least one this vertex touches otherTr, then return none
					return;
				}
			}
		}
	}
// if this has only three points, then otherTr does not intersect 
	if(thisTr.size() == 3) {
		bool otherTotOutside = true;
		unsigned int i = 0;
		while (i < otherTr.size()) {
			if (otherTr[i].outside == false && otherTr[i].boundary == false) {
				otherTotOutside = false;
				break;
			}
			i++;
		}
		if (otherTotOutside) {
// otherTr is totaly outside this, return other
			if (thisTr[0].outside == true || thisTr[1].outside == true || thisTr[2].outside == true) {
					rColTri.push_back(new CGALTriangle2d(*this));
					return;
			} else {
// otherTr encompases this
				return;
			}
		} else {
// is otherTr totaly inside this
			if (otherTr[0].outside == false && otherTr[0].outside == false && otherTr[0].outside == false
				&& otherTr[0].boundary == false && otherTr[1].boundary == false && otherTr[2].boundary == false) {
					Point_2_CGAL arrT[3];
					arrT[0] = this->_point[0];
					arrT[1] = this->_point[1];
					arrT[2] = this->_point[2];
					Point_2_CGAL arrR[3];
					arrR[0] = other._point[0];
					arrR[1] = other._point[1];
					arrR[2] = other._point[2];
					TriangulateHole(arrT, arrR, rColTri);
					return;
			}
		}
	}

////////////// trivial case

// indexes of current point
	unsigned int indThis;
	unsigned int indOther = 0;
	bool processingThis;
// vector of subtracted polygones
	vector<Point_2_CGAL> resultPolyPoints;
// list of indexes of subtracted polygon - used to detect when the polygone is closed
	list<int> resultPolyInd;
	while (true) {
// find the first point that is outside other
		indThis = 0;
		while (indThis < thisTr.size() && (thisTr[indThis].outside == false || thisTr[indThis].processed == true)) {
			indThis++;
		}
		if (indThis == thisTr.size()) {
			break;
		}
// start with processing of this points
		processingThis = true;
		resultPolyPoints.clear();
		resultPolyInd.clear();
		while ( true ) {
			if (processingThis) {
// check if new point is the last of already inserted
// if so stop and form a new polygon
				if (!resultPolyInd.empty() && (indThis == resultPolyInd.front())) {
					break;
				}
// if not, add the point to the polygon
				resultPolyPoints.push_back(thisTr[indThis].point);
				resultPolyInd.push_back(indThis);
				thisTr[indThis].processed = true;
// is added point the intersection?
				if (thisTr[indThis].cross != -1) {
// if so, switch processing to other points
					processingThis = false;
					indOther = thisTr[indThis].cross;
					if( indOther == 0 ) {
						indOther = otherTr.size() - 1;
					} else {
						--indOther;
					}
// index of next point is indOther
				} else {
// if added point is not intersectio, increment indThis and iterate
					if (indThis == thisTr.size() - 1) {
						indThis = 0;
					} else {
						++indThis;
					}
				}
			} else {
				resultPolyPoints.push_back(otherTr[indOther].point);
				resultPolyInd.push_back(-1);
				otherTr[indOther].processed = true;
// is added point the intersection?
				if (otherTr[indOther].cross != -1) {
// if it is switch processing to this points
					processingThis = true;
					indThis = otherTr[indOther].cross;
					thisTr[indThis].processed = true;
					if (indThis == thisTr.size() -1) {
						indThis = 0;
					} else {
						++indThis;
					}
					thisTr[indThis].processed = true;
				} else {
// if added point is not intersection decrement indOther 
					if (indOther == 0) {
						indOther = otherTr.size() - 1;
					} else {
						--indOther;
					}
				}
			}
		} 

//cout << endl << endl;
//cout << "resultPoly: vertex coords----------------" << endl;
//int ind = 0;
//for (auto it : resultPolyPoints) {
//	cout << ind++ << " " << setw(10) << setprecision(5) << it.x()
//		<< " " << setw(10) << setprecision(5) << it.y() << endl;
//}

		if (resultPolyPoints.size() == 3) {
			rColTri.push_back(new CGALTriangle2d(resultPolyPoints[0], 
				resultPolyPoints[1], resultPolyPoints[2]));
		} else if (resultPolyPoints.size() == 4){
			TriangulateFourConvexPoints(resultPolyPoints[0],
				resultPolyPoints[1], resultPolyPoints[2], resultPolyPoints[3], rColTri);
		} else if (resultPolyPoints.size() == 5) {
			TriangulateFiveConcavePoints(resultPolyPoints, rColTri);
		} else if (resultPolyPoints.size() == 6) {
			TriangulateSixConcavePoints(resultPolyPoints, rColTri);
		} else if (resultPolyPoints.size() == 7) {
			TriangulateSevenConcavePoints(resultPolyPoints, rColTri);
		} else {
// error!
			++lNoExceptionsSubtract;
		}
	}
}


// function that triangulates four point polygon into two triangles, using Delaunay triangulation
// pre:		four points (rA, rB, rC i rD) positively oriented
// post:	two Dealunay triangles in vector rColTri
void CGALTriangle2d::TriangulateFourConvexPoints(const Point_2_CGAL& rA, const Point_2_CGAL& rB, 
	const Point_2_CGAL& rC, const Point_2_CGAL& rD, vector<CGALTriangle2d*>& rColTri) {
	if (IsDealunay(rA, rB, rC, rD)){
		rColTri.push_back(new CGALTriangle2d(rA, rB, rC));
		rColTri.push_back(new CGALTriangle2d(rA, rC, rD));
	} else  {
		rColTri.push_back(new CGALTriangle2d(rA, rB, rD));
		rColTri.push_back(new CGALTriangle2d(rB, rC, rD));
	}
}

// opis:	triangluates convex polygon with n points, with positive orientation
// prima:	vector arrPoints with at least five points
// vraæa:	vector rColTri with resulting triangles
void CGALTriangle2d::TriangulateMoreThanFourConvexPoints(const vector<Point_2_CGAL>& arrPoints, vector<CGALTriangle2d*>& rColTri) {
	int i, j, iX, iTemp;
	int i0, i1, i2;
// which point belongs to which edge - 0 is start, 1 is end, 2 is left, 3 is right
	int ep[3][4];
// which point belong to which triangle
	int tp[4][3];
// flag indicating the end of check
	bool bDTFinished = false;
// array belonging to the flipped edge - it contains all edges of its triangles
	int efliped[5][2];
// array belonging to the edge scheduled for checking
	int echecked[5][2];
// overlap of fliped and checked
	int iefc;

	int iN = arrPoints.size() - 3;
// setting initial values for edges
	for (i = 0; i < iN; i++) {
		ep[i][0] = 0;
		ep[i][1] = 2 + i;
		ep[i][2] = 1 + i;
		ep[i][3] = 3 + i;
	}
	do {
		for (iX = 0; iX < iN; iX++) {
			if (!IsDealunay(arrPoints[ep[iX][0]], arrPoints[ep[iX][2]], arrPoints[ep[iX][1]], arrPoints[ep[iX][3]])) {
				break;
			}
		}
// if checking is over, and all is OK
		if (iX == iN) {
			bDTFinished = true;
		} else {
// flips the edge, then reviews and sets all others
// creates the edges of two adjacent triangles
			efliped[0][0] = min(ep[iX][0], ep[iX][1]); efliped[0][1] = max(ep[iX][0], ep[iX][1]);
			efliped[1][0] = min(ep[iX][0], ep[iX][2]); efliped[1][1] = max(ep[iX][0], ep[iX][2]);
			efliped[2][0] = min(ep[iX][0], ep[iX][3]); efliped[2][1] = max(ep[iX][0], ep[iX][3]);
			efliped[3][0] = min(ep[iX][1], ep[iX][2]); efliped[3][1] = max(ep[iX][1], ep[iX][2]);
			efliped[4][0] = min(ep[iX][1], ep[iX][3]); efliped[4][1] = max(ep[iX][1], ep[iX][3]);
// flipping the edge iX
			iTemp = ep[iX][0];
			ep[iX][0] = ep[iX][2];
			ep[iX][2] = ep[iX][1];
			ep[iX][1] = ep[iX][3];
			ep[iX][3] = iTemp;
// check other edges
			for (int i = 0; i < iN; i++) {
				if (i != iX) {
// is this neigbouring one
// creates the edges of two adjacent triangles
					echecked[0][0] = min(ep[i][0], ep[i][1]); echecked[0][1] = max(ep[i][0], ep[i][1]);
					echecked[1][0] = min(ep[i][0], ep[i][2]); echecked[1][1] = max(ep[i][0], ep[i][2]);
					echecked[2][0] = min(ep[i][0], ep[i][3]); echecked[2][1] = max(ep[i][0], ep[i][3]);
					echecked[3][0] = min(ep[i][1], ep[i][2]); echecked[3][1] = max(ep[i][1], ep[i][2]);
					echecked[4][0] = min(ep[i][1], ep[i][3]); echecked[4][1] = max(ep[i][1], ep[i][3]);
					iefc = 0;
					for (int j = 0; j < 5 && iefc < 2; j++) {
						for (int k = 0; k < 5 && iefc < 2; k++) {
							if (efliped[j][0] == echecked[k][0] && efliped[j][1] == echecked[k][1]) {
								iefc++;
							}
						}
					}
					if (iefc > 1) {
						if ((ep[i][2] == ep[iX][2]) ||
							(ep[i][2] == ep[iX][3])) {
							if ((ep[iX][0] != ep[i][0]) &&
								(ep[iX][0] != ep[i][1])) {
								ep[i][2] = ep[iX][0];
							} else if ((ep[iX][1] != ep[i][0]) &&
								(ep[iX][1] != ep[i][1])) {
								ep[i][2] = ep[iX][1];
							}
						} else if ((ep[i][3] == ep[iX][2]) ||
							(ep[i][3] == ep[iX][3])) {
							if ((ep[iX][0] != ep[i][0]) &&
								(ep[iX][0] != ep[i][1])) {
								ep[i][3] = ep[iX][0];
							} else if ((ep[iX][1] != ep[i][0]) &&
								(ep[iX][1] != ep[i][1])) {
								ep[i][3] = ep[iX][1];
							}
						}
					}
				}
			}
		}
	} while (!bDTFinished);
// triangulation done, creating triangles
	iX = 0;
	for (i = 0; i < iN; i++) {
// creating left triangle, with indexes sorted accending
		if (ep[i][0] < min(ep[i][1], ep[i][2])) {
			i0 = ep[i][0];
			if (ep[i][1] < ep[i][2]) {
				i1 = ep[i][1];
				i2 = ep[i][2];
			} else {
				i1 = ep[i][2];
				i2 = ep[i][1];
			}
		} else if (ep[i][1] < min(ep[i][0], ep[i][2])) {
			i0 = ep[i][1];
			if (ep[i][0] < ep[i][2]) {
				i1 = ep[i][0];
				i2 = ep[i][2];
			} else {
				i1 = ep[i][2];
				i2 = ep[i][0];
			}
		} else {
			i0 = ep[i][2];
			if (ep[i][0] < ep[i][1]) {
				i1 = ep[i][0];
				i2 = ep[i][1];
			} else {
				i1 = ep[i][1];
				i2 = ep[i][0];
			}
		}
// add it if it is not already there
		for (j = 0; j < iX; j++) {
			if (tp[j][0] == i0 && tp[j][1] == i1 && tp[j][2] == i2) {
				break;
			}
		}
		if (j == iX) {
			tp[iX][0] = i0;
			tp[iX][1] = i1;
			tp[iX][2] = i2;
			iX++;
		}
// creating left triangle, with indexes sorted accending
		if (ep[i][0] < min(ep[i][1], ep[i][3])) {
			i0 = ep[i][0];
			if (ep[i][1] < ep[i][3]) {
				i1 = ep[i][1];
				i2 = ep[i][3];
			} else {
				i1 = ep[i][3];
				i2 = ep[i][1];
			}
		} else if (ep[i][1] < min(ep[i][0], ep[i][3])) {
			i0 = ep[i][1];
			if (ep[i][0] < ep[i][3]) {
				i1 = ep[i][0];
				i2 = ep[i][3];
			} else {
				i1 = ep[i][3];
				i2 = ep[i][0];
			}
		} else {
			i0 = ep[i][3];
			if (ep[i][0] < ep[i][1]) {
				i1 = ep[i][0];
				i2 = ep[i][1];
			} else {
				i1 = ep[i][1];
				i2 = ep[i][0];
			}
		}
// add it if it is not already there
		for (j = 0; j < iX; j++) {
			if (tp[j][0] == i0 && tp[j][1] == i1 && tp[j][2] == i2) {
				break;
			}
		}
		if (j == iX) {
			tp[iX][0] = i0;
			tp[iX][1] = i1;
			tp[iX][2] = i2;
			iX++;
		}
	}
	for (i = 0; i < iX; i++) {
		rColTri.push_back(new CGALTriangle2d(arrPoints[tp[i][0]], arrPoints[tp[i][1]], arrPoints[tp[i][2]]));
	}
}


// function triangulates concave polygon with 5 points (not checking for Dealunay criterion)
// pre:		vector arrP with 5 points
// post:	vector rColTri with 3 triangles
void CGALTriangle2d::TriangulateFiveConcavePoints(const vector<Point_2_CGAL>& arrP, vector<CGALTriangle2d*>& rColTri) {
// find which point has angle larger than 180 degrees 
// that point is pivot - vertex of other triangle penetrating the firs one
// from that point create edges to #vertex+2, #vertex+3 i #vertex+4

	int i, pi, si;

// find pivot
	for (i = 0; i < 5; i++) {
		pi = (i - 1 + 5) % 5;
		si = (i + 1) % 5;
		if (Line_2_CGAL(arrP[i], arrP[pi]).has_on_positive_side(arrP[si])) {
			break;
		}
	}
	if (i == 5) {
// pivot is point that is not on the line 
		int j;
		for (j = 0; j < 5; ++j) {
// make line from #j+1 to #j+2
// if #j is on the line it is not pivot
			Line_2_CGAL l = Line_2_CGAL( arrP[(j+1)%5], arrP[(j+2)%5] );
			if (!l.has_on(arrP[j]) && l.has_on(arrP[(j + 3) % 5])) {
				break;
			}
		}
		if (j == 5) {
			++lNoExceptionsSubtract;
// error!
		}
		rColTri.push_back(new CGALTriangle2d(arrP[(j +4) % 5], arrP[j], arrP[(j + 3) % 5]));
		rColTri.push_back(new CGALTriangle2d(arrP[(j + 3) % 5], arrP[j], arrP[(j + 2) % 5]));
		rColTri.push_back(new CGALTriangle2d(arrP[(j + 2) % 5], arrP[j], arrP[(j + 1) % 5]));
	} else {
		pi = (i + 1) % 5;	si = (i + 2) % 5;
		rColTri.push_back(new CGALTriangle2d(arrP[i], arrP[pi], arrP[si]));
		pi = (i + 2) % 5;	si = (i + 3) % 5;
		rColTri.push_back(new CGALTriangle2d(arrP[i], arrP[pi], arrP[si]));
		pi = (i + 3) % 5;	si = (i + 4) % 5;
		rColTri.push_back(new CGALTriangle2d(arrP[i], arrP[pi], arrP[si]));
	}

}

// function triangulates concave polygon with 6 points (not checking for Dealunay criterion)
// pre:		vector arrP with 6 points
// post:	vector rColTri with 4 triangles
void CGALTriangle2d::TriangulateSixConcavePoints(const vector<Point_2_CGAL>& arrP, vector<CGALTriangle2d*>& rColTri) {
// find which point has angle larger than 180 degrees 
// that point is pivot - vertex of other triangle penetrating the firs one
// from that point create edges to #vertex+2, #vertex+3 i #vertex+4

	int i, j, pi, si;

// checking for boundary case
	for (i = 0; i < 5; i++) {
		for (j = i+1; j < 6; j++) {
			if (arrP[i] == arrP[j]) {
				pi = (i + 1) % 6;	si = (i - 1 + 6) % 6;
				rColTri.push_back(new CGALTriangle2d(arrP[i], arrP[pi], arrP[si]));
				if (Line_2_CGAL(arrP[(i + 1) % 6], arrP[(i + 2) % 6]).has_on_positive_side(arrP[(j + 1) % 6])) {
					rColTri.push_back(new CGALTriangle2d(arrP[(i + 1) % 6], arrP[(i - 2 + 6) % 6], arrP[(i - 1 + 6) % 6]));
					rColTri.push_back(new CGALTriangle2d(arrP[(i + 2) % 6], arrP[(j + 1) % 6], arrP[(i + 1) % 6]));
				} else {
					rColTri.push_back(new CGALTriangle2d(arrP[(i+1)%6], arrP[(i+2)%6], arrP[(i-1+6)%6]));
					rColTri.push_back(new CGALTriangle2d(arrP[(i + 2) % 6], arrP[(j + 1) % 6], arrP[(j + 2) % 6]));
				}
				pi = (j + 1) % 6;	si = (j - 1 + 6) % 6;
				rColTri.push_back(new CGALTriangle2d(arrP[j], arrP[pi], arrP[si]));
				return;
			}
		}
	}

// find pivot
	for (i = 0; i < 6; i++) {
		pi = (i - 1 + 6) % 6;
		si = (i + 1) % 6;
		if (Line_2_CGAL(arrP[i], arrP[pi]).has_on_positive_side(arrP[si])) {
			break;
		}
	}
	if (i == 6) {
		++lNoExceptionsSubtract;
// error!
	}
	pi = (i + 1) % 6;	si = (i + 2) % 6;
	rColTri.push_back(new CGALTriangle2d(arrP[i], arrP[pi], arrP[si]));
	pi = (i + 2) % 6;	si = (i + 3) % 6;
	rColTri.push_back(new CGALTriangle2d(arrP[i], arrP[pi], arrP[si]));
	pi = (i + 3) % 6;	si = (i + 4) % 6;
	rColTri.push_back(new CGALTriangle2d(arrP[i], arrP[pi], arrP[si]));
	pi = (i + 4) % 6;	si = (i + 5) % 6;
	rColTri.push_back(new CGALTriangle2d(arrP[i], arrP[pi], arrP[si]));

}

// function triangulates concave polygon with 7 points (not checking for Dealunay criterion)
// pre:		vector arrP with 7 points
// post:	vector rColTri with 4 triangles
void CGALTriangle2d::TriangulateSevenConcavePoints(const vector<Point_2_CGAL>& arrP, vector<CGALTriangle2d*>& rColTri) {
// find which point has angle larger than 180 degrees 
// that point is pivot - vertex of other triangle penetrating the firs one
// from that point create edges to #vertex+2, #vertex+3 i #vertex+4

	int i, pi, si;
	int j, pj, sj;

// find pivot
	for (i = 0; i < 7; i++) {
		pi = (i - 1 + 7) % 7;
		si = (i + 1) % 7;
		if (Line_2_CGAL(arrP[i], arrP[pi]).has_on_positive_side(arrP[si])) {
			j = (i + 1) % 7;
			pj = (j - 1 + 7) % 7;
			sj = (j + 1) % 7;
			if (Line_2_CGAL(arrP[j], arrP[pj]).has_on_positive_side(arrP[sj])) {
				break;
			}
		}
	}
	if (i == 7) {
		++lNoExceptionsSubtract;
// error!
	}
	pj = (j + 1) % 7;	sj = (j + 2) % 7;
	rColTri.push_back(new CGALTriangle2d(arrP[j], arrP[pj], arrP[sj]));
	pj = (j + 2) % 7;	sj = (j + 3) % 7;
	rColTri.push_back(new CGALTriangle2d(arrP[j], arrP[pj], arrP[sj]));	
	if (Line_2_CGAL(arrP[i], arrP[j]).has_on_positive_side(arrP[sj])) {
		rColTri.push_back(new CGALTriangle2d(arrP[j], arrP[sj], arrP[i]));
		pi = (i + 4) % 7;	si = (i + 5) % 7;
		rColTri.push_back(new CGALTriangle2d(arrP[i], arrP[pi], arrP[si]));
	} else {
		pj = (j + 3) % 7;	sj = (j + 4) % 7;
		rColTri.push_back(new CGALTriangle2d(arrP[j], arrP[pj], arrP[sj]));
		rColTri.push_back(new CGALTriangle2d(arrP[j], arrP[sj], arrP[i]));
	}
	pi = (i + 5) % 7;	si = (i + 6) % 7;
	rColTri.push_back(new CGALTriangle2d(arrP[i], arrP[pi], arrP[si]));
}

// function checks if four points are triangulated according to Delaunay
// pre:		4 points (rA, rB, rC i rD), positive order
// post:	returns true if triangulation using edge AC is a Delaunay one
bool CGALTriangle2d::IsDealunay(const Point_2_CGAL& rA, const Point_2_CGAL& rB, const Point_2_CGAL& rC, const Point_2_CGAL& rD) {
// matrix coefficients
	double f11, f12, f13;
	double f21, f22, f23;
	double f31, f32, f33;

	f11 = CGAL::to_double(rA.x()) - CGAL::to_double(rD.x());
	f21 = CGAL::to_double(rB.x()) - CGAL::to_double(rD.x());
	f31 = CGAL::to_double(rC.x()) - CGAL::to_double(rD.x());
	f12 = CGAL::to_double(rA.y()) - CGAL::to_double(rD.y());
	f22 = CGAL::to_double(rB.y()) - CGAL::to_double(rD.y());
	f32 = CGAL::to_double(rC.y()) - CGAL::to_double(rD.y());
	f13 = pow(f11, 2) + pow(f12, 2);
	f23 = pow(f21, 2) + pow(f22, 2);
	f33 = pow(f31, 2) + pow(f32, 2);

// if the determinant is positive then triangulating edge should be BD, else AC
	if (f11*f22*f33 - f11*f23*f32 - f12*f21*f33 + f12*f23*f31 + f13*f21*f32 - f13*f22*f31 > 0.0) {
		return false;
	} else {
		return true;
	}
}

// function triangulates concave polygon with a hole (not checking for Dealunay criterion)
// pre:		two triangles, the first one is outer boundary, the second one is inner boundary
// post:	vector rColTri with triangles
void CGALTriangle2d::TriangulateHole(Point_2_CGAL arrThis[], Point_2_CGAL arrRTri[], vector<CGALTriangle2d*>& rColTri) {
	double dist, temp;
	int rTriZero;

// which other vertex is the closest one to this[0]
	dist = CGAL::to_double(squared_distance(arrThis[0], arrRTri[0]));
	rTriZero = 0;
	temp = CGAL::to_double(squared_distance(arrThis[0], arrRTri[1]));
	if ( temp < dist) {
		dist = temp;
		rTriZero = 1;
	}
	temp = CGAL::to_double(squared_distance(arrThis[0], arrRTri[2]));
	if( temp < dist) {
		rTriZero = 2;
	}

// triangulate three four points polygons 
	int t0, t1, r0, r1;

	t0 = 0; t1 = 1;
	r0 = rTriZero; r1 = (rTriZero + 1 ) % 3;
	if (Line_2_CGAL(arrRTri[r1], arrRTri[r0]).has_on_positive_side(arrThis[t0])) {
		rColTri.push_back(new CGALTriangle2d(arrThis[t0], arrThis[t1], arrRTri[r1]));
		rColTri.push_back(new CGALTriangle2d(arrThis[t0], arrRTri[r1], arrRTri[r0]));
	} else {
		rColTri.push_back(new CGALTriangle2d(arrThis[t0], arrThis[t1], arrRTri[r0]));
		rColTri.push_back(new CGALTriangle2d(arrThis[t1], arrRTri[r1], arrRTri[r0]));
	}

	t0 = 1; t1 = 2;
	r0 = (rTriZero + 1) % 3; r1 = (rTriZero + 2) % 3;
	if (Line_2_CGAL(arrRTri[r1], arrRTri[r0]).has_on_positive_side(arrThis[t0])) {
		rColTri.push_back(new CGALTriangle2d(arrThis[t0], arrThis[t1], arrRTri[r1]));
		rColTri.push_back(new CGALTriangle2d(arrThis[t0], arrRTri[r1], arrRTri[r0]));
	} else {
		rColTri.push_back(new CGALTriangle2d(arrThis[t0], arrThis[t1], arrRTri[r0]));
		rColTri.push_back(new CGALTriangle2d(arrThis[t1], arrRTri[r1], arrRTri[r0]));
	}

	t0 = 2; t1 = 0;
	r0 = ( rTriZero + 2 ) % 3; r1 = rTriZero;
	if (Line_2_CGAL(arrRTri[r1], arrRTri[r0]).has_on_positive_side(arrThis[t0])) {
		rColTri.push_back(new CGALTriangle2d(arrThis[t0], arrThis[t1], arrRTri[r1]));
		rColTri.push_back(new CGALTriangle2d(arrThis[t0], arrRTri[r1], arrRTri[r0]));
	} else {
		rColTri.push_back(new CGALTriangle2d(arrThis[t0], arrThis[t1], arrRTri[r0]));
		rColTri.push_back(new CGALTriangle2d(arrThis[t1], arrRTri[r1], arrRTri[r0]));
	}
}