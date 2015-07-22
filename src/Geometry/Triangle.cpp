/************************
 //  \file   Triangle.cpp
 //  \brief
 // 
 //  \author czhou
 //  \date   3 oct. 2014 
 ***********************/

#include "Triangle.h"
#include "Vector.h"
#include "Line.h"
#include "Relation.h"
#include "../Utility/ArrayList.h"
#include <iostream>

namespace Larus {

Triangle2D::Triangle2D(const Point2D& a, const Point2D& b, const Point2D& c) :
		array_3<Point2D>(a, b, c) {
}

void Triangle2D::show() const {
	std::cout << "Triangle2D \n";
	std::cout << "Point 0 (" << elems[0].x << ", " << elems[0].y << ")\n";
	std::cout << "Point 1 (" << elems[1].x << ", " << elems[1].y << ")\n";
	std::cout << "Point 2 (" << elems[2].x << ", " << elems[2].y << ")\n";
}

bool Triangle2D::isExist() const {
	if (elems[0] == elems[1] && elems[1] == elems[2] && elems[2] == elems[0]) {
		return false;
	} else {
		return true;
	}
}

bool Triangle2D::isColine() const {
	if (0.0 == CROSS(elems[0], elems[1], elems[2])) {
		return true;
	} else {
		return false;
	}
}

//=====================================

bool Triangle3D::isIEqual(idx_t i) const {
	if (elems[0][i] == elems[1][i] && elems[1][i] == elems[2][i]
			&& elems[2][i] == elems[0][i]) {
		return true;
	} else {
		return false;
	}
}

Triangle3D::Triangle3D(const Point3D& a, const Point3D& b, const Point3D& c) :
		array_3<Point3D>(a, b, c) {
}

void Triangle3D::show() const {
	std::cout << "Triangle3D \n";
	std::cout << "Point 0 (" << elems[0].x << ", " << elems[0].y << ", "
			<< elems[0].z << ")\n";
	std::cout << "Point 1 (" << elems[1].x << ", " << elems[1].y << ", "
			<< elems[1].z << ")\n";
	std::cout << "Point 2 (" << elems[2].x << ", " << elems[2].y << ", "
			<< elems[2].z << ")\n";
}

void Triangle3D::Transfer(Float dx, Float dy, Float dz) {
	for (int i = 0; i < 3; i++) {
		elems[i].x += dx;
		elems[i].y += dy;
		elems[i].z += dz;
	}
}

Triangle2D Triangle3D::projectXY() const {
	Point2D a(elems[0].x, elems[0].y);
	Point2D b(elems[1].x, elems[1].y);
	Point2D c(elems[2].x, elems[2].y);
	return Triangle2D(a, b, c);
}
Triangle2D Triangle3D::projectYZ() const {
	Point2D a(elems[0].y, elems[0].z);
	Point2D b(elems[1].y, elems[1].z);
	Point2D c(elems[2].y, elems[2].z);
	return Triangle2D(a, b, c);
}
Triangle2D Triangle3D::projectZX() const {
	Point2D a(elems[0].z, elems[0].x);
	Point2D b(elems[1].z, elems[1].x);
	Point2D c(elems[2].z, elems[2].x);
	return Triangle2D(a, b, c);
}

Triangle2D Triangle3D::project() const {
	if ((!isXEqual()) && (!isYEqual())) {
		return projectXY();
	}
	if ((!isYEqual()) && (!isZEqual())) {
		return projectYZ();
	}
	if ((!isZEqual()) && (!isXEqual())) {
		return projectZX();
	}
	ASSERT_MSG(false, "triangle project fail");
	return Triangle2D();
}

bool Triangle3D::isExist() const {
	if (elems[0] == elems[1] && elems[1] == elems[2] && elems[2] == elems[0]) {
		return false;
	} else {
		return true;
	}
}

bool Triangle3D::isColine() const {
	Triangle2D xy = this->projectXY();
	Triangle2D yz = this->projectYZ();
	Triangle2D zx = this->projectZX();
	if (xy.isColine() && yz.isColine() && zx.isColine()) {
		return true;
	} else {
		return false;
	}
}

int Triangle3D::onWhichSide(const Point3D& p) const {
	Float res = CROSS(elems[0], elems[1], p, elems[2]);
	if (res == 0.0) {
		return 0;
	} else if (res > 0.0) {
		return 1;
	} else {
		return -1;
	}
}

Plane Triangle3D::getPlane() const {
	Vector3D vab, vac;
	for (idx_t i = 0; i < 3; i++) {
		vab[i] = elems[1][i] - elems[0][i];
		vac[i] = elems[2][i] - elems[0][i];
	}
	Vector3D nout = cross(vab, vac);
	Float D = (nout[0] * elems[0].x + nout[1] * elems[0].y
			+ nout[2] * elems[0].z);
	return Plane(nout[0], nout[1], nout[2], D);
}

void Triangle3D::reverse() {
	Point3D tmp(elems[0]);
	elems[0] = elems[2];
	elems[2] = tmp;
}

bool Triangle3D::isXEqual() const {
	return isIEqual(0);
}
bool Triangle3D::isYEqual() const {
	return isIEqual(1);
}
bool Triangle3D::isZEqual() const {
	return isIEqual(2);
}

Segment3D Triangle3D::getSegment(idx_t i) const {
	ASSERT(i < 3);
	return Segment3D(elems[i], elems[(i + 1) % 3]);
}

//Triangle instersection is base on Devillers and Guigue Algorithm

//step 1;
void T1_instersect_PI2(const Triangle3D& t1, const Triangle3D& t2,
		arrayList_int& res) {
	for (idx_t i = 0; i < 3; i++) {
		res[i] = t2.onWhichSide(t1[i]);
	}
}
//step 2 analyse cases
int analyse_cases(arrayList_int& arr) {
	//case 1 same sign
	if ((arr[0] == 1 && arr[1] == 1 && arr[2] == 1)
			|| (arr[0] == -1 && arr[1] == -1 && arr[2] == -1)) {
		return -1;
	}
	//case 2 all equal to zero, co-plane
	int num_zero = arr.countEq(0);
	if (num_zero == 3) {
		return -2;  //==>>> 2D
	}
	//case 3 diffrent signs
	//3-2 2 zero
	if (num_zero == 2) {  //two points on PI
		return 2;  //==>>> 2D
	}
	//3-1 1 zero
	if (num_zero == 1) { //only one point on PI
		if (arr.countEq(1) == 2 || arr.countEq(-1) == 2) {
			return -3; //other two are at same side  //==>>> 2D
		}
		//other two are at different side
	}
	return 1;  //general case 3D
}
Triangle3D permutation(const Triangle3D& t, idx_t i) {
	return Triangle3D(t[i], t[(i + 1) % 3], t[(i + 2) % 3]);
}
void reverse(Triangle3D& t, idx_t i) {
	t.swap(((i + 4) % 3), ((i + 2) % 3));
}
//make sure the intersect type is TRI_GEN
void cal_intersect_general(const Triangle3D& T1,
		const Triangle3D& T2, Segment3D& seg){
	const Point_3D<Float>& p1 = T1[0];
	const Point_3D<Float>& q1 = T1[1];
	const Point_3D<Float>& r1 = T1[2];
	const Point_3D<Float>& p2 = T2[0];
	const Point_3D<Float>& q2 = T2[1];
	const Point_3D<Float>& r2 = T2[2];
	Point3D ps, pe;
	Float cjk=CROSS(p1,r1,q2,p2);
	Float cil=CROSS(p1,q1,r2,p2);
	Line3D line;
	Plane plane;
	if(cjk>=0){ // get k
		line=T2.getSegment(0).getLine();
		plane=T1.getPlane();
	}else{  //get j
		line=T1.getSegment(0).getLine();
		plane=T2.getPlane();
	}
	ps=intersect(line,plane);
	if(cil>=0){ // get i
		pe=intersect(T1.getSegment(0).getLine(),T2.getPlane());
	}else{     //get l
		pe=intersect(T2.getSegment(2).getLine(),T1.getPlane());
	}

	seg.reconstruct(ps,pe);
}

//case 7
TRI_INTERSECT_TYPE general_case(const Triangle3D& t1, const Triangle3D& t2,
		arrayList_int& arr1, arrayList_int& arr2, Segment3D& seg) {
	int count_1 = 0, count_2 = 0;
	idx_t idx1p = 0, idx1m = 0;
	idx_t idx2p = 0, idx2m = 0;
	for (idx_t i = 0; i < 3; i++) {
		if (arr1[i] == 1) {
			count_1 += 1;
			idx1p = i;
		} else {  //zero and -1 are -1
			idx1m = i;
		}
		if (arr2[i] == 1) {
			count_2 += 1;
			idx2p = i;
		} else {  //zero and -1 are -1
			idx2m = i;
		}
	}
	Triangle3D t1tmp;
	Triangle3D t2tmp(t2);
	if (count_1 == 1) {
		t1tmp = permutation(t1, idx1p);
	} else { // has 2 one
		t1tmp = permutation(t1, idx1m);
		if (count_2 == 1) {
			reverse(t2tmp, idx2p);
		} else {
			reverse(t2tmp, idx2m);
		}
	}
	if (count_2 == 1) {
		t2tmp = permutation(t2, idx2p);
	} else { // has 2 one
		t2tmp = permutation(t2, idx2m);
		if (count_1 == 1) {
			reverse(t1tmp, idx1p);
		} else {
			reverse(t1tmp, idx1m);
		}
	}
	//t1tmp.show();
	//t2tmp.show();
	const Point_3D<Float>& p1 = t1tmp[0];
	const Point_3D<Float>& q1 = t1tmp[1];
	const Point_3D<Float>& r1 = t1tmp[2];
	const Point_3D<Float>& p2 = t2tmp[0];
	const Point_3D<Float>& q2 = t2tmp[1];
	const Point_3D<Float>& r2 = t2tmp[2];
	//decision tree
	//i k
	Float cik= CROSS(p1, q1, q2, p2);
	if(cik >= 0){
		return TRI_NO_INTERSECT;
	}else{ //cik<0 i>k
		Float cjl= CROSS(p1, r1, r2, p2);
		if(cjl<=0){
			return TRI_NO_INTERSECT;
		}else{ //cjl>0 j<l
			cal_intersect_general(t1tmp,t2tmp,seg);
			return TRI_GEN ;
		}
	}
}


//main
TRI_INTERSECT_TYPE get_intersect_type(const Triangle3D& T1,
		const Triangle3D& T2, Segment3D& seg) {
	//T1 intersect PI2
	arrayList_int acase1(3);
	T1_instersect_PI2(T1, T2, acase1);
	for (idx_t i = 0; i < 3; i++) {
		std::cout << "arr1 " << acase1[i] << std::endl;
	}
	int ct1 = analyse_cases(acase1);
	if (ct1 < 0) {
		return TRI_NO_INTERSECT;
	} else if (ct1 == 2) { // 2D condition;
		//==>>> 2D
	} else {  //Additional tests
		arrayList_int acase2(3);
		T1_instersect_PI2(T2, T1, acase2);
		for (idx_t i = 0; i < 3; i++) {
			std::cout << "arr2 " << acase2[i] << std::endl;
		}
		int ct2 = analyse_cases(acase2);
		if (ct2 ==2 ) {  // 2D condition;

		} else {  // 3D condition
			return general_case(T1, T2, acase1, acase2, seg);
		}
	}
	return TRI_NO_INTERSECT;
}

bool isIntersect(const Triangle3D& T1, const Triangle3D& T2){
	Segment3D seg;
	TRI_INTERSECT_TYPE type=get_intersect_type(T1,T2,seg);
	if(type==TRI_NO_INTERSECT){
		return false;
	}else{
		return true;
	}
}

bool cal_intersect(const Triangle3D& T1,const Triangle3D& T2, Segment3D&seg){
	TRI_INTERSECT_TYPE type=get_intersect_type(T1,T2,seg);
	if(type==TRI_NO_INTERSECT){
		return false;
	}else{
		return true;
	}
}

}  //end of namespace
