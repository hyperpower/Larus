/************************
 //  \file   Polygon.cpp
 //  \brief
 // 
 //  \author czhou
 //  \date   24 juil. 2014 
 ***********************/

#include "Polygon.h"
#include "Segment.h"
#include <iostream>

namespace Larus {

Polygon::Polygon() {
}

Polygon::Polygon(const Polygon &a) {
	this->AP = a.AP;
}

Polygon::Polygon(const ArrayP &a) {
	assert(a.size() >= 3);
	AP = a;
	trimSamePoints();
	assert(isSimple(AP));
}

bool Polygon::isSimple(const ArrayP &ap) {
	int nap = ap.size();
	//assert(nap >= 3);
	if (nap == 3) {
		return true;
	}
	int i = 0;
	for (int j = i + 2; j < nap - 1; j++) {
		if (isIntersect(ap[i], ap[i + 1], ap[j], ap[j + 1])) {
			return false;
		}
	}
	for (i = 1; i < nap - 2; i++) {
		for (int j = i + 2; j < nap; j++) {
			if (isIntersect(ap[i], ap[i + 1], ap[j], ap[(j + 1) % nap])) {
				return false;
			}
		}
	}
	return true;
}

void Polygon::trimSamePoints() {
	for (int i = 0; i < AP.size() - 1; i++) {
		if (AP[i] == AP[i + 1]) {
			AP.erase(i);
			i--;
		}
	}
	if (AP[0] == AP[AP.size() - 1]) {
		AP.pop_back();
	}
}

Polygon& Polygon::operator=(const Polygon &a) {
	if (this == &a) {
		return *this;
	} else {
		this->AP = a.AP;
	}
	return *this;
}

bool Polygon::isExist() const {
	if (AP.size() == 0) {
		return false;
	} else {
		return true;
	}
}

Float Polygon::calPerimeter() const{
	assert(isExist());
	Float res=0;
	for (int i = 1; i < AP.size() - 1; i++) {
		res+=calDistance(AP[i], AP[i+1]);
	}
	res+=calDistance(AP[AP.size()-1],AP[0]);
	return res;
}

Float Polygon::area() {
	if (!isExist()) {
		return 0.0;
	}
	Float s = 0.00;
	for (int i = 1; i < AP.size() - 1; i++) {
		s = s + CROSS(AP[i + 1], AP[i], AP[0]); // det to cro
	}
	return abs(s) / 2.00;
}

void Polygon::show() {
	for (int i = 0; i < AP.size(); i++) {
		std::cout << "> AP[ " << i << " ]=( " << AP[i].x << " , " << AP[i].y
				<< " )" << std::endl;
	}
}

void Polygon::reverse() {
	if (!isExist()) {
		return;
	} else {
		//ArrayP tmp = AP;
		//for (int i = 0, j = AP.size() - 1; i < AP.size(); i++, j--) {
		//	AP[i] = tmp[j];
		//}
		AP.reverse();
	}
}

arrayListT<Segment2D> Polygon::toArraySegment() const {
	int n = AP.size();
	arrayListT<Segment2D> as(n);
	for (int i = 0; i < n - 1; i++) {
		as[i].reconstruct(AP[i], AP[i + 1]);
	}
	as[n - 1].reconstruct(AP[n - 1], AP[0]);
	return as;
}

Point2D Polygon::getVertex(int i) const {
	assert(i < AP.size());
	return AP[i];
}

Segment2D Polygon::getSegment(int i) const {
	return Segment2D(AP[i], AP[(i + 1) % AP.size()]);
}

Float Polygon::getMaxX() const {
	assert(AP.size() > 0);
	Float max = AP[0].x;
	for (int i = 1; i < AP.size(); i++) {
		if (AP[i].x > max) {
			max = AP[i].x;
		}
	}
	return max;
}
Float Polygon::getMinX() const {
	assert(AP.size() > 0);
	Float min = AP[0].x;
	for (int i = 1; i < AP.size(); i++) {
		if (AP[i].x < min) {
			min = AP[i].x;
		}
	}
	return min;

}

Float Polygon::getMaxY() const {
	assert(AP.size() > 0);
	Float max = AP[0].y;
	for (int i = 1; i < AP.size(); i++) {
		if (AP[i].x > max) {
			max = AP[i].y;
		}
	}
	return max;
}
Float Polygon::getMinY() const {
	assert(AP.size() > 0);
	Float min = AP[0].y;
	for (int i = 1; i < AP.size(); i++) {
		if (AP[i].x < min) {
			min = AP[i].y;
		}
	}
	return min;

}


Float Polygon::WindingNum(const Point2D& ref, const Point2D& vi,
		 const Point2D& vip) const{
	Point2D refh(ref.x + 1.0, ref.y);
	Float wn = 0;
	Float a = CROSS(refh, vi, ref);
	Float b = CROSS(refh, vip, ref);
	if (a == 0 && b == 0) {
		return wn;
	}
	if (a * b < 0) {
		//vi vi+1 crosses the x
		Float c = CROSS(vip, ref, vi);
		if ((c > 0 && a<0)||(c < 0 && a>0)) {
			//vi vi+1 crosses the positive x
			if (a < 0) {
				wn++;
			} else {
				wn--;
			}
		}

	} else if (a == 0 && (vi.x > ref.x)) {
		if (b > 0) {
			wn = wn + 0.5;
		} else {
			wn = wn - 0.5;
		}
	} else if (b == 0 && (vip.x > ref.x)) {
		if (a < 0) {
			wn = wn + 0.5;
		} else {
			wn = wn - 0.5;
		}
	}
	return wn;
}

Float Polygon::calWindingNumber(const Point2D& ref) const {
	Float wn = 0;    // the  winding number counter
	// loop through all edges of the polygon
	for (int i = 0; i < AP.size() - 1; i++) {   // edge from V[i] to  V[i+1]
		wn+= WindingNum (ref, AP[i] , AP[i+1]);
	}
	return wn+=WindingNum (ref,AP[AP.size()-1],AP[0]);
}

bool Polygon::isOutsideofPolygon(const Point2D& ref) const {
	return (0 == calWindingNumber(ref)) ? true : false;
}

//===============================================
PolygonPlane3D::PolygonPlane3D(const ListT<Point3D>& list){
	if(list.empty()){
		return;
	}
	for(ListT<Point3D>::const_iterator iter=list.begin(); iter!=list.end(); iter++){
		this->push_back((*iter));
	}
}

}

