/************************
 //  \file   Polygon.h
 //  \brief
 // 
 //  \author czhou
 //  \date   24 juil. 2014 
 ***********************/
#ifndef POLYGON_H_
#define POLYGON_H_

#include <assert.h>
#include "../TypeDef.h"
#include "../Utility/ArrayList.h"
#include "../Utility/List.h"
#include "Point.h"
#include "Segment.h"

namespace Larus {

class Polygon: public ObjectBase {
private:
	typedef arrayListT<Point2D> ArrayP;
	bool isSimple(const ArrayP &ap);
	void trimSamePoints();
	Float WindingNum(const Point2D& ref, const Point2D& vi,
			const Point2D& vip) const;
	ArrayP AP;
public:
	Polygon(); //
	Polygon(const Polygon &a); //
	Polygon(const ArrayP &a);

	Polygon& operator=(const Polygon &a);
	Float area();
	bool isExist() const; //
	void show();
	void reverse();

	arrayListT<Segment2D> toArraySegment() const;

	inline int getNumVertexs() const; //
	inline int getNumSegments() const;  //
	Point2D getVertex(int i) const;  //
	Segment2D getSegment(int i) const;

	Float getMaxX() const; //
	Float getMinX() const; //
	Float getMaxY() const; //
	Float getMinY() const; //

	Float calPerimeter() const;

	Float calWindingNumber(const Point2D& ref) const;
	bool isOutsideofPolygon(const Point2D&) const;

};

inline int Polygon::getNumVertexs() const {
	return AP.size();
}

inline int Polygon::getNumSegments() const {
	return AP.size();
}
//==========================================
template<typename POINT>
class PolygonT: public ListT<POINT> {
public:
	PolygonT() {
	}
	//PolygonT(const Polygon &a);
	//Polygon& operator=(const Polygon &a);

	bool isExist() const;
	inline int getNumVertexs() const {
		return this->size();
	}
	inline int getNumSegments() const {
		return this->size();
	}
	POINT getVertex(int i) const {
		return this->get(i);
	}
	typename POINT::value_type getMaxX() const;
	typename POINT::value_type getMinX() const;
	typename POINT::value_type getMaxY() const;
	typename POINT::value_type getMinY() const;
};

template<typename POINT>
typename POINT::value_type PolygonT<POINT>::getMaxX() const {
	ASSERT(this->size() > 0);
	typename POINT::value_type max = this->at(0).x;
	for (typename PolygonT<POINT>::iterator iter = this->begin();
			iter != this->end(); iter++) {
		if ((*iter).x > max) {
			max = (*iter).x;
		}
	}
	return max;
}
template<typename POINT>
typename POINT::value_type PolygonT<POINT>::getMinX() const {
	ASSERT(this->size() > 0);
	typename POINT::value_type min = this->at(0).x;
	for (typename PolygonT<POINT>::iterator iter = this->begin();
			iter != this->end(); iter++) {
		if ((*iter).x < min) {
			min = (*iter).x;
		}
	}
	return min;
}
template<typename POINT>
typename POINT::value_type PolygonT<POINT>::getMaxY() const {
	ASSERT(this->size() > 0);
	typename POINT::value_type max = this->at(0).y;
	for (typename PolygonT<POINT>::iterator iter = this->begin();
			iter != this->end(); iter++) {
		if ((*iter).y > max) {
			max = (*iter).y;
		}
	}
	return max;
}
template<typename POINT>
typename POINT::value_type PolygonT<POINT>::getMinY() const {
	ASSERT(this->size() > 0);
	typename POINT::value_type min = this->at(0).y;
	for (typename PolygonT<POINT>::iterator iter = this->begin();
			iter != this->end(); iter++) {
		if ((*iter).y < min) {
			min = (*iter).y;
		}
	}
	return min;
}

//==========================================
class PolygonPlane3D : public PolygonT<Point3D>{
public:
	PolygonPlane3D() {
	}
	PolygonPlane3D(const ListT<Point3D>&);
	//PolygonT(const Polygon &a);
	//Polygon& operator=(const Polygon &a);
};

}

#endif /* POLYGON_H_ */
