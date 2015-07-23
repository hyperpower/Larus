/************************
 //  \file   Triangle.h
 //  \brief
 // 
 //  \author czhou
 //  \date   3 oct. 2014 
 ***********************/
#ifndef TRIANGLE_H_
#define TRIANGLE_H_

#include "../TypeDef.h"
#include "GeoDef.h"
#include "Point.h"
#include "Plane.h"
#include "Segment.h"

namespace Larus {

enum TRI_INTERSECT_TYPE {
	TRI_NO_INTERSECT = 1 << 0,
	TRI_INTERSECT = 1 << 1,
	TRI_GEN       = 1 << 2
};

class Triangle2D: public array_3<Point2D> {
public:
	Triangle2D() {
	}
	Triangle2D(const Point2D& a, const Point2D& b, const Point2D& c);
	//Triangle2D(const Triangle2D& a);
	//Triangle2D& operator=(const Triangle2D &a);

	inline Point2D get(int i) const;
	void show() const;
	//void scale(Float xfactor, Float yfactor);
	bool isColine() const;
	bool isExist() const;
};

inline Point2D Triangle2D::get(int i) const {
	assert(i >= 0 && i <= 2);
	return elems[i];
}

class Triangle3D: public array_3<Point3D> {
protected:
	bool isIEqual(LarusDef::size_type i) const;
public:
	Triangle3D() {
	}
	Triangle3D(const Point3D& a, const Point3D& b, const Point3D& c);
	//Triangle3D(Point3D a, Point3D b, Point3D c);
	//Vector2D(const Vector2D& a);

	inline Point3D get(idx_t i) const;
	Segment3D getSegment(idx_t i) const;
	//inline Point3D& getRef(int i);
	void show() const;
	//void scale(Float xfactor, Float yfactor, Float zfactor);

	//Float area() const;
	Plane getPlane() const;

	void Transfer(Float dx, Float dy, Float dz);

	Triangle2D projectXY() const;
	Triangle2D projectYZ() const;
	Triangle2D projectZX() const;
	Triangle2D project() const;

	bool isExist() const;
	bool isColine() const;
	bool isXEqual() const;
	bool isYEqual() const;
	bool isZEqual() const;
	int onWhichSide(const Point3D&) const;

	void reverse();

};

bool cal_intersect(const Triangle3D&, const Triangle3D&, Segment3D&);

inline Point3D Triangle3D::get(idx_t i) const {
	assert(i >= 0 && i <= 2);
	return elems[i];
}


}

#endif /* TRIANGLE_H_ */
