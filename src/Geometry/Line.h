/************************
 //  \file   Line.h
 //  \brief
 // 
 //  \author zhou
 //  \date   23 mai 2014 
 ***********************/
#ifndef LINE_H_
#define LINE_H_

#include "Point.h"

#include "../TypeDef.h"
#include "GeoDef.h"
#include "Vector.h"
#include "../Utility/Pair.h"

namespace Larus {

//==========================================

class Line: public array_3<Float> {
	//The Line function defined as Ax+By=C
public:
	Line() {
	}
	Line(const Float& a, const Float& b, const Float& c) :
			array_3<Float>(a, b, c) {
	}
	Line(Float ax, Float ay, Float bx, Float by);
	Line(const Point2D &a, const Point2D &b);
	void reconstruct(Float a, Float b, Float c);
	inline Float getA() const;
	inline Float getB() const;
	inline Float getC() const;
	Float calX(Float y) const;
	Float calY(Float x) const;
	Float getSlope() const;
	Float getInterseptX() const;
	Float getInterseptY() const;
	Float getNormX() const;
	Float getNormY() const;
	//Vector2D getNormVector() const;
	bool isExist() const;
	bool isVertical() const;
	bool isHorizontal() const;
	void show() const;
};

inline Float Line::getA() const {
	return this->elems[0];
}

inline Float Line::getB() const {
	return this->elems[1];
}

inline Float Line::getC() const {
	return this->elems[2];
}

class Line3D: public Pair<Point3D, Vector3D> {
public:
	//typedef line_3d_tag self_tag;
	typedef Line3D self_type;

	Line3D() {
	}

	Line3D(const Point3D& p, const Vector3D& v) {
		this->first = p;
		this->second = v;
	}
	Line3D(const Point3D &a, const Point3D &b) {
		this->first = a;
		this->second = Vector3D(a, b);
	}

	void show() const;
};

}

#endif /* LINE_H_ */
