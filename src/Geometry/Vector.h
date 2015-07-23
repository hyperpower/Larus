/*
 * Vector.h
 *
 *  Created on: Dec 24, 2014
 *      Author: zhou
 */

#ifndef GEOMETRY_VECTOR_H_
#define GEOMETRY_VECTOR_H_

#include "../TypeDef.h"
#include "Point.h"

namespace Larus {

template<typename TYPE>
class Vector_2D: public Point_2D<TYPE> {
public:
	Vector_2D() :
			Point_2D<TYPE>() {
	}
	Vector_2D(const TYPE& a, const TYPE& b) :
			Point_2D<TYPE>(a, b) {
	}

	bool isExist() const {
		return (this->x == 0 && this->y == 0) ? false : true;
	}
};

template<typename TYPE>
class Vector_3D: public Point_3D<TYPE> {
public:
	Vector_3D() :
			Point_3D<TYPE>() {
	}
	Vector_3D(const TYPE& a, const TYPE& b, const TYPE& c) :
			Point_3D<TYPE>(a, b, c) {
	}
	Vector_3D(const Point_3D<TYPE>& a, const Point_3D<TYPE>& b):
	        Point_3D<TYPE>(b.x-a.x, b.y-a.y, b.z-a.z){
	}
	bool isExist() const {
		return (this->x == 0 && this->y == 0 && this->z == 0) ? false : true;
	}
};

template<typename TYPE>
Vector_3D<TYPE> cross(const Vector_3D<TYPE>& a, const Vector_3D<TYPE>& b) {
	Float i =  a[1] * b[2] - b[1] * a[2];
	Float j = -a[0] * b[2] + b[0] * a[2];
	Float k =  a[0] * b[1] - b[0] * a[1];
	return Vector_3D<TYPE>(i, j, k);
}

//====================
typedef Vector_2D<Float> Vector2D;
typedef Vector_3D<Float> Vector3D;

inline void normalize(Vector2D& v) {
	Float l = sqrt(Float(v.x * v.x + v.y * v.y));
	v.x = v.x / l;
	v.y = v.y / l;
}

inline void normalize(Vector3D& v) {
	Float l = sqrt(Float(v.x * v.x + v.y * v.y + v.z * v.z));
	v.x = v.x / l;
	v.y = v.y / l;
	v.z = v.z / l;
}

}
#endif /* GEOMETRY_VECTOR_H_ */
