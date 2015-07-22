/************************
 //  \file   Point.cpp
 //  \brief
 // 
 //  \author zhou
 //  \date   23 mai 2014 
 ***********************/
#include "../TypeDef.h"
#include "Point.h"
#include <iostream>
#include "Predicates.h"

namespace Larus {

//===============================================
// Calculate the distance between 2 Points
// for Point2D
// @param   p1 Point1
// @param   p2 Point1
// @return     Distance
//-----------------------------------------------
Float calDistance(const Point2D &p1, const Point2D &p2) {
	return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

//===============================================
// Dot multiply (sp-op).(ep-op)
// for Point2D
// @param    p1 Point1
// @param    p2 Point1
// @return      the resualt of dot multiply
//-----------------------------------------------
Float dot(const Point2D &sp, const Point2D &ep, const Point2D &op) {
	return ((sp.x - op.x) * (ep.x - op.x) + (sp.y - op.y) * (ep.y - op.y));
}

bool isEqual(const Point2D &p1, const Point2D &p2) {
	return (isEqual(p1.x, p2.x) && isEqual(p1.y, p2.y));
}

Float cro(const Point2D &sp, const Point2D &ep, const Point2D &op) {
	return ((sp.x - op.x) * (ep.y - op.y) - (ep.x - op.x) * (sp.y - op.y));
}

Float crofast(const Point2D &sp, const Point2D &ep, const Point2D &op) {
	REAL asp[2];
	asp[0] = sp.x;
	asp[1] = sp.y;
	REAL aep[2];
	aep[0] = ep.x;
	aep[1] = ep.y;
	REAL aop[2];
	aop[0] = op.x;
	aop[1] = op.y;
	return orient2dfast(asp, aep, aop);
}

Float croexact(const Point2D &sp, const Point2D &ep, const Point2D &op) {
	REAL asp[2];
	asp[0] = sp.x;
	asp[1] = sp.y;
	REAL aep[2];
	aep[0] = ep.x;
	aep[1] = ep.y;
	REAL aop[2];
	aop[0] = op.x;
	aop[1] = op.y;
	return orient2dexact(asp, aep, aop);
}

Float croslow(const Point2D &sp, const Point2D &ep, const Point2D &op) {
	REAL asp[2];
	asp[0] = sp.x;
	asp[1] = sp.y;
	REAL aep[2];
	aep[0] = ep.x;
	aep[1] = ep.y;
	REAL aop[2];
	aop[0] = op.x;
	aop[1] = op.y;
	return orient2dslow(asp, aep, aop);
}

Float croadapt(const Point2D &sp, const Point2D &ep, const Point2D &op) {
	REAL asp[2];
	asp[0] = sp.x;
	asp[1] = sp.y;
	REAL aep[2];
	aep[0] = ep.x;
	aep[1] = ep.y;
	REAL aop[2];
	aop[0] = op.x;
	aop[1] = op.y;
	return orient2d(asp, aep, aop);
}

//==================================================================
//Point3d===========================================================

template<typename TYPE>
Point_2D<TYPE> Point_3D<TYPE>::projectXY() const {
	return Point_2D<TYPE>(this->x, this->y);
}
template<typename TYPE>
Point_2D<TYPE> Point_3D<TYPE>::projectYZ() const {
	return Point_2D<TYPE>(this->y, this->z);
}
template<typename TYPE>
Point_2D<TYPE> Point_3D<TYPE>::projectZX() const {
	return Point_2D<TYPE>(this->z, this->x);
}

Float cro(const Point3D &v1, const Point3D &v2, const Point3D &v3,
		const Point3D &v4) {
	Float a[3][3];
	for (int i = 0; i != 3; ++i) {
		a[0][i] = v1[i] - v4[i];
		a[1][i] = v2[i] - v4[i];
		a[2][i] = v3[i] - v4[i];
	}

	return a[0][0] * a[1][1] * a[2][2] + a[0][1] * a[1][2] * a[2][0]
			+ a[0][2] * a[1][0] * a[2][1] - a[0][2] * a[1][1] * a[2][0]
			- a[0][1] * a[1][0] * a[2][2] - a[0][0] * a[1][2] * a[2][1];
}

Float croadapt(const Point3D &v1, const Point3D &v2, const Point3D &v3,
		const Point3D &v4) {
	REAL a[3], b[3], c[3], d[3];
	for (int i = 0; i != 3; ++i) {
		a[i] = v1[i];
		b[i] = v2[i];
		c[i] = v3[i];
		d[i] = v4[i];
	}
	return orient3d(a, b, c, d);
}

// this is the end of Point3D
//----------------------------

//===============================================
// Calculate the distance between 2 Points
// for Point3D
// @param   p1 Point1
// @param   p2 Point1
// @return     Distance
//-----------------------------------------------
Float calDistance(const Point3D &p1, const Point3D &p2) {
	return sqrt(
			(p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y)
					+ (p1.z - p2.z) * (p1.z - p2.z));
}

Point2D calMid(const Point2D &p1, const Point2D &p2){
	return Point2D((p1.x + p2.x)*0.5, (p1.y + p2.y)*0.5);
}
Point3D calMid(const Point3D &p1, const Point3D &p2){
	return Point3D((p1.x + p2.x)*0.5, (p1.y + p2.y)*0.5, (p1.z + p2.z)*0.5);
}

Float CROSS(const Point_2D<Float> &v1, const Point_2D<Float> &v2, const Point_2D<Float> &v3){
	return croadapt(v1,v2,v3);
}

Float CROSS(const Point_3D<Float> &v1, const Point_3D<Float> &v2, const Point_3D<Float> &v3,
		const Point_3D<Float> &v4){
	return croadapt(v1,v2,v3,v4);
}

}
