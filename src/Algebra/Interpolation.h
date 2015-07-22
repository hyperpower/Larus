/************************
 //  \file   Interpolation.h
 //  \brief
 // 
 //  \author czhou
 //  \date   2 févr. 2015 
 ***********************/
#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

#include "../TypeDef.h"

#include "../Utility/ArrayList.h"
#include "../Geometry/Point.h"
#include "../Geometry/GeoDef.h"
#include "Arithmetic.h"
#include "Expression.h"
#include <iostream>

namespace Larus {
//Assume that p,p0 and p1 is colinear
//p0 is original point, p1 is on positive direction
//The function decide p is on positive or on negtive position
//This function will not check whether these points are on the same line
template<typename POINT>
inline Float point_sign_distance( //
		const POINT& p, //
		const POINT& p0, //
		const POINT& p1 //
		) {
	Float d =  calDistance(p0, p1); //the distance between p0 and p1
	Float d0 = calDistance(p0, p);  //the distance between p0 and p
	Float d1 = calDistance(p1, p);  //the distance between p1 and p
	int max, mid, small;
	sort(d, d0, d1, Comp_great, max, mid, small);
	//d or d0 is largest one
	//p is on positive side of p0, p1
	return (max == 0 || max == 1) ? d0 : -d0;
}

//Linear interpolation===========================
template<typename TYPE>
inline Float linear_interpolation( //
		const TYPE& x, //
		const TYPE& x0, const TYPE& y0, //
		const TYPE& x1, const TYPE& y1  //
		) {
	return (y1 - y0) / (x1 - x0) * (x - x0) + y0;
}
//interploat on center
//      d0       d1
// |---------*-------|
// v0       res      v1
template<typename TYPE>
inline Float linear_interpolation_center( //
		const TYPE& d0, const TYPE& v0, //
		const TYPE& d1, const TYPE& v1  //
		) {
	return (v1 - v0) / (d0 + d1) * d0 + v0;
}

template<typename TYPE>
inline Float linear_weight_interpolation( //
		const TYPE& w0, const TYPE& y0, //
		const TYPE& w1, const TYPE& y1  //
		) {
	return ((y0 / w0) + (y1 / w1)) / ((1.0 / w0) + (1.0 / w1));
}

//Linear gradient===========================
//x1 greater than x0
template<typename TYPE>
inline Float linear_gradient( //
		const TYPE& x0, const TYPE& y0, //
		const TYPE& x1, const TYPE& y1  //
		) {
	return (y1 - y0) / (x1 - x0);
}
//Assume that p,p0 and p1 is colinear
//This function will not check whether these points are on the same line
template<typename TYPE>
inline Float linear_interpolation( //
		const Point2D& p, //
		const Point2D& p0, const TYPE& y0, //
		const Point2D& p1, const TYPE& y1  //
		) {
	Float x1 = calDistance(p0, p1); //the distance between p0 and p1;
	//d or d0 is largest one
	//p is on positive side of p0, p1
	Float x = point_sign_distance(p, p0, p1);
	return linear_interpolation(x, 0.0, y0, x1, y1);
}
//Assume that p,p0 and p1 is colinear
//This function will not check whether these points are on the same line
template<typename TYPE>
inline Float linear_interpolation( //
		const Point3D& p, //
		const Point3D& p0, const TYPE& y0, //
		const Point3D& p1, const TYPE& y1  //
		) {
	Float x1 = calDistance(p0, p1); //the distance between p0 and p1;
	//d or d0 is largest one
	//p is on positive side of p0, p1
	Float x = point_sign_distance(p, p0, p1);
	return linear_interpolation(x, 0.0, y0, x1, y1);
}

Expression linear_interpolation_expression( //
		const Float& x,  //
		const Float& x0, const Expression& y0, //
		const Float& x1, const Expression& y1 //
		);

Expression linear_gradient_expression( //
		const Float& x0, const Expression& y0, //
		const Float& x1, const Expression& y1 //
		);

Expression linear_weight_interpolation_expression( //
		const Float& x0, const Expression& y0, //
		const Float& x1, const Expression& y1 //
		);

Expression linear_interpolation_expression( //
		const Point2D& p,  //
		const Point2D& p0, const Expression& y0, //
		const Point2D& p1, const Expression& y1 //
		);

Expression linear_interpolation_expression( //
		const Point3D& p,  //
		const Point3D& p0, const Expression& y0, //
		const Point3D& p1, const Expression& y1 //
		);

//Second order interpolation ====================
//===============================================
template<typename TYPE>
inline Float second_order_interpolation( //
		const TYPE& x, //
		const TYPE& x1, const TYPE& y1, //
		const TYPE& x2, const TYPE& y2, //
		const TYPE& x3, const TYPE& y3  //
		) {
	return y1 * (x - x2) * (x - x3) / (x1 - x2) / (x1 - x3)
			+ y2 * (x - x1) * (x - x3) / (x2 - x1) / (x2 - x3)
			+ y3 * (x - x1) * (x - x2) / (x3 - x1) / (x3 - x2);
}

template<typename TYPE>
inline Float second_order_gradient( //
		const TYPE& x, //
		const TYPE& x1, const TYPE& y1, //
		const TYPE& x2, const TYPE& y2, //
		const TYPE& x3, const TYPE& y3  //
		) {
	return y1 * (2.0 * x - x2 - x3) / (x1 - x2) / (x1 - x3)
			+ y2 * (2.0 * x - x1 - x3) / (x2 - x1) / (x2 - x3)
			+ y3 * (2.0 * x - x1 - x2) / (x3 - x1) / (x3 - x2);
}

template<typename TYPE>
inline Float second_order_interpolation( //
		const Point2D& p, //
		const Point2D& p1, const TYPE& y1, //
		const Point2D& p2, const TYPE& y2, //
		const Point2D& p3, const TYPE& y3  //
		) {
	Float d2 = calDistance(p1, p2); //the distance between p0 and p1
	Float d3 = point_sign_distance(p3, p1, p2);
	Float dp = point_sign_distance(p, p1, p2);
	return second_order_interpolation(dp, 0.0, y1, d2, y2, d3, y3);
}

Expression second_order_interpolation_expression( //
		const Point2D&,  //
		const Point2D&, const Expression&, //
		const Point2D&, const Expression&, //
		const Point2D&, const Expression&   //
		);

Expression second_order_interpolation_expression( //
		const Float& x,  //
		const Float& x1, const Expression& y1, //
		const Float& x2, const Expression& y2, //
		const Float& x3, const Expression& y3  //
		);

Expression second_order_gradient_expression( //
		const Float& x,  //
		const Float& x1, const Expression& y1, //
		const Float& x2, const Expression& y2, //
		const Float& x3, const Expression& y3  //
		);

//Bilinear interpolation=========================
template<typename TYPE>
inline Float bilinear_interpolation(
//
		const TYPE& x, const TYPE& y, //
		const TYPE& x1, const TYPE& y1, //
		const TYPE& x2, const TYPE& y2, //
		const TYPE& q11, const TYPE& q12, //
		const TYPE& q21, const TYPE& q22) {
	return (q11 * (y2 - y) * (x2 - x) + q21 * (y2 - y) * (x - x1)
			+ q12 * (y - y1) * (x2 - x) + q22 * (y - y1) * (x - x1)) / (x2 - x1)
			/ (y2 - y1);
}
Float bilinear_interpolation(const Point2D& p, const Point2D& p11,
		const Point2D& p22, const Float& q11, const Float& q12,
		const Float& q21, const Float& q22);

//Trilinear interpolation========================
template<typename TYPE>
inline Float trilinear_interpolation(
//
		const TYPE& x, const TYPE& y, const TYPE& z, //
		const TYPE& x1, const TYPE& y1, const TYPE& z1, //
		const TYPE& x2, const TYPE& y2, const TYPE& z2, //
		const TYPE& q111, const TYPE& q121, //
		const TYPE& q211, const TYPE& q221, //
		const TYPE& q112, const TYPE& q122, //
		const TYPE& q212, const TYPE& q222) {
	return (q111 * (x2 - x) * (y2 - y) * (z2 - z)
			+ q211 * (x - x1) * (y2 - y) * (z2 - z)
			+ q121 * (x2 - x) * (y - y1) * (z2 - z)
			+ q221 * (x - x1) * (y - y1) * (z2 - z)
			+ q112 * (x2 - x) * (y2 - y) * (z - z1)
			+ q212 * (x - x1) * (y2 - y) * (z - z1)
			+ q122 * (x2 - x) * (y - y1) * (z - z1)
			+ q222 * (x - x1) * (y - y1) * (z - z1)) / (x2 - x1) / (y2 - y1)
			/ (z2 - z1);
}

Float trilinear_interpolation(
//
		const Point3D& p, //
		const Point3D& p11, //
		const Point3D& p22, //
		const Float& q111, const Float& q121, //
		const Float& q211, const Float& q221, //
		const Float& q112, const Float& q122, //
		const Float& q212, const Float& q222);

//template function =============================
template<typename POINT>
inline Float _space_linear_interpolation(point_2d_tag, const POINT& p, //
		const POINT& p11, //
		const POINT& p22, //
		const Float& q111, const Float& q121, //
		const Float& q211, const Float& q221, //
		const Float& q112 = 0, const Float& q122 = 0, //
		const Float& q212 = 0, const Float& q222 = 0) {
	return bilinear_interpolation(p, //
			p11, //
			p22, //
			q111, q121, //
			q211, q221);
}

template<typename POINT>
inline Float _space_linear_interpolation(point_3d_tag, const POINT& p, //
		const POINT& p11, //
		const POINT& p22, //
		const Float& q111, const Float& q121, //
		const Float& q211, const Float& q221, //
		const Float& q112, const Float& q122, //
		const Float& q212, const Float& q222) {
	return trilinear_interpolation(p, //
			p11, //
			p22, //
			q111, q121, //
			q211, q221, //
			q112, q122, //
			q212, q222);
}

template<typename POINT>
Float space_linear_interpolation(const POINT& p, //
		const POINT& p11, //
		const POINT& p22, //
		const Float& q111, const Float& q121, //
		const Float& q211, const Float& q221, //
		const Float& q112 = 0, const Float& q122 = 0, //
		const Float& q212 = 0, const Float& q222 = 0) {
	typedef typename geometry_traits<POINT>::self_tag _sg;
	return _space_linear_interpolation(_sg(), p, //
			p11, //
			p22, //
			q111, q121, //
			q211, q221, //
			q112, q122, //
			q212, q222);
}

//polynomial interpolation ======================
//Lagrange polynomial interpolation =============
Float interpolation_Lagrange(Float x, arrayList& arrx, arrayList& arry);

}

#endif /* INTERPOLATION_H_ */
