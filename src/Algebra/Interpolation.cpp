/************************
 //  \file   Interpolation.cpp
 //  \brief
 // 
 //  \author czhou
 //  \date   2 févr. 2015 
 ***********************/

#include "Interpolation.h"

namespace Larus {

Float bilinear_interpolation(const Point2D& p,  //point unknow
		const Point2D& p11, //p mm
		const Point2D& p22, //p pp
		const Float& q11, const Float& q12,  //the 4 value
		const Float& q21, const Float& q22) {
	ASSERT(p11 != p22);
	return bilinear_interpolation(p.x, p.y, p11.x, p11.y, p22.x, p22.y, q11,
			q12, q21, q22);
}

Float trilinear_interpolation(
//
		const Point3D& p, //
		const Point3D& p111, //
		const Point3D& p222, //
		const Float& q111, const Float& q121, //
		const Float& q211, const Float& q221, //
		const Float& q112, const Float& q122, //
		const Float& q212, const Float& q222) {
	ASSERT(p111 != p222);
	return trilinear_interpolation(p.x, p.y, p.z, p111.x, p111.y, p111.z,
			p222.x, p222.y, p222.z, q111, q121, q211, q221, q112, q122, q212,
			q222);
}

Float interpolation_Lagrange(Float x, arrayList& arrx, arrayList& arry) {
	if (arrx.size() != arry.size()) {
		std::cerr << " !> Array Length is not equal. " << arrx.size() << "!="
				<< arrx.size() << "\n";
		return -1;
	}
	Float lag = 0.0;
	for (int i = 0; i < arrx.size(); i++) {
		Float ji = 1.0;
		for (int j = 0; j < arrx.size(); j++) {
			if (i != j)
				ji = ji * ((x - arrx[j]) / (arrx[i] - arrx[j])); //Lagrange basis polynomials
		}
		lag = lag + ji * arry[i];                     //interpolation polynomial
	}
	return lag;
}

Expression linear_interpolation_expression( //
		const Float& x,  //
		const Float& x0, const Expression& y0, //
		const Float& x1, const Expression& y1 //
		) {
	Expression res = y1 - y0;
	res.times((x - x0) / (x1 - x0));
	res.plus(y0);
	return res;
}

Expression linear_gradient_expression( //
		const Float& x0, const Expression& y0, //
		const Float& x1, const Expression& y1 //
		) {
	Expression res = y1 - y0;
	res.times(1.0 / (x1 - x0));
	return res;
}

Expression linear_weight_interpolation_expression( //
		const Float& w0, const Expression& y0, //
		const Float& w1, const Expression& y1 //
		){
	Float d =((1.0 / w0) + (1.0 / w1));
	Expression res(y0);
	res.times(1.0/w0/d);
	Expression res1(y1);
	res1.times(1.0/w1/d);
	res.plus(res1);
	return res;
}

Expression linear_interpolation_expression( //
		const Point2D& p,  //
		const Point2D& p0, const Expression& y0, //
		const Point2D& p1, const Expression& y1 //
		) {
	Float x1 = calDistance(p0, p1); //the distance between p0 and p1;
	//d or d0 is largest one
	//p is on positive side of p0, p1
	Float x = point_sign_distance(p, p0, p1);
	//linear_interpolation(x, 0.0, y0, x1, y1);
	Expression res = y1 - y0;
	res.times(x / x1);
	res.plus(y0);
	return res;
}

Expression linear_interpolation_expression( //
		const Point3D& p,  //
		const Point3D& p0, const Expression& y0, //
		const Point3D& p1, const Expression& y1 //
		) {
	Float x1 = calDistance(p0, p1); //the distance between p0 and p1;
	//d or d0 is largest one
	//p is on positive side of p0, p1
	Float x = point_sign_distance(p, p0, p1);
	//linear_interpolation(x, 0.0, y0, x1, y1);
	Expression res = y1 - y0;
	res.times(x / x1);
	res.plus(y0);
	return res;
}

Expression second_order_interpolation_expression( //
		const Point2D& p,  //
		const Point2D& p1, const Expression& y1, //
		const Point2D& p2, const Expression& y2, //
		const Point2D& p3, const Expression& y3  //
		) {
	Float x2 = calDistance(p1, p2); //the distance between p0 and p1
	Float x3 = point_sign_distance(p3, p1, p2);
	Float x = point_sign_distance(p, p1, p2);
	Expression res(y1);
	res.times((x - x2) * (x - x3) / x2 / x3);
	Expression res2(y2);
	res2.times((x) * (x - x3) / (x2) / (x2 - x3));
	Expression res3(y3);
	res3.times((x) * (x - x2) / (x3) / (x3 - x2));
	res.plus(res2);
	res.plus(res3);
	return res;
}

Expression second_order_interpolation_expression( //
		const Float& x,  //
		const Float& x1, const Expression& y1, //
		const Float& x2, const Expression& y2, //
		const Float& x3, const Expression& y3  //
		) {
	Expression res(y1);
	res.times((x - x2) * (x - x3) / (x1 - x2) / (x1 - x3));
	Expression res2(y2);
	res2.times((x - x1) * (x - x3) / (x2 - x1) / (x2 - x3));
	Expression res3(y3);
	res3.times((x - x1) * (x - x2) / (x3 - x1) / (x3 - x2));
	res.plus(res2);
	res.plus(res3);
	return res;
}

Expression second_order_gradient_expression( //
		const Float& x,  //
		const Float& x1, const Expression& y1, //
		const Float& x2, const Expression& y2, //
		const Float& x3, const Expression& y3  //
		) {
	Expression res(y1);
	res.times((2.0 * x - x2 - x3) / (x1 - x2) / (x1 - x3));
	Expression res2(y2);
	res2.times((2.0 * x - x1 - x3) / (x2 - x1) / (x2 - x3));
	Expression res3(y3);
	res3.times((2.0 * x - x1 - x2) / (x3 - x1) / (x3 - x2));
	res.plus(res2);
	res.plus(res3);
	return res;
}

}

