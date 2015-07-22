/************************
 //  \file   test_expression.h
 //  \brief
 // 
 //  \author czhou
 //  \date   3 févr. 2015 
 ***********************/
#ifndef TEST_EXPRESSION_H_
#define TEST_EXPRESSION_H_

#include "../src/IO/IO_gnuplot.h"
#include "../src/IO/IO_vtk.h"
#include "../src/Geometry/Triangle.h"
#include "../src/Geometry/Line.h"
#include "../src/Utility/Array.h"
#include "../src/Geometry/Relation.h"
#include "../src/Geometry/Plane.h"
#include "../src/Grid/SPTreeNode.h"
#include "../src/Grid/SPTree.h"
#include "../src/Grid/Cell.h"
#include "../src/Algebra/Space.h"
#include "../src/Calculation/VOF.h"
#include "../src/Algebra/Arithmetic.h"
#include "../src/Algebra/Cube.h"
#include "../src/TypeDef.h"
#include "../src/Algebra/Interpolation.h"
#include "../src/Algebra/Expression.h"

using namespace std;
namespace Larus {
void test_expression() {
	Expression exp;
	exp.Insert(ExpTerm(1, 0.3));
	//exp.Insert(ExpTerm(2, 0.3));
	//exp.Insert(ExpTerm(3, 0.4));
	//exp.Insert(ExpTerm(4, 0.5));

	Expression exp2;
	exp2.Insert(ExpTerm(5, 1.0));
	//exp2.Insert(ExpTerm(3, 0.0));
	//exp2.Insert(ExpTerm(4, 0.5));
	//exp2.Insert(ExpTerm(2, 1.0));
	Expression exp3;
	exp3.Insert(ExpTerm(7, 2.0));
	exp.show();
	cout << "=================\n";
	exp2.show();
	exp2.PrintTree();

	Point2D p(1.3, 0);
	Point2D p0(1, 0);
	Point2D p1(2, 0);
	Point2D p2(3.1, 0);
	Expression expres;
	expres = second_order_interpolation_expression(
			1,         //
			0.1, exp,  //
			0.2, exp2, //
			0.3, exp3
			);
	cout << "res================\n";
	expres.show();
	expres.PrintTree();
	cout << second_order_interpolation(p,  //
			p0, 0.3, //
			p1, 1.0,
			p2, 2.0);

	//cout<<"zero = "<<exp2.count_zero()<<endl;
	//exp2.trim_zero();
	//exp2.show();
	//exp2.PrintTree();
}
}

#endif /* TEST_EXPRESSION_H_ */
