/*
 * test_triangle.h
 *
 *  Created on: Dec 25, 2014
 *      Author: zhou
 */

#ifndef TEST_TRIANGLE_H_
#define TEST_TRIANGLE_H_

#include "../src/IO/IO_gnuplot.h"
#include "../src/IO/IO_vtk.h"
#include "../src/Geometry/Triangle.h"
#include "../src/Geometry/Line.h"
#include "../src/Utility/Array.h"
#include "../src/Geometry/Relation.h"
#include <list>
#include <math.h>

using namespace std;
namespace Larus {
void test_triangle() {
	cout << "========test triangle===========" << endl;
	//Point3D a(1, 2, 1);
	//Point3D b(1, 1, 1);
	//Point3D c(2, 1.3,1);
	Point3D a(0, 0, 0);
	Point3D b(1, 0, 0);
	Point3D c(0, 1, 0);
	Triangle3D T1(a, b, c);
	T1.show();
	drawtofile_vtk("tri1.vtk", T1);
	//Point3D a2(1.5, 1.1, 0);
	//Point3D b2(1.5, 2.1, 0);
	//Point3D c2(1.3, 1.5, 1.5);
	Point3D a2(0.5,-0.3, -1);
	Point3D b2(0.5, 0.9, -1);
	Point3D c2(0.5, 0.48, 1);
	Triangle3D T2(a2, b2, c2);
	T2.show();
	drawtofile_vtk("tri2.vtk", T2);
	//cout<<isIntersect(T1,T2)<<endl;
	Segment3D seg;
	cal_intersect(T1, T2, seg);
	seg.show();
	drawtofile_vtu("seg.vtu", seg);
}
}

#endif /* TEST_TRIANGLE_H_ */
