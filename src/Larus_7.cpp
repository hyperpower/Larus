//============================================================================
// Name        : Larus0_1.cpp
// Author      : Zhou Chengsi
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "TypeDef.h"
#include "Utility/Array.h"
#include "Utility/AVLTree.h"
#include "Utility/BigInt.h"
#include "Utility/Clock.h"
#include "Geometry/Point.h"
#include "Geometry/Line.h"
#include "Geometry/Segment.h"
#include "Geometry/Polygon.h"
#include "../test/test_triangle.h"
#include "../test/test_octree.h"
#include "../test/TestCal.h"
#include "../test/test_expression.h"


#include "../test/test_blas.h"
#include "../test/test_interpolate.h"
#include "../test/test_solver.h"
#include "../test/test_taylor2.h"
#include "../test/test_poisson.h"
#include "../test/test_advection.h"
#include "../test/test_opencl.h"
//#include "../test/test_python.h"
//#include "../test/test_opencl.h"
//#include "../test/test_ts_heap.h"
//#include "OpenCL/cl_utility.h"


using namespace std;
using namespace Larus;
using namespace cl;
//using namespace LarusTS;

int main() {
	//processMain();
	//test_quadtree();
	//test_nei_3D();
	//test_vof();
	//test_interploate_value_one_node();
	//neck_group();
	//pressure_contour_line();
	//test_gauss_e();
	//cl::work_flow_control();
	//test_sp1();
	//show_opencl_info();
	//vector_dot();
	run_test_contraction();
	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!
	return 0;
}
