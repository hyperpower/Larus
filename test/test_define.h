/*
 * test_solver.h
 *
 *  Created on: May 3, 2015
 *      Author: zhou
 */

#ifndef TEST_SOLVER_H_
#define TEST_SOLVER_H_

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
#include "../src/Algebra/Solver_matrix.h"

namespace Larus {

using namespace std;

#ifdef __APPLE__
const string workdir_out = "/Users/zhou/Workspace/Larus_7/out/";
const string workdir = "/Users/zhou/Workspace/Larus_7/";
const string fn = "re-2-13.00.txt";
#else
const string workdir = "/home/czhou/Dropbox/ProgramWorkbench/Larus9/";
const string fn = "re-2-13.00.txt";
const string workdir_out = "/home/czhou/Dropbox/ProgramWorkbench/Larus9/out/";
#endif


} //end namespace

#endif /* TEST_SOLVER_H_ */
