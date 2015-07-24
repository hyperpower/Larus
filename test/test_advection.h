#ifndef _TEST_ADVECTION_H_
#define _TEST_ADVECTION_H_

#include "../src/IO/IO_gnuplot.h"
#include "../src/IO/IO_vtk.h"
#include "../src/IO/mmio.h"
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
#include "../src/Calculation/Boundary.h"
#include "../src/Calculation/Poisson.h"
#include "../src/Calculation/Adaptive.h"
#include "../src/Calculation/Advection.h"

namespace Larus {

void test_advection_uni() {
	int L = 1;
	int N = 1;
	int max_level = 4;
	Forest2D forest(L, N, -0.5, -0.5, 1, max_level);
	for (int i = 0; i < forest.size(); i++) {
		pQuadTree tree = forest.getpTree_1d(i);
		if (tree != NULL_PTR) {
			tree->CreatFullTree();
		}
	}
	forest.ConnectTrees();
	forest.show_info();
	//forest.draw_to_gnuplot_leaf("tree.txt");
	resize_array_on_center_leaf(forest, 6);
	set_index_on_center_leaf(forest, Idx_IDX);

	arrayList_st arridx(5);  //data index
	arridx.assign_forward(1, 1);
	arrayList arrval(5);     //data plus
	arrval[0] = 1;
	arrval[1] = 0;
	arrval[2] = 0;
	arrval[3] = 0;
	arrval[4] = 0;
	plus_scalar_on_leaf( // 2D Forest
			forest,      // pQuadTree
			arridx,      // data index
			arrval       // data plus
			);
	BCManager<Dimension_2D> bcm(&forest);
	bcm.new_ghost_nodes();

	Advection_Eq<Dimension_2D> ae(&forest, &bcm, 1, 2, 3, 4);



	//gnuplot_show(lr);

	cout << "End of test =========\n";
}

}

#endif
