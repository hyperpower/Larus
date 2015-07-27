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

int fun_bc(const BoundaryCondition<Dimension_2D>& bc, pQTNode pn) {
	if (bc.direction == SPD_IM) {
		pn->data->aCenterData[bc.value_idx] = 1;
	}
	if (bc.direction == SPD_IP) {
		pn->data->aCenterData[bc.value_idx] = 0;
	}
	if (bc.direction == SPD_JM) {
		pn->data->aCenterData[bc.value_idx] = 0;
	}
	if (bc.direction == SPD_JP) {
		pn->data->aCenterData[bc.value_idx] = 0;
	}
	return 1;
}

int fun_bc_advance(const BoundaryCondition<Dimension_2D>& bc, pQTNode pn) {
	if (bc.direction == SPD_IM) {
		pn->data->aCenterData[bc.value_idx] = 0;
	}
	if (bc.direction == SPD_IP) {
		pn->data->aCenterData[bc.value_idx] = 0;
	}
	if (bc.direction == SPD_JM) {
		pn->data->aCenterData[bc.value_idx] = 0;
	}
	if (bc.direction == SPD_JP) {
		pn->data->aCenterData[bc.value_idx] = 0;
	}
	return 1;
}

Float set_circle(Float x, Float y, Float z) {
	Float cx = -0.3;
	Float cy = -0.3;
	Float r = 0.1;
	if ((x - cx) * (x - cx) + (y - cy) * (y - cy) < r * r) {
		return 1;
	} else {
		return 0;
	}
}

Float set_half1half0(Float x, Float y, Float z) {
	if (x < y) {
		return 1;
	} else {
		return 0;
	}
}

void test_advection_uni() {
	int L = 1;
	int N = 1;
	int max_level = 6;
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
	resize_array_on_center_leaf(forest, 7);
	set_index_on_center_leaf(forest, Idx_IDX);

	arrayList_st arridx(6);  //data index
	arridx.assign_forward(1, 1);
	arrayList arrval(6);     //data plus
	arrval[0] = 0;    // phi_idx;
	arrval[1] = 0;    // phi_idx;
	arrval[2] = 0.005; // dt_idx;
	arrval[3] = 1;   // u_idx;
	arrval[4] = 1;   // v_idx;
	arrval[5] = 0;   // w_idx;
	plus_scalar_on_leaf( // 2D Forest
			forest,      // pQuadTree
			arridx,      // data index
			arrval       // data plus
			);

	BCManager<Dimension_2D> bcm(&forest);
	bcm.new_ghost_nodes();

	Advection_Eq<Dimension_2D> ae(&forest, &bcm, 1, 2, 3, 4, 5);
	set_scalar_on_leaf_by_function( // 2D Forest
			forest,        //pQuadTree
			ae.phi_idx,    //data index
			set_half1half0 //data plus
			);
	gnuplot_show_as_surface((*ae.pforest), ae.phi_idx);

	for (int i = 0; i < bcm.pforest->size(); ++i) {
		pQuadTree pt = bcm.pforest->getpTree_1d(i);
		if (pt != NULL_PTR) {
			for (int ii = 4; ii <= 7; ii++) {
				if (pt->getNeighborpTree(toDirection(ii)) == NULL) {
					BoundaryCondition<Dimension_2D> bc;
					bc.direction = toDirection(ii);
					bc.tree_idx = i;
					bc.value_idx = ae.phi_idx;
					bc.pfun = fun_bc;
					bcm.add_BC(bc);
				}
			}
		}
	}
	bcm.set_bc();
	//ae.advance(30);
	cout<< "slover  ="<<ae.slove(1e-6)<<endl;
	//gnuplot_show_as_contour(ae);
	int n = 1000;
	arrayList lx(n);
	arrayList lv(n);
	Float y = 0;
	Float xm = -0.5;
	Float xp = 0.5;
	for (int i = 0; i < n; i++) {
		arrayList_st arridx(1);  //idx
		arrayList arrres(1);
		arridx[0] = ae.phi_idx;
		Float dx = (xp - xm) / n;
		Point2D p(xm + i * dx, y);
		interpolate((*ae.pforest), p, arridx, arrres);
		lx[i] = p.x;
		lv[i] = arrres[0];
	}
	Gnuplot gp("lines");
	gp.plot_2(lx,lv," with lines lw 2");

	gnuplot_show_as_surface((*ae.pforest), ae.phi_idx);

	cout << "End of test =========\n";
}

void test_advection_uni_advance() {
	int L = 1;
	int N = 1;
	int max_level = 6;
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
	resize_array_on_center_leaf(forest, 7);
	set_index_on_center_leaf(forest, Idx_IDX);

	arrayList_st arridx(6);  //data index
	arridx.assign_forward(1, 1);
	arrayList arrval(6);     //data plus
	arrval[0] = 0;    // phi_idx;
	arrval[1] = 0;    // phi_idx;
	arrval[2] = 0.005; // dt_idx;
	arrval[3] = 1;   // u_idx;
	arrval[4] = 1;   // v_idx;
	arrval[5] = 0;   // w_idx;
	plus_scalar_on_leaf( // 2D Forest
			forest,      // pQuadTree
			arridx,      // data index
			arrval       // data plus
			);

	BCManager<Dimension_2D> bcm(&forest);
	bcm.new_ghost_nodes();

	Advection_Eq<Dimension_2D> ae(&forest, &bcm, 1, 2, 3, 4, 5);
	set_scalar_on_leaf_by_function( // 2D Forest
			forest, //pQuadTree
			ae.phi_idx, //data index
			set_circle //data plus
			);
	//gnuplot_show_as_contour(ae);
	gnuplot_show_as_surface((*ae.pforest), ae.phi_idx);

	for (int i = 0; i < bcm.pforest->size(); ++i) {
		pQuadTree pt = bcm.pforest->getpTree_1d(i);
		if (pt != NULL_PTR) {
			for (int ii = 4; ii <= 7; ii++) {
				if (pt->getNeighborpTree(toDirection(ii)) == NULL) {
					BoundaryCondition<Dimension_2D> bc;
					bc.direction = toDirection(ii);
					bc.tree_idx = i;
					bc.value_idx = ae.phi_idx;
					bc.pfun = fun_bc_advance;
					bcm.add_BC(bc);
				}
			}
		}
	}
	bcm.set_bc();
	ae.advance(30);
	//ae.slove(1e-6);
	//gnuplot_show_as_contour(ae);
	gnuplot_show_as_surface((*ae.pforest), ae.phi_idx);

	cout << "End of test =========\n";
}

}

#endif
