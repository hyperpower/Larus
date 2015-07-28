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
#include <math.h>

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
int fun_bc2(const BoundaryCondition<Dimension_2D>& bc, pQTNode pn) {
	if (bc.direction == SPD_IM) {
		Float y = pn->cell->get(CSAxis_Y, eCPL_M);
		Float val = 0.3414;
		if (y <= val) {
			refcVal(pn, bc.value_idx) = sin(
					LarusDef::PI / 2
							* MAX(1 - (ABS(y - val / 2)) / (val / 2), 0.0));
		} else {
			refcVal(pn, bc.value_idx) = 0;
		}
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
int fun_bc3(const BoundaryCondition<Dimension_2D>& bc, pQTNode pn) {
	if (bc.direction == SPD_IM) {
		Float y = pn->cell->get(CSAxis_Y, eCPL_M);
		Float val = 0.3;
		if (y <= val) {
			refcVal(pn, bc.value_idx) = 1;
		} else {
			refcVal(pn, bc.value_idx) = 0;
		}
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

int fun_bc4(const BoundaryCondition<Dimension_2D>& bc, pQTNode pn) {
	if (bc.direction == SPD_IM) {
		refcVal(pn, bc.value_idx) = 0;
	}
	if (bc.direction == SPD_IP) {
		pn->data->aCenterData[bc.value_idx] = 0;
	}
	if (bc.direction == SPD_JM) {
		Float x = pn->cell->get(CSAxis_X, eCPL_M);
		Float val = -0.5;
		if (val < x && x < 0.0) {
			refcVal(pn, bc.value_idx) = 1;
		} else {
			refcVal(pn, bc.value_idx) = 0;
		}

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

Float set_case4u(Float x, Float y, Float z) {
	return 2.0*y*(1.0-x*x);
}

Float set_case4v(Float x, Float y, Float z) {
	return -2.0*x*(1.0-y*y);
}
// the four test base on paper
// M.S. Darwish, F. Moukalled , TVD schemes for unstructured grids
void test_advection_uni(int si) {
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
	resize_array_on_center_leaf(forest, 8);
	set_index_on_center_leaf(forest, Idx_IDX);

	arrayList_st arridx(7);  //data index
	arridx.assign_forward(1, 1);
	arrayList arrval(7);     //data plus
	arrval[0] = 0;    // phi_idx;
	arrval[1] = 0;    // phi_idx;
	arrval[2] = 0;    // phi_idx;
	arrval[3] = 0.005; // dt_idx;
	arrval[4] = 1;   // u_idx;
	arrval[5] = 1;   // v_idx;
	arrval[6] = 0;   // w_idx;
	plus_scalar_on_leaf( // 2D Forest
			forest,      // pQuadTree
			arridx,      // data index
			arrval       // data plus
			);

	BCManager<Dimension_2D> bcm(&forest);
	bcm.new_ghost_nodes();

	Advection_Eq<Dimension_2D> ae(&forest, &bcm, si, 1, 2, 3, 4, 5, 6);
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
	cout << "slover  =" << ae.slove(1e-6) << endl;
//gnuplot_show_as_contour(ae);
	ListT<pQTNode> listnode;     //as output
	getListpNode_leaf_on_line(listnode,           //as output
			*(ae.pforest),                  // Forest
			0.0, CSAxis_Y);
	ListT<Float> lx;
	ListT<Float> lv;
	for (ListT<pQTNode>::const_iterator iter = listnode.begin();
			iter != listnode.end(); iter++) {
		lx.push_back((*iter)->cell->get(CSAxis_X, eCPL_C));
		lv.push_back(getcVal((*iter), ae.phi_idx));
	}
//Gnuplot gp("lines");
//gp.plot_2(lx, lv, " with lines lw 2");

	gnuplot_show_as_surface((*ae.pforest), ae.phi_idx);

//cout << "End of test =========\n";
}

void test_advection_uni2(int si) {
	int L = 1;
	int N = 1;
	int max_level = 6;
	Forest2D forest(L, N, 0.0, 0.0, 1, max_level);
	for (int i = 0; i < forest.size(); i++) {
		pQuadTree tree = forest.getpTree_1d(i);
		if (tree != NULL_PTR) {
			tree->CreatFullTree();
		}
	}
	forest.ConnectTrees();
	forest.show_info();
//forest.draw_to_gnuplot_leaf("tree.txt");
	resize_array_on_center_leaf(forest, 8);
	set_index_on_center_leaf(forest, Idx_IDX);

	arrayList_st arridx(7);  //data index
	arridx.assign_forward(1, 1);
	arrayList arrval(7);     //data plus
	arrval[0] = 0;    // phi_idx;
	arrval[1] = 0;    // phi_idx;
	arrval[2] = 0;    // phi_idx;
	arrval[3] = 0.005; // dt_idx;
	arrval[4] = 1;   // u_idx;
	arrval[5] = 1;   // v_idx;
	arrval[6] = 0;   // w_idx;
	plus_scalar_on_leaf( // 2D Forest
			forest,      // pQuadTree
			arridx,      // data index
			arrval       // data plus
			);

	BCManager<Dimension_2D> bcm(&forest);
	bcm.new_ghost_nodes();

	Advection_Eq<Dimension_2D> ae(&forest, &bcm, si, 1, 2, 3, 4, 5, 6);
	set_scalar_on_leaf_by_function( // 2D Forest
			forest,        //pQuadTree
			ae.phi_idx,    //data index
			set_circle //data plus
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
					bc.pfun = fun_bc2;
					bcm.add_BC(bc);
				}
			}
		}
	}
	bcm.set_bc();
//ae.advance(30);
	cout << "slover  =" << ae.slove(1e-6) << endl;
//gnuplot_show_as_contour(ae);
	ListT<pQTNode> listnode;     //as output
	getListpNode_leaf_on_line(listnode,           //as output
			*(ae.pforest),                  // Forest
			0.0, CSAxis_Y);
	ListT<Float> lx;
	ListT<Float> lv;
	for (ListT<pQTNode>::const_iterator iter = listnode.begin();
			iter != listnode.end(); iter++) {
		lx.push_back((*iter)->cell->get(CSAxis_X, eCPL_C));
		lv.push_back(getcVal((*iter), ae.phi_idx));
	}
//Gnuplot gp("lines");
//gp.plot_2(lx, lv, " with lines lw 2");

	gnuplot_show_as_surface((*ae.pforest), ae.phi_idx);

//cout << "End of test =========\n";
}

void test_advection_uni3(int si) {
	int L = 1;
	int N = 1;
	int max_level = 6;
	Forest2D forest(L, N, 0.0, 0.0, 1, max_level);
	for (int i = 0; i < forest.size(); i++) {
		pQuadTree tree = forest.getpTree_1d(i);
		if (tree != NULL_PTR) {
			tree->CreatFullTree();
		}
	}
	forest.ConnectTrees();
	forest.show_info();
//forest.draw_to_gnuplot_leaf("tree.txt");
	resize_array_on_center_leaf(forest, 8);
	set_index_on_center_leaf(forest, Idx_IDX);

	arrayList_st arridx(7);  //data index
	arridx.assign_forward(1, 1);
	arrayList arrval(7);     //data plus
	arrval[0] = 0;    // phi_idx;
	arrval[1] = 0;    // phi_idx;
	arrval[2] = 0;    // phi_idx;
	arrval[3] = 0.005; // dt_idx;
	arrval[4] = 1;   // u_idx;
	arrval[5] = 1;   // v_idx;
	arrval[6] = 0;   // w_idx;
	plus_scalar_on_leaf( // 2D Forest
			forest,      // pQuadTree
			arridx,      // data index
			arrval       // data plus
			);

	BCManager<Dimension_2D> bcm(&forest);
	bcm.new_ghost_nodes();

	Advection_Eq<Dimension_2D> ae(&forest, &bcm, si, 1, 2, 3, 4, 5, 6);
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
					bc.pfun = fun_bc3;
					bcm.add_BC(bc);
				}
			}
		}
	}
	bcm.set_bc();
//ae.advance(30);
	cout << "slover  =" << ae.slove(1e-6) << endl;
//gnuplot_show_as_contour(ae);
	ListT<pQTNode> listnode;     //as output
	getListpNode_leaf_on_line(listnode,           //as output
			*(ae.pforest),                  // Forest
			0.0, CSAxis_Y);
	ListT<Float> lx;
	ListT<Float> lv;
	for (ListT<pQTNode>::const_iterator iter = listnode.begin();
			iter != listnode.end(); iter++) {
		lx.push_back((*iter)->cell->get(CSAxis_X, eCPL_C));
		lv.push_back(getcVal((*iter), ae.phi_idx));
	}
//Gnuplot gp("lines");
//gp.plot_2(lx, lv, " with lines lw 2");

	gnuplot_show_as_surface((*ae.pforest), ae.phi_idx);

//cout << "End of test =========\n";
}

void test_advection_uni4(int si) {
	int L = 2;
	int N = 1;
	int max_level = 6;
	Forest2D forest(L, N, -1.0, 0.0, 1, max_level);
	for (int i = 0; i < forest.size(); i++) {
		pQuadTree tree = forest.getpTree_1d(i);
		if (tree != NULL_PTR) {
			tree->CreatFullTree();
		}
	}
	forest.ConnectTrees();
	forest.show_info();
//forest.draw_to_gnuplot_leaf("tree.txt");
	resize_array_on_center_leaf(forest, 8);
	set_index_on_center_leaf(forest, Idx_IDX);

	arrayList_st arridx(7);  //data index
	arridx.assign_forward(1, 1);
	arrayList arrval(7);     //data plus
	arrval[0] = 0;    // phi_idx;
	arrval[1] = 0;    // phi_idx;
	arrval[2] = 0;    // phi_idx;
	arrval[3] = 0.005; // dt_idx;
	arrval[4] = 1;   // u_idx;
	arrval[5] = 1;   // v_idx;
	arrval[6] = 0;   // w_idx;
	plus_scalar_on_leaf( // 2D Forest
			forest,      // pQuadTree
			arridx,      // data index
			arrval       // data plus
			);

	BCManager<Dimension_2D> bcm(&forest);
	bcm.new_ghost_nodes();

	Advection_Eq<Dimension_2D> ae(&forest, &bcm, si, 1, 2, 3, 4, 5, 6);
	set_scalar_on_leaf_by_function( // 2D Forest
			forest,        //pQuadTree
			ae.u_idx,    //data index
			set_case4u //data plus
			);
	set_scalar_on_leaf_by_function( // 2D Forest
			forest,        //pQuadTree
			ae.v_idx,    //data index
			set_case4v //data plus
			);
	//gnuplot_show_as_surface((*ae.pforest), ae.phi_idx);
	gnuplot_show_as_surface((*ae.pforest), ae.v_idx);
	for (int i = 0; i < bcm.pforest->size(); ++i) {
		pQuadTree pt = bcm.pforest->getpTree_1d(i);
		if (pt != NULL_PTR) {
			for (int ii = 4; ii <= 7; ii++) {
				if (pt->getNeighborpTree(toDirection(ii)) == NULL) {
					BoundaryCondition<Dimension_2D> bc;
					bc.direction = toDirection(ii);
					bc.tree_idx = i;
					bc.value_idx = ae.phi_idx;
					bc.pfun = fun_bc4;
					bcm.add_BC(bc);
				}
			}
		}
	}
	bcm.set_bc();
//ae.advance(30);
	cout << "slover  =" << ae.slove(1e-6) << endl;
//gnuplot_show_as_contour(ae);
	ListT<pQTNode> listnode;     //as output
	getListpNode_leaf_on_line(listnode,           //as output
			*(ae.pforest),                  // Forest
			0.0, CSAxis_Y);
	ListT<Float> lx;
	ListT<Float> lv;
	for (ListT<pQTNode>::const_iterator iter = listnode.begin();
			iter != listnode.end(); iter++) {
		lx.push_back((*iter)->cell->get(CSAxis_X, eCPL_C));
		lv.push_back(getcVal((*iter), ae.phi_idx));
	}
//Gnuplot gp("lines");
//gp.plot_2(lx, lv, " with lines lw 2");

	gnuplot_show_as_surface((*ae.pforest), ae.phi_idx);

//cout << "End of test =========\n";
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
	resize_array_on_center_leaf(forest, 8);
	set_index_on_center_leaf(forest, Idx_IDX);

	arrayList_st arridx(7);  //data index
	arridx.assign_forward(1, 1);
	arrayList arrval(7);     //data plus
	arrval[0] = 0;    // phi_idx;
	arrval[1] = 0;    // phi_idx;
	arrval[2] = 0;    // phi_idx;
	arrval[3] = 0.005; // dt_idx;
	arrval[4] = 1;   // u_idx;
	arrval[5] = 1;   // v_idx;
	arrval[6] = 0;   // w_idx;
	plus_scalar_on_leaf( // 2D Forest
			forest,      // pQuadTree
			arridx,      // data index
			arrval       // data plus
			);

	BCManager<Dimension_2D> bcm(&forest);
	bcm.new_ghost_nodes();

	Advection_Eq<Dimension_2D> ae(&forest, &bcm, 1, 1, 2, 3, 4, 5, 6);
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
