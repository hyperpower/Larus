/************************
 //  \file   test_poisson.h
 //  \brief
 //
 //  \author czhou
 //  \date   23 févr. 2015
 ***********************/

#ifndef _TEST_POISSON_H_
#define _TEST_POISSON_H_

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

namespace Larus {

int Fun_bc(const BoundaryCondition<Dimension_2D>& bc, pQTNode pn) {
	if (bc.direction == SPD_IM || bc.direction == SPD_IP) {
		pn->data->aCenterData[bc.value_idx] = 0;
	}
	if (bc.direction == SPD_JM) {
		//if(bc.tree_idx == 2){
		pn->data->aCenterData[bc.value_idx] = 0;
		//}else{
		//	pn->data->aCenterData[bc.value_idx] = 10;
		//}
	}
	if (bc.direction == SPD_JP) {
		pn->data->aCenterData[bc.value_idx] = 0;
	}
	return 1;
}

Float f_fun(Float x, Float y, Float z) {
	return 1;
	Float xc = 1;
	Float yc = 1;
	Float r = 0.5;
	if ((x - xc) * (x - xc) + (y - yc) * (y - yc) < r * r) {
		return 10;
	} else {
		return 0;
	}
}

void test_poisson_uni() {
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
	forest.draw_to_gnuplot_leaf("tree.txt");
	resize_array_on_center_leaf(forest, 4);
	set_index_on_center_leaf(forest, Idx_IDX);

	arrayList_st arridx(3);  //data index
	arridx.assign_forward(1, 1);
	arrayList arrval(3);     //data plus
	arrval[0] = 1;
	arrval[1] = 0;
	arrval[2] = 0;
	plus_scalar_on_leaf( // 2D Forest
			forest,      // pQuadTree
			arridx,      // data index
			arrval       // data plus
			);
	BCManager<Dimension_2D> bcm(&forest);
	bcm.new_ghost_nodes();

	Poisson_Eq<Dimension_2D> pe(&forest, &bcm, 1, 2, 3);

	for (int i = 0; i < bcm.pforest->size(); ++i) {
		pQuadTree pt = bcm.pforest->getpTree_1d(i);
		if (pt != NULL_PTR) {
			for (int ii = 4; ii <= 7; ii++) {
				if (pt->getNeighborpTree(toDirection(ii)) == NULL) {
					BoundaryCondition<Dimension_2D> bc;
					bc.direction = toDirection(ii);
					bc.tree_idx = i;
					bc.value_idx = pe.phi_idx;
					bc.pfun = Fun_bc;
					bcm.add_BC(bc);
				}
			}
		}
	}
	bcm.set_bc();
	//bcm.draw_ghost_node("gh.txt");
	Float tol = 1e-6;
	pe.set_f_term(f_fun);
	pe.pforest->show_as_contour(3);
	pe.solve(tol);
	pe.show();
	//bcm.show_boundary_contour(3);

	//gnuplot_show(lr);

	cout << "End of test =========\n";
}

Float pfun_circle(Float x, Float y, Float z) {
	return (x ) * (x ) + (y - 0) * (y - 0);
}

void test_poissons_adp() {
	int L = 3;
	int N = 3;
	int max_level = 5;
	Forest2D forest(L, N, 0.0, 0.0, 1, max_level);
	ListT<Point2D> lp;
	initial_adaptation_le(forest, pfun_circle, 0.09, 3, 5);
	//for (int i = 0; i < forest.size(); i++) {
	//	pQuadTree tree = forest.getpTree_1d(i);
	//	if (tree != NULL_PTR) {
	//		tree->CreatFullTree();
	//	}
	//}
	forest.ConnectTrees();
	forest.show_info();
	//forest.show();
	forest.draw_to_gnuplot_leaf("tree.txt");
	resize_array_on_center_leaf(forest, 4);
	set_index_on_center_leaf(forest, Idx_IDX);

	arrayList_st arridx(3);  //data index
	arridx.assign_forward(1, 1);
	arrayList arrval(3);     //data plus
	arrval[0] = 1;
	arrval[1] = 0;
	arrval[2] = 0;
	plus_scalar_on_leaf( // 2D Forest
			forest,      // pQuadTree
			arridx,      // data index
			arrval       // data plus
			);
	BCManager<Dimension_2D> bcm(&forest);
	bcm.new_ghost_nodes();

	Poisson_Eq<Dimension_2D> pe(&forest, &bcm, 1, 2, 3);

	for (int i = 0; i < bcm.pforest->size(); ++i) {
		pQuadTree pt = bcm.pforest->getpTree_1d(i);
		if (pt != NULL_PTR) {
			for (int ii = 4; ii <= 7; ii++) {
				if (pt->getNeighborpTree(toDirection(ii)) == NULL) {
					BoundaryCondition<Dimension_2D> bc;
					bc.direction = toDirection(ii);
					bc.tree_idx = i;
					bc.value_idx = pe.phi_idx;
					bc.pfun = Fun_bc;
					bcm.add_BC(bc);
				}
			}
		}
	}
	bcm.set_bc();
	//bcm.show_boundary_contour(2);
	bcm.draw_ghost_node("gh.txt");
	MatrixSCR<Float> mat;
	arrayListV<Float> b;
	Float tol = 1e-10;
	pe.set_f_term(f_fun);
	Float t = get_wall_time();
	pe.solve(tol);
	cout.precision(7);
	cout << "time  " << (get_wall_time() - t) << endl;
	pe.show();

	cout << "End of test =========\n";
}

Float f_fun_1(Float x, Float y, Float z) {
	Float k = 3;
	Float l = 3;
	Float pi = LarusDef::PI;
	//return 0;
	return -pi * pi * (k * k + l * l) * sin(pi * k * x) * sin(pi * l * y);
}
Float exact_1(Float x, Float y, Float z) {
	Float k = 3;
	Float l = 3;
	Float pi = LarusDef::PI;
	return sin(pi * k * x) * sin(pi * l * y);
}

int f_bc_1(const BoundaryCondition<Dimension_2D>& bc, pQTNode pn) {
	pn->data->aCenterData[bc.value_idx] = exact_1(pn->cell->getCPX(),
			pn->cell->getCPY(), 0);
	return 1;
}

void test_poissons_1(int level, Float& norm1, Float& norm2, Float& normi) {
	cout << "Poissons slover test " << level << endl;
	int L = 1;
	int N = 1;
	int max_level = level;
	Forest2D forest(L, N, -0.5, -0.5, 1, max_level);
	ListT<Point2D> lp;
	//initial_adaptation_eq(forest, pfun_circle, 1, 3, 5);
	for (int i = 0; i < forest.size(); i++) {
		pQuadTree tree = forest.getpTree_1d(i);
		if (tree != NULL_PTR) {
			tree->CreatFullTree();
		}
	}
	forest.ConnectTrees();
	//forest.show_info();
	//forest.show();
	//forest.draw_to_gnuplot_leaf("tree.txt");
	resize_array_on_center_leaf(forest, 5);
	set_index_on_center_leaf(forest, Idx_IDX);

	arrayList_st arridx(4);  //data index
	arridx.assign_forward(1, 1);
	arrayList arrval(4);     //data plus
	arrval[0] = 1;
	arrval[1] = 0;
	arrval[2] = 0;
	arrval[3] = 0;  //error
	plus_scalar_on_leaf( // 2D Forest
			forest,      // pQuadTree
			arridx,      // data index
			arrval       // data plus
			);
	BCManager<Dimension_2D> bcm(&forest);
	bcm.new_ghost_nodes();

	Poisson_Eq<Dimension_2D> pe(&forest, &bcm, 1, 2, 3);

	for (int i = 0; i < bcm.pforest->size(); ++i) {
		pQuadTree pt = bcm.pforest->getpTree_1d(i);
		if (pt != NULL_PTR) {
			for (int ii = 4; ii <= 7; ii++) {
				if (pt->getNeighborpTree(toDirection(ii)) == NULL) {
					BoundaryCondition<Dimension_2D> bc;
					bc.direction = toDirection(ii);
					bc.tree_idx = i;
					bc.value_idx = pe.phi_idx;
					bc.pfun = f_bc_1;
					bcm.add_BC(bc);
				}
			}
		}
	}
	bcm.set_bc();
	//bcm.show_boundary_contour(2);
	//bcm.draw_ghost_node("gh.txt");
	MatrixSCR<Float> mat;
	arrayListV<Float> b;
	Float tol = 1e-10;
	pe.set_f_term(f_fun_1);
	Float t = get_wall_time();
	pe.solve(tol);
	cout.precision(7);
	cout << "time  " << (get_wall_time() - t) << endl;
	// set exact error and calculate error
	norm1 = 0;
	Float ssum = 0;
	Float svol = 0;
	for (Forest2D::iterator iter = pe.pforest->begin();
			iter != pe.pforest->end(); ++iter) {
		Float x = iter->cell->getCPX();
		Float y = iter->cell->getCPY();
		iter->data->aCenterData[3] = exact_1(x, y, 0);
		// error = res - exact;
		Float& error = iter->data->aCenterData[4];
		Float& res = iter->data->aCenterData[2];
		Float& exact = iter->data->aCenterData[3];
		Float volume = iter->cell->volume();
		svol += volume;
		error = res - exact;
		norm1 += (ABS(error) * volume);
		ssum += (error * error) * volume;
		if (iter == pe.pforest->begin()) {
			normi = ABS(error);
		} else {
			if (normi < ABS(error)) {
				normi = ABS(error);
			}
		}
	}
	norm1 = norm1 / svol;
	norm2 = sqrt(ssum) / svol;
	//set_scalar_on_leaf_by_function(*pe.pforest, 3, exact_1);

	//pe.show();
	//pe.pforest->show_as_contour(3);
}

void test_poissons_1_adp(int level, Float& norm1, Float& norm2, Float& normi) {
	cout << "Poissons slover test " << level << endl;
	int L = 2;
	int N = 2;
	int max_level = level;
	Forest2D forest(L, N, -0.5, -0.5, 0.5, max_level);
	ListT<Point2D> lp;
	initial_adaptation_le(forest, pfun_circle, 0.09, level - 1, level);
	forest.ConnectTrees();
	//forest.show_info();
	//forest.show();
	//forest.draw_to_gnuplot_leaf("tree.txt");
	resize_array_on_center_leaf(forest, 5);
	set_index_on_center_leaf(forest, Idx_IDX);

	arrayList_st arridx(4);  //data index
	arridx.assign_forward(1, 1);
	arrayList arrval(4);     //data plus
	arrval[0] = 1;
	arrval[1] = 0;
	arrval[2] = 0;
	arrval[3] = 0;  //error
	plus_scalar_on_leaf( // 2D Forest
			forest,      // pQuadTree
			arridx,      // data index
			arrval       // data plus
			);
	BCManager<Dimension_2D> bcm(&forest);
	bcm.new_ghost_nodes();

	Poisson_Eq<Dimension_2D> pe(&forest, &bcm, 1, 2, 3);

	for (int i = 0; i < bcm.pforest->size(); ++i) {
		pQuadTree pt = bcm.pforest->getpTree_1d(i);
		if (pt != NULL_PTR) {
			for (int ii = 4; ii <= 7; ii++) {
				if (pt->getNeighborpTree(toDirection(ii)) == NULL) {
					BoundaryCondition<Dimension_2D> bc;
					bc.direction = toDirection(ii);
					bc.tree_idx = i;
					bc.value_idx = pe.phi_idx;
					bc.pfun = f_bc_1;
					bcm.add_BC(bc);
				}
			}
		}
	}
	bcm.set_bc();
	//bcm.show_boundary_contour(2);
	//bcm.draw_ghost_node("gh.txt");
	MatrixSCR<Float> mat;
	arrayListV<Float> b;
	Float tol = 1e-7;
	pe.set_f_term(f_fun_1);
	Float t = get_wall_time();
	pe.solve(tol);
	cout.precision(7);
	cout << "time  " << (get_wall_time() - t) << endl;
	// set exact error and calculate error
	norm1 = 0;
	Float ssum = 0;
	Float svol = 0;
	for (Forest2D::iterator iter = pe.pforest->begin();
			iter != pe.pforest->end(); ++iter) {
		Float x = iter->cell->getCPX();
		Float y = iter->cell->getCPY();
		iter->data->aCenterData[3] = exact_1(x, y, 0);
		// error = res - exact;
		Float& error = iter->data->aCenterData[4];
		Float& res = iter->data->aCenterData[2];
		Float& exact = iter->data->aCenterData[3];
		Float volume = iter->cell->volume();
		svol += volume;
		error = res - exact;
		norm1 += (ABS(error) * volume);
		ssum += (error * error) * volume;
		if (iter == pe.pforest->begin()) {
			normi = ABS(error);
		} else {
			if (normi < ABS(error)) {
				normi = ABS(error);
			}
		}
	}
	norm1 = norm1 / svol;
	norm2 = sqrt(ssum) / svol;
	//set_scalar_on_leaf_by_function(*pe.pforest, 3, exact_1);

	//pe.show();
	//pe.pforest->show_as_contour(4);
}

void test_case_sigle(){
	int level=5;
	Float norm1, norm2, normi;
	test_poissons_1_adp(level, norm1, norm2,
					normi);
	cout<<"Level = "<<level<<endl;
	cout<<"norm 1  "<<norm1<<endl;
	cout<<"norm 2  "<<norm2<<endl;
	cout<<"norm inf"<<normi<<endl;
}

void test_case_order() {
	int n  = 5;
	int bl = 3;
	arrayList_st arr_level(n);
	arrayList arr_norm1(n);
	arrayList arr_norm2(n);
	arrayList arr_normi(n);
	for (int i = 0; i < n; i++) {
		arr_level[i] = bl + i;
	}
	// put the function here
	for (int i = 0; i < n; i++) {
		test_poissons_1_adp(arr_level[i], arr_norm1[i], arr_norm2[i],
				arr_normi[i]);
	}
	//output all result
	cout << "error output =======================" << endl;
	cout << "  i         n1         n2         ni" << endl;
	for (int ii = 0; ii < n; ii++) {
		cout << std::scientific;
		cout << bl + ii << "  ";
		cout << arr_norm1[ii] << "  ";
		cout << arr_norm2[ii] << "  ";
		cout << arr_normi[ii] << endl;
	}
	cout << "order output =======================" << endl;
	cout << "  i         n1o         n2o         nio" << endl;
	Float sumn1 = 0.0;
	Float sumn2 = 0.0;
	Float sumni = 0.0;
	for (int ii = 1; ii < n; ii++) {
		cout << " " << bl + ii - 1 << "-" << bl + ii << "  ";
		cout << cal_order(arr_norm1[ii - 1], arr_norm1[ii]) << "  ";
		cout << cal_order(arr_norm2[ii - 1], arr_norm2[ii]) << "  ";
		cout << cal_order(arr_normi[ii - 1], arr_normi[ii]) << "  " << endl;
		sumn1 += cal_order(arr_norm1[ii - 1], arr_norm1[ii]);
		sumn2 += cal_order(arr_norm2[ii - 1], arr_norm2[ii]);
		sumni += cal_order(arr_normi[ii - 1], arr_normi[ii]);
	}
	cout << " ================================ " << endl;
	cout << "aver  ";
	cout << sumn1 / (n - 1) << "  ";
	cout << sumn2 / (n - 1) << "  ";
	cout << sumni / (n - 1) << "  " << endl;
}

}
#endif
