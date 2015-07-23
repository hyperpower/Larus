/*
 * test_interprolate.h
 *
 *  Created on: May 1, 2015
 *      Author: zhou
 */

#ifndef TEST_INTERPOLATE_H_
#define TEST_INTERPOLATE_H_

#include <string>
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
#include "test_define.h"
#include "../src/Calculation/Exp.h"
namespace Larus {

Float set_sin(Float x, Float y) {
	//center
	Float ox = 0.5;
	Float oy = 0.5;
	Float dis = calDistance(ox, oy, x, y);
	return sin(3 * dis);
}

inline Float cal_order(Float ec, Float ef) {
	return log(ABS(ec) / ABS(ef)) / log(2);
}

//this function test the order of interpolate function
//
//
void test_interploate_2() {
	//test interpolate property in 2d
	//build a forest
	int nbe = 3;
	int nin = 7;
	arrayList ainorm1(nin);
	arrayList ainorm2(nin);
	arrayList ainormi(nin);
	for (int ii = nbe; ii < nbe + nin; ii++) {
		int maxlevel = ii;
		Forest2D forest(1, 1, 0.0, 0.0, 1, maxlevel);
		pQuadTree pt = forest.getFirstpeTree();
		pt->CreatFullTree();
		forest.ConnectTrees();
		//forest.show_info();
		//forest.draw_to_gnuplot_leaf(workdir_out + "f_grid.txt");
		//new values
		new_array_on_center_leaf(pt, 2); //2D
		FILE* file = open_file(workdir_out + "sc.txt", 1);
		for (Forest2D::iterator iter = forest.begin(); iter != forest.end();
				++iter) {
			//set values ==========================
			Float x = iter->cell->getCenterPoint().x;
			Float y = iter->cell->getCenterPoint().y;
			iter->data->aCenterData[0] = set_sin(x, y);
			//output scatter data =================
			fprintf(file, "%lf %lf ", x, y);
			fprintf(file, "%lf \n", iter->data->aCenterData[0]);
		}
		fclose(file);

		//initial a line ============================
		Point2D sp(0.1, 0.1);
		Point2D ep(0.9, 0.9);
		int n = 100;
		ListT<Point2D> lp;
		generate_ListPoint2D(sp, ep, n, lp);

		arrayList ainter(lp.size());
		arrayList aorig(lp.size());
		arrayList_st arridx(1);
		arridx[0] = 0;
		arrayList arrres(1);
		int i = 0;
		for (ListT<Point2D>::iterator iter = lp.begin(); iter != lp.end();
				++iter) {
			//iter->show();
			interpolate_1order(forest, (*iter), arridx, arrres);
			ainter[i] = arrres[0];
			aorig[i] = set_sin(iter->x, iter->y);
			i++;
		}
		arrayList arrerr = ainter - aorig;
		//output results ============================
		std::stringstream ss;
		ss << ii;
		file = open_file(workdir_out + "res-" + ss.str() + ".txt", 1);
		for (int i = 0; i < aorig.size(); ++i) {
			//output scatter data =================
			fprintf(file, "%d ", i);
			fprintf(file, "%lf %lf %lf \n", aorig[i], ainter[i], arrerr[i]);
		}
		fclose(file);
		ainorm1[ii - nbe] = nrm1(arrerr);
		ainorm2[ii - nbe] = nrm2(arrerr);
		ainormi[ii - nbe] = nrminf(arrerr);
	}
	//output all result
	cout << "error output =======================" << endl;
	cout << "  i         n1         n2         ni" << endl;
	for (int ii = 0; ii < nin; ii++) {
		cout << std::scientific;
		cout << nbe + ii << "  ";
		cout << ainorm1[ii] << "  ";
		cout << ainorm2[ii] << "  ";
		cout << ainormi[ii] << endl;
	}
	cout << "order output =======================" << endl;
	cout << "  i         n1o         n2o         nio" << endl;
	Float sumn1 = 0.0;
	Float sumn2 = 0.0;
	Float sumni = 0.0;
	for (int ii = 1; ii < nin; ii++) {
		cout << " " << nbe + ii - 1 << "-" << nbe + ii << "  ";
		cout << cal_order(ainorm1[ii - 1], ainorm1[ii]) << "  ";
		cout << cal_order(ainorm2[ii - 1], ainorm2[ii]) << "  ";
		cout << cal_order(ainormi[ii - 1], ainormi[ii]) << "  " << endl;
		sumn1 += cal_order(ainorm1[ii - 1], ainorm1[ii]);
		sumn2 += cal_order(ainorm2[ii - 1], ainorm2[ii]);
		sumni += cal_order(ainormi[ii - 1], ainormi[ii]);
	}
	cout << "      ================================" << endl;
	cout << "aver  ";
	cout << sumn1 / (nin - 1) << "  ";
	cout << sumn2 / (nin - 1) << "  ";
	cout << sumni / (nin - 1) << "  " << endl;
}

void test_interpolate_exp() {
	int maxlevel = 5;
	Forest2D forest(1, 1, 0.0, 0.0, 1, maxlevel);
	pQuadTree pt = forest.getFirstpeTree();
	pt->CreatFullTree();
	forest.ConnectTrees();
	//forest.show_info();
	//forest.draw_to_gnuplot_leaf(workdir_out + "f_grid.txt");
	//new values
	new_array_on_center_leaf(pt, 11); //2D
	set_index_on_center_leaf(forest, Idx_IDX);
	for (Forest2D::iterator iter = forest.begin(); iter != forest.end();
			++iter) {
		//set values ==========================
		Float x = iter->cell->getCenterPoint().x;
		Float y = iter->cell->getCenterPoint().y;
		iter->data->aCenterData[0] = set_sin(x, y);
	}
	show_gnuplot_as_contour(forest, 0);

	//initial a point
	Point2D sp(0.13, 0.1);
	Forest2D::iterator it_p = forest.end();
	for (Forest2D::iterator iter = forest.begin(); iter != forest.end();
			++iter) {
		if (iter->cell->isInOnCell(sp)) {
			cout << "node found";
			it_p = iter;
			break;
		}
	}
	Float x = it_p->cell->getCenterPoint(CSAxis_X);
	Float y = it_p->cell->getCenterPoint(CSAxis_Y);
	cout << " x " << x;
	cout << " y " << y << endl;
	Float dis = it_p->cell->getDx() / 3.0;
	Point2D calp(x, y + dis);
	//
	arrayList_st arridx(1);
	arridx[0] = 0;
	arrayList arrres(1);
	interpolate_1order(forest, sp, arridx, arrres);
	cout << "Exact " << set_sin(sp.x, sp.y) << endl;
	cout << "res    " << arrres[0] << endl;
	//
	//Expression exp;
	//interpolate_expression_on_axis_2f( // 2D QuadTree Node
	//		it_p.get_pointer(),  //node
	//		CSAxis_X, //axix
	//		dis,   //distance to center of pn with sign
	//		exp   //Expression
	//		);
	//exp.show();
	Matrix3x3 m;
	m.show();
}

} //end namespace

#endif /* TEST_INTERPROLATE_H_ */
