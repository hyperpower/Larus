/************************
 //  \file   test_octree.h
 //  \brief
 // 
 //  \author czhou
 //  \date   15 janv. 2015 
 ***********************/
#ifndef TEST_OCTREE_H_
#define TEST_OCTREE_H_

#include "../src/IO/IO_gnuplot.h"
#include "../src/IO/IO_vtk.h"
#include "../src/Geometry/Triangle.h"
#include "../src/Geometry/Line.h"
#include "../src/Geometry/Polygon.h"
#include "../src/Utility/Array.h"
#include "../src/Geometry/Relation.h"
#include "../src/Geometry/Plane.h"
#include "../src/Grid/SPTreeNode.h"
#include "../src/Grid/SPTree.h"
#include "../src/Grid/Cell.h"
#include "../src/Grid/Forest.h"
#include "../src/Algebra/Space.h"
#include "../src/Calculation/VOF.h"
#include "../src/Calculation/Levelset.h"
#include "../src/Calculation/Scalar.h"
#include "../src/Algebra/Arithmetic.h"
#include "../src/Algebra/Cube.h"
#include "../src/TypeDef.h"
#include "../src/Algebra/Interpolation.h"
#include "../src/Algebra/Expression.h"
#include "../src/Grid/Dimension.h"
#include "../src/IO/IO_gerris.h"
#include <list>
#include <math.h>
#include <sstream>

using namespace std;
namespace Larus {

void test_nei_3D() {
	Point3D p1(0, 0, 0);
	Point3D p2(1, 1, 1);
	Cell3D c(p1, p2);
	OcTree tree(c, 4);
	tree.CreatFullTree();
	tree.show_info();
	tree.draw_to_vtk_leaf("octree.vtk");
	Point3D pf(0.3, 0.45, 0.58);
	pOCNode ori = tree.getpNode(pf);
	cout << ori->getIdx() << endl;
	//ori->draw_self("ori.vtk", 1);

	//pOCNode neimm = tree.getCorNeighbor(ori, SPD_MM_ZX);
	//pOCNode neimp = tree.getCorNeighbor(ori, SPD_MP_ZX);
	//pOCNode neipm = tree.getCorNeighbor(ori, SPD_PM_ZX);
	//pOCNode neipp = tree.getCorNeighbor(ori, SPD_PP_ZX);

	//pOCNode neippp = tree.getCorNeighbor_XYZ(ori, SPD_PMP);
	//cout<<nei->cell->volume()<<endl;
	//neimm->draw_self("neimm.vtk", 1);
	//neimp->draw_self("neimp.vtk", 1);
	//neipm->draw_self("neipm.vtk", 1);
	//neipp->draw_self("neipp.vtk", 1);
	//neippp->draw_self("neippp.vtk", 1);
}
void test_octree() {
	Point3D p1(0, 0, 0);
	Point3D p2(1, 1, 1);
	Cell3D c(p1, p2);
	OCNode node(NULL_PTR, 1, 2, 0, c);
	node.creatChild_all(4);
	cout << "idx                " << node.getEIdx() << endl;
	cout << "has child          " << node.hasChild() << endl;
	cout << "Height             " << node.Height() << endl;
	cout << "isadj              " << node.isAdjacent(SPB_JP) << endl;
	//node.draw_self("node.vtk", 1);

	OcTree oct(c, 3);
	oct.CreatFullTree();
	oct.show_info();
	oct.draw_to_vtk_leaf("oct.vtk");

	//Space====================
	SpaceT<Float, 3> sp(2, 5, 4);
	sp(1, 2, 3) = 100;
	sp.set(150, 1, 2, 2);
	cout << sp(1, 2, 2) << endl;
	//Float* f=sp.getpValue(1,2,2);
}

void test_forest() {
	Forest<Dimension_3D> forest(5, 6, 0.0, 0.0, 1.0, 3, 7, 0.0);
	//Forest<OcTree,3> forest(5, 6, 0.0, 0.0, 1.0, 3, 7, 0.0);
	cout << forest.get_dim() << endl;
	cout << forest.size() << endl;
	cout << forest.getNumEnableTree() << endl;
	forest.show_info();
	////forest.getpTree(2,2,0)->getpRootCell()->show();
	//forest.draw_to_vtk("forest.vtk");
}

void test_io_gerris() {
	//Point3D p(-0.5, -0.5, -0.5);
	//GerrisFileAdp<Dimension_3D> gfa("re-3-7.20.txt", 4, 7, p , 1);
	//gfa.show_file_info();
	//gfa.forest.show_info();
	////gfa.forest.getpTree(0,0,0)->show_info();
	//gfa.forest.check_leaf_with_data();
	//Point3D pp(0.1, 0.25, 0.34);
	////gfa.forest.getpTree(3,0,0)->draw_to_vtk_scalars("t00.vtk","vof_data",5);
	//gfa.forest.getpTree(0,0,0)->Find(pp)->data->aCenterData.show();
	//gfa.forest.draw_to_vtk_leaf_vectors("forest.vtk","v_data",2,3,4);

	Point3D p1(0, 0, 0);
	Point3D p2(1, 1, 1);
	Cell3D c(p1, p2);
	Point3D p(0.75, 0.2, 0.3);
	Point3D pp(0.1, 0.75, 0.34);
	OcTree oct(c, 5);
	arrayList a(3);
	a.assign(2);
	arrayList b(3);
	b.assign(3);
	pCellData3D pcd = new CellData3D();
	pcd->aCenterData = a;
	pCellData3D pcd2 = new CellData3D();
	pcd2->aCenterData = b;
	//_insert_partial_with_data(Dimension_3D(), &oct, p, 2, pcd);
	//_insert_partial_with_data(Dimension_3D(), &oct, pp, 2, pcd2);
	oct.show_info();
	//oct.draw_to_vtk_all("tree.vtk");
	//oct.draw_to_vtk_vectors("tree.vtk", "tdata", 0, 1, 2);
}

void test_vof() {
	//Cell 3D
	Float x1 = 1;
	Float x2 = 2;
	Float x3 = 3;
	Float m1 = -0.31;
	Float m2 = -0.5;
	Float m3 = -0.1;
	Point3D p1(0, 0, 0);
	Point3D p2(x1, x2, x3);
	Cell3D c(p1, p2);
	c.draw_to_vtk("cell.vtk");
	Point3D v1(MIN(0.0, m1), MIN(0.0, m2), MIN(0.0, m3));
	Point3D v2(MAX(0.0, m1), MAX(0.0, m2), MAX(0.0, m3));
	Cell3D v(v1, v2);
	v.draw_to_vtk("cellv.vtk");

	int n = 100;
	Float maxalpha = ABS(m1) * x1 + ABS(m2) * x2 + ABS(m3) * x3;
	for (int i = 0; i < n - 2; i++) {
		Plane p(m1, m2, m3, (i + 1) * (-maxalpha / n));
		cout << "test case ================== " << i << endl;
		cout << "Volume " << calFractionVolume(p, x1, x2, x3) << endl;
		if (i == 34) {
			cout << "stop" << endl;
		}
		p.show();
		ListT<Point3D> listp = calInterctPoints(p, x1, x2, x3);
		ostringstream sstream;
		sstream.precision(0);
		sstream << i;
		string filename = "planes/plane-" + sstream.str() + ".vtk";
		if (listp.size() > 0) {
			drawtofile_vtk(filename, listp);
		}
	}
}

void test_equation() {
	Float a = 1;
	Float b = 0;
	Float c = -27;
	Float d = 54;
	Float res1 = 0;
	Float res2 = 0;
	Float res3 = 0;
	int numr = solve_cubic_equation(a, b, c, d, res1, res2, res3);
	cout << "numr   " << numr << endl;
	cout << "x1     " << res1 << endl;
	cout << "x2     " << res2 << endl;
	cout << "x3     " << res3 << endl;

	Float x1 = 1;
	Float x2 = 1;
	//Float x3 = 1;
	Point2D p1(0, 0);
	Point2D p2(x1, x2);
	Float q1 = 0;
	Float q2 = 1;
	Float q3 = 1;
	Float q4 = 1;
	Point2D p(0.5, 0.5);
	cout << space_linear_interpolation(p, p1, p2, q1, q2, q3, q4) << endl;
}

void test_scalar() {
	Point3D p1(0, 0, 0);
	Point3D p2(1, 1, 1);
	Point3D pc(0.5, 0.5, 0.5);
	Cell3D c(p1, p2);
	OcTree oct(c, 5);
	oct.CreatFullTree();
	oct.show_info();
	//Float a = 3.0;
	new_array_on_center_leaf(&oct, 3);
	//set_const_on_center_leaf(&oct, a, 0);
	set_sphere_ls(&oct, pc, 0.1, 0);
	draw_vtu_point_leaf_scalars("tree.vtu", &oct, "levelset", 0);
}

void test_refine() {
	Point3D p1(0, 0, 0);
	Point3D p2(1, 1, 1);
	Point3D pc(0.43, 0.41, 0.41);
	Cell3D c(p1, p2);
	OcTree oct(c, 5);
	set_sphere_ls_initial_refine(&oct, pc, 0.1, 0);
	new_array_on_center_leaf(&oct, 3);
	set_sphere_ls(&oct, pc, 0.1, 0);
	oct.show_info();
	oct.draw_to_vtk_leaf("tree.vtk");
	draw_vtu_point_leaf_scalars("tree.vtu", &oct, "levelset", 0);
}

void test_io_gerris_vof() {
	Point3D p(-1, -1, -1);
	GerrisFileAdp<Dimension_3D> gfa("re-3-10.60.txt", 5, 8, p, 2);
	gfa.show_file_info();
	gfa.forest.show_info();
	//gfa.forest.getpTree(0,0,0)->show_info();
	gfa.forest.check_leaf_with_data();
	//gfa.forest.draw_to_vtk_leaf_vectors("forest.vtk","v_data",2,3,4);
	//gfa.forest.draw_to_vtk_leaf_scalars("forest.vtk", "vof_data", 5);
	ListT<pOCNode> listvof;
	gfa.forest.getpTree(2, 0, 0)->show_info();
	//gfa.forest.getpTree(2,0,0)->toList_leaf(listvof);
	//===
	//int ac=0;
	//for(ListT<pOCNode>::iterator iter=listvof.begin(); iter!=listvof.end(); iter++){
	//	if(isInRange(0.0,(*iter)->data->aCenterData[5],1.0,Range_oo)){
	//		ac++;
	//	}
	//}
	//cout<<ac<<endl;
	//===
	getListpNode_leaf_center_data_in_range(	//
			listvof, //as output
			gfa.forest, // ptree
			5, //data index
			0.0, //min value
			1.0, //max value
			Range_oo);
	cout << listvof.size() << endl;
	draw_vtu_vof("vof.vtu", gfa.forest, 5, 6, 7, 8, 9);
	//Point3D pp(0.1, 0.25, 0.34);
	//gfa.forest.getpTree(3,0,0)->draw_to_vtk_scalars("t00.vtk","vof_data",5);
	//gfa.forest.getpTree(0,0,0)->Find(pp)->data->aCenterData.show();
	//gfa.forest.draw_to_vtk_leaf_vectors("forest.vtk","v_data",2,3,4);
}
//=======================================================
//Process Taylor 3D =====================================

string getfilename(int group, Float i) {
	stringstream ss;
	//ss << group << "-" << i << "-resEnd.txt";
	ss << "re-" << group;
	ss.precision(2);
	ss.setf(std::ios::fixed, std::ios::floatfield);
	ss << "-" << i << ".txt";
	string str = ss.str();
	return str;
}



void process_output_vof(GerrisFileAdp<Dimension_3D>& gfa, arrayList& arg) {
	//Float group = arg[0];
	Float i = arg[1];
	stringstream ss;
	//ss << group << "-" << i << "-resEnd.txt";
	ss << "/home/czhou/Gerris/Taylor-s/5-3d/vof/";
	ss << "vof-";
	ss.precision(2);
	ss.setf(std::ios::fixed, std::ios::floatfield);
	ss << i << ".vtu";
	string str = ss.str();
	draw_vtu_vof(str, gfa.forest, 5, 6, 7, 8, 9);
}

void test_Main_octree() {
	string filedic = "/home/czhou/Gerris/Taylor-s/5-3d/result/";
	ListT<Float> listfilm;
	Float group = 3;
	for (Float i = 5; i <= 17.6; i = i + 0.2) {
		string filename = getfilename(group, i);
		string strss = filedic + filename;
		//name============================
		Point3D p(-1, -1, -1);
		GerrisFileAdp<Dimension_3D> gfa(strss, 5, 8, p, 2);
		//gfa.show_file_info();
		//cout << "=========================\n";
		//gfa.forest.show_info();
		//cout << "End read=========================\n";
		//prepare arg
		arrayList arg(4);
		arg[0] = group;
		arg[1] = i;

		process_output_vof(gfa, arg);
		cout<<i<<endl;
		//listfilm.push_back(process_filmthickness(gfa, arg));
	}

}

}

#endif /* TEST_OCTREE_H_ */
