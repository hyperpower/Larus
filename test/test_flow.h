#ifndef _TEST_FLOW_H_
#define _TEST_FLOW_H_

#include <stdio.h>

#include "../src/Algebra/Arithmetic.h"
#include "../src/IO/IO_gerris.h"

#include "../src/IO/IO_gnuplot.h"
#include "../src/IO/Gnuplot.h"
#include "../src/IO/IO.h"
#include "../src/Geometry/Triangle.h"

#include "../src/Algebra/Expression.h"
#include "../src/Algebra/Sort.h"
#include "../src/Grid/Forest.h"

#include "../src/Calculation/Scalar.h"
#include "../src/Calculation/VOF.h"
#include "../src/Calculation/Velocity.h"

#include "../src/Utility/Array.h"
#include "../src/Utility/ArrayList.h"

#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>

#include "test_gerris.h"
#include "test_taylor2.h"

using namespace std;

namespace Larus {

/*
 *  the function use show the difference of the mesh level
 */
inline void test_flow_mesh() {
	//
	// 1. get the gerris file
	//
	string dir = "/home/czhou/Gerris/Taylor/ffr_flow2/result1/";
	int group = 0;
	Float time = 30.00;
	GerrisFileAdp<Dimension_2D>* gfa;
	int lmin = 4;
	int lmax = 7;
	new_gerris_file(gfa, dir, group, time, lmin, lmax);
	// ------------------------------------------
	gfa->show_file_info();
	cout << "num leaf " << gfa->forest.count_leaf() << endl;
	stringstream ss;
	string fn = setfilename("_bdl_", group, time, "");
	ss << "./" << fn;
	string dfn = ss.str();
	gnuplot_datafile_different_level(dfn, gfa->forest);
	// output surface ----------------------------
	ListT<Segment2D> lsur;
	get_surface(gfa, lsur);
	ss.str("");
	fn = setfilename("_sur_", group, time, "");
	ss << "./" << fn;
	dfn = ss.str();
	drawtofile_gnuplot(dfn, lsur, 1);

	//
	// delete
	//
	if (gfa != NULL) {
		delete gfa;
	}
}

void output_basic(const string& filename, Forest2D& f, int outputmode) {
	//expand arg
	ListT<Segment2D> l_surface;
	Float head_x;
	Float tail_x;
	Float max_y;
	Float tail_xc;
	Float u;
	Float v;
	Float x;
	Float y;
	//
	//surface=========================
	getListSegment( //
			l_surface, //list Segment
			f, //forest
			Gerris_T, //data color index
			Gerris_T_x, //data index
			Gerris_T_y, //data index
			Gerris_T_alpha  //data index
			);
	//get taylor bubble head====================
	head_x = get_head_x(l_surface);

	//get taylor bubble tail====================
	tail_x = get_tail_x(l_surface);

	//get max y=================================
	max_y = get_surface_max_y(l_surface);

	//get bubble tail center ===================
	Float xr = tail_x + (head_x - tail_x) / 2.0;
	ListT<Segment2D>::const_iterator iters = l_surface.begin();
	for (ListT<Segment2D>::const_iterator iter = l_surface.begin();
			iter != l_surface.end(); ++iter) {
		Float tmpy = MIN(iter->PSY(), iter->PEY());
		Float tmpx = MIN(iter->PSX(), iter->PEX());
		if (tmpy <= MIN(iters->PSY(), iters->PEY()) && tmpx < xr) {
			iters = iter;
		}
	}
	tail_xc = MIN(iters->PSX(), iters->PEX());
	//
	ListT<pQTNode> listpn;
	getListpNode_leaf_center_data_in_range(listpn, f, Gerris_T, 0.0, 1.0,
			Range_co);
	//get bubble velocity=======================
	get_bubble_velocity(listpn, u, v);
	get_center_of_bubble(listpn, x, y);

	FILE* pf = open_file(filename, outputmode);
	fprintf(pf, "%f ", head_x);
	fprintf(pf, "%f ", tail_x);
	fprintf(pf, "%f ", tail_xc);
	fprintf(pf, "%f ", max_y);
	fprintf(pf, "%f ", u);
	fprintf(pf, "%f ", v);
	fprintf(pf, "%f ", x);
	fprintf(pf, "%f ", y);
	fprintf(pf, "\n");
	fclose(pf);
}

/*
 * test flow velcity
 */
inline void test_flow_one() {
	//
	// 1. get the gerris file
	//
	string in_dir = "/home/czhou/Copy/Paper/slug/data/mesh/";
	string out_dir = "/home/czhou/Copy/Paper/slug/data/mesh/";
	int group = 4;
	Float time = 30.00;
	GerrisFileAdp<Dimension_2D>* gfa;
	int lmin = 5;
	int lmax = 8;
	new_gerris_file(gfa, in_dir, "re-5-8.txt", lmin, lmax);
	// output gird different level----------------
	gfa->show_file_info();
	cout << "num leaf " << gfa->forest.count_leaf() << endl;
	stringstream ss;
	string fn = setfilename("_dfl_", group, time, "");
	ss << out_dir << fn;
	string dfn = ss.str();
	gnuplot_datafile_different_level(dfn, gfa->forest);
	// output surface ----------------------------
	ListT<Segment2D> lsur;
	get_surface(gfa, lsur);
	ss.str("");
	fn = setfilename("_sur_", group, time, "");
	ss << out_dir << fn;
	dfn = ss.str();
	drawtofile_gnuplot(dfn, lsur, 1);
	// output velocity on min level --------------
	arrayList_st arridx(2);
	arridx[0] = Gerris_U;
	arridx[1] = Gerris_V;
	ListT<Pair<Point2D, arrayList> > arrres;
	get_average_value_on_level(gfa->forest, lmin, arridx, arrres);
	cout << " num points :" << arrres.size() << std::endl;
	dfn = set_dir_and_filename(out_dir, "_veo_", group, time, "");
	drawtofile_gnuplot(dfn, arrres, 1);
	// output grid -------------------------------
	dfn = set_dir_and_filename(out_dir, "_gri_", group, time, "");
	drawtofile_gnuplot(dfn, gfa->forest, 1);
	// output basic ------------------------------
	dfn = set_dir_and_filename(out_dir, "_bas_", group, time, "");
	output_basic(dfn, gfa->forest, 1);
	/**
	 * delete
	 */
	if (gfa != NULL) {
		delete gfa;
	}
}

}

#endif

