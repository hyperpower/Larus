/************************
 //  \file   ProcessTaylor.cpp
 //  \brief
 //
 //  \author czhou
 //  \date   10 déc. 2014
 ***********************/

#ifndef _TEST_TAYLOR2_H_
#define _TEST_TAYLOR2_H_
#include <stdio.h>
#include "TestCal.h"
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

using namespace std;

namespace Larus {

const string workdir_taylor =
		"/home/czhou/Gerris/Taylor_data/straight/paper-araujo/";

#ifdef __APPLE__
const string strDIR = "";
#else
const string strDIR =
"/home/czhou/Gerris/Taylor_data/Expansion/exp-re31eo100r1.4L4/";
#endif

void cal_surface(GerrisFileAdp<Dimension_2D>& gfa, ListT<Segment2D>& list) {
	getListSegment( //
			list, //list Segment
			gfa.forest, //forest
			Gerris_T, //data color index
			Gerris_T_x, //data index
			Gerris_T_y, //data index
			Gerris_T_alpha  //data index
			);
}

Float get_head_x(const ListT<Segment2D>& sg) {
	assert(sg.size() > 0);
	Float res = sg.begin()->PSX();
	for (ListT<Segment2D>::const_iterator iter = sg.begin(); iter != sg.end();
			++iter) {
		Float tmp = MAX(iter->PSX(), iter->PEX());
		if (tmp > res) {
			res = tmp;
		}
	}
	return res;
}

Float get_tail_x(const ListT<Segment2D>& sg) {
	assert(sg.size() > 0);
	Float res = sg.begin()->PSX();
	for (ListT<Segment2D>::const_iterator iter = sg.begin(); iter != sg.end();
			++iter) {
		Float tmp = MIN(iter->PSX(), iter->PEX());
		if (tmp < res) {
			res = tmp;
		}
	}
	return res;
}

Float get_tail_xc(const ListT<Segment2D>& sg) {
	assert(sg.size() > 0);
	Float head_x = get_head_x(sg);
	Float tail_x = get_tail_x(sg);
	Float xr = tail_x + (head_x - tail_x) / 2.0;
	ListT<Segment2D>::const_iterator iters = sg.begin();
	for (ListT<Segment2D>::const_iterator iter = sg.begin(); iter != sg.end();
			++iter) {
		Float tmpy = MIN(iter->PSY(), iter->PEY());
		Float tmpx = MIN(iter->PSX(), iter->PEX());
		if (tmpy <= MIN(iters->PSY(), iters->PEY()) && tmpx < xr) {
			iters = iter;
		}
	}
	return MIN(iters->PSX(), iters->PEX());
}

Float get_surface_max_y(const ListT<Segment2D>& sg) {
	assert(sg.size() > 0);
	Float res = sg.begin()->PSY();
	for (ListT<Segment2D>::const_iterator iter = sg.begin(); iter != sg.end();
			++iter) {
		Float tmp = MAX(iter->PSY(), iter->PEY());
		if (tmp > res) {
			res = tmp;
		}
	}
	return res;
}

Segment2D get_surface_max_segment(const ListT<Segment2D>& sg) {
	assert(sg.size() > 0);
	Float res = sg.begin()->PSY();
	Segment2D resseg;
	for (ListT<Segment2D>::const_iterator iter = sg.begin(); iter != sg.end();
			++iter) {
		Float tmp = MAX(iter->PSY(), iter->PEY());
		if (tmp > res) {
			res = tmp;
			resseg = (*iter);
		}
	}
	return resseg;
}

void get_bubble_velocity(const ListT<pQTNode>& lnode, Float& u, Float& v) {
	Float area = 0;
	Float sumveou = 0;
	Float sumveov = 0;
	for (ListT<pQTNode>::const_iterator pn = lnode.begin(); pn != lnode.end();
			++pn) {
		sumveou += (*pn)->cell->area() * (*pn)->cell->getMM().y
				* (1 - (*pn)->data->aCenterData[Gerris_T])
				* (*pn)->data->aCenterData[Gerris_U];
		//mirro=====================
		sumveou += (*pn)->cell->area() * (*pn)->cell->getMM().y
				* (1 - (*pn)->data->aCenterData[Gerris_T])
				* (*pn)->data->aCenterData[Gerris_U];

		sumveov += (*pn)->cell->getMM().y * (*pn)->cell->area()
				* (1 - (*pn)->data->aCenterData[Gerris_T])
				* (*pn)->data->aCenterData[Gerris_V];
		//mirro=====================
		sumveov += -(*pn)->cell->getMM().y * (*pn)->cell->area()
				* (1 - (*pn)->data->aCenterData[Gerris_T])
				* (*pn)->data->aCenterData[Gerris_V];

		area += 2 * (*pn)->cell->getMM().y * (*pn)->cell->area()
				* (1 - (*pn)->data->aCenterData[Gerris_T]);
	}
	u = sumveou / area;
	v = sumveov / area;
}

void get_center_of_bubble(ListT<pQTNode>& listpn, Float& x, Float& y) {
	Float area = 0;
	Float sumveox = 0;
	Float sumveoy = 0;
	for (ListT<pQTNode>::const_iterator pn = listpn.begin(); pn != listpn.end();
			++pn) {
		sumveox += (*pn)->cell->area() * (*pn)->cell->getMM().y
				* (1 - (*pn)->data->aCenterData[Gerris_T])
				* (*pn)->cell->getCenterPoint().x;
		//mirro=====================
		sumveox += (*pn)->cell->area() * (*pn)->cell->getMM().y
				* (1 - (*pn)->data->aCenterData[Gerris_T])
				* (*pn)->cell->getCenterPoint().x;

		sumveoy += (*pn)->cell->getMM().y * (*pn)->cell->area()
				* (1 - (*pn)->data->aCenterData[Gerris_T])
				* (*pn)->cell->getCenterPoint().y;
		//mirro=====================
		sumveoy += -(*pn)->cell->getMM().y * (*pn)->cell->area()
				* (1 - (*pn)->data->aCenterData[Gerris_T])
				* (*pn)->cell->getCenterPoint().y;

		area += 2 * (*pn)->cell->getMM().y * (*pn)->cell->area()
				* (1 - (*pn)->data->aCenterData[Gerris_T]);
	}
	x = sumveox / area;
	y = sumveoy / area;
}

void save_list_seg(const ListT<Segment2D>& ls, string filename) {
	Matrix m(ls.size() * 2, 2);
	int i = 0;

}

string setfilename(string prefix, int group, Float i, string end) {
	stringstream ss;
	//ss << group << "-" << i << "-resEnd.txt";
	ss << prefix << group;
	ss.precision(2);
	ss.setf(std::ios::fixed, std::ios::floatfield);
	ss << "-" << i << end << ".txt";
	string str = ss.str();
	return str;
}

string setfilename(string prefix, Float i, string end) {
	stringstream ss;
	//ss << group << "-" << i << "-resEnd.txt";
	ss << prefix;
	ss.precision(2);
	ss.setf(std::ios::fixed, std::ios::floatfield);
	ss << "" << i << end << ".txt";
	string str = ss.str();
	return str;
}

string get_gerrisfilename(int group, Float i) {
	stringstream ss;
	//ss << group << "-" << i << "-resEnd.txt";
	ss << "re-" << group;
	ss.precision(2);
	ss.setf(std::ios::fixed, std::ios::floatfield);
	ss << "-" << i << ".txt";
	string str = ss.str();
	return str;
}

string get_gerrisfilename(Float i) {
	stringstream ss;
	//ss << group << "-" << i << "-resEnd.txt";
	ss << "re";
	ss.precision(2);
	ss.setf(std::ios::fixed, std::ios::floatfield);
	ss << "-" << i << ".txt";
	string str = ss.str();
	return str;
}

bool compare_seg(const Segment2D& s1, const Segment2D& s2) {
	if (s1.getCP().x < s2.getCP().x) {
		return true;
	} else if (s1.getCP().x == s2.getCP().x && s1.getCP().y > s2.getCP().y) {
		return true;
	} else {
		return false;
	}
}

void make_ListPoint2D(ListT<Point2D>& list, Point2D& pstr, Point2D& pend,
		int n) {
	Float dx = ABS(pend.x - pstr.x);
	Float dy = ABS(pend.y - pstr.y);
	Float ddx = dx / n;
	Float ddy = dy / n;
	for (int i = 1; i < n - 1; i++) {
		Point2D p(pstr.x + ddx * i, pstr.y + ddy * i);
		list.push_back(p);
	}
}
//this function will clear list
//        | zmax
//        |
//     ___| zsig
//    | v
//    |
//    |     zmin
//    r
void generate_solid_sig(ListT<Point2D>& listp, //list point
		Float r,               //location of x
		Float v,               //sig length
		Float zmin,            //
		Float zsig,            //
		Float zmax)            //
		{
	listp.clear();
	Point2D bp(r, zmin);
	Point2D sbp(r, zsig);
	generate_ListPoint2D(bp, sbp, 500, listp);
	Point2D sep(r + v, zsig);
	generate_ListPoint2D(sbp, sep, 50, listp);
	Point2D ep(r + v, zmax);
	generate_ListPoint2D(sep, ep, 500, listp);
}

//divide surface into two part;
void seperate_surface(const Float& tail_x,    //
		const Float& tail_xc,   //
		const ListT<Segment2D>& l_surface, //
		ListT<Segment2D>& b_sur,    //body surface
		ListT<Segment2D>& t_sur) {  //body surface
	ListT<Segment2D> p_surface;
	Point2D pp(tail_x, 0);
	for (ListT<Segment2D>::const_iterator iter = l_surface.begin();
			iter != l_surface.end(); iter++) {
		if (iter->PEX() <= tail_xc && iter->PSX() <= tail_xc) {
			if (iter->onWhichside3(pp) == -1) {
				t_sur.push_back((*iter));
			} else {
				b_sur.push_back((*iter));
			}
		} else {
			b_sur.push_back((*iter));
		}
	}
}

void new_gerris_file(GerrisFileAdp<Dimension_2D>*& gfa, string dir, int group,
		Float time) {
	//string filename = get_gerrisfilename(group, time);
	string filename = get_gerrisfilename(time);
	string strss = dir + filename;
	Point2D op(-0.25, 0.0);  //intial point
	int lmin = 5;            //level min
	int lmax = 8;            //level MAX
	Float boxl = 0.5;        //box lenght
	gfa = new GerrisFileAdp<Dimension_2D>(strss, lmin, lmax, op, boxl);
	gfa->forest.ConnectTrees();
	//
	//gfa->show_file_info();
	//gfa->forest.show_info();
}

void basic_sig(Forest2D& f, arrayListT<utPointer>& arg) //0->8
		{
	//expand arg
	int b_idx = 0;
	ListT<Segment2D>& l_surface = (*CAST(ListT<Segment2D>*, arg[b_idx]));
	Float& head_x = (*CAST(Float*, arg[b_idx + 1]));
	Float& tail_x = (*CAST(Float*, arg[b_idx + 2]));
	Float& max_y = (*CAST(Float*, arg[b_idx + 3]));
	Float& tail_xc = (*CAST(Float*, arg[b_idx + 4]));
	Float& u = (*CAST(Float*, arg[b_idx + 5]));
	Float& v = (*CAST(Float*, arg[b_idx + 6]));
	Float& x = (*CAST(Float*, arg[b_idx + 7]));
	Float& y = (*CAST(Float*, arg[b_idx + 8]));
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
}

void neck_sig(Forest2D& f, arrayListT<utPointer>& arg) {
	//expand arg
	ListT<Segment2D>& l_surface = (*CAST(ListT<Segment2D>*, arg[0]));
	Float& tail_x = (*CAST(Float*, arg[2]));
	Float& tail_xc = (*CAST(Float*, arg[4]));
	ListT<Segment2D>& b_sur = (*CAST(ListT<Segment2D>*, arg[9]));
	ListT<Segment2D>& t_sur = (*CAST(ListT<Segment2D>*, arg[10]));
	Point2D& p_neck = (*CAST(Point2D*, arg[11]));
	//
	seperate_surface(tail_x, tail_xc, l_surface, b_sur, t_sur);
	drawtofile_gnuplot("body.txt", b_sur, 1);
	drawtofile_gnuplot("tail.txt", t_sur, 1);

	//cut tail ================================
	ListT<Segment2D> lcal;
	for (ListT<Segment2D>::iterator iter = b_sur.begin(); iter != b_sur.end();
			++iter) {
		if (iter->is_gt_x(tail_xc + 0.1)) {
			lcal.push_back(*iter);
		}
	}

	arrayListT<Segment2D> arr_seg(lcal.size());
	int i = 0;
	for (ListT<Segment2D>::iterator iter = lcal.begin(); iter != lcal.end();
			++iter, ++i) {
		arr_seg[i] = (*iter);
	}
	//sort ========================
	Sort_Bubble(arr_seg, compare_seg);
	//getSlop =====================
	arrayListT<Point2D> arr_slop(arr_seg.size());
	for (arrayListT<Segment2D>::size_type i = 0; i < arr_slop.size(); ++i) {
		arr_slop[i].x = arr_seg[i].PC().x;
		arr_slop[i].y = arr_seg[i].getSlope();
	}
	//drawtofile_gnuplot("slop.txt", arr_slop, 1);

	//monotonic
	p_neck = Point2D(0.0, 0.0);
	int s = sign(arr_slop[0].y);
	int countch = 0;
	for (arrayListT<Segment2D>::size_type i = 1; i < arr_slop.size(); ++i) {
		//cout<< arr_slop[i].y <<endl;
		if (sign(arr_slop[i].y) != s) {
			countch++;
			if (countch == 2) { //the second point
				p_neck = arr_seg[i].PC();
			}
		}
		s = sign(arr_slop[i].y);
	}
}

void arg_output(const arrayListT<utPointer>& arg) {
	//expand arg
	//ListT<Segment2D>& l_surface = (*CAST(ListT<Segment2D>*, arg[0]));
	Float& head_x = (*CAST(Float*, arg[1]));
	Float& tail_x = (*CAST(Float*, arg[2]));
	Float& max_y = (*CAST(Float*, arg[3]));
	Float& tail_xc = (*CAST(Float*, arg[4]));
	Float& u = (*CAST(Float*, arg[5]));
	Float& v = (*CAST(Float*, arg[6]));
	Float& x = (*CAST(Float*, arg[7]));
	Float& y = (*CAST(Float*, arg[8]));
	//
	printf("head x  :%10.3lf\n", head_x);
	printf("tail x  :%10.3lf\n", tail_x);
	printf("tail xc :%10.3lf\n", tail_xc);
	printf("max y   :%10.3lf\n", max_y);
	printf("veo     :%10.3lf  %10.3lf\n", u, v);
	printf("c o mass:%10.3lf  %10.3lf\n", x, y);
}

void out_put_time_idx(string filename, Float t_str, Float t_end, Float dt) {
	int ic = 0;
	for (Float t = t_str; t <= t_end; t += dt) {
		ic++;
	}
	arrayList arrt(ic);
	arrayList arri(ic);
	int i = 0;
	for (Float t = t_str; t <= t_end; t += dt) {
		arri[i] = i;
		arrt[i] = t;
		i++;
	}
	drawtofile_gnuplot(filename, arri, arrt, 1);
}

void run_test() {
	//known value
	const Float sigx = 5.75;

	arrayListT<utPointer> arg(12);

	//loop for group ===============================
	string dir_ori = strDIR + "ori/";
	string dir_out = strDIR + "out_each_c/";
	string dir_outg = strDIR + "out_group_c/";

	int group = 3;
	Float t_str = 0;
	Float t_end = 25;
	Float dt = 0.2;
	out_put_time_idx(dir_outg + "time_idx.txt", t_str, t_end, dt);
	int ic = 0;

	for (Float t = t_str; t <= t_end; t += dt) {
		//arg set =====================
		ListT<Segment2D> l_surface;
		Float head_x;
		Float tail_x;
		Float max_y;
		Float tail_xc;
		Float u;
		Float v;
		Float x;
		Float y;
		ListT<Segment2D> b_sur;
		ListT<Segment2D> t_sur;
		Point2D p_neck;
		// =============================
		arg[0] = &l_surface;
		arg[1] = &head_x;
		arg[2] = &tail_x;
		arg[3] = &max_y;
		arg[4] = &tail_xc;
		arg[5] = &u;
		arg[6] = &v;
		arg[7] = &x;
		arg[8] = &y;
		arg[9] = &b_sur;  //
		arg[10] = &t_sur;
		arg[11] = &p_neck;
		//new forest ==================
		GerrisFileAdp<Dimension_2D>* gfa = NULL_PTR;
		new_gerris_file(gfa, dir_ori, group, t);
		basic_sig(gfa->forest, arg);

		//out_put_basic
		arrayList arrb(8);
		arrb[0] = head_x;
		arrb[1] = tail_x;
		arrb[2] = max_y;
		arrb[3] = tail_xc;
		arrb[4] = u;
		arrb[5] = v;
		arrb[6] = x;
		arrb[7] = y;
		arrayList arrid(8);
		arrid.assign_forward(0, 1);
		drawtofile_gnuplot(  //
				dir_out + setfilename("_bas_", t, ""), arrid, arrb, 1);
		//

		//out_put_surface
		drawtofile_gnuplot(  //
				dir_out + setfilename("_sur_", t, ""), l_surface, 1);
		//

		neck_sig(gfa->forest, arg);
		//out_put_neck_point
		drawtofile_gnuplot(  //
				dir_out + setfilename("_neckp_", t, ""), p_neck, 1);
		//p_neck.show();

		//

		delete gfa;
		ic++;
		cout << "File  ---> ";
		cout << setw(3) << ic << " " << setw(5) << setprecision(2) << t << "  "
				<< setw(10) << setprecision(5) << head_x << endl;
	}

	ListT<Point2D> listsig;
	generate_solid_sig(listsig, 0.5, 0.05, -0.25, sigx, 15);
	drawtofile_gnuplot(dir_outg + "line_sig_1.1.txt", listsig, 1);
	generate_solid_sig(listsig, 0.5, 0.10, -0.25, sigx, 15);
	drawtofile_gnuplot(dir_outg + "line_sig_1.2.txt", listsig, 1);
	generate_solid_sig(listsig, 0.5, 0.15, -0.25, sigx, 15);
	drawtofile_gnuplot(dir_outg + "line_sig_1.3.txt", listsig, 1);
	generate_solid_sig(listsig, 0.5, 0.20, -0.25, sigx, 15);
	drawtofile_gnuplot(dir_outg + "line_sig_1.4.txt", listsig, 1);
	cout<<"Finish run =====================\n"<<endl;
}
void run_test_contraction() {
	string dir = "/home/czhou/Gerris/Taylor_data/Contraction/re135eo150r0.90L4/";
	//known value
	const Float sigx = 10.5;

	arrayListT<utPointer> arg(12);

	//loop for group ===============================
	string dir_ori = dir + "ori/";
	string dir_out = dir + "out_each_c/";
	string dir_outg = dir + "out_group_c/";

	int group = 3;
	Float t_str = 0;
	Float t_end = 30;
	Float dt = 0.1;
	out_put_time_idx(dir_outg + "time_idx.txt", t_str, t_end, dt);
	int ic = 0;

	for (Float t = t_str; t <= t_end; t += dt) {
		//arg set =====================
		ListT<Segment2D> l_surface;
		Float head_x;
		Float tail_x;
		Float max_y;
		Float tail_xc;
		Float u;
		Float v;
		Float x;
		Float y;
		ListT<Segment2D> b_sur;
		ListT<Segment2D> t_sur;
		Point2D p_neck;
		// =============================
		arg[0] = &l_surface;
		arg[1] = &head_x;
		arg[2] = &tail_x;
		arg[3] = &max_y;
		arg[4] = &tail_xc;
		arg[5] = &u;
		arg[6] = &v;
		arg[7] = &x;
		arg[8] = &y;
		arg[9] = &b_sur;  //
		arg[10] = &t_sur;
		arg[11] = &p_neck;
		//new forest ==================
		GerrisFileAdp<Dimension_2D>* gfa = NULL_PTR;
		new_gerris_file(gfa, dir_ori, group, t);
		basic_sig(gfa->forest, arg);

		//out_put_basic
		arrayList arrb(8);
		arrb[0] = head_x;
		arrb[1] = tail_x;
		arrb[2] = max_y;
		arrb[3] = tail_xc;
		arrb[4] = u;
		arrb[5] = v;
		arrb[6] = x;
		arrb[7] = y;
		arrayList arrid(8);
		arrid.assign_forward(0, 1);
		drawtofile_gnuplot(  //
				dir_out + setfilename("_bas_", t, ""), arrid, arrb, 1);
		//

		//out_put_surface
		drawtofile_gnuplot(  //
				dir_out + setfilename("_sur_", t, ""), l_surface, 1);
		//

		//neck_sig(gfa->forest, arg);
		//out_put_neck_point
		//drawtofile_gnuplot(  //
		//		dir_out + setfilename("_neckp_", t, ""), p_neck, 1);
		//p_neck.show();

		//

		delete gfa;
		ic++;
		cout << "File  ---> ";
		cout << setw(3) << ic << " " << setw(5) << setprecision(2) << t << "  "
				<< setw(10) << setprecision(5) << head_x << endl;
	}

	cout<<"Finish run =====================\n"<<endl;
}
void pressure_on_centerline() {
	const string dir =
			"/Volumes/Untitled/Taylor_data/Expansion/exp-re100eo100r1.2L4/";
	string dir_ori = dir + "ori/";
	string dir_out = dir + "out_each_c/";
	ListT<Float> listfilm;

	Float sigratio = 1.2;
	Float sigx = 5.75;
	Float yo = 0.5;
	Float yexp = yo * sigratio;

	int group = 3;
	Float t_str = 0;
	Float t_end = 25;
	Float dt = 0.2;

	//loop ------------------------
	for (Float time = t_str; time <= t_end; time += dt) {
		//new forest ------------------
		GerrisFileAdp<Dimension_2D>* gfa = NULL_PTR;
		new_gerris_file(gfa, dir_ori, group, time);
		cout << "time  " << time << " \n";
		cout << "Finsh load tree =====       \n\n";
		//shift pressure --------------
		Forest2D::iterator itlast = gfa->forest.last();
		Float p0 = itlast->data->aCenterData[Gerris_P];
		Float xc0 = itlast->cell->getCenterPoint().x;
		arrayList_st arridx(1);
		arridx[0] = Gerris_P;
		arrayList arrval(1);
		arrval[0] = -p0;
		plus_scalar_on_leaf( // 2D Forest
				gfa->forest, //pQuadTree
				arridx, //data index
				arrval //data plus
				);
		//static pressure revise=======
		for (Forest2D::iterator iter = gfa->forest.begin();
				iter != gfa->forest.end(); ++iter) {
			Float sv = xc0 - iter->cell->getCenterPoint().x;
			iter->data->aCenterData[Gerris_P] -= sv;
		}

		ListT<Segment2D> l_surface;
		cal_surface((*gfa), l_surface);
		//list center of film ==============
		Float head_x = get_head_x(l_surface);
		Float tail_x = get_tail_x(l_surface);
		ListT<Point2D> lcp;
		int numtap = 50;
		Float dx = (head_x - tail_x) / numtap;
		for (int i = 0; i < numtap + 1; i++) {
			Segment2D seg(Point2D(head_x - i * dx, 2),
					Point2D(head_x - i * dx, -0.5));
			int numc = 0;
			Segment2D seginter(Point2D(0, -0.5), Point2D(1, -0.5));
			for (ListT<Segment2D>::iterator iter = l_surface.begin();
					iter != l_surface.end(); ++iter) {
				if (isIntersect(seg, (*iter))) {
					numc++;
					if (iter->getCP().y > seginter.getCP().y) {
						seginter = (*iter);
					}
				}
			}
			if (numc > 0) {
				Float b = seg.PEX() > sigx ? yexp : yo;
				Point2D p(seg.PEX(),
						(b - seginter.getCP().y) / 2.0 + seginter.getCP().y);
				lcp.push_back(p);
			}
		}
		//center stream line
		//drawtofile_gnuplot("cpl.txt", lcp, 1);
		//drawtofile_gnuplot("bubble.txt", l_surface, 1);

		ListT<Pair<int, arrayList> > lres;
		interpolate(	// 2D Forest
				gfa->forest,	//Forest
				lcp,	//point
				arridx,	//data index
				lres	//data res if int=-1 result is wrong
				);
		//output =====
		FILE* data = open_file(dir_out + setfilename("_clp_", time, ""), 1);
		ListT<Point2D>::iterator iterp = lcp.begin();
		for (ListT<Pair<int, arrayList> >::iterator iter = lres.begin();
				iter != lres.end(); ++iter) {
			fprintf(data, "%f %f ", iterp->x, iterp->y);
			fprintf(data, "%d ", iter->first);
			for (arrayList::size_type in = 0; in < iter->second.size(); in++) {
				if (iter->first == 1) {
					fprintf(data, "%f ", iter->second[in]);
				} else {
					fprintf(data, "nan ");
				}
			}
			fprintf(data, "\n");
			++iterp;
		}
		fclose(data);
		delete gfa;
	}
}

void pressure_contour_line() {
	const string dir =
			"/Volumes/Untitled/Taylor_data/Expansion/exp-re100eo100r1.2L4/";
	string dir_ori = dir + "ori/";
	string dir_out = dir + "out_each_c/";
	ListT<Float> listfilm;

	Float sigratio = 1.2;
	Float sigx = 5.75;
	Float yo = 0.5;
	Float yexp = yo * sigratio;

	int group = 3;
	Float t_str = 0;
	Float t_end = 25;
	Float dt = 0.2;

	//loop ------------------------
	for (Float time = t_str; time <= t_end; time += dt) {
		//new forest ------------------
		//Float time =10.0;
		GerrisFileAdp<Dimension_2D>* gfa = NULL_PTR;
		new_gerris_file(gfa, dir_ori, group, time);
		cout << "time  " << time << " \n";
		cout << "Finsh load tree =====       \n\n";
		Forest2D::iterator itlast = gfa->forest.last();
		Float p0 = itlast->data->aCenterData[Gerris_P];
		Float xc0 = itlast->cell->getCenterPoint().x;
		cout << "p0   " << p0 << endl;

		//advance(it,300);
		//Float p1= it->data->aCenterData[Gerris_P];
		//Float xc1= it->cell->getCenterPoint().x;
		//cout<<"p1   "<<p1<<endl;
		//cout<<"dp   "<<p1 - p0<<endl;
		//cout<<"dx   "<<xc1 - xc0<<endl;
		//shift pressure ==============
		arrayList_st arridx(1);
		arridx[0] = Gerris_P;
		arrayList arrval(1);
		arrval[0] = -p0;
		plus_scalar_on_leaf(	// 2D Forest
				gfa->forest,	//pQuadTree
				arridx,	//data index
				arrval	//data plus
				);
		//static pressure revise=======
		for (Forest2D::iterator iter = gfa->forest.begin();
				iter != gfa->forest.end(); ++iter) {
			Float sv = xc0 - iter->cell->getCenterPoint().x;
			iter->data->aCenterData[Gerris_P] -= sv;
		}
		Float max_p = get_max_value(gfa->forest, Gerris_P);
		Float min_p = get_min_value(gfa->forest, Gerris_P);
		cout.precision(5);
		cout << "max p  " << max_p << endl;
		cout << "min p  " << min_p << endl;
		ListT<Segment2D> l_surface;
		cal_surface((*gfa), l_surface);
		Float head_x = get_head_x(l_surface);
		Float tail_x = get_tail_x(l_surface);
		cout << "head x " << head_x << endl;
		cout << "tail x " << tail_x << endl;
		//draw pressure con============
		draw_gnuplot_as_contour(	// 2D QuadTree
				dir_out + setfilename("_pcon_", time, ""),//filename
				gfa->forest,	//tree
				1,	//mode
				Gerris_P	//idx x
				);
		//draw pressure line ===============
		int nl = 10;
		arrayList alevel(nl);
		//cal level
		Float d = (max_p - min_p) / (nl + 1);
		for (int i = 0; i < nl; i++) {
			alevel[i] = min_p + (i + 1) * d;
		}
		ListT<pQTNode> lnode;
		getListpNode_leaf_center_data_in_range(lnode, gfa->forest, Gerris_T, 0,
				1, Range_oc);
		draw_gnuplot_as_contour_line(dir_out + setfilename("_pconl_", time, ""), lnode, 1, Gerris_P,
				alevel);
		delete gfa;
	}
}

}

#endif
