#include "test_gerris.h"

namespace Larus {

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

string setfilename(string prefix, int group, Float i, string end) {
	stringstream ss;
	//ss << group << "-" << i << "-resEnd.txt";
	ss << prefix << group;
	ss.precision(2);
	ss.setf(std::ios::fixed, std::ios::floatfield);
	ss << "_" << i << end << ".txt";
	string str = ss.str();
	return str;
}

string set_dir_and_filename(string dir, string prefix, int group, Float i,
		string end) {
	stringstream ss;
	//ss << group << "-" << i << "-resEnd.txt";
	ss << dir << prefix << group;
	ss.precision(2);
	ss.setf(std::ios::fixed, std::ios::floatfield);
	ss << "_" << i << end << ".txt";
	string str = ss.str();
	return str;
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

void new_gerris_file(GerrisFileAdp<Dimension_2D>*& gfa, string dir, int group,
		Float time, int lmin, int lmax) {
	string filename = get_gerrisfilename(group, time);
	string strss = dir + filename;
	Point2D op(-0.25, 0.0);  //intial point
	Float boxl = 0.5;        //box lenght
	gfa = new GerrisFileAdp<Dimension_2D>(strss, lmin, lmax, op, boxl);
	gfa->forest.ConnectTrees();
	//
	//gfa->show_file_info();
	//gfa->forest.show_info();
}

void new_gerris_file( //
		GerrisFileAdp<Dimension_2D>*& gfa, //
		string dir, string name, //
		int lmin, //
		int lmax) {
	string strss = dir + name;
	Point2D op(-0.25, 0.0);  //intial point
	Float boxl = 0.5;        //box lenght
	gfa = new GerrisFileAdp<Dimension_2D>(strss, lmin, lmax, op, boxl);
	gfa->forest.ConnectTrees();
	//
	//gfa->show_file_info();
	//gfa->forest.show_info();
}

void get_surface(GerrisFileAdp<Dimension_2D>* gfa, ListT<Segment2D>& lseg) {
	getListSegment( //
			lseg, //list Segment
			gfa->forest, //forest
			Gerris_T, //data color index
			Gerris_T_x, //data index
			Gerris_T_y, //data index
			Gerris_T_alpha  //data index
			);
}

}  // end of name space
