#ifndef _TEST_GERRIS_H_
#define _TEST_GERRIS_H_

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

using namespace std;

namespace Larus {
// get file name
string get_gerrisfilename(int group, Float i);
string get_gerrisfilename(Float i);
string setfilename(string prefix, int group, Float i, string end);
string set_dir_and_filename(string dir, string prefix, int group, Float i,
		string end);
/*
 *  new gerris file
 */
void new_gerris_file(
		GerrisFileAdp<Dimension_2D>*& gfa,
        string dir,
		int group,
		Float time);
void new_gerris_file(
		GerrisFileAdp<Dimension_2D>*& gfa,
        string dir,
		int group,
		Float time,
		int lmin,
		int lmax);
void new_gerris_file( //
		GerrisFileAdp<Dimension_2D>*& gfa, //
		string dir, string name, //
		int lmin, //
		int lmax);
void get_surface(GerrisFileAdp<Dimension_2D>* gfa,
		ListT<Segment2D>& lseg);
}


#endif



