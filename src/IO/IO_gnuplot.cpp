/************************
 //  \file   GeoIO_gnuplot.cpp
 //  \brief
 // 
 //  \author czhou
 //  \date   23 juil. 2014 
 ***********************/

#include "IO_gnuplot.h"

#include "IO.h"

#include <assert.h>
#include <math.h>
#include <string.h>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>

namespace Larus {

using namespace std;

void drawtofile_gnuplot(string filename, Segment2D& s, int mode) {
	FILE* data = open_file(filename, mode);
	fprintf(data, "%f %f \n", s.PSX(), s.PSY());
	fprintf(data, "%f %f \n\n", s.PEX(), s.PEY());
	fclose(data);
}

int readfile_gnuplot_countline(string filename) {
	FILE *fs;
	int nline = 0;
	fs = fopen(filename.c_str(), "r");
	if (fs != NULL) {
		while (!feof(fs)) {
			char mid = fgetc(fs);
			if (mid == '\n')
				nline++;
		}
	} else {
		cerr << ">! Can't find file " << filename << endl;
	}
	fclose(fs);
	return nline;
}

void FILE_Point2D(FILE *&data, const Point2D& p) {
	fprintf(data, "%f %f\n", p.x, p.y);
}

void drawtofile_gnuplot(string filename, const Point2D& p, int mode) {
	FILE* data = open_file(filename, mode);
	FILE_Point2D(data, p);
	fclose(data);
}

void FILE_Arrow2D(FILE *&data, const Arrow2D& p) {
	fprintf(data, "%f %f %f %f\n", p[0], p[1], p[2], p[3]);
}

void drawtofile_gnuplot(string filename, const Arrow2D& p, int mode) {
	FILE* data = open_file(filename, mode);
	FILE_Arrow2D(data, p);
	fclose(data);
}

void drawtofile_gnuplot(string filename, Polygon& p, int mode) {
	FILE* data = open_file(filename, mode);
	if (p.getNumVertexs() == 0) {
		std::cerr
				<< " >! Function drawtofile_gnuplot(string, Polygon, int)\n   Polygon is Empty! \n";
		fclose(data);
		return;
	}
	for (int i = 0; i < p.getNumVertexs() - 1; i++) {
		fprintf(data, "%f %f \n", p.getVertex(i).x, p.getVertex(i).y);
		fprintf(data, "%f %f \n\n", p.getVertex(i + 1).x, p.getVertex(i + 1).y);
	}
	fprintf(data, "%f %f \n", p.getVertex(p.getNumVertexs() - 1).x,
			p.getVertex(p.getNumVertexs() - 1).y);
	fprintf(data, "%f %f \n\n", p.getVertex(0).x, p.getVertex(0).y);
	fclose(data);
}

void drawtofile_gnuplot(string filename, const ListT<Point2D>& listp,
		int mode) {
	FILE* data = open_file(filename, mode);
	for (ListT<Point2D>::const_iterator it = listp.begin(); it != listp.end();
			it++) {
		FILE_Point2D(data, (*it));
	}
	fclose(data);
}

void drawtofile_gnuplot(string filename, const arrayListT<Point2D>& arrp,
		int mode) {
	FILE* data = open_file(filename, mode);
	for (arrayListT<Point2D>::size_type i = 0; i < arrp.size(); ++i) {
		FILE_Point2D(data, arrp[i]);
	}
	fclose(data);
}

void drawtofile_gnuplot(string filename, const Forest2D& f, int mode) {
	FILE* data = open_file(filename, mode);
	for (Forest2D::const_iterator it = f.begin(); it != f.end(); ++it) {
		gnuplot_file_data(data, *(it->cell));
	}
	fclose(data);
}

void drawtofile_gnuplot(string filename, const ListT<Segment2D>& listp,
		int mode) {
	FILE* data = open_file(filename, mode);
	for (ListT<Segment2D>::const_iterator iter = listp.begin();
			iter != listp.end(); iter++) {
		fprintf(data, "%f %f \n", (*iter).PS().x, (*iter).PS().y);
		fprintf(data, "%f %f \n", (*iter).PE().x, (*iter).PE().y);
		fprintf(data, "\n");
	}
	fclose(data);
}

void drawtofile_gnuplot(string filename, const arrayListT<Segment2D>& arrs,
		int mode) {
	FILE* data = open_file(filename, mode);
	for (arrayListT<Segment2D>::size_type i = 0; i < arrs.size(); ++i) {
		fprintf(data, "%f %f %d\n", arrs[i].PSX(), arrs[i].PSY(), i);
		fprintf(data, "%f %f %d\n", arrs[i].PEX(), arrs[i].PEY(), i);
		fprintf(data, "\n");
	}
	fclose(data);
}

void drawtofile_gnuplot(string filename, const arrayList& arr, int mode) {
	FILE* data = open_file(filename, mode);
	for (int i = 0; i < arr.size(); i++) {
		fprintf(data, "%f \n", arr[i]);
	}
	fclose(data);
}

void drawtofile_gnuplot(string filename, const arrayList& arr1,
		const arrayList& arr2, int mode) {
	ASSERT(arr1.size() == arr2.size());
	FILE* data = open_file(filename, mode);
	for (int i = 0; i < arr1.size(); i++) {
		fprintf(data, "%f %f\n", arr1[i], arr2[i]);
	}
	fclose(data);
}

void drawtofile_gnuplot(string filename, const ListT<Arrow2D>& listp,
		int mode) {
	FILE* data = open_file(filename, mode);
	for (ListT<Arrow2D>::const_iterator it = listp.begin(); it != listp.end();
			it++) {
		FILE_Arrow2D(data, (*it));
	}
	fclose(data);
}

void drawtofile_gnuplot(string filename,
		const ListT<Pair<Point2D, arrayList> >& listp, int mode) {
	FILE* data = open_file(filename, mode);
	for (ListT<Pair<Point2D, arrayList> >::const_iterator it = listp.begin();
			it != listp.end(); it++) {
		// out put point
		fprintf(data, "%f %f", it->first.x, it->first.y);
		for (arrayList::const_iterator ita = it->second.begin();
				ita != it->second.end(); ita++) {
			fprintf(data, " %f", (*ita));
		}
		fprintf(data, "\n");
	}
	fclose(data);
}

void generate_ListPoint2D(const Point2D& b, const Point2D& e, int num,
		ListT<Point2D>& listp) {
	Float ddx = e.x - b.x;
	Float ddy = e.y - b.y;
	Float dx = ddx / num;
	Float dy = ddy / num;
	listp.push_back(b);
	for (int i = 0; i < num; ++i) {
		Point2D ptmp(b.x + dx * (i + 1), b.y + dy * (i + 1));
		listp.push_back(ptmp);
	}
	listp.push_back(e);
}

void generate_ListPoint2D_circle(const Point2D& c, Float r, int num,
		ListT<Point2D>& listp) {
	assert(num > 0);
	Float dsita = 2 * LarusDef::PI / num;
	for (int i = 0; i <= num; ++i) {
		Float sita = i * dsita;
		Float x = c.x + r * cos(sita);
		Float y = c.y + r * sin(sita);
		listp.push_back(Point2D(x, y));
	}
}

void gnuplot_file_data(FILE* data, const Cell2D& cell) {
	fprintf(data, "%f %f\n", cell.get(CSAxis_X, eCPL_M),
			cell.get(CSAxis_Y, eCPL_M));
	fprintf(data, "%f %f\n", cell.get(CSAxis_X, eCPL_P),
			cell.get(CSAxis_Y, eCPL_M));
	fprintf(data, "%f %f\n", cell.get(CSAxis_X, eCPL_P),
			cell.get(CSAxis_Y, eCPL_P));
	fprintf(data, "%f %f\n", cell.get(CSAxis_X, eCPL_M),
			cell.get(CSAxis_Y, eCPL_P));
	fprintf(data, "%f %f\n", cell.get(CSAxis_X, eCPL_M),
			cell.get(CSAxis_Y, eCPL_M));
}

void gnuplot_file_data(FILE* data, const QTNodeFace& face) {
	switch (face.direction) {
	case SPD_IM: {
		fprintf(data, "%f %f\n", face.pnode->cell->get(CSAxis_X, eCPL_M),
				face.pnode->cell->get(CSAxis_Y, eCPL_M));
		fprintf(data, "%f %f\n", face.pnode->cell->get(CSAxis_X, eCPL_M),
				face.pnode->cell->get(CSAxis_Y, eCPL_P));
		break;
	}
	case SPD_IP: {
		fprintf(data, "%f %f\n", face.pnode->cell->get(CSAxis_X, eCPL_P),
				face.pnode->cell->get(CSAxis_Y, eCPL_M));
		fprintf(data, "%f %f\n", face.pnode->cell->get(CSAxis_X, eCPL_P),
				face.pnode->cell->get(CSAxis_Y, eCPL_P));
		break;
	}
	case SPD_JP: {
		fprintf(data, "%f %f\n", face.pnode->cell->get(CSAxis_X, eCPL_M),
				face.pnode->cell->get(CSAxis_Y, eCPL_P));
		fprintf(data, "%f %f\n", face.pnode->cell->get(CSAxis_X, eCPL_P),
				face.pnode->cell->get(CSAxis_Y, eCPL_P));
		break;
	}
	case SPD_JM: {
		fprintf(data, "%f %f\n", face.pnode->cell->get(CSAxis_X, eCPL_M),
				face.pnode->cell->get(CSAxis_Y, eCPL_M));
		fprintf(data, "%f %f\n", face.pnode->cell->get(CSAxis_X, eCPL_P),
				face.pnode->cell->get(CSAxis_Y, eCPL_M));
		break;
	}
	default: {
		ASSERT_MSG(false, " !> Error face direction (gnuplot_inline_data)\n");
	}
	}
}

}
