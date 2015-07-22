/************************
 //  \file   IO_gerris.h
 //  \brief
 // 
 //  \author czhou
 //  \date   5 févr. 2015 
 ***********************/
#ifndef IO_GERRIS_H_
#define IO_GERRIS_H_

#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include "../Grid/SPTree.h"
#include "../TypeDef.h"
#include "../Algebra/Matrix.h"
#include "../Utility/List.h"
#include "../Grid/Forest.h"

using namespace std;

namespace Larus {
const int Gerris_P = 0;
const int Gerris_Pmac = 1;
const int Gerris_U = 2;
const int Gerris_V = 3;
const int Gerris_T = 4;
const int Gerris_T_x = 5;
const int Gerris_T_y = 6;
const int Gerris_T_alpha = 7;
const int Gerris_T1 = 8;
const int Gerris_K = 9;
const int Gerris_kmax = 10;

int compareLineX(const vector<Float> &a, const vector<Float> &b);

int compareLineY(const vector<Float> &a, const vector<Float> &b);

int cal_level(Float x, Float b, Float L = 1, int ls = 4, int le = 7);

template<typename DIM_TYPE>
class GerrisFileAdp {
public:
	typedef typename DIM_TYPE::Tree Tree;
	typedef typename DIM_TYPE::pTree pTree;
	typedef typename DIM_TYPE::Cell Cell;
	typedef typename DIM_TYPE::pCell pCell;
	typedef typename DIM_TYPE::Node Node;
	typedef typename DIM_TYPE::pNode pNode;
	typedef typename DIM_TYPE::CellData CellData;
	typedef typename DIM_TYPE::pCellData pCellData;
	typedef typename DIM_TYPE::Point Point;
	typedef typename DIM_TYPE::pPoint pPoint;

	typedef typename DIM_TYPE::size_type size_type;
private:
	//data ======================================
	int m_col;
	int m_row;
	string m_filename;
	ifstream m_fs;
	Float LengthofGerrisbox;
	//
	int countCol();
	int countRow();
	void creatForest();
	void buildTree();
	//void creatListCell();
	//void creatListTree();
public:
	//data ======================================
	Forest<DIM_TYPE> forest;
	vector<vector<Float> > vv;
	Point op;
	int stalevel;
	int endlevel;

	GerrisFileAdp(string, int, int, Point&, Float = 1);
	~GerrisFileAdp();
	//File=======================================
	int getCol();
	int getRow();
	string getFilename();
	Float getColmin(int j);
	Float getColmax(int j);
	Float getXmin();
	Float getXmax();
	Float getYmin();
	Float getYmax();
	Float getZmin();
	Float getZmax();
	void show_file_info();
	//=======================================

	int calNumBoxX();
	int calNumBoxY();
	int calNumBoxZ();
	int totalNumBox();
	int countPointsinCell(Cell& c);
	int countLevelPointsinCell(Cell& c, int level);
	//IO=====================================
	void output_original_file(string filename, int mode);

};

//Class GerrisFileAdp
template<typename DIM_TYPE>
int GerrisFileAdp<DIM_TYPE>::countRow() {
	string tmp, firstline;
	m_fs.clear();
	m_fs.seekg(0, ios::beg);
	getline(m_fs, firstline);
	int n = 0;
	if (m_fs.is_open()) {
		while (getline(m_fs, tmp)) {
			n++;
		}
		return n;
	} else {
		return -1;
	}
}
template<typename DIM_TYPE>
int GerrisFileAdp<DIM_TYPE>::countCol() {
	string firstline;
	m_fs.clear();
	m_fs.seekg(0, ios::beg);
	getline(m_fs, firstline);
	int col = 0;
	for (string::iterator it = firstline.begin(); it != firstline.end(); ++it) {
		if (*it == ':') {
			col++;
		}
	}
	return col;
}
template<typename DIM_TYPE>
GerrisFileAdp<DIM_TYPE>::GerrisFileAdp(string filename, int tl, int dl,
		GerrisFileAdp<DIM_TYPE>::Point& ori, Float lgbox) {
	LengthofGerrisbox = lgbox;
	stalevel = tl;
	endlevel = dl;
	m_filename = filename;
	string firstline;
	op = ori;
	m_fs.open(filename.c_str(), ios::in);
	if (m_fs.is_open()) {
		m_col = countCol();
		m_row = countRow();
	} else {
		cerr << ">! Can't find file \n" << filename << endl;
		exit(-1);
	}
	m_fs.clear();
	m_fs.seekg(0, ios::beg);
	getline(m_fs, firstline); //ignore the first line
	for (int i = 0; i < m_row; i++) {
		vector<Float> tmv;
		for (int c = 0; c < m_col; c++) {
			Float tmpv;
			m_fs >> tmpv;
			tmv.push_back(tmpv);
		}
		vv.push_back(tmv);
	}

	//creat======================================
	creatForest();
	buildTree();
	//this->forest.ConnectTrees();
}
template<typename DIM_TYPE>
void GerrisFileAdp<DIM_TYPE>::show_file_info() {
	std::cout << "Gerris File Info: \n";
	std::cout << " > " << m_filename << endl;
	std::cout << "Dim:    " << DIM_TYPE::DIM << "\n";
	std::cout.flags(std::ios::right);
	std::cout << "Row:    " << getRow() << endl;
	std::cout << "Col:    " << getCol() << endl;
	std::cout << "Xrange: " << getXmin() << " => " << getXmax() << endl;
	std::cout << "Yrange: " << getYmin() << " => " << getYmax() << endl;
	if (DIM_TYPE::DIM == 3) {
		std::cout << "Zrange: " << getZmin() << " => " << getZmax() << endl;
	}
	std::cout << "BoxLen: " << LengthofGerrisbox << endl;
	std::cout << "Cell nX " << calNumBoxX() << endl;
	std::cout << "Cell nY " << calNumBoxY() << endl;
	if (DIM_TYPE::DIM == 3) {
		std::cout << "Cell nZ " << calNumBoxZ() << endl;
	}
	std::cout << "Level   " << stalevel << " => " << endlevel << endl;
}
template<typename DIM_TYPE>
void GerrisFileAdp<DIM_TYPE>::output_original_file(string filename, int mode) {
	FILE *data = NULL;
	if (mode == 1) {
		data = fopen(filename.c_str(), "w"); //write
	} else if (mode == 2) {
		data = fopen(filename.c_str(), "a"); //append
		if (data == NULL) {
			std::cerr << "!> Can't find file. " << filename.c_str() << " \n";
			exit(0);
		}
	}
	//output root
	for (int i = 0; i < m_row; i++) {
		for (int c = 0; c < m_col; c++) {
			//std::cout<<std::scientific<<vv[i][c]<<" ";
			fprintf(data, "%e ", vv[i][c]);
		}
		fprintf(data, "\n");
		//std::cout<<"\n";
	}

}
template<typename DIM_TYPE>
int GerrisFileAdp<DIM_TYPE>::calNumBoxX() {
	Float max = getXmax();
	for (int i = 1; i < 1000; i++) {
		if (max < (op.x + i * LengthofGerrisbox)) {
			return i;
			break;
		}
	}
	return 0;
}
template<typename DIM_TYPE>
int GerrisFileAdp<DIM_TYPE>::calNumBoxY() {
	Float max = getYmax();
	for (int i = 1; i < 1000; i++) {
		if (max < (op.y + i * LengthofGerrisbox)) {
			return i;
			break;
		}
	}
	return 0;
}
template<typename DIM_TYPE>
int GerrisFileAdp<DIM_TYPE>::calNumBoxZ() {
	if (DIM_TYPE::DIM == 2) {
		return 0;
	}
	Float max = getZmax();
	for (int i = 1; i < 1000; i++) {
		if (max < (op[2] + i * LengthofGerrisbox)) {
			return i;
			break;
		}
	}
	return 0;
}
template<typename DIM_TYPE>
GerrisFileAdp<DIM_TYPE>::~GerrisFileAdp() {
	m_fs.close();
}
template<typename DIM_TYPE>
int GerrisFileAdp<DIM_TYPE>::getCol() {
	return m_col;
}
template<typename DIM_TYPE>
int GerrisFileAdp<DIM_TYPE>::getRow() {
	return m_row;
}
template<typename DIM_TYPE>
int GerrisFileAdp<DIM_TYPE>::totalNumBox() {
	return forest.getNumEnableTree();
}
template<typename DIM_TYPE>
Float GerrisFileAdp<DIM_TYPE>::getColmin(int j) {
	Float res = vv[0][j];
	for (unsigned int i = 1; i < vv.size(); i++) {
		if (vv[i][j] < res) {
			res = vv[i][j];
		}
	}
	return res;
}
template<typename DIM_TYPE>
Float GerrisFileAdp<DIM_TYPE>::getColmax(int j) {
	Float res = vv[0][j];
	for (unsigned int i = 1; i < vv.size(); i++) {
		if (vv[i][j] > res) {
			res = vv[i][j];
		}
	}
	return res;
}

template<typename DIM_TYPE>
void GerrisFileAdp<DIM_TYPE>::creatForest() {
	int nx = calNumBoxX();
	int ny = calNumBoxY();
	int nz = DIM_TYPE::DIM == 3 ? calNumBoxZ() : 1;
	Point minp = op;
	forest.reconstruct(nx, ny, nz);
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			for (int k = 0; k < nz; k++) {
				Point mm(minp[0] + i * LengthofGerrisbox,
						minp[1] + j * LengthofGerrisbox,
						(DIM_TYPE::DIM == 3 ?
								(minp[2] + k * LengthofGerrisbox) : 0));
				Point pp(minp[0] + (i + 1) * LengthofGerrisbox,
						minp[1] + (j + 1) * LengthofGerrisbox,
						(DIM_TYPE::DIM == 3 ?
								(minp[2] + (k + 1) * LengthofGerrisbox) : 0));
				Cell c(mm, pp);
				pTree p = new Tree(c, endlevel);
				forest.setpTree(p, i, j, k);
			}
		}
	}
}

//call back function

void _t_condition_point_at_which_child(arrayList& arrt, pQTNode pnode,
		utPointer p);

void _t_condition_point_at_which_child(arrayList& arrt, pOCNode pnode,
		utPointer p);

void _visit_insert_level_copydata(pOCNode pn, utPointer p);
void _visit_insert_level_copydata(pQTNode pn, utPointer p);

template<typename DIM_TYPE>
void _insert_partial_with_data(DIM_TYPE type, typename DIM_TYPE::pTree pt,
		typename DIM_TYPE::Point& point, int level,
		typename DIM_TYPE::pCellData pcd) {
	arrayListT<utPointer> vp(4);
	vp[0] = &point;
	vp[1] = pt;
	vp[2] = &level;
	vp[3] = pcd;
	pt->Conditional_Traversal(_t_condition_point_at_which_child,
			_visit_insert_level_copydata, &vp);
}

template<typename DIM_TYPE>
void GerrisFileAdp<DIM_TYPE>::buildTree() {
	for (int ii = 0; ii < forest.iLen(); ii++) {
		for (int jj = 0; jj < forest.jLen(); jj++) {
			for (int kk = 0; kk < (DIM_TYPE::DIM == 3 ? forest.kLen() : 1);
					kk++) {
				pCell c = forest.getpTree(ii, jj, kk)->getpRootNode()->cell;
				for (int i = 0; i < m_row; i++) {
					Point p(vv[i][0], vv[i][1], vv[i][2]);
					if (c->isInCell(p)) {
						//IN cell, then creat a tree
						int level = cal_level(vv[i][0],
								c->getPoint(eCPL_M, eCPL_M, eCPL_M).x,
								c->getDx(), stalevel, endlevel);
						arrayList atmp(m_col - 3);
						for (int j = 0; j < atmp.size(); j++) {
							atmp[j] = vv[i][3 + j];
						}
						pCellData pc = new CellData();   //new!!!!!!!!!!!!!!!!
						pc->aCenterData = atmp;
						_insert_partial_with_data(DIM_TYPE(),
								forest.getpTree(ii, jj, kk), p, level, pc);
					}
				}
				if (forest.getpTree(ii, jj, kk)->isNewTree()) {
					forest.set_attribution(ATT_DISABLE, ii, jj, kk);
				}
			}
		}
	}
}
template<typename DIM_TYPE>
Float GerrisFileAdp<DIM_TYPE>::getXmin() {
	return getColmin(0);
}
template<typename DIM_TYPE>
Float GerrisFileAdp<DIM_TYPE>::getXmax() {
	return getColmax(0);
}
template<typename DIM_TYPE>
Float GerrisFileAdp<DIM_TYPE>::getYmin() {
	return getColmin(1);
}
template<typename DIM_TYPE>
Float GerrisFileAdp<DIM_TYPE>::getYmax() {
	return getColmax(1);
}
template<typename DIM_TYPE>
Float GerrisFileAdp<DIM_TYPE>::getZmin() {
	return getColmin(2);
}
template<typename DIM_TYPE>
Float GerrisFileAdp<DIM_TYPE>::getZmax() {
	return getColmax(2);
}

template<typename DIM_TYPE>
string GerrisFileAdp<DIM_TYPE>::getFilename() {
	return m_filename;
}

//count points in cell
template<typename DIM_TYPE>
int GerrisFileAdp<DIM_TYPE>::countPointsinCell(
		GerrisFileAdp<DIM_TYPE>::Cell& c) {
	int res = 0;
	for (int i = 0; i < m_row; i++) {
		if (c.isInCell(Point2D(vv[i][0], vv[i][1]))) {
			//IN cell
			res++;
		}
	}
	return res;
}
template<typename DIM_TYPE>
int GerrisFileAdp<DIM_TYPE>::countLevelPointsinCell(
		GerrisFileAdp<DIM_TYPE>::Cell& c, int level) {
	int res = 0;
	for (int i = 0; i < m_row; i++) {
		if (c.isInCell(Point2D(vv[i][0], vv[i][1]))) {
			//IN cell
			if (level
					== cal_level(vv[i][0], c.getMM().x, c.getDx(), level,
							level)) {
				res++;
			}
		}
	}
	return res;
}

}

#endif /* IO_GERRIS_H_ */
