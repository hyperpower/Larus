/************************
 //  \file   VOF.cpp
 //  \brief
 // 
 //  \author czhou
 //  \date   27 janv. 2015 
 ***********************/

#include "VOF.h"
#include "../Utility/Pair.h"

#include "../Algebra/Arithmetic.h"

#include "../IO/IO.h"
#include "../IO/IO_vtk.h"
#include <iostream>

namespace Larus {

//===============================================

//===============================================

const int CASE_3D_8[2][2][2] = { //
		//
				{ { 0, 4 }, { 3, 7 } }, //
				{ { 1, 5 }, { 2, 6 } } };

int whichcase_3D_8(Float a, Float b, Float c) {
	int idx_a = (a >= 0) ? 0 : 1;
	int idx_b = (b >= 0) ? 0 : 1;
	int idx_c = (c >= 0) ? 0 : 1;
	return CASE_3D_8[idx_a][idx_b][idx_c];
}

int whichcase_3D_8(const Plane& p) {
	Float a = p.getA();
	Float b = p.getB();
	Float c = p.getC();
	return whichcase_3D_8(a, b, c);
}

const int COE_CAL_ALPHA[8][3] = { //
		//
				{ 0, 0, 0 }, //
				{ 1, 0, 0 }, //
				{ 1, 1, 0 }, //
				{ 0, 1, 0 }, //
				{ 0, 0, 1 }, //
				{ 1, 0, 1 }, //
				{ 1, 1, 1 }, //
				{ 0, 1, 1 }, //
		};

Float newAlpha_Forward(int ca, const Plane& p, Float x1, Float x2, Float x3) {
	return p.getD() - COE_CAL_ALPHA[ca][0] * p[0] * x1
			- COE_CAL_ALPHA[ca][1] * p[1] * x2
			- COE_CAL_ALPHA[ca][2] * p[2] * x3;
}

bool Comp_le_c(const Pair<Float, Float>& a, const Pair<Float, Float>& b) {
	return a.first * a.second <= b.first * b.second;
}
void coe_premutation( //
		const Pair<Float, Float>& o1, //
		const Pair<Float, Float>& o2, //
		const Pair<Float, Float>& o3, //
		int& p1, // out
		int& p2, // out
		int& p3) { //out
	sort(o1, o2, o3, Comp_le_c, p1, p2, p3);
}

Float forward_problem_3D(Float m1, Float m2, Float m3, Float a, Float c1,
		Float c2, Float c3) {
	Float V = 0;
	Float m12 = m1 * c1 + m2 * c2;
	Float max_a = m1 * c1 + m2 * c2 + m3 * c3;
	//Forward problem
	if (0 <= a && a < m1) { //case 1
		V = a * a * a / (6 * m1 * m2 * m3);
	} else if (m1 <= a && a < m2) { //case 2
		V = (a * a * c1 - a * m1 * c1 * c1) / (2 * m2 * m3)
				+ (m1 * m1 * c1 * c1 * c1) / (6 * m2 * m3);
	} else if (m2 <= a && a < MIN(m12, m3 * c3)) { //case 3
		V = (a * a * (3 * m1 * c1 + 3 * m2 * c2 - a)
				+ m1 * m1 * (m1 * c1 * c1 * c1 - 3 * a * c1 * c1)
				+ m2 * m2 * (m2 * c2 * c2 * c2 - 3 * a * c2 * c2))
				/ (6 * m1 * m2 * m3);
	} else if (MIN(m12, m3 * c3) <= a && a < 0.5 * max_a) {
		if (m3 * c3 < m12) {  //case 4
			V = (a * a * (3 * m12 + 3 * m3 * c3 - 2 * a)
					+ m1 * m1 * (m1 * c1 * c1 * c1 - 3 * a * c1 * c1)
					+ m2 * m2 * (m2 * c2 * c2 * c2 - 3 * a * c2 * c2)
					+ m3 * m3 * (m3 * c3 * c3 * c3 - 3 * a * c3 * c3))
					/ (6 * m1 * m2 * m3);
		} else { //case 5
			V = c1 * c2 * (2 * a - m12) / (2 * m3);
		}
	}
	return V;
}

Float inverse_problem_3D(Float m1, Float m2, Float m3, Float V, Float c1,
		Float c2, Float c3) {
	Float a = 0;
	Float mc12 = m1 * c1 + m2 * c2;
	Float V1 = m1 * m1 * c1 * c1 * c1 / (MAX(6.0 * m2 * m3, SMALL));
	Float V2 = V1 + (m2 * c1 - m1 * c1 * c1) / (2.0 * m3);
	Float V31 = (m3 * m3 * (3 * m1 * c1 + 3 * m2 * c2 - m3)
			+ m1 * m1 * (m1 * c1 * c1 * c1 - 3 * m3 * c1 * c1)
			+ m2 * m2 * (m2 * c2 * c2 * c2 - 3 * m3 * c2 * c2))
			/ (6 * m1 * m2 * m3);
	Float V32 = c1 * c2 * (m1 * c1 + m2 * c2) / (2 * m3);
	Float max_v = c1 * c2 * c3;
	//Forward problem
	if (0 <= V && V < V1) { //case 1
		a = pow(6.0 * m1 * m2 * m3 * V, 1.0 / 3.0);
	} else if (V1 <= V && V < V2) { //case 2
		a = (m1 * c1 * c1
				+ sqrt(
						8.0 * m2 * m3 * c1 * V
								- m1 * m1 * c1 * c1 * c1 * c1 / 3.0)) / 2.0
				* c1;
	} else if (V2 <= V && V < MIN(V31, V32)) { //case 3
		Float ta = -1;
		Float tb = 3.0 * mc12;
		Float tc = -3.0 * (m1 * m1 * c1 * c1 + m2 * m2 * c2 * c2);
		Float td = m1 * m1 * m1 * c1 * c1 * c1 + m2 * m2 * m2 * c2 * c2 * c2
				- 6.0 * m1 * m2 * m3 * V;
		Float res0 = 0;
		Float res2 = 0;
		solve_cubic_equation(ta, tb, tc, td, res0, a, res2);
	} else if (MIN(V31, V32) <= V && V < 0.5 * max_v) {
		if (V31 == MIN(V31, V32)) {  //case 4
			Float ta = -2.0;
			Float tb = 3.0 * (mc12 + m3 * c3);
			Float tc =
					-3.0
							* (m1 * m1 * c1 * c1 + m2 * m2 * c2 * c2
									+ m3 * m3 * c3 * c3);
			Float td = m1 * m1 * m1 * c1 * c1 * c1 + m2 * m2 * m2 * c2 * c2 * c2
					+ m3 * m3 * m3 * c3 * c3 * c3 - 6.0 * m1 * m2 * m3 * V;
			Float res0 = 0;
			Float res2 = 0;
			solve_cubic_equation(ta, tb, tc, td, res0, a, res2);
		} else { //case 5
			a = (mc12 + 2.0 * m3 * V) / 2 / c1 / c2;
		}
	}
	return a;
}

int middle_problem_3D(ListT<Point3D>& list, Float m1, Float m2, Float m3,
		Float a, Float c1, Float c2, Float c3) {
	int V = 0;
	Float m12 = m1 * c1 + m2 * c2;
	Float max_a = m1 * c1 + m2 * c2 + m3 * c3;
	//Forward problem
	list.clear();
	if (0 <= a && a < m1 * c1) { //case 1  ==>>> 3 points
		list.push_back(Point3D(a / m1, 0, 0)); //I
		list.push_back(Point3D(0, a / m2, 0)); //J
		list.push_back(Point3D(0, 0, a / m3)); //K
		V = 1;
	} else if (m1 * c1 <= a && a < m2 * c2) { //case 2  ==>>> 4 points
		list.push_back(Point3D(0, 0, a / m3));              //K
		list.push_back(Point3D(c1, 0, (a - m1 * c1) / m3)); //P
		list.push_back(Point3D(c1, (a - m1 * c1) / m2, 0)); //Q
		list.push_back(Point3D(0, a / m2, 0));              //J
		V = 2;
	} else if (m2 * c2 <= a && a < MIN(m12, m3 * c3)) { //case 3  ==>>> 5 points
		list.push_back(Point3D(0, 0, a / m3));              //K
		list.push_back(Point3D(c1, 0, (a - m1 * c1) / m3)); //P
		list.push_back(Point3D(c1, (a - m1 * c1) / m2, 0)); //Q
		list.push_back(Point3D((a - m2 * c2) / m1, c2, 0)); //R
		list.push_back(Point3D(0, c2, (a - m2 * c2) / m3)); //S
		V = 3;
	} else if (MIN(m12, m3 * c3) <= a && a <= 0.5 * max_a) {
		if (m3 * c3 < m12) {  //case 4  ==>>> 6 points
			list.push_back(Point3D((a - m3 * c3) / m1, 0, c3)); //T
			list.push_back(Point3D(c1, 0, (a - m1 * c1) / m3)); //P
			list.push_back(Point3D(c1, (a - m1 * c1) / m2, 0)); //Q
			list.push_back(Point3D((a - m2 * c2) / m1, c2, 0)); //R
			list.push_back(Point3D(0, c2, (a - m2 * c2) / m3)); //S
			list.push_back(Point3D(0, (a - m3 * c3) / m2, c3)); //U
			V = 4;
		} else { //case 5
			list.push_back(Point3D(0, 0, a / m3));              //K
			list.push_back(Point3D(c1, 0, (a - m1 * c1) / m3)); //P
			list.push_back(Point3D(c1, c2, (a - m1 * c1 - m2 * c2) / m3)); //W
			list.push_back(Point3D(0, c2, (a - m2 * c2) / m3)); //S
			V = 5;
		}
	}
	return V;
}

Float calFractionVolume(const Plane& p, Float x1, Float x2, Float x3) {
	array_3<Pair<Float, Float> > ori;
	ori[0].reconstruct(abs(p[0]), x1);
	ori[1].reconstruct(abs(p[1]), x2);
	ori[2].reconstruct(abs(p[2]), x3);
	array_3<int> tidx;
	coe_premutation(ori[0], ori[1], ori[2], tidx[0], tidx[1], tidx[2]);

	int ca = whichcase_3D_8(p);
	std::cout << "co case " << ca << std::endl;
	Float alpha = newAlpha_Forward(ca, p, x1, x2, x3);

	Float& m1 = ori[tidx[0]].first;
	Float& m2 = ori[tidx[1]].first;
	Float& m3 = ori[tidx[2]].first;
	Float& c1 = ori[tidx[0]].second;
	Float& c2 = ori[tidx[1]].second;
	Float& c3 = ori[tidx[2]].second;

	Float maxalpha = m1 * c1 + m2 * c2 + m3 * c3;
//The first special case
	if (alpha < 0) {
		return 0.0;
	}
	Float v_box = c1 * c2 * c3;
//The second special case
	if (alpha >= maxalpha) {
		return v_box;
	}
//normal case
//This step produce error;
	if (isZero(m1)) {
		m1 = SMALL;
	}
	if (isZero(m2)) {
		m2 = SMALL;
	}
	if (isZero(m3)) {
		m3 = SMALL;
	}
	if (alpha < maxalpha * 0.5) {
		return forward_problem_3D(m1, m2, m3, alpha, c1, c2, c3);
	} else {
		return -forward_problem_3D(m1, m2, m3, maxalpha - alpha, c1, c2, c3)
				+ v_box;
	}
}

void point_symmetry(Point3D& p, const Point3D& ps) {
	p.x = 2 * ps.x - p.x;
	p.y = 2 * ps.y - p.y;
	p.z = 2 * ps.z - p.z;
}

void point_symmetry(ListT<Point3D>& list, const Point3D& p) {
	for (ListT<Point3D>::iterator iter = list.begin(); iter != list.end();
			iter++) {
		point_symmetry((*iter), p);
	}
}

void point_transfer(int ca, Point3D& p, Float x1, Float x2, Float x3) {
	p.x = (COE_CAL_ALPHA[ca][0] == 1) ? (-p.x + x1) : p.x;
	p.y = (COE_CAL_ALPHA[ca][1] == 1) ? (-p.y + x2) : p.y;
	p.z = (COE_CAL_ALPHA[ca][2] == 1) ? (-p.z + x3) : p.z;
}

void point_transfer(int ca, ListT<Point3D>& list, Float x1, Float x2,
		Float x3) {
	for (ListT<Point3D>::iterator iter = list.begin(); iter != list.end();
			iter++) {
		point_transfer(ca, (*iter), x1, x2, x3);
	}
}

void shift_to_ori(ListT<Point3D>& list, const array_3<int>& tidx) {
	for (ListT<Point3D>::iterator iter = list.begin(); iter != list.end();
			iter++) {
		Point3D tmp;
		for (Point3D::size_type i = 0; i < tmp.size(); i++) {
			tmp[tidx[i]] = (*iter)[i];
		}
		for (Point3D::size_type i = 0; i < tmp.size(); i++) {
			(*iter)[i] = tmp[i];
		}
	}
}

ListT<Point3D> calInterctPoints(const Plane &p, Float x1, Float x2, Float x3) {
	array_3<Pair<Float, Float> > ori;
	ori[0].reconstruct(abs(p[0]), x1);
	ori[1].reconstruct(abs(p[1]), x2);
	ori[2].reconstruct(abs(p[2]), x3);
	array_3<int> tidx;
	coe_premutation(ori[0], ori[1], ori[2], tidx[0], tidx[1], tidx[2]);

	int ca = whichcase_3D_8(p);
	//std::cout << "co case " << ca << std::endl;
	Float alpha = newAlpha_Forward(ca, p, x1, x2, x3);

	Float& m1 = ori[tidx[0]].first;
	Float& m2 = ori[tidx[1]].first;
	Float& m3 = ori[tidx[2]].first;
	Float& c1 = ori[tidx[0]].second;
	Float& c2 = ori[tidx[1]].second;
	Float& c3 = ori[tidx[2]].second;

	//std::cout << "m1     " << "2     " << "3     " << std::endl;
	//std::cout << m1 << "  " << m2 << "  " << m3 << "  " << std::endl;
	//std::cout << c1 << "  " << c2 << "  " << c3 << "  " << std::endl;
	//std::cout << tidx[0] << "  " << tidx[1] << "  " << tidx[2] << "  "
	//		<< std::endl;
	//std::cout << "alpha " << alpha << std::endl;
	Float maxalpha = m1 * c1 + m2 * c2 + m3 * c3;
	//std::cout << " alphamax " << maxalpha << std::endl;
	ListT<Point3D> listres;
//The first special case
	if (alpha <= 0 || alpha >= maxalpha) {
		return listres;
	}
//normal case
//This step produce error;
	if (isZero(m1)) {
		m1 = SMALL;
	}
	if (isZero(m2)) {
		m2 = SMALL;
	}
	if (isZero(m3)) {
		m3 = SMALL;
	}
	//int casenum = -1;
	if (alpha < maxalpha * 0.5) {
		middle_problem_3D(listres, m1, m2, m3, alpha, c1, c2, c3);
	} else {
		middle_problem_3D(listres, m1, m2, m3, maxalpha - alpha, c1, c2, c3);
		point_symmetry(listres, Point3D(0.5 * c1, 0.5 * c2, 0.5 * c3));
	}
	shift_to_ori(listres, tidx);
	point_transfer(ca, listres, x1, x2, x3);
	//std::cout << "case " << casenum << std::endl;
	return listres;
}

void transfer_ListP(ListT<Point3D> &list, Float x1, Float x2, Float x3) {
	for (ListT<Point3D>::iterator iter = list.begin(); iter != list.end();
			iter++) {
		iter->transfer(x1, x2, x3);
	}
}

void scale_ListP(ListT<Point3D> &list, Float x1, Float x2, Float x3) {
	for (ListT<Point3D>::iterator iter = list.begin(); iter != list.end();
			iter++) {
		iter->scale(x1, x2, x3);
	}
}

//===============================================
//draw===========================================
void _draw_vtu_vof(std::string filename,  //filename
		ListT<pOCNode> list,      //tree
		LarusDef::size_type idc, //data color index
		LarusDef::size_type idx, //data index
		LarusDef::size_type idy, //data index
		LarusDef::size_type idz, //data index
		LarusDef::size_type ida  //data index
		) {
	LarusDef::size_type num_node = list.size();
	LarusDef::size_type num_points = 0;
	arrayListT<ListT<Point3D> > llp(num_node);
	LarusDef::size_type arridx = 0;
	for (ListT<pOCNode>::iterator iter = list.begin(); iter != list.end();
			iter++, arridx++) {
		Plane p((*iter)->data->aCenterData[idx],
				(*iter)->data->aCenterData[idy],
				(*iter)->data->aCenterData[idz],
				(*iter)->data->aCenterData[ida]);
		llp[arridx] = calInterctPoints(p, 1.0, 1.0, 1.0);
		scale_ListP(llp[arridx], (*iter)->cell->getDx() / 1.0,
				(*iter)->cell->getDy() / 1.0, (*iter)->cell->getDz() / 1.0);
		transfer_ListP(llp[arridx], (*iter)->cell->getMMM().x,
				(*iter)->cell->getMMM().y, (*iter)->cell->getMMM().z);
		num_points += llp[arridx].size();
	}
	//================================
	std::ofstream fs;
	open_file(fs, filename, 1);

	vtu_unstructured_grid_file_head(fs);
	fs << "<Piece NumberOfPoints=\"" << num_points << "\" NumberOfCells=\""
			<< num_node << "\">";
	fs << "<Points>";
	fs
			<< "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";
	for (LarusDef::size_type i = 0; i < llp.size(); i++) {
		drawtofile_vtu(fs, llp[i]);
	}
	fs << "</DataArray>";
	fs << "</Points>";
	fs << "<Cells>";
	fs << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
	for (LarusDef::size_type i = 0; i < num_points; i++) {
		fs << i << " ";
	}
	fs << "</DataArray>";
	fs << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">";
	LarusDef::size_type offset = 0;
	for (LarusDef::size_type i = 0; i < llp.size(); i++) {
		offset += llp[i].size();
		fs << offset << " ";
	}
	fs << "</DataArray>";
	fs << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">";
	for (LarusDef::size_type i = 0; i < llp.size(); i++) {
		fs << 7 << " ";
	}
	fs << "</DataArray>";
	fs << "</Cells>";
	fs << "</Piece>";
	vtu_unstructured_grid_file_end(fs);
	fs.close();
}
void draw_vtu_vof( // draw vof 3D==============
		std::string filename,  //filename
		pOCTree oct,      //tree
		LarusDef::size_type idc, //data color index
		LarusDef::size_type idx, //data index
		LarusDef::size_type idy, //data index
		LarusDef::size_type idz, //data index
		LarusDef::size_type ida  //data index
		) {
	//================================
	ListT<pOCNode> list;
	getListpNode_leaf_center_data_in_range(  //
			list, //as output
			oct, // ptree
			idc, //data index
			0.0, //min value
			1.0, //max value
			Range_oo);
	_draw_vtu_vof(filename, list, idc, idx, idy, idz, ida);
}

void draw_vtu_vof( // draw vof 3D==============
		std::string filename,  //filename
		Forest3D& forest,      //tree
		LarusDef::size_type idc, //data color index
		LarusDef::size_type idx, //data index
		LarusDef::size_type idy, //data index
		LarusDef::size_type idz, //data index
		LarusDef::size_type ida  //data index
		) {
	//================================
	ListT<pOCNode> list;
	getListpNode_leaf_center_data_in_range(  //
			list, //as output
			forest, // forest
			idc, //data index
			0.0, //min value
			1.0, //max value
			Range_oo);
	_draw_vtu_vof(filename, list, idc, idx, idy, idz, ida);
}

//Function for 2D================================
//===============================================

const int CASE_2D_8[2][2][2] = { //
		//
				{ { 2, 1 }, { 7, 8 } }, //
				{ { 3, 4 }, { 6, 5 } }  //
		};
/**
 * \brief   Calculate the direction of the line normal vetor
 *          <pre>
 *                          ^y
 *                +         |         +
 *                  +    3  |  2    +
 *                    +     |     +
 *                      +   |   +
 *                  4     + | +      1
 *          --+------+------+------+------+--> x
 *                        + | +
 *                  5   +   |   +    8
 *                    +     |     +
 *                  +     6 | 7     +
 *                +         |         +
 *          </pre>
 * \param   Float a a in ax+by=C
 * \param   Float b b in ax+by=C
 * \return  the case
 */
inline int whichCase8(Float a, Float b) {
	return CASE_2D_8[a >= 0 ? 0 : 1][b >= 0 ? 0 : 1][abs(b) >= abs(a) ? 0 : 1];
}

const int CASE_2D_4[2][2] = { { 1, 4 }, { 2, 3 } };
/**
 * \brief   Calculate the direction of the line normal vetor
 *          <pre>
 *                          ^y
 *                          |
 *                          |
 *                    2     |     1
 *                          |
 *                          |
 *          --+------+------+------+------+--> x
 *                          |
 *                          |
 *                    3     |     4
 *                          |
 *                          |
 *          </pre>
 * \param   Float a a in (a,b)
 * \param   Float b b in (a,b)
 * \return  the case
 */
inline int whichCase4(Float a, Float b) {
	return CASE_2D_4[a >= 0 ? 0 : 1][b >= 0 ? 0 : 1];
}

/**
 * \brief   known a,b in ax+by=alpha and C, calculate alpha \n
 *          no matter what a and b are, they will be change to abs(a) and abs(b)
 *          return alpha, ax+by=alpha, a>b>0;
 * \param   Float a a in ax+by=alpha
 * \param   Float b b in ax+by=alpha
 * \param   Float C the color function
 * \return  alpha
 */
Float calAlpha(Float a, Float b, Float C) {
	Float c1, c2, alpha;
	Float absa = (a < 0) ? (-a) : a;
	Float absb = (b < 0) ? (-b) : b;
	Float m, n;
	n = (absa >= absb) ? (absa) : absb;
	m = (absa <= absb) ? (absa) : absb;
	c1 = m / 2 / n;
	c2 = 1 - c1;
	if (C >= 0 && C <= c1) {
		alpha = sqrt(2 * C * m * n);
	} else if (C > c1 && C < c2) {
		alpha = (2 * C * n + m) / 2;
	} else { //(C>=c2 && C<=1)
		alpha = m + n - sqrt(2 * (1 - C) * m * n);
	}
	return alpha;
}

/**
 * \brief   known a,b in ax+by=alpha and C, calculate alpha \n
 *
 *          return alpha, ax+by=alpha;
 * \param   Float a a in ax+by=alpha
 * \param   Float b b in ax+by=alpha
 * \param   Float C the color function
 * \return  Line
 */
Line calLine(Float a, Float b, Float C) {
	if (a == 0.0 && b == 0.0) {
		a = 2 * SMALL;
	} else if (isEqual(a, 0.0) && isEqual(b, 0.0)) {
		a = 2 * SMALL * a / abs(a);
		b = 2 * SMALL * b / abs(b);
	}
	Line res;
	Float alpha;
	Float absa = abs(a);
	Float absb = abs(b);
	Float m, n;

	n = (absa >= absb) ? (absa) : absb;
	m = (absa <= absb) ? (absa) : absb;
	alpha = calAlpha(a, b, C);
	switch (whichCase8(a, b)) {
	case 1: {
		res.reconstruct(n, m, alpha);
		break;
	}
	case 2: {
		res.reconstruct(m, n, alpha);
		break;
	}
	case 3: {
		res.reconstruct(-m, n, alpha - m);
		break;
	}
	case 4: {
		res.reconstruct(-n, m, alpha - n);
		break;
	}
	case 5: {
		res.reconstruct(-n, -m, alpha - m - n);
		break;
	}
	case 6: {
		res.reconstruct(-m, -n, alpha - n - m);
		break;
	}
	case 7: {
		res.reconstruct(m, -n, alpha - n);
		break;
	}
	case 8: {
		res.reconstruct(n, -m, alpha - m); //
		break;
	}
	}
	return res;
}

/**
 * \brief   known a Line, calculate the fraction area in a rectangle \n
 *          <pre>
 *            y
 *            ^
 *            +
 *            |+
 *            | +
 *         -  +--+---+------+
 *         |  |***+         |
 *        c2  |****+Line    |
 *         |  |*****+       |
 *         -  +------+------+------+--> x
 *            0,0     +
 *            |------c1-----|
 *          </pre>
 * \param   Line
 * \param   Float the lenght x
 * \param   Float the lenght y
 * \return  Fractional area
 */
Float calFractionArea(const Line &l, Float c1, Float c2) {
	Float m1 = l.getA();
	Float m2 = l.getB();
	Float alpha = l.getC();

	Float am1 = abs(m1);
	Float am2 = abs(m2);
	Float aalpha;
	int lcase = whichCase4(m1, m2);
	switch (lcase) {
	case 1: {
		aalpha = alpha;
		break;
	}
	case 2: {
		aalpha = alpha - m1 * c1;
		break;
	}
	case 3: {
		aalpha = alpha - m1 * c1 - m2 * c2;
		break;
	}
	case 4: {
		aalpha = alpha - m2 * c2;
		break;
	}
	}
	assert(!(isZero(am1) && isZero(am2)));
	//This step produce error;
	if (isZero(am1)) {
		am1 = SMALL;
	}
	if (isZero(am2)) {
		am2 = SMALL;
	}
	//The first special case
	if (aalpha < 0) {
		return 0.0;
	}
	//The second special case
	if ((aalpha - am1 * c1) / am2 > c2 && (aalpha - am2 * c2) / am1 > c1) {
		return c1 * c2;
	}
	//The normal case
	Float r1, r2;
	if (StepFun(aalpha - c1 * am1) == 1) {
		r1 = ((aalpha - c1 * am1) / aalpha) * ((aalpha - c1 * am1) / aalpha);
	} else {
		r1 = 0;
	}
	if (StepFun(aalpha - c2 * am2) == 1) {
		r2 = ((aalpha - c2 * am2) / aalpha) * ((aalpha - c2 * am2) / aalpha);
	} else {
		r2 = 0;
	}
	return 0.5 * aalpha * aalpha / am1 / am2 * (1.0 - r1 - r2);
}

/**
 * \brief   known a Line, calculate the fraction area in a rectangle \n
 *          Original point could be set \n
 *          <pre>
 *            y
 *            ^
 *            +
 *            |+
 *            | +
 *         -  +--+---+------+
 *         |  |***+         |
 *        c2  |****+Line    |
 *         |  |*****+       |
 *         -  +------+------+------+--> x
 *          xo,yo     +
 *            |------c1-----|
 *          </pre>
 * \param   Line
 * \param   Float xo Orignal point x
 * \param   Float yo Orignal point y
 * \param   Float the lenght x
 * \param   Float the lenght y
 * \return  Fractional area
 */
Float calFractionArea(const Line &l, Float xo, Float yo, Float c1, Float c2) {
	assert(l.isExist());
	Float m = l.getA();
	Float n = l.getB();
	Float alpha = l.getC();
	Line tmpl(m, n, alpha - m * xo - n * yo);
	return calFractionArea(tmpl, c1, c2);
}

/**
 * \brief   known a Line, calculate the Intersect Point in a rectangle \n
 *
 * \param   Line
 * \param   Float the lenght x
 * \param   Float the lenght y
 * \return  Fractional area
 */
Segment2D calInterctPoints(const Line &l, Float c1, Float c2) {
	Float m1 = l.getA();
	Float m2 = l.getB();
	Float alpha = l.getC();

	Float am1, am2, aalpha;
	int lcase = whichCase4(m1, m2);
	switch (lcase) {
	case 1: {
		am1 = m1;
		am2 = m2;
		aalpha = alpha;
		break;
	}
	case 2: {
		am1 = -m1;
		am2 = m2;
		aalpha = alpha - m1 * c1;
		break;
	}
	case 3: {
		am1 = -m1;
		am2 = -m2;
		aalpha = alpha - m1 * c1 - m2 * c2;
		break;
	}
	case 4: {
		am1 = m1;
		am2 = -m2;
		aalpha = alpha - m2 * c2;
		break;
	}
	}

	assert(!(isZero(am1) && isZero(am2)));
	//This step produce error;
	if (isZero(am1)) {
		am1 = SMALL;
	}
	if (isZero(am2)) {
		am2 = SMALL;
	}

	//The first special case
	if (aalpha < 0) {
		return Segment2D();
	}
	//The second special case
	if ((aalpha - am1 * c1) / am2 > c2 && (aalpha - am2 * c2) / am1 > c1) {
		return Segment2D();
	}

	//The normal case
	Point2D p1, p2;
	if (StepFun(aalpha - c1 * am1) == 1) {
		//Point F
		p1.x = c1;
		p1.y = (aalpha - am1 * c1) / am2;
	} else {
		//Point E
		p1.x = aalpha / am1;
		p1.y = 0.0;
	}
	if (StepFun(aalpha - c2 * am2) == 1) {
		//Point G
		p2.x = (aalpha - am2 * c2) / am1;
		p2.y = c2;
	} else {
		//Point E
		p2.x = 0.0;
		p2.y = aalpha / am2;
	}
	if (p1 == p2) {
		std::cerr << " >! Intersect Point Equal. \n"
				<< "(VOF.h -->> Segment2D calInterctPoints(const Line2D &l, Float c1, Float c2)"
				<< std::endl;
		return Segment2D();
	} else {
		//The region on the leftside of segment equal ONE
		switch (lcase) {
		case 1: {
			return Segment2D(p2, p1);
			break;
		}
		case 2: {
			p1.x = -p1.x + c1;
			p2.x = -p2.x + c1;
			return Segment2D(p1, p2);
			break;
		}
		case 3: {
			p1.x = -p1.x + c1;
			p2.x = -p2.x + c1;
			p1.y = -p1.y + c2;
			p2.y = -p2.y + c2;
			return Segment2D(p2, p1);
			break;
		}
		case 4: {
			p1.y = -p1.y + c2;
			p2.y = -p2.y + c2;
			return Segment2D(p1, p2);
			break;
		}
		default: {
			std::cerr << "Return case error. "
					<< "(VOF.h -->> Segment2D calInterctPoints(const Line2D &l, Float c1, Float c2)"
					<< std::endl;

		}
		}
		return Segment2D();
	}
}

/**
 * \brief   Calculate the segment of gerris, this function will clear the list
 *
 * \param   List of segments
 *
 * \return  void
 */
void getListSegment( //
		ListT<Segment2D>& lseg, //list Segment
		Forest2D& forest,        //forest
		LarusDef::size_type idc, //data color index
		LarusDef::size_type idx, //data index
		LarusDef::size_type idy, //data index
		LarusDef::size_type ida  //data index
		) {
	lseg.clear();
	ListT<pQTNode> lnode;
	getListpNode_leaf_center_data_in_range(  //
			lnode, //as output
			forest, // forest
			idc,    //data index
			0.0, //min value
			1.0, //max value
			Range_oo);
	for (ListT<pQTNode>::iterator iter = lnode.begin(); iter != lnode.end();
			++iter) {
		Line l((*iter)->data->aCenterData[idx], (*iter)->data->aCenterData[idy],
				(*iter)->data->aCenterData[ida]);
		Segment2D seg = calInterctPoints(l, 1, 1);
		seg.scale((*iter)->cell->getDx(), (*iter)->cell->getDy());
		seg.transfer((*iter)->cell->getMM().x, (*iter)->cell->getMM().y);
		lseg.push_back(seg);
	}
}

Float calInterfaceLength(const ListT<Segment2D> &sg) {
	assert(sg.size() > 0);
	Float pe = 0.0;
	for (ListT<Segment2D>::const_iterator iter = sg.begin(); iter != sg.end(); ++iter) {
		Float a = iter->getLength();
		pe += a;
	}
	return pe;
}

} //end of namespace
