/************************
 //  \file   Cell.cpp
 //  \brief
 // 
 //  \author czhou
 //  \date   4 juin 2014 
 ***********************/

#include "Cell.h"
#include "../Algebra/Arithmetic.h"
#include <iostream>
#include <fstream>
#include <stdio.h>

namespace Larus {

typedef Point2D (Cell2D::*pf_getPoint2D)() const;

pf_getPoint2D A_PF_GETPOINT2D[3][3] = { { &Cell2D::getMM, &Cell2D::getMC,
		&Cell2D::getMP }, { &Cell2D::getCM, &Cell2D::getCC, &Cell2D::getCP }, {
		&Cell2D::getPM, &Cell2D::getPC, &Cell2D::getPP } };

const CellPointLocation A_CELL3D_ORDER_VTK[8][3] = {
//VTK_VOXEL (=11)
		{ eCPL_M, eCPL_M, eCPL_M }, //
		{ eCPL_P, eCPL_M, eCPL_M }, //
		{ eCPL_M, eCPL_P, eCPL_M }, //
		{ eCPL_P, eCPL_P, eCPL_M }, //
		{ eCPL_M, eCPL_M, eCPL_P }, //
		{ eCPL_P, eCPL_M, eCPL_P }, //
		{ eCPL_M, eCPL_P, eCPL_P }, //
		{ eCPL_P, eCPL_P, eCPL_P } };

const CellPointLocation A_CELL2D_ORDER_VTK[4][2] = {
//VTK_PIXEL (=8)
		{ eCPL_M, eCPL_M }, //
		{ eCPL_P, eCPL_M }, //
		{ eCPL_M, eCPL_P }, //
		{ eCPL_P, eCPL_P } //
};

const CellPointLocation VERTEX_ORD[8][3] = {
//          x        y      z
		{ eCPL_P, eCPL_M, eCPL_M }, //
		{ eCPL_M, eCPL_M, eCPL_M }, //
		{ eCPL_P, eCPL_P, eCPL_M }, //
		{ eCPL_P, eCPL_M, eCPL_M }, //
		{ eCPL_P, eCPL_M, eCPL_P }, //
		{ eCPL_M, eCPL_M, eCPL_P }, //
		{ eCPL_P, eCPL_P, eCPL_P }, //
		{ eCPL_P, eCPL_M, eCPL_P } };
//===========================================

//begin Cell2D===============================
Cell2D::Cell2D() {
}

Cell2D::Cell2D(const Point2D& mm, const Point2D& pp) {
	assert(mm.x < pp.x && mm.y < pp.y);
	this->elems[0] = mm;
	this->elems[1] = pp;
}

Cell2D::Cell2D(Float mx, Float my, Float px, Float py, Float zo, Float zoo) {
	//zo and zoo just make compiler happy!
	assert(mx < px && my < py);
	this->elems[0].x = mx;
	this->elems[0].y = my;
	this->elems[1].x = px;
	this->elems[1].y = py;
}

Point2D Cell2D::getCenterPoint() const {
	return Point2D(
			this->elems[0].x + (this->elems[1].x - this->elems[0].x) * 0.5,
			this->elems[0].y + (this->elems[1].y - this->elems[0].y) * 0.5);
}
Float Cell2D::getCenterPoint(const CSAxis& aix) const {
	ASSERT(aix == CSAxis_X || aix == CSAxis_Y);
	return this->elems[0][aix]
			+ (this->elems[1][aix] - this->elems[0][aix]) * 0.5;

}
Point2D Cell2D::getPC() const {
	return Point2D(this->elems[1].x,
			this->elems[0].y + (this->elems[1].y - this->elems[0].y) * 0.5);
}
Point2D Cell2D::getMC() const {
	return Point2D(this->elems[0].x,
			this->elems[0].y + (this->elems[1].y - this->elems[0].y) * 0.5);
}
Point2D Cell2D::getCP() const {
	return Point2D(
			this->elems[0].x + (this->elems[1].x - this->elems[0].x) * 0.5,
			this->elems[1].y);
}

Point2D Cell2D::getCM() const {
	return Point2D(
			this->elems[0].x + (this->elems[1].x - this->elems[0].x) * 0.5,
			this->elems[0].y);
}
Point2D Cell2D::getCC() const {
	return Point2D(
			this->elems[0].x + (this->elems[1].x - this->elems[0].x) * 0.5,
			this->elems[0].y + (this->elems[1].y - this->elems[0].y) * 0.5);
}
Point2D Cell2D::getPP() const {
	return this->elems[1];
}
Point2D Cell2D::getMP() const {
	return Point2D(this->elems[0].x, this->elems[1].y);
}
Point2D Cell2D::getPM() const {
	return Point2D(this->elems[1].x, this->elems[0].y);
}
Point2D Cell2D::getMM() const {
	return this->elems[0];
}

Point2D Cell2D::getPoint(CellPointLocation x, CellPointLocation y,
		CellPointLocation z) const {
	return (this->*A_PF_GETPOINT2D[x][y])();
}

Float Cell2D::getD(CSAxis axi) const {
	return this->elems[1][axi] - this->elems[0][axi];
}

Float Cell2D::get(CSAxis axi, CellPointLocation l) const {
	assert(l != ErrCellPointLocation);
	Float xv = 0;
	switch (l) {
	case eCPL_M: {
		xv = this->elems[0][axi];
		break;
	}
	case eCPL_C: {
		xv = this->elems[0][axi] + getD(axi) * 0.5;
		break;
	}
	case eCPL_P: {
		xv = this->elems[1][axi];
		break;
	}
	default: {
		break;
	}
	}
	return xv;
}
Float Cell2D::getX(CellPointLocation l) const {
	return get(CSAxis_X, l);
}
Float Cell2D::getY(CellPointLocation l) const {
	return get(CSAxis_Y, l);
}
Float Cell2D::getZ(CellPointLocation l) const {
	return 0.0;
}

Float Cell2D::area() const {
	return getDx() * getDy();
}

Float Cell2D::face(const CSAxis& axis) const {
	ASSERT(axis != ErrCSAxis && axis != CSAxis_Z);
	if (axis == CSAxis_X) {
		return this->elems[1][CSAxis_Y] - this->elems[0][CSAxis_Y];
	}
	if (axis == CSAxis_Y) {
		return this->elems[1][CSAxis_X] - this->elems[0][CSAxis_X];
	}
	return 0.0;
}

Polygon Cell2D::toPolygon() const {
	arrayListT<Point2D> ap(4);
	ap[0] = getMM();
	ap[1] = getPM();
	ap[2] = getPP();
	ap[3] = getMP();
	return Polygon(ap);
}

void Cell2D::show() const {
	std::cout << std::scientific << getMP().x << "____ " << getCP().x << "____ "
			<< this->elems[1].x << std::endl;
	std::cout << std::scientific << getMP().y << "     " << getCP().y << "     "
			<< this->elems[1].y << std::endl;
	std::cout << "     |                               |     " << std::endl;
	std::cout << "     |                               |     " << std::endl;
	std::cout << getMC().x << "    " << getCenterPoint().x << "    "
			<< getPC().x << std::endl;
	std::cout << getMC().y << "    " << getCenterPoint().y << "    "
			<< getPC().y << std::endl;
	std::cout << "     |                               |     " << std::endl;
	std::cout << "     |                               |     " << std::endl;
	std::cout << this->elems[0].x << "____ " << getCM().x << "____ "
			<< getPM().x << std::endl;
	std::cout << this->elems[0].y << "     " << getCM().y << "     "
			<< getPM().y << std::endl;
}

//If the point is on the boundary of Cell --> false
bool Cell2D::isInCell(const Point2D &pt) const {
	return (isInRange(this->elems[0].x, pt.x, this->elems[1].x, Range_oo)
			&& isInRange(this->elems[0].y, pt.y, this->elems[1].y, Range_oo));
}

bool Cell2D::isInCell(Point2D &pt) {
	return (isInRange(this->elems[0].x, pt.x, this->elems[1].x, Range_oo)
			&& isInRange(this->elems[0].y, pt.y, this->elems[1].y, Range_oo));
}

bool Cell2D::isInOnCell(const Point2D &pt) const { //
	return (isInRange(this->elems[0].x, pt.x, this->elems[1].x, Range_cc)
			&& isInRange(this->elems[0].y, pt.y, this->elems[1].y, Range_cc));
}
bool Cell2D::isInOnCell(Point2D &pt) {
	return (isInRange(this->elems[0].x, pt.x, this->elems[1].x, Range_cc)
			&& isInRange(this->elems[0].y, pt.y, this->elems[1].y, Range_cc));
}

void Cell2D::output_vertex_in_vtk_order(FILE *f) const {
	for (int i = 0; i < NUM_VERTEXES; i++) {
		fprintf(f, "%f %f %f \n",
				getPoint(A_CELL2D_ORDER_VTK[i][0], A_CELL2D_ORDER_VTK[i][1]).x,
				getPoint(A_CELL2D_ORDER_VTK[i][0], A_CELL2D_ORDER_VTK[i][1]).y,
				0.0);
	}
}
void Cell2D::output_vertex_in_vtk_order(std::ofstream& fs) const {
	for (int i = 0; i < NUM_VERTEXES; i++) {
		fs << getPoint(A_CELL2D_ORDER_VTK[i][0], A_CELL2D_ORDER_VTK[i][1]).x
				<< " "
				<< getPoint(A_CELL2D_ORDER_VTK[i][0], A_CELL2D_ORDER_VTK[i][1]).y
				<< " " << 0.0 << " ";
	}
}

//this is the end of class Cell2D==============
//=============================================

Cell3D::Cell3D() {
}

Cell3D::Cell3D(const Point3D& mmm, const Point3D& ppp) :
		array_2<Point3D>(mmm, ppp) {
	assert(mmm.x < ppp.x && mmm.y < ppp.y && mmm.z < ppp.z);
	//this->elems[0] = mmm;
	//this->elems[1] = ppp;
}
Cell3D::Cell3D(Float mx, Float my, Float mz, Float px, Float py, Float pz) {
	assert(mx < px && my < py && mz < pz);
	this->elems[0] = Point3D(mx, my, mz);
	this->elems[1] = Point3D(px, py, pz);
}

Float Cell3D::getX(CellPointLocation l) const {
	assert(l != ErrCellPointLocation);
	Float xv = 0;
	switch (l) {
	case eCPL_M: {
		xv = this->elems[0].x;
		break;
	}
	case eCPL_C: {
		xv = this->elems[0].x + getDx() * 0.5;
		break;
	}
	case eCPL_P: {
		xv = this->elems[1].x;
		break;
	}
	default: {
		break;
	}
	}
	return xv;
}

Float Cell3D::getY(CellPointLocation l) const {
	assert(l != ErrCellPointLocation);
	Float yv = 0;
	switch (l) {
	case eCPL_M: {
		yv = this->elems[0].y;
		break;
	}
	case eCPL_C: {
		yv = this->elems[0].y + getDy() * 0.5;
		break;
	}
	case eCPL_P: {
		yv = this->elems[1].y;
		break;
	}
	default: {
		break;
	}
	}
	return yv;
}

Float Cell3D::getZ(CellPointLocation l) const {
	Float zv = 0;
	switch (l) {
	case eCPL_M: {
		zv = this->elems[0].z;
		break;
	}
	case eCPL_C: {
		zv = this->elems[0].z + getDz() * 0.5;
		break;
	}
	case eCPL_P: {
		zv = this->elems[1].z;
		break;
	}
	default: {
		break;
	}
	}
	return zv;
}

Point3D& Cell3D::getMMM() {
	return this->elems[0];
}
const Point3D& Cell3D::getMMM() const {
	return this->elems[0];
}

Float Cell3D::face(const CSAxis& axis) const {
	ASSERT(axis != ErrCSAxis);
	if (axis == CSAxis_X) {
		return getDy() * getDz();
	}
	if (axis == CSAxis_Y) {
		return getDx() * getDz();
	}
	if (axis == CSAxis_Z) {
		return getDy() * getDx();
	}
	return 0.0; //make complier happy
}

Point3D Cell3D::getPoint(CellPointLocation x, CellPointLocation y,
		CellPointLocation z) const {
	return Point3D(getX(x), getY(y), getZ(z));
}

Point3D Cell3D::getVertex(LarusDef::size_type idx) const {
	ASSERT(idx >= 0 && idx < 8);
	return getPoint(VERTEX_ORD[idx][0], VERTEX_ORD[idx][1], VERTEX_ORD[idx][2]);

}
Point3D Cell3D::getCenterPoint() const {
	return getPoint(eCPL_C, eCPL_C, eCPL_C);
}

void Cell3D::show() const {
	std::cout << std::scientific;
	std::cout << "Cell3D \n";
	std::cout << "Z = M \n";
	std::cout << "Zm = " << this->elems[0].z << " Zp = " << this->elems[1].z
			<< "\n";
	std::cout << getPoint(eCPL_M, eCPL_P, eCPL_M).x << "____ "
			<< getPoint(eCPL_C, eCPL_P, eCPL_M).x << "____ " << this->elems[1].x
			<< std::endl;
	std::cout << getPoint(eCPL_M, eCPL_P, eCPL_M).y << "     "
			<< getPoint(eCPL_C, eCPL_P, eCPL_M).y << "     " << this->elems[1].y
			<< std::endl;
	std::cout << "     |                               |     " << std::endl;
	std::cout << "     |                               |     " << std::endl;
	std::cout << getPoint(eCPL_M, eCPL_C, eCPL_M).x << "     "
			<< getPoint(eCPL_C, eCPL_C, eCPL_M).x << "     "
			<< getPoint(eCPL_P, eCPL_C, eCPL_M).x << std::endl;
	std::cout << getPoint(eCPL_M, eCPL_C, eCPL_M).y << "     "
			<< getPoint(eCPL_C, eCPL_C, eCPL_M).y << "     "
			<< getPoint(eCPL_P, eCPL_C, eCPL_M).y << std::endl;
	std::cout << "     |                               |     " << std::endl;
	std::cout << "     |                               |     " << std::endl;
	std::cout << this->elems[0].x << "____ "
			<< getPoint(eCPL_C, eCPL_M, eCPL_M).x << "____ "
			<< getPoint(eCPL_P, eCPL_M, eCPL_M).x << std::endl;
	std::cout << this->elems[0].y << "     "
			<< getPoint(eCPL_C, eCPL_M, eCPL_M).y << "     "
			<< getPoint(eCPL_P, eCPL_M, eCPL_M).y << std::endl;
// pppp========================
}

void Cell3D::output_vertex_in_vtk_order(FILE *f) const {
	for (int i = 0; i < NUM_VERTEXES; i++) {
		fprintf(f, "%f %f %f \n",
				getPoint(A_CELL3D_ORDER_VTK[i][0], A_CELL3D_ORDER_VTK[i][1],
						A_CELL3D_ORDER_VTK[i][2]).x,
				getPoint(A_CELL3D_ORDER_VTK[i][0], A_CELL3D_ORDER_VTK[i][1],
						A_CELL3D_ORDER_VTK[i][2]).y,
				getPoint(A_CELL3D_ORDER_VTK[i][0], A_CELL3D_ORDER_VTK[i][1],
						A_CELL3D_ORDER_VTK[i][2]).z);
	}
}

void Cell3D::output_vertex_in_vtk_order(std::ofstream& fs) const {
	for (int i = 0; i < NUM_VERTEXES; i++) {
		fs
				<< getPoint(A_CELL3D_ORDER_VTK[i][0], A_CELL3D_ORDER_VTK[i][1],
						A_CELL3D_ORDER_VTK[i][2]).x << " "
				<< getPoint(A_CELL3D_ORDER_VTK[i][0], A_CELL3D_ORDER_VTK[i][1],
						A_CELL3D_ORDER_VTK[i][2]).y << " "
				<< getPoint(A_CELL3D_ORDER_VTK[i][0], A_CELL3D_ORDER_VTK[i][1],
						A_CELL3D_ORDER_VTK[i][2]).z << " ";
	}
}

void Cell3D::draw_to_vtk(std::string filename) const {
	FILE *data;
	data = fopen(filename.c_str(), "w");
	fprintf(data, "# vtk DataFile Version 3.0\n");
	fprintf(data, "Gird output\n");
	fprintf(data, "ASCII\n");
	fprintf(data, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(data, "POINTS %d float\n", 8);
	output_vertex_in_vtk_order(data);
	fprintf(data, "\n");
	fprintf(data, "CELLS %d %d \n", 1, 9);
	fprintf(data, "%d ", 8);
	for (int i = 0; i < 8; i++) {
		fprintf(data, "%d ", i);
	}
	fprintf(data, "\n\n");
	fprintf(data, "CELL_TYPES %d\n", 1);
	fprintf(data, "%d \n", 11);

	fclose(data);
}

//If the point is on the boundary of Cell --> false
bool Cell3D::isInCell(const Point3D &pt) const { // equal to Segment.h isInBox()
	return (((this->elems[0].x < pt.x) && (pt.x < this->elems[1].x))
			|| ((this->elems[1].x < pt.x) && (pt.x < this->elems[0].x)))
			&& (((this->elems[0].y < pt.y) && (pt.y < this->elems[1].y))
					|| ((this->elems[1].y < pt.y) && (pt.y < this->elems[0].y)))
			&& (((this->elems[0].z < pt.z) && (pt.z < this->elems[1].z))
					|| ((this->elems[1].z < pt.z) && (pt.z < this->elems[0].z)));
}

bool Cell3D::isInOnCell(const Point3D &pt) const { //
	return (isInRange(this->elems[0].x, pt.x, this->elems[1].x, Range_cc)
			&& isInRange(this->elems[0].y, pt.y, this->elems[1].y, Range_cc)
			&& isInRange(this->elems[0].z, pt.z, this->elems[1].z, Range_cc));
}
bool Cell3D::isInOnCell(Point3D &pt) {
	return (isInRange(this->elems[0].x, pt.x, this->elems[1].x, Range_cc)
			&& isInRange(this->elems[0].y, pt.y, this->elems[1].y, Range_cc)
			&& isInRange(this->elems[0].z, pt.z, this->elems[1].z, Range_cc));
}

} //end of namespace
