/************************
 //  \file   Cell.h
 //  \brief
 // 
 //  \author czhou
 //  \date   4 juin 2014 
 ***********************/
#ifndef CELL_H_
#define CELL_H_

#include "../TypeDef.h"
#include "../Geometry/Point.h"
#include "../Geometry/Polygon.h"
#include <string>
#include <stdio.h>
#include <fstream>

namespace Larus {

enum CellFace {
	ErrCellFace = -1,
	CF_IP = 2,
	CF_IM = 0,
	CF_JP = 1,
	CF_JM = 3,
	CF_KP = 4,
	CF_KM = 5,
	CF_E = 2,
	CF_W = 0,
	CF_N = 1,
	CF_S = 3,
	CF_F = 4,
	CF_B = 5,
};

enum CellPointLocation {
	ErrCellPointLocation = -1, eCPL_M = 0, eCPL_C = 1, eCPL_P = 2,
};

enum Cell2DDirection {
	ErrCell2DDirection = -1,
	CD_IP = 2,
	CD_IM = 0,
	CD_JP = 1,
	CD_JM = 3,
	CD_E = 2,
	CD_W = 0,
	CD_N = 1,
	CD_S = 3,
	CD_MP = 4,
	CD_NW = 4,
	CD_PP = 5,
	CD_NE = 5,
	CD_PM = 6,
	CD_SE = 6,
	CD_MM = 7,
	CD_SW = 7,
};


const static CellPointLocation CELL_POINT_ORDER[8][3] = {
		{eCPL_M,eCPL_P,eCPL_M},   //
		{eCPL_M,eCPL_M,eCPL_M},   //
		{eCPL_P,eCPL_P,eCPL_M},   //
		{eCPL_P,eCPL_M,eCPL_M},   //
		{eCPL_M,eCPL_P,eCPL_P},   //
		{eCPL_M,eCPL_M,eCPL_P},   //
		{eCPL_P,eCPL_P,eCPL_P},   //
		{eCPL_P,eCPL_M,eCPL_P},   //
};

class Cell2D: public array_2<Point2D> {
public:
	typedef Point2D Point;
	static const int Dim = 2;
	static const int NUM_VERTEXES = 4;  // make complier happy;
	static const int NUM_FACES = 4;
	Cell2D();
	Cell2D(const Point2D& mm, const Point2D& pp);
	Cell2D(Float mx, Float my, Float px, Float py, Float = 0.0, Float = 0.0);
	Point2D getCenterPoint() const;
	Float   getCenterPoint(const CSAxis&) const;
	Point2D getPC() const;
	Point2D getMC() const;
	Point2D getCP() const;
	Point2D getCM() const;
	Point2D getCC() const;
	Point2D getPP() const;
	Point2D getMP() const;
	Point2D getPM() const;
	Point2D getMM() const;
	Point2D getPoint(CellPointLocation x, CellPointLocation y,
			CellPointLocation = eCPL_M) const;

	inline Float getCPX() const;
	inline Float getCPY() const;
	inline Float getCPZ() const;

	Float getD(CSAxis axi) const;
	Float get(CSAxis axi, CellPointLocation l) const;
	Float getX(CellPointLocation l) const;
	Float getY(CellPointLocation l) const;
	Float getZ(CellPointLocation l) const;
	inline Float getDx() const;
	inline Float getDy() const;
	inline Float getDz() const;
	inline Float gethDx() const;
	inline Float gethDy() const;
	inline Float gethDz() const;
	inline Float getDvertex() const;
	Float area() const;
	inline Float volume() const;
	Float face(const CSAxis&) const;
	Polygon toPolygon() const;
	void show() const;
	bool isInCell(const Point2D &pt) const; // equal to Segment.h isInBox()
	bool isInCell(Point2D &pt);
	bool isInOnCell(const Point2D &pt) const; //
	bool isInOnCell(Point2D &pt);

	inline void transfer(Float dx, Float dy);

	void output_vertex_in_vtk_order(FILE *f) const;
	void output_vertex_in_vtk_order(std::ofstream& f) const;
	~Cell2D() {
	}
};

inline Float Cell2D::getDx() const {
	return this->elems[1].x - this->elems[0].x;
}

inline Float Cell2D::getDy() const {
	return this->elems[1].y - this->elems[0].y;
}

inline Float Cell2D::getDz() const {
	return 0.0;
}

inline Float Cell2D::gethDx() const {
	return getDx() * 0.5;
}
inline Float Cell2D::gethDy() const {
	return getDy() * 0.5;
}
inline Float Cell2D::gethDz() const {
	return 0.0;
}
inline Float Cell2D::getCPX() const {
	return this->elems[0].x + (this->elems[1].x - this->elems[0].x) * 0.5;
}
inline Float Cell2D::getCPY() const {
	return this->elems[0].y + (this->elems[1].y - this->elems[0].y) * 0.5;
}
inline Float Cell2D::getCPZ() const {
	return 0;
}

inline Float Cell2D::volume() const {
	return getDx() * getDy();
}
inline Float Cell2D::getDvertex() const {
	return sqrt(0.25 * getDx() * getDx() + 0.25 * getDy() * getDy());
}

inline void Cell2D::transfer(Float dx, Float dy) {
	this->elems[0][0] += dx;
	this->elems[0][1] += dy;
	this->elems[1][0] += dx;
	this->elems[1][1] += dy;
}

typedef Cell2D* pCell2D;

//                 v0______________________v2
//                  /|                    /|
//                 / |                   / |
//                /  |                  /  |
//               /___|_________________/   |
//            v4|    |                 |v6 |
//              |    |                 |   |
//              |    |                 |   |
//              |    |                 |   |
//              |    |_________________|___|
//              |   / v1               |   /v3
//              |  /                   |  /
//              | /                    | /
//              |/_____________________|/
//            v5                        v7

class Cell3D: public array_2<Point3D> {
public:
	typedef Point3D Point;
	static const int Dim = 3;
	static const int NUM_VERTEXES = 8;  //make complier happy;
	static const int NUM_FACES = 6;
	Cell3D();
	Cell3D(const Point3D& mmm, const Point3D& ppp);
	Cell3D(Float mx, Float my, Float mz, Float px, Float py, Float pz);
	Point3D getC() const;
	// ^
	// y
	// P MP CP PP
	// C MC CC PC
	// M MM CM PM
	//   M  C  P x>
	Float getX(CellPointLocation l) const;
	Float getY(CellPointLocation l) const;
	Float getZ(CellPointLocation l) const;
	Point3D& getMMM();
	const Point3D& getMMM() const;
	Point3D getPoint(CellPointLocation x, CellPointLocation y,
			CellPointLocation z) const;
	Point3D getVertex(LarusDef::size_type) const;
	Point3D getCenterPoint() const;
	inline Float getDx() const;
	inline Float getDy() const;
	inline Float getDz() const;
	inline Float gethDx() const;
	inline Float gethDy() const;
	inline Float gethDz() const;
	inline Float getDiagnalLength() const;
	inline Float gethDiagnalLength() const;

	inline Float getCPX() const {
		return this->elems[0].x + (this->elems[1].x - this->elems[0].x) * 0.5;
	}
	inline Float getCPY() const {
		return this->elems[0].y + (this->elems[1].y - this->elems[0].y) * 0.5;
	}
	inline Float getCPZ() const {
		return this->elems[0].z + (this->elems[1].z - this->elems[0].z) * 0.5;
	}
	inline Float volume() const;
	Float face(const CSAxis&) const;
	void show() const;
	bool isInCell(const Point3D &pt) const; // equal to Segment.h isInBox()
	bool isInOnCell(const Point3D &pt) const; //
	bool isInOnCell(Point3D &pt);
	void output_vertex_in_vtk_order(FILE *f) const;
	void output_vertex_in_vtk_order(std::ofstream& f) const;
	void draw_to_vtk(std::string filename) const;
	~Cell3D() {
	}
};

inline Float Cell3D::getDx() const {
	return this->elems[1].x - this->elems[0].x;
}
inline Float Cell3D::getDy() const {
	return this->elems[1].y - this->elems[0].y;
}
inline Float Cell3D::getDz() const {
	return this->elems[1].z - this->elems[0].z;
}
inline Float Cell3D::getDiagnalLength() const {
	return sqrt(getDx() * getDy() * getDz());
}
inline Float Cell3D::volume() const {
	return getDx() * getDy() * getDz();
}

inline Float Cell3D::gethDx() const {
	return getDx() * 0.5;
}
inline Float Cell3D::gethDy() const {
	return getDy() * 0.5;
}
inline Float Cell3D::gethDz() const {
	return getDz() * 0.5;
}
inline Float Cell3D::gethDiagnalLength() const {
	return getDiagnalLength() * 0.5;
}

} //this is the end of namespace

#endif /* CELL_H_ */
