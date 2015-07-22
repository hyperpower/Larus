/************************
 //  \file   Segment.h
 //  \brief
 //
 //  \author zhou
 //  \date   23 mai 2014
 ***********************/
#ifndef SEGMENT_H_
#define SEGMENT_H_

#include "../TypeDef.h"
//#include "../Array.h"
//#include "../ArrayT.h"
#include "Line.h"
#include "Point.h"

namespace Larus
{

enum IntersectType
{
	NO_INTERSECT = 0x10000,
	INTERSECT = 0x1,
	START_1 = 0x10,
	END_1 = 0x20,
	START_2 = 0x100,
	END_2 = 0x200,
};

void parseIntersectType(int inter);

template<typename TYPE>
class Segment_2D: public array_2<Point_2D<TYPE> >
{
private:
	Line calLine(const Point_2D<TYPE> &, const Point_2D<TYPE> &) const;
	Float calLength(const Point_2D<TYPE> &a, const Point_2D<TYPE> &b) const;
public:
	typedef Point_2D<TYPE> value_type_point;
	typedef Point_2D<TYPE>& reference_point;
	typedef const Point_2D<TYPE>& const_reference_point;

	Segment_2D()
	{
	}
	Segment_2D(const Point_2D<TYPE>&, const Point_2D<TYPE>&);
	Segment_2D(TYPE ax, TYPE ay, TYPE bx, TYPE by);
	void reconstruct(const Point_2D<TYPE>&, const Point_2D<TYPE>&);
	void reconstruct(TYPE ax, TYPE ay, TYPE bx, TYPE by);

	bool operator==(const Segment_2D&) const;

	Line getLine() const;

	reference_point PS();
	const_reference_point PS() const;
	reference_point PE();
	const_reference_point PE() const;

	value_type_point PC() const;

	TYPE PSX() const;
	TYPE PEX() const;
	TYPE PSY() const;
	TYPE PEY() const;
	Float getLength() const;
	Float getSlope() const;
	value_type_point getCP() const;

	void scale(TYPE xfactor, TYPE yfactor);
	void transfer(TYPE dx, TYPE dy);
	bool isExist() const;
	void show() const;

	bool isVertical() const;
	bool isHorizontal() const;

	Float cross(const Point_2D<TYPE> &pt) const;
	bool isInbox(const Point_2D<TYPE> &pt) const;

	int onWhichside3(const Point_2D<TYPE> &pt) const;

	//compare ==================================
	bool is_gt_x(const TYPE& v) const;    //>=
	bool is_gt_y(const TYPE& v) const;    //>=
	bool is_ge_x(const TYPE& v) const;    //>=
	bool is_ge_y(const TYPE& v) const;    //>=

	bool is_lt_x(const TYPE& v) const;    //<
	bool is_lt_y(const TYPE& v) const;    //<
	bool is_le_x(const TYPE& v) const;    //<=
	bool is_le_y(const TYPE& v) const;    //<=
};
//=============
template<typename TYPE>
bool isBoxCross(  //
		const Segment_2D<TYPE>&, //
		const Segment_2D<TYPE>&);
template<typename TYPE>
int IntersectType( //
		const Segment_2D<TYPE> &s1, //
		const Segment_2D<TYPE> &s2);
template<typename TYPE>
bool isIntersect(   //
		const Segment_2D<TYPE> &s1, //
		const Segment_2D<TYPE> &s2);
template<typename TYPE>
bool isIntersect( //
		const Point_2D<TYPE> &s1s, //
		const Point_2D<TYPE> &s1e, //
		const Point_2D<TYPE> &s2s, //
		const Point_2D<TYPE> &s2e); //
template<typename TYPE>
Point2D calIntersect(const Segment_2D<TYPE> &s1, const Segment_2D<TYPE> &s2);
//=============
typedef Segment_2D<Float> Segment2D;

//=============
template<typename TYPE>
Segment_2D<TYPE>::Segment_2D(const Point_2D<TYPE>& a, const Point_2D<TYPE>& b) :
		array_2<Point_2D<TYPE> >(a, b)
{
	assert(a != b);
	if (a == b) {
		this->elems[0].x = 0.0;
		this->elems[0].y = 0.0;
		this->elems[1].x = 0.0;
		this->elems[1].y = 0.0;
	}
}

template<typename TYPE>
Segment_2D<TYPE>::Segment_2D(TYPE ax, TYPE ay, TYPE bx, TYPE by)
{
	Point_2D<TYPE> a(ax, ay);
	Point_2D<TYPE> b(bx, by);
	assert(a != b);
	this->elems[0] = a;
	this->elems[1] = b;
}

template<typename TYPE>
void Segment_2D<TYPE>::reconstruct(const Point_2D<TYPE>& a,
		const Point_2D<TYPE>& b)
{
	assert(a != b);
	this->elems[0] = a;
	this->elems[1] = b;
}

template<typename TYPE>
void Segment_2D<TYPE>::reconstruct(TYPE ax, TYPE ay, TYPE bx, TYPE by)
{
	Point_2D<TYPE> a(ax, ay);
	Point_2D<TYPE> b(bx, by);
	assert(a != b);
	this->elems[0] = a;
	this->elems[1] = b;
}

template<typename TYPE>
Line Segment_2D<TYPE>::calLine(const Point_2D<TYPE> &ps,
		const Point_2D<TYPE> &pe) const
{
	assert(ps != pe);
	Float a;
	Float b;
	Float c;
	a = -(pe.y - ps.y);
	b = pe.x - ps.x;
	c = pe.y * (pe.x - ps.x) - pe.x * (pe.y - ps.y);
	return Line(a, b, c);
}
template<typename TYPE>
bool Segment_2D<TYPE>::operator==(const Segment_2D<TYPE>& s) const
{
	return (s[0] == this->elems[0] && s[1] == this->elems[1]) ? true : false;
}

template<typename TYPE>
Line Segment_2D<TYPE>::getLine() const
{
	return calLine(this->elems[0], this->elems[1]);
}
template<typename TYPE>
Point_2D<TYPE>& Segment_2D<TYPE>::PS()
{
	return this->elems[0];
}
template<typename TYPE>
const Point_2D<TYPE>& Segment_2D<TYPE>::PS() const
{
	return this->elems[0];
}
template<typename TYPE>
Point_2D<TYPE>& Segment_2D<TYPE>::PE()
{
	return this->elems[1];
}
template<typename TYPE>
const Point_2D<TYPE>& Segment_2D<TYPE>::PE() const
{
	return this->elems[1];
}
template<typename TYPE>
typename Segment_2D<TYPE>::value_type_point Segment_2D<TYPE>::PC() const
{
	return value_type_point((PEX() + PSX()) *0.5, (PEY() + PSY()) * 0.5);
}

template<typename TYPE>
TYPE Segment_2D<TYPE>::PSX() const
{
	return this->elems[0].x;
}
template<typename TYPE>
TYPE Segment_2D<TYPE>::PEX() const
{
	return this->elems[1].x;
}
template<typename TYPE>
TYPE Segment_2D<TYPE>::PSY() const
{
	return this->elems[0].y;
}
template<typename TYPE>
TYPE Segment_2D<TYPE>::PEY() const
{
	return this->elems[1].y;
}
template<typename TYPE>
Float Segment_2D<TYPE>::calLength(const Point_2D<TYPE> &a,
		const Point_2D<TYPE> &b) const
{
	Float len = 0.0;
	len = sqrt(Float((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y)));
	return len;
}

template<typename TYPE>
Float Segment_2D<TYPE>::getLength() const
{
	return calLength(this->elems[0], this->elems[1]);
}

template<typename TYPE>
Float Segment_2D<TYPE>::getSlope() const
{
	return (PE().y - PS().y) / (PE().x - PS().x + SMALL);
}
template<typename TYPE>
typename Segment_2D<TYPE>::value_type_point Segment_2D<TYPE>::getCP() const
{
	return value_type_point((PEX() + PSX()) / 2, (PEY() + PSY()) / 2);
}
template<typename TYPE>
void Segment_2D<TYPE>::scale(TYPE xfactor, TYPE yfactor)
{
	this->elems[0].x = this->elems[0].x * xfactor;
	this->elems[0].y = this->elems[0].y * yfactor;
	this->elems[1].x = this->elems[1].x * xfactor;
	this->elems[1].y = this->elems[1].y * yfactor;
	if (this->elems[0] == this->elems[1]) {
		this->elems[0].x = 0.0;
		this->elems[0].y = 0.0;
		this->elems[1].x = 0.0;
		this->elems[1].y = 0.0;
	}
}
template<typename TYPE>
void Segment_2D<TYPE>::transfer(TYPE dx, TYPE dy)
{
	if (isExist()) {
		this->elems[0].x = this->elems[0].x + dx;
		this->elems[0].y = this->elems[0].y + dy;
		this->elems[1].x = this->elems[1].x + dx;
		this->elems[1].y = this->elems[1].y + dy;
	}
}
template<typename TYPE>
bool Segment_2D<TYPE>::isExist() const
{
	if (this->elems[0].x == 0.0 && this->elems[0].y == 0.0
			&& this->elems[1].x == 0.0 && this->elems[1].y == 0.0) {
		return false;
	} else {
		return true;
	}
}
template<typename TYPE>
void Segment_2D<TYPE>::show() const
{
	std::cout.precision(4);
	std::cout << "( " << PS().x << ", " << PS().y << " )--->(" << PE().x << ", "
			<< PE().y << " ) \n";
}
template<typename TYPE>
bool Segment_2D<TYPE>::isVertical() const
{
	ASSERT(isExist());
	return PS().x == PE().x;
}
template<typename TYPE>
bool Segment_2D<TYPE>::isHorizontal() const
{
	ASSERT(isExist());
	return PS().y == PE().y;
}

template<typename TYPE>
bool Segment_2D<TYPE>::isInbox(const Point_2D<TYPE> &pt) const
{
	ASSERT(isExist());
	if (isHorizontal()) {
		return (((PS().x <= pt.x) && (pt.x <= PE().x))
				|| ((PE().x <= pt.x) && (pt.x <= PS().x)));
	}
	if (isVertical()) {
		return (((PS().y <= pt.y) && (pt.y <= PE().y))
				|| ((PE().y <= pt.y) && (pt.y <= PS().y)));
	}
	return (((PS().x <= pt.x) && (pt.x <= PE().x))
			|| ((PE().x <= pt.x) && (pt.x <= PS().x)))
			&& (((PS().y <= pt.y) && (pt.y <= PE().y))
					|| ((PE().y <= pt.y) && (pt.y <= PS().y)));
}

template<typename TYPE>
Float Segment_2D<TYPE>::cross(const Point_2D<TYPE> &pt) const
{
	return CROSS(PE(), pt, PS());
}

template<typename TYPE>
int Segment_2D<TYPE>::onWhichside3(const Point_2D<TYPE> &pt) const
{
	Float rcro = cross(pt);
	if (rcro == 0.0) {
		return 0;
	} else if (rcro < 0) {
		return -1;
	} else {
		return 1;
	}
}
template<typename TYPE>
bool Segment_2D<TYPE>::is_gt_x(const TYPE& v) const
{
	return (this->PEX() > v && this->PSX() > v);
}
template<typename TYPE>
bool Segment_2D<TYPE>::is_gt_y(const TYPE& v) const
{
	return (this->PEY() > v && this->PSY() > v);
}

template<typename TYPE>
bool Segment_2D<TYPE>::is_ge_x(const TYPE& v) const
{
	return (this->PEX() >= v && this->PSX() >= v);
}
template<typename TYPE>
bool Segment_2D<TYPE>::is_ge_y(const TYPE& v) const
{
	return (this->PEY() >= v && this->PSY() >= v);
}

//==============================================
template<typename TYPE>
bool isBoxCross(const Segment_2D<TYPE>& s1, const Segment_2D<TYPE>& s2)
{
	return s1.isInbox(s2.PS()) || s1.isInbox(s2.PE()) || s2.isInbox(s1.PS())
			|| s2.isInbox(s1.PE());
}

template<typename TYPE>
int IntersectType(const Segment_2D<TYPE> &s1, const Segment_2D<TYPE> &s2)
{
	if (!isBoxCross(s1, s2)) {
		return NO_INTERSECT;
	} else {
		//step 1
		int s12s = s1.onWhichside3(s2.PS());
		int s12e = s1.onWhichside3(s2.PE());
		if (s12s == s12e) { //ignore the both equal to 0, overlap is not intersect
			return NO_INTERSECT;
		}
		int s21s = s2.onWhichside3(s1.PS());
		int s21e = s2.onWhichside3(s1.PE());
		if (s21s == s21e) { //ignore the both equal to 0, overlap is not intersect
			return NO_INTERSECT;
		}
		if ((s12s + s12e) == 0 && (s21s + s21e) == 0) {
			return INTERSECT;
		}
		int res = INTERSECT;
		if (s12s == 0)
			res = res | START_2;
		if (s12e == 0)
			res = res | END_2;
		if (s21s == 0)
			res = res | START_1;
		if (s21e == 0)
			res = res | END_1;
		return res;
	}
}

template<typename TYPE>
bool isIntersect(const Segment_2D<TYPE> &s1, const Segment_2D<TYPE> &s2)
{
	int type = IntersectType(s1, s2);
	return (type | INTERSECT) == type ? true : false;
}

template<typename TYPE>
bool isIntersect(const Point_2D<TYPE> &s1s, const Point_2D<TYPE> &s1e,
		const Point_2D<TYPE> &s2s, const Point_2D<TYPE> &s2e)
{
	Segment2D s1(s1s, s1e);
	Segment2D s2(s2s, s2e);
	return isIntersect(s1, s2);
}

//end 2D========================================

//Segment3D=====================================
template<typename TYPE>
class Segment_3D: public array_2<Point_3D<TYPE> >
{
private:
	Float calLength(const Point_3D<TYPE> &a, const Point_3D<TYPE> &b) const;
public:
	typedef Point_3D<TYPE> value_type_point;
	typedef Point_3D<TYPE>& reference_point;
	typedef const Point_3D<TYPE>& const_reference_point;

	Segment_3D()
	{
	}
	Segment_3D(const Point_3D<TYPE>&, const Point_3D<TYPE>&);
	Segment_3D(TYPE ax, TYPE ay, TYPE az, TYPE bx, TYPE by, TYPE bz);
	void reconstruct(const Point_3D<TYPE>&, const Point_3D<TYPE>&);
	void reconstruct(TYPE ax, TYPE ay, TYPE az, TYPE bx, TYPE by, TYPE bz);

	//Line getLine() const;

	reference_point PS();
	const_reference_point PS() const;
	reference_point PE();
	const_reference_point PE() const;

	Line3D getLine() const;

	Float getLength() const;
	//Float getSlope() const;

	void scale(const TYPE&, const TYPE&, const TYPE&);
	void transfer(const TYPE&, const TYPE&, const TYPE&);
	bool isExist() const;
	void show() const;

	bool isXY() const;
	bool isYZ() const;
	bool isZX() const;

	//bool isVertical() const;
	//bool isHorizontal() const;

	//int cross(const Point_3D<TYPE> &pt) const;
	bool isInbox(const Point_3D<TYPE> &pt) const;

	//int onWhichside3(const Point_3D<TYPE> &pt) const;
};

//=============
typedef Segment_3D<Float> Segment3D;

template<typename TYPE>
Segment_3D<TYPE>::Segment_3D(const Point_3D<TYPE>& a, const Point_3D<TYPE>& b) :
		array_2<Point_3D<TYPE> >(a, b)
{
	assert(a != b);
	if (a == b) {
		this->elems[0] = Point_3D<TYPE>();
		this->elems[1] = Point_3D<TYPE>();
	}
}
template<typename TYPE>
Segment_3D<TYPE>::Segment_3D(TYPE ax, TYPE ay, TYPE az, TYPE bx, TYPE by,
		TYPE bz)
{
	Point_3D<TYPE> a(ax, ay, az);
	Point_3D<TYPE> b(bx, by, bz);
	assert(a != b);
	this->elems[0] = a;
	this->elems[1] = b;
}
template<typename TYPE>
void Segment_3D<TYPE>::reconstruct(const Point_3D<TYPE>& a,
		const Point_3D<TYPE>& b)
{
	assert(a != b);
	this->elems[0] = a;
	this->elems[1] = b;
}
template<typename TYPE>
void Segment_3D<TYPE>::reconstruct(TYPE ax, TYPE ay, TYPE az, TYPE bx, TYPE by,
		TYPE bz)
{
	Point_3D<TYPE> a(ax, ay, az);
	Point_3D<TYPE> b(bx, by, bz);
	assert(a != b);
	this->elems[0] = a;
	this->elems[1] = b;
}
template<typename TYPE>
Point_3D<TYPE>& Segment_3D<TYPE>::PS()
{
	return this->elems[0];
}
template<typename TYPE>
const Point_3D<TYPE>& Segment_3D<TYPE>::PS() const
{
	return this->elems[0];
}
template<typename TYPE>
Point_3D<TYPE>& Segment_3D<TYPE>::PE()
{
	return this->elems[1];
}
template<typename TYPE>
const Point_3D<TYPE>& Segment_3D<TYPE>::PE() const
{
	return this->elems[1];
}

template<typename TYPE>
Float Segment_3D<TYPE>::calLength(const Point_3D<TYPE> &a,
		const Point_3D<TYPE> &b) const
{
	Float len = 0.0;
	len = sqrt(
			Float(
					(a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y)
							+ (a.z - b.z) * (a.z - b.z)));
	return len;
}

template<typename TYPE>
Float Segment_3D<TYPE>::getLength() const
{
	return calLength(this->elems[0], this->elems[1]);
}
template<typename TYPE>
void Segment_3D<TYPE>::scale(const TYPE& xfactor, const TYPE& yfactor,
		const TYPE& zfactor)
{
	this->elems[0].x = this->elems[0].x * xfactor;
	this->elems[0].y = this->elems[0].y * yfactor;
	this->elems[0].z = this->elems[0].z * zfactor;
	this->elems[1].x = this->elems[1].x * xfactor;
	this->elems[1].y = this->elems[1].y * yfactor;
	this->elems[1].z = this->elems[1].z * zfactor;
	if (this->elems[0] == this->elems[1]) {
		this->assign(Point_3D<TYPE>());
	}
}
template<typename TYPE>
void Segment_3D<TYPE>::transfer(const TYPE& dx, const TYPE& dy, const TYPE& dz)
{
	if (!isExist()) {
		return;
	}
	PE().transfer(dx, dy, dz);
	PS().transfer(dx, dy, dz);
}
template<typename TYPE>
bool Segment_3D<TYPE>::isExist() const
{
	Point_3D<TYPE> p;
	if (PS() == p && PE() == p) {
		return false;
	} else {
		return true;
	}
}
template<typename TYPE>
void Segment_3D<TYPE>::show() const
{
	std::cout << "( " << PS().x << ", " << PS().y << ", " << PS().z << " )--->("
			<< PE().x << ", " << PE().y << ", " << PE().z << " ) \n";
}

template<typename TYPE>
bool Segment_3D<TYPE>::isXY() const
{
	if (PS().z == PE().z) {
		return true;
	} else {
		return false;
	}
}
template<typename TYPE>
bool Segment_3D<TYPE>::isYZ() const
{
	if (PS().x == PE().x) {
		return true;
	} else {
		return false;
	}
}
template<typename TYPE>
bool Segment_3D<TYPE>::isZX() const
{
	if (PS().y == PE().y) {
		return true;
	} else {
		return false;
	}
}
template<typename TYPE>
bool Segment_3D<TYPE>::isInbox(const Point_3D<TYPE> &pt) const
{
	assert(isExist());
	if (isInRange_cc(PS().x, pt.x, PE().x) && isInRange_cc(PS().y, pt.y, PE().y)
			&& isInRange_cc(PS().z, pt.z, PE().z)) {
		return true;
	} else {
		return false;
	}
}

template<typename TYPE>
Line3D Segment_3D<TYPE>::getLine() const
{
	//TYPE==Float
	return Line3D(this->elems[0], this->elems[1]);
}

}	//end namespace

#endif /* SEGMENT_H_ */
