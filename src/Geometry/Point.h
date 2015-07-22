/************************
 //  \file   Point.h
 //  \brief
 // 
 //  \author zhou
 //  \date   19 janv. 2014 
 ***********************/
#ifndef POINT_H_
#define POINT_H_

//#define CIR_PREDICATE

#include "../TypeDef.h"
#include "../Utility/Array.h"
#include "GeoDef.h"
#include <iostream>
#include <iomanip>


#include <math.h>

namespace Larus
{

//====================================================
//2D template ========================================
template<typename TYPE>
class Point_2D: public ObjectBase
{
public:
	typedef point_2d_tag self_tag;
	typedef LarusDef::size_type size_type;
	typedef TYPE value_type;
	typedef TYPE& reference;
	typedef const TYPE& const_reference;

	TYPE x;
	TYPE y;
	//constructor
	Point_2D()
	{
		x = 0;    //make sure initial to zero
		y = 0;
	}
	Point_2D(
			const TYPE& a,
			const TYPE& b,
			const TYPE& c = 0)
	{
		x = a;
		y = b;
	}
	Point_2D(
			const Point_2D& a)
	{
		x = a.x;
		y = a.y;
	}
	Point_2D& operator=(
			const Point_2D& rhs)
	{
		if (&rhs == this)
		{
			return *this;
		}
		else
		{
			x = rhs.x;
			y = rhs.y;
		}
		return *this;
	}
	const_reference operator[](
			size_type index) const
	{  //overload []
		ASSERT(index >= 0 && index < 2);
		return index == 0 ? x : y;
	}
	reference operator[](
			size_type index)
	{
		ASSERT(index >= 0 && index < 2);
		return index == 0 ? x : y;
	}
	void reconstruct(
			const TYPE& a,
			const TYPE& b)
	{
		x = a;
		y = b;
	}
	bool operator==(
			const Point_2D &a) const
	{
		return (x == a.x && y == a.y) ? true : false;
	}
	bool operator!=(
			const Point_2D &a) const
	{
		return !((x == a.x && y == a.y) ? true : false);
	}
	void transfer(
			const TYPE&dx,
			const TYPE&dy)
	{
		x = x + dx;
		y = y + dy;
	}
	void scale(
			const TYPE&dx,
			const TYPE&dy)
	{
		x = x * dx;
		y = y * dy;
	}
	inline size_type size() const
	{
		return 2;
	}
	void show() const
	{
		std::cout.precision(5);
		std::cout << "( "<< std::setw(10) <<  x << " , " << std::setw(10) << y << " )\n";
	}
}
;

typedef Point_2D<Float> Point2D;
//====================================================

//===============================================
template<typename TYPE>
class Point_3D: public ObjectBase
{
public:
	typedef point_3d_tag self_tag;
	typedef LarusDef::size_type size_type;
	typedef TYPE value_type;
	typedef TYPE& reference;
	typedef const TYPE& const_reference;

	TYPE x;  //
	TYPE y;  //
	TYPE z;  //
	//constructor
	Point_3D()
	{
		x = 0;   //make sure initial to zero
		y = 0;
		z = 0;
	}

	Point_3D(
			const TYPE& a,
			const TYPE& b,
			const TYPE& c)
	{
		x = a;
		y = b;
		z = c;
	}
	Point_3D(
			const Point_3D& a)
	{
		x = a.x;
		y = a.y;
		z = a.z;
	}

	Point_3D& operator=(
			const Point_3D& rhs)
	{
		if (&rhs == this)
		{
			return *this;
		}
		else
		{
			x = rhs.x;
			y = rhs.y;
			z = rhs.z;
		}
		return *this;
	}
	const_reference operator[](
			size_type index) const
	{  //overload []
		ASSERT(index >= 0 && index < 3);
		return index == 0 ? x : (index == 1 ? y : z);
	}
	reference operator[](
			size_type index)
	{
		ASSERT(index >= 0 && index < 3);
		return index == 0 ? x : (index == 1 ? y : z);
	}
	void reconstruct(
			const TYPE& a,
			const TYPE& b,
			const TYPE& c)
	{
		x = a;
		y = b;
		z = c;
	}

	bool operator==(
			const Point_3D &a) const
	{
		return (x == a.x && y == a.y && z == a.z) ? true : false;
	}
	bool operator!=(
			const Point_3D &a) const
	{
		return !((x == a.x && y == a.y && z == a.z) ? true : false);
	}
	void show() const
	{
		std::cout << std::scientific << "( " << x << " , " << y << " , " << z
				<< " )\n";
	}
	void transfer(
			const TYPE&dx,
			const TYPE&dy,
			const TYPE&dz)
	{
		x = x + dx;
		y = y + dy;
		z = z + dz;
	}
	void scale(
			const TYPE&dx,
			const TYPE&dy,
			const TYPE&dz)
	{
		x = x * dx;
		y = y * dy;
		z = z * dz;
	}
	inline size_type size() const
	{
		return 3;
	}

	Point_2D<TYPE> projectXY() const;
	Point_2D<TYPE> projectYZ() const;
	Point_2D<TYPE> projectZX() const;
};

typedef Point_3D<Float> Point3D;
//function point3D===============================
template <class TYPE>
inline TYPE calDistance(
		const TYPE &x1,
		const TYPE &y1,
		const TYPE &x2,
		const TYPE &y2){
	return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}


Float calDistance(
		const Point3D &p1,
		const Point3D &p2);
Float calDistance(
		const Point2D &p1,
		const Point2D &p2);



Point2D calMid(
		const Point2D &p1,
		const Point2D &p2);
Point3D calMid(
		const Point3D &p1,
		const Point3D &p2);
bool isEqual(
		const Point2D &p1,
		const Point2D &p2);

//cross<TYPE>
Float CROSS(
		const Point_2D<Float> &v1,
		const Point_2D<Float> &v2,
		const Point_2D<Float> &v3);
Float cross(
		const Point_2D<int> &v1,
		const Point_2D<int> &v2,
		const Point_2D<int> &v3);

Float CROSS(
		const Point_3D<Float> &v1,
		const Point_3D<Float> &v2,
		const Point_3D<Float> &v3,
		const Point_3D<Float> &v4);
Float cross(
		const Point_3D<int> &v1,
		const Point_3D<int> &v2,
		const Point_3D<int> &v3,
		const Point_3D<int> &v4);
//Point T ====================================
template<typename TYPE, int DIM>
class PointT: public array<TYPE, DIM>
{
public:
	//typedef point__tag self_tag;
	typedef LarusDef::size_type size_type;
	typedef TYPE value_type;
	typedef TYPE& reference;
	typedef const TYPE& const_reference;

	//constructor
	PointT() :
			array<TYPE, DIM>()
	{
	}

	PointT(
			const TYPE& a,
			const TYPE& b,
			const TYPE& c = 0) :
			array<TYPE, DIM>()
	{
		this->elems[0] = a;
		this->elems[1] = b;
		if (DIM == 3)
		{
			this->elems[2] = c;
		}
	}

	const_reference x() const
	{
		return this->elems[0];
	}

	reference x()
	{
		return this->elems[0];
	}

	const_reference y() const
	{
		return this->elems[1];
	}

	reference y()
	{
		return this->elems[1];
	}

	const_reference z() const
	{
		ASSERT(DIM == 3);
		return this->elems[2];
	}

	reference z()
	{
		ASSERT(DIM == 3);
		return this->elems[2];
	}

	void reconstruct(
			const TYPE& a,
			const TYPE& b,
			const TYPE& c = 0)
	{
		this->elems[0] = a;
		this->elems[1] = b;
		if (DIM == 3)
		{
			this->elems[2] = c;
		}
	}

	bool operator==(
			const PointT<TYPE, DIM> &a) const
	{
		if (DIM == 2)
		{
			return (this->elems[0] == a[0] && this->elems[1] == a[1]) ?
					true : false;
		}
		else
		{
			return (this->elems[0] == a[0] && this->elems[1] == a[1]
					&& this->elems[2] == a[2]) ? true : false;
		}
	}
	bool operator!=(
			const PointT<TYPE, DIM> &a) const
	{
		if (DIM == 2)
		{
			return !(
					(this->elems[0] == a[0] && this->elems[1] == a[1]) ?
							true : false);
		}
		else
		{
			return !(
					(this->elems[0] == a[0] && this->elems[1] == a[1]
							&& this->elems[2] == a[2]) ? true : false);
		}
	}
	void show() const
	{
		std::cout << std::scientific << "( " << this->elems[0] << " , "
				<< this->elems[1];
		if (DIM == 3)
		{
			std::cout << " , " << this->elems[2] << " )\n";
		}
		else
		{
			std::cout << " )\n";
		}
	}
	template<typename T>
	void transfer(
			const T&dx,
			const T&dy,
			const T&dz)
	{
		this->elems[0] = this->elems[0] + TYPE(dx);
		this->elems[1] = this->elems[1] + TYPE(dy);
		if (DIM == 3)
		{
			this->elems[2] = this->elems[2] + TYPE(dz);
		}
	}
	template<typename T>
	void scale(
			const T&dx,
			const T&dy,
			const T&dz)
	{
		this->elems[0] = this->elems[0] * TYPE(dx);
		this->elems[1] = this->elems[1] * TYPE(dy);
		if (DIM == 3)
		{
			this->elems[2] = this->elems[2] * TYPE(dz);
		}
	}
	inline size_type size() const
	{
		return size_type(DIM);
	}
};

} //end namespace

#endif /* POINT_H_ */
