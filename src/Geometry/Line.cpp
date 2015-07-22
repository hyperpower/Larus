/************************
 //  \file   Line.cpp
 //  \brief
 // 
 //  \author zhou
 //  \date   23 mai 2014 
 ***********************/

#include "Line.h"
#include <assert.h>
#include <iostream>

namespace Larus {

Line::Line(const Point2D &a, const Point2D &b) {
	assert(a != b);
	if (isEqual(a.x, b.x)) {
		this->elems[0] = 1;
		this->elems[1] = 0;
		this->elems[2] = a.x;
	} else if (isEqual(a.y, b.y)) {
		this->elems[0] = 0;
		this->elems[1] = 1;
		this->elems[2] = a.y;
	} else {
		this->elems[0] = 1.0 / (a.x - b.x);
		this->elems[1] = -1.0 / (a.y - b.y);
		this->elems[2] = b.x / (a.x - b.x) - b.y / (a.y - b.y);
	}
}

Line::Line(Float ax, Float ay, Float bx, Float by) {
	//assert(!isEqual(ax, bx) || !isEqual(ay,by));
	Point2D p1(ax, ay), p2(bx, by);
	if (isEqual(p1.x, p2.x)) {
		this->elems[0] = 1;
		this->elems[1] = 0;
		this->elems[2] = p1.x;
	} else if (isEqual(p1.y, p2.y)) {
		this->elems[0] = 0;
		this->elems[1] = 1;
		this->elems[2] = p1.y;
	} else {
		this->elems[0] = 1.0 / (p1.x - p2.x);
		this->elems[1] = -1.0 / (p1.y - p2.y);
		this->elems[2] = p2.x / (p1.x - p2.x) - p2.y / (p1.y - p2.y);
	}
}

void Line::reconstruct(Float a, Float b, Float c) {
	if (a == 0.0 && b == 0.0) {
		a = 1.1 * SMALL;
	} else if (isZero(a) && isZero(b)) {
		a = 1.1 * SMALL * a / abs(a);
		b = 1.1 * SMALL * b / abs(b);
	}
	assert(!isZero(a) || !isZero(b));
	this->elems[0] = a;
	this->elems[1] = b;
	this->elems[2] = c;
}

Float Line::calX(Float y) const {
	if (this->elems[0] == 0.0) {
		std::cerr << "*> warning: Function Line::calX(Float ) \n";
		std::cerr << "*> Line Ax+By=C, Coefficient A=0 \n";
		std::cerr << "*> Automactically change to A= " << SMALL << std::endl;
		return (this->elems[2] - this->elems[1] * y) / SMALL;
	} else {
		return (this->elems[2] - this->elems[1] * y) / this->elems[0];
	}
}

Float Line::calY(Float x) const {
	if (this->elems[1] == 0.0) {
		std::cerr << "*> warning: Function Line::calY(Float ) \n";
		std::cerr << "*> Line Ax+By=C, Coefficient B=0 \n";
		std::cerr << "*> Automactically change to A= " << SMALL << std::endl;
		return (this->elems[2] - this->elems[0] * x) / SMALL;
	} else {
		return (this->elems[2] - this->elems[0] * x) / this->elems[1];
	}
}

Float Line::getSlope() const {
	return -this->elems[0] /(this->elems[1]+SMALL);
}

Float Line::getInterseptX() const {
	return this->elems[2] / (this->elems[0]+SMALL);
}

Float Line::getInterseptY() const {
	return this->elems[2] / (this->elems[1]+SMALL);
}

Float Line::getNormX() const {
	return this->elems[0];
}

Float Line::getNormY() const {
	return this->elems[1];
}

bool Line::isExist() const {
	if (this->elems[0] != 0.0 || this->elems[1] != 0.0) {
		return true;
	} else {
		return false;
	}
}

bool Line::isVertical() const {
	if (isZero(this->elems[1])) {
		return true;
	} else {
		return false;
	}
}

bool Line::isHorizontal() const {
	if (isZero(this->elems[0])) {
		return true;
	} else {
		return false;
	}
}

void Line::show() const {
	std::cout << this->elems[0] << " X + " << this->elems[1] << " Y= " << this->elems[2] << std::endl;
}
//Line3D========================================
void Line3D::show() const{
	std::cout<< " X = " << first[0] << " + " << second[0]<< " * t " << std::endl;
	std::cout<< " Y = " << first[1] << " + " << second[1]<< " * t " << std::endl;
	std::cout<< " Z = " << first[2] << " + " << second[2]<< " * t " << std::endl;
}

}
