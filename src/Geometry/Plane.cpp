/************************
 //  \file   Plane.cpp
 //  \brief
 // 
 //  \author czhou
 //  \date   3 oct. 2014 
 ***********************/

#include "Plane.h"
#include "GeoDef.h"
#include <assert.h>
#include <iostream>

namespace Larus {
Plane::Plane() {
	this->assign(0);
}
Plane::Plane(Float a, Float b, Float c, Float d) : array_4<Float>(a,b,c,d){
	assert(!(a == 0.0 && b == 0.0 && c == 0.0));
}
Plane::Plane(const Point3D & p, const Vector3D& n) {
	assert(n.isExist());
	this->elems[0] = n.x;
	this->elems[1] = n.y;
	this->elems[2] = n.z;
	this->elems[3] = (this->elems[0] * p.x + this->elems[1] * p.y + this->elems[2] * p.z);
}

void Plane::reconstruct(Float a, Float b, Float c, Float alpha) {
	assert(!(a == 0.0 && b == 0.0 && c == 0.0));
	this->elems[0] = a;
	this->elems[1] = b;
	this->elems[2] = c;
	this->elems[3] = alpha;
}

Float Plane::calX(Float y, Float z) const {
	if (this->elems[0] == 0.0) {
		return (this->elems[3] - this->elems[1] * y - this->elems[2] * z) / SMALL;
	} else {
		return (this->elems[3] - this->elems[1] * y - this->elems[2] * z) / this->elems[0];
	}
}
Float Plane::calY(Float x, Float z) const {
	if (this->elems[1] == 0.0) {
		return (this->elems[3] - this->elems[0] * x - this->elems[2] * z) / SMALL;
	} else {
		return (this->elems[3] - this->elems[0] * x - this->elems[2] * z) / this->elems[1];
	}
}
Float Plane::calZ(Float x, Float y) const {
	if (this->elems[2] == 0.0) {
		return (this->elems[3] - this->elems[0] * x - this->elems[1] * y) / SMALL;
	} else {
		return (this->elems[3] - this->elems[0] * x - this->elems[1] * y) / this->elems[2];
	}
}

Float Plane::getIntersept(CSAxis cd) const {
	switch (cd) {
	case CSAxis_X: {
		return this->calX(0, 0);
		break;
	}
	case CSAxis_Y: {
		return this->calY(0, 0);
		break;
	}
	case CSAxis_Z: {
		return this->calZ(0, 0);
		break;
	}
	default:
		return 0;
	}
}

Float Plane::getNorm(CSAxis cd) const {
	switch (cd) {
	case CSAxis_X: {
		return this->getA();
		break;
	}
	case CSAxis_Y: {
		return this->getB();;
		break;
	}
	case CSAxis_Z: {
		return this->getC();;
		break;
	}
	default:
		return 0;
	}
}

bool Plane::isExist() const {
	if (this->elems[0] == 0 && this->elems[1] == 0 && this->elems[2] == 0) {
		return false;
	} else {
		return true;
	}
}
bool Plane::isXY() const {
	if (isZero(this->elems[2])) {
		return true;
	} else {
		return false;
	}
}
bool Plane::isYZ() const {
	if (isZero(this->elems[0])) {
		return true;
	} else {
		return false;
	}
}
bool Plane::isZX() const {
	if (isZero(this->elems[1])) {
		return true;
	} else {
		return false;
	}
}

void Plane::show() const{
	std::cout << this->elems[0] << " X + " << this->elems[1] << " Y + " << this->elems[2] <<" Z = "<<this->elems[3]<< std::endl;
}

} //end of name space
