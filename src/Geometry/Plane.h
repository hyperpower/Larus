/************************
 //  \file   Plane.h
 //  \brief
 // 
 //  \author czhou
 //  \date   3 oct. 2014 
 ***********************/
#ifndef PLANE_H_
#define PLANE_H_

#include "Point.h"
#include "Vector.h"
#include "GeoDef.h"
#include "../TypeDef.h"


namespace Larus{

class Plane : public array_4<Float>{
	//The Plane function defined as Ax+By+Cz=D
	//
public:
	//typedef plane_tag self_tag;
	typedef Plane self_type;
	Plane();
	Plane(Float, Float, Float, Float);
	Plane(const Point3D &, const Vector3D &);
	~Plane(){};
	void reconstruct(Float, Float, Float, Float);
	inline Float getA() const;
	inline Float getB() const;
	inline Float getC() const;
	inline Float getD() const;
	Float calX(Float y,Float z) const;
	Float calY(Float x,Float z) const;
	Float calZ(Float x,Float y) const;

	Float getIntersept(CSAxis) const;
	Float getNorm(CSAxis) const;
	//Float getNormY() const;

	bool isExist() const;
	bool isXY() const;
	bool isYZ() const;
	bool isZX() const;

	void show() const;

};

inline Float Plane::getA() const{
	return this->elems[0];
}
inline Float Plane::getB() const{
	return this->elems[1];
}
inline Float Plane::getC() const{
	return this->elems[2];
}
inline Float Plane::getD() const{
	return this->elems[3];
}


}



#endif /* PLANE_H_ */
