/************************
 //  \file   Sphere.h
 //  \brief
 // 
 //  \author czhou
 //  \date   16 févr. 2015 
 ***********************/
#ifndef SPHERE_H_
#define SPHERE_H_

#include "../TypeDef.h"
#include "../Utility/Array.h"
#include "GeoDef.h"
#include <iostream>

#include <math.h>

namespace Larus {

template<typename TYPE>
class Sphere: public ObjectBase {
protected:
	Point_3D<TYPE> _cp;
	TYPE r;
public:
	Sphere();
	Sphere(const Point_3D<TYPE>& p, const TYPE& r);
	Sphere(const Sphere&);

};

}

#endif /* SPHERE_H_ */
