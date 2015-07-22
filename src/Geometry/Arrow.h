/*
 * arrow.h
 *
 *  Created on: Dec 24, 2014
 *      Author: zhou
 */

#ifndef GEOMETRY_ARROW_H_
#define GEOMETRY_ARROW_H_

#include "../TypeDef.h"
#include "../Utility/Array.h"
#include "Point.h"
#include "Vector.h"

namespace Larus {
template<typename TYPE>
class arrow_2D: public array_4<TYPE> {
public:
	//constructor
	arrow_2D() :
			array_4<TYPE>() {
	}
	arrow_2D(const TYPE& x, const TYPE& y, const TYPE& dx, const TYPE& dy) :
			array_4<TYPE>(x, y, dx, dy) {
	}
	arrow_2D(const Point_2D<TYPE>& p, const Vector_2D<TYPE>& v) :
			array_4<TYPE>(p.x, p.y, v.x, v.y) {
	}
};

typedef arrow_2D<Float> Arrow2D;

}

#endif /* GEOMETRY_ARROW_H_ */
