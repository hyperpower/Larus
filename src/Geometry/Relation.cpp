//
//  Relation.cpp
//  LarusX_0_1
//
//  Created by zhou on 27/12/14.
//  Copyright (c) 2014 zhou. All rights reserved.
//

#include "../TypeDef.h"
#include "GeoDef.h"
#include "Line.h"
#include "Plane.h"

namespace Larus {
Point3D intersect(const Line3D& l, const Plane& p) {
	const Float& A = p[0];
	const Float& B = p[1];
	const Float& C = p[2];
	const Float& D = p[3];
	Float t;
	Float s = (D - A * l.first[0] - B * l.first[1] - C * l.first[2]);
	Float d = (A * l.second[0] + B * l.second[1] + C * l.second[2]);
	if (d == 0) {
		t = s / SMALL;
	} else {
		t = s / d;
	}
	return Point3D(l.first[0] + l.second[0] * t, l.first[1] + l.second[1] * t,
			l.first[2] + l.second[2] * t);
}

bool is_in(const Point3D& p, const Point3D& pc, const Float& r){
	if(calDistance(p,pc)-r<0){
		return true;
	}else{
		return false;
	}
}


}
