//
//  Relation.h
//  LarusX_0_1
//
//  Created by zhou on 27/12/14.
//  Copyright (c) 2014 zhou. All rights reserved.
//

#ifndef Relation_h
#define Relation_h

#include "../TypeDef.h"
#include "GeoDef.h"
#include "Line.h"
#include "Plane.h"

namespace Larus {
    //interset===================================
    Point3D intersect(const Line3D& l, const Plane& p);
    

    //Predicate==================================
    //                 point     (sphere: center point, radius       )
    bool is_in(const Point3D& p, const Point3D& pc, const Float& r);



   /* template <class A, class B>
    void
    _intersect(const A& geo_a,
               line_3d_tag,
               const B& geo_b,
               plane_tag)
    {
        __intersect(geo_a,geo_b);
    }


    template <class res, class A, class B>
    res interset(const A& geo_a, const B& geo_b){
    	typedef typename geometry_traits<A>::self_tag a_tag;
    	typedef typename geometry_traits<B>::self_tag b_tag;
    	return static_cast<res>(_intersect(geo_a, a_tag(), geo_b, b_tag()));
    	printf("ashdfasdfa \n");
    } */

}

#endif
