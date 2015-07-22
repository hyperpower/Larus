/*
 * GeoDef.h
 *
 *  Created on: Dec 23, 2014
 *      Author: zhou
 */

#ifndef GEOMETRY_GEODEF_H_
#define GEOMETRY_GEODEF_H_

#include "../TypeDef.h"

namespace Larus {

enum CSAxis {
	ErrCSAxis = -1, CSAxis_X = 0, CSAxis_Y = 1, CSAxis_Z = 2,
};

enum CSDirection {
	ErrCSDirection = -1,
	CSDirection_Xm = -2,
	CSDirection_Xp = 0,
	CSDirection_Ym = -3,
	CSDirection_Yp = 1,
	CSDirection_Zm = -4,
	CSDirection_Zp = 2,
};

    
    //Define the Geometry tags=====================
    struct point_2d_tag {};
    struct point_3d_tag {};
    //struct line_2d_tag {};
    //struct line_3d_tag {};
    //struct segment_2d_tag {};
    //struct segment_3d_tag {};
    //struct triangle_2d_tag {};
    //struct triangle_3d_tag {};
    //struct plane_tag {};
    
    template <class Geo>
    struct geometry_traits {
        typedef typename Geo::self_tag self_tag;
    };

}



#endif /* GEOMETRY_GEODEF_H_ */
