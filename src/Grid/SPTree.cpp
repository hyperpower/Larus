/************************
 //  \file   Octree.cpp
 //  \brief
 // 
 //  \author zhou
 //  \date   3 sept. 2014 
 ***********************/

#include "Cell.h"
#include "SPTree.h"
#include "../TypeDef.h"

#include <stdio.h>
#include <iostream>
#include <math.h>

namespace Larus {
//QuadTree=======================================
//===============================================
void visit_draw_gnuplot(pQTNode pn, utPointer p) {
	if (condition_is_leaf(pn)) {
		FILE* file = CAST(FILE*, p);
		draw_gnuplot_boundary(pn,file,SP_XY);
	}
}
//-----------------------------------------------
//-----------------------------------------------

//OcTree=======================================
//===============================================


} //this is the end of namespace
