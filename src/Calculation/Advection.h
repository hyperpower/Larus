/*
 * Advection.h
 *
 *  Created on: Mar 11, 2015
 *      Author: zhou
 */

#ifndef _ADVECTION_H_
#define _ADVECTION_H_

#include "../TypeDef.h"
#include "CalDef.h"
#include "../Grid/SPTree.h"
#include "../Grid/SPTreeNode.h"
#include "../Grid/Forest.h"
#include "../Algebra/Arithmetic.h"
#include "../Algebra/Expression.h"

namespace Larus {

//k-scheme ======================================
//===============================================
Float k_scheme_CDS(Float c, Float dc, Float gp, Float gm);    //central difference
Float k_scheme_QUICK(Float c, Float dc, Float gp, Float gm);  //Quadratic-upwind
Float k_scheme_CUI(Float c, Float dc, Float gp, Float gm);    //cubic-upwind
Float k_scheme_Fromm(Float c, Float dc, Float gp, Float gm);  //Fromm's scheme
Float k_scheme_SOU(Float c, Float dc, Float gp, Float gm);    //second-order upwind

//Traversal ListFace=============================
void traversal_to_ListFace(pQuadTree ptree, ListT<QTNodeFace>& listf);
//draw ==========================================
void draw_gnuplot_ListFace(std::string filename,  //filename
		ListT<QTNodeFace>& listf, //list
		int mode              //mode
		);

} //end of namespace



#endif /* CALULATION_ADVECTION_H_ */
