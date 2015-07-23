/*
 * Adaptive.h
 *
 *  Created on: Jul 6, 2015
 *      Author: zhou
 */

#ifndef ADAPTIVE_H_
#define ADAPTIVE_H_

#include "../TypeDef.h"
#include "CalDef.h"
#include "../Grid/SPTree.h"
#include "../Grid/SPTreeNode.h"
#include "../Grid/Forest.h"
#include "../Algebra/Arithmetic.h"
#include "../Algebra/Expression.h"

namespace Larus {

// initial adaptation ==========================
// ===========================================
void initial_adaptation( //
		Forest2D& forest,      // Forest 2D
		Point2D& point,  // A Point2D
		int minlevel,          // minlevel
		int maxlevel);         // maxlevel

void initial_adaptation( //
		Forest2D& forest,      // Forest 2D
		ListT<Point2D>& lp,    // A Point2D
		int minlevel,          // minlevel
		int maxlevel);         // maxlevel

typedef Float (*pFun)(Float, Float, Float);

void initial_adaptation_eq( //
		Forest2D& forest,      // Forest 2D
		pFun pf,               // A Point2D
		Float value,           // Threshold value
		int minlevel,          // minlevel
		int maxlevel);         // maxlevel
void initial_adaptation_le( //
		Forest2D& forest,      // Forest 2D
		pFun pf,               // A Point2D
		Float value,           // Threshold value
		int minlevel,          // minlevel
		int maxlevel) ;         // maxlevel
// Dynamic adaptive =========================
// ==========================================




}



#endif /* CALCULATION_ADAPTIVE_H_ */
