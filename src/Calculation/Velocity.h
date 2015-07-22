/************************
 //  \file   Velocity.h
 //  \brief
 // 
 //  \author czhou
 //  \date   13 mars 2015 
 ***********************/
#ifndef VELOCITY_H_
#define VELOCITY_H_

#include "Scalar.h"
#include "../TypeDef.h"
#include "CalDef.h"
#include "../Grid/SPTree.h"
#include "../Grid/SPTreeNode.h"
#include "../Grid/Forest.h"
#include "../Algebra/Arithmetic.h"
#include "../Algebra/Expression.h"

namespace Larus {

void cal_stream_line(  //Forest 2D
		Point2D& bp,  //the start point
		Forest2D& forest,  //forest
		Float dt,     //dt
		Float nump,   //number of points
		ListT<Point2D>& list,   //result list
		int idxu,   //idx u
		int idxv);  //idx v

}

#endif /* VELOCITY_H_ */
