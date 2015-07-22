/************************
 //  \file   Velocity.cpp
 //  \brief
 // 
 //  \author czhou
 //  \date   13 mars 2015 
 ***********************/

#include "Velocity.h"

namespace Larus {


void cal_stream_line(  //Forest 2D
		Point2D& bp,  //the start point
		Forest2D& forest,  //forest
		Float dt,     //dt
		Float nump,   //number of points
		ListT<Point2D>& list,   //result list
		LarusDef::size_type idxu,   //idx u
		LarusDef::size_type idxv) {  //idx v
	pQuadTree pt = forest.getpTree(bp);
	_IF_TRUE_RETRUN(NULL_PTR==pt);
	_IF_TRUE_RETRUN(NULL_PTR==pt->Find(bp));
	list.push_back(bp);
	Point2D op, np;
	op = bp;
	arrayList_st arridx(2);
	arridx[0] = idxu;
	arridx[1] = idxv;
	arrayList arrres(2);
	for (int i = 0; i < nump; i++) {
		interpolate( // 2D QuadTree
				pt,  //pQuadTree
				op,  //point
				arridx,  //data index
				arrres   //data res
				);
		np.x = op.x + arrres[0] * dt;
		np.y = op.y + arrres[1] * dt;
		list.push_back(np);
		op = np;
		pt = forest.getpTree(op);
		if (NULL_PTR==pt) {
			i = nump;
		}
	}
}

}
