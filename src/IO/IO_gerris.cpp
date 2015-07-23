/************************
 //  \file   IO_gerris.cpp
 //  \brief
 // 
 //  \author zhou
 //  \date   22 août 2014 
 ***********************/

#include "IO_gerris.h"
#include "../Algebra/Arithmetic.h"
//#include "../Grid/QuadTreeNode.h"
#include "../Grid/Cell.h"
#include "../Utility/Array.h"
//#include "../Grid/TreeFactory.h"
//#include "../Calculation/VOF.h"

using namespace std;

namespace Larus {

int compareLineX(const vector<Float> &a, const vector<Float> &b) {
	if (*(a.begin()) < *(b.begin()))
		return 1;
	else
		return 0;
}

int compareLineY(const vector<Float> &a, const vector<Float> &b) {
	if (*(a.begin() + 1) < *(b.begin() + 1))
		return 1;
	else
		return 0;
}

int cal_level(Float x, Float b, Float L, int ls, int le) {
	arrayList arr(le - ls + 1);
	for (int i = ls; i <= le; i++) {
		Float level = (x - b) / L * pow(2, i) - 0.5;
		arr[i - ls] = ABS(level - round(level));
		if (ABS(level - round(level)) < 0.01) {
			return i;
		}
	}
	return arr.findMinIdx() + ls;
}

void _t_condition_point_at_which_child(arrayList& arrt, pOCNode pnode,
		utPointer p) {
	if (pnode != NULL_PTR && p != NULL_PTR) {
		utPointer ap = CAST(arrayListT<utPointer>*,p)->at(0);
		arrt[pnode->whichChild((*CAST(Point3D*, ap)))] = 1;
	}
}

void _t_condition_point_at_which_child(arrayList& arrt, pQTNode pnode,
		utPointer p) {
	if (pnode != NULL_PTR && p != NULL_PTR) {
		utPointer ap = CAST(arrayListT<utPointer>*,p)->at(0);
		arrt[pnode->whichChild((*CAST(Point2D*, ap)))] = 1;
	}
}

void _visit_insert_level_copydata(pQTNode pn, utPointer p) {
	utPointer ap2 = CAST(arrayListT<utPointer>*,p)->at(2);
	int level = *(CAST(int*, ap2));
	if (pn->getLevel() < level) {
		utPointer ap1 = CAST(arrayListT<utPointer>*,p)->at(1);
		utPointer ap0 = CAST(arrayListT<utPointer>*,p)->at(0);
		array<bool, QTNode::NUM_CELLS> ctt;
		int idx_child = point_at_which_child(pn, *(CAST(Point2D*, ap0)));
		for (int i = 0; i < pn->NUM_CELLS; i++) {
			if (i == idx_child) {
				ctt[i] = true;
			} else {
				ctt[i] = false;
			}
		}
		CAST(pQuadTree, ap1)->Creatchildren_partial(pn, ctt[0], ctt[1], ctt[2],
				ctt[3]);
	} else if (pn->getLevel() == level) {
		utPointer ap3 = CAST(arrayListT<utPointer>*,p)->at(3);
		pn->data = CAST(pCellData2D, ap3);
	}
}

void _visit_insert_level_copydata(pOCNode pn, utPointer p) {
	// if (condition_is_leaf(pn)) {
	utPointer ap2 = CAST(arrayListT<utPointer>*,p)->at(2);
	int level = *(CAST(int*, ap2));
	if (pn->getLevel() < level) {
		utPointer ap1 = CAST(arrayListT<utPointer>*,p)->at(1);
		utPointer ap0 = CAST(arrayListT<utPointer>*,p)->at(0);
		array<bool, OCNode::NUM_CELLS> ctt;
		int idx_child = point_at_which_child(pn, *(CAST(Point3D*, ap0)));
		for (int i = 0; i < pn->NUM_CELLS; i++) {
			if (i == idx_child) {
				ctt[i] = true;
			} else {
				ctt[i] = false;
			}
		}
		CAST(pOCTree, ap1)->Creatchildren_partial(pn, ctt[0], ctt[1], ctt[2],
				ctt[3], ctt[4], ctt[5], ctt[6], ctt[7]);
	} else if (pn->getLevel() == level) {
		utPointer ap3 = CAST(arrayListT<utPointer>*,p)->at(3);
		pn->data = CAST(pCellData3D, ap3);
	}
	// }
}

}

