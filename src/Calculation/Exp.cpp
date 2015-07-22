/*
 * Exp.cpp
 *
 *  Created on: Jul 22, 2015
 *      Author: zhou
 */

#include "Exp.h"


namespace Larus{
void interpolate_expression_on_axis_2f( // 2D QuadTree Node
		pQTNode pc,       //node
		CSAxis axis,      //axix
		Float dis,        //distance to center of pn with sign
		Exp2D& exp   //Expression
		) {
	ASSERT(pc != NULL_PTR);
	ASSERT(T_dis_check(pc, axis, dis));
	ASSERT(axis != CSAxis_Z && axis != ErrCSAxis);
	//special case: there is no neighbor node
	if (dis == 0) {
		T_exp_1node(pc, exp);
		return;
	}
	// axis and dis to direction
	SPDirection dir =
			(axis == CSAxis_X) ?
					(dis < 0 ? SPD_IM : SPD_IP) : (dis < 0 ? SPD_JM : SPD_JP);
	pQTNode pn = NULL_PTR;
	if (0 == T_find_neighbor_2(pc, dir, pn)) {
		T_exp_1node(pc, exp);
		return;
	} else { //pn!=NULL
		_exp_interplolate_2node(pc,  //node
				pn, axis,                 //axis
				dis,                 //distance to center of pn
				exp);
	}
}

void interpolate_expression_on_axis_2f( // 2D QuadTree Node
		pQTNode pc,       //node
		CSAxis axis,      //axix
		Float dis,        //distance to center of pn with sign
		Exp2D& exp,        //Expression
		BCManager<Dimension_2D>& bcm
		) {
	ASSERT(pc != NULL_PTR);
	ASSERT(T_dis_check(pc, axis, dis));
	ASSERT(axis != CSAxis_Z && axis != ErrCSAxis);
	//special case: there is no neighbor node
	if (dis == 0) {
		T_exp_1node(pc, exp);
		return;
	}
	// axis and dis to direction
	SPDirection dir =
			(axis == CSAxis_X) ?
					(dis < 0 ? SPD_IM : SPD_IP) : (dis < 0 ? SPD_JM : SPD_JP);
	pQTNode pn = NULL_PTR;
	if (0 == T_find_neighbor_2(pc, dir, pn)) {
		pn = bcm.find_ghost(pc->data->aCenterData[Idx_IDX], dir);
		if (pn == NULL_PTR){
			T_exp_1node(pc, exp);
			return;
		}
	}
	//pn!=NULL
	_exp_interplolate_2node(pc,  //node
			pn, axis,                 //axis
			dis,                 //distance to center of pn
			exp);
}

void interpolate_expression_on_axis_3fb( // 2D QuadTree Node
		pQTNode pc,  //node
		CSAxis axis, //axix
		Float dis,   //distance to center of pn with sign
		Exp2D& exp   //Expression
		) {
	ASSERT(pc != NULL_PTR);
	ASSERT(T_dis_check(pc, axis, dis));
	ASSERT(axis != CSAxis_Z && axis != ErrCSAxis);
	//special case: there is no neighbor node
	if (dis == 0) {
		T_exp_1node(pc, exp);
		return;
	}
	// axis and dis to direction
	SPDirection dirf =
			(axis == CSAxis_X) ?
					(dis < 0 ? SPD_IM : SPD_IP) : (dis < 0 ? SPD_JM : SPD_JP);
	pQTNode pf = NULL_PTR;
	pQTNode pb = NULL_PTR;
	int re_find = T_find_neighbor_3fb(pc, dirf, pf, pb);
	if (0 == re_find) {
		T_exp_1node(pc, exp);
		return;
	} else if (ABS(re_find) == 1) { //
		_exp_interplolate_2node(pc,         //node
				(re_find == -1) ? pf : pb,  //node neighbor
				axis,                       //axis
				dis,                        //distance to center of pn
				exp);
		return;
	} else {
		_exp_interplolate_3node(pb, pc, pf, axis, dis, exp);
		return;
	}
}

void interpolate_expression_on_axis_3fb( // 2D QuadTree Node
		pQTNode pc,  //node
		CSAxis axis, //axix
		Float dis,   //distance to center of pn with sign
		Exp2D& exp,   //Expression
		BCManager<Dimension_2D>& bcm) {
	ASSERT(pc != NULL_PTR);
	ASSERT(T_dis_check(pc, axis, dis));
	ASSERT(axis != CSAxis_Z && axis != ErrCSAxis);
	//special case: there is no neighbor node
	if (dis == 0) {
		T_exp_1node(pc, exp);
		return;
	}
	// axis and dis to direction
	SPDirection dirf =
			(axis == CSAxis_X) ?
					(dis < 0 ? SPD_IM : SPD_IP) : (dis < 0 ? SPD_JM : SPD_JP);
	pQTNode pf = NULL_PTR;
	pQTNode pb = NULL_PTR;
	int re_find = T_find_neighbor_3fb(pc, dirf, pf, pb);
	switch (re_find) {
	case 0: {
		pf = bcm.find_ghost(getIDX(pc), dirf);
		pb = bcm.find_ghost(getIDX(pc), oppositeDirection(dirf));
		break;
	}
	case -1: {
		pb = bcm.find_ghost(getIDX(pc), oppositeDirection(dirf));
		break;
	}
	case 1: {
		pf = bcm.find_ghost(getIDX(pc), dirf);
		break;
	}
	default: {
		break;
	}
	}
	if (re_find != 3) {
		if (pf == NULL_PTR && pb == NULL_PTR) {
			re_find = 0;
		}
		if (pf == NULL_PTR) {
			re_find = 1;
		}
		if (pb == NULL_PTR) {
			re_find = -1;
		}
	}
	if (0 == re_find) {
		T_exp_1node(pc, exp);
		return;
	} else if (ABS(re_find) == 1) { //
		_exp_interplolate_2node(pc,         //node
				(re_find == -1) ? pf : pb,  //node neighbor
				axis,                       //axis
				dis,                        //distance to center of pn
				exp);
		return;
	} else {
		_exp_interplolate_3node(pb, pc, pf, axis, dis, exp);
		return;
	}
}

}


