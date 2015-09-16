/************************
 //  \file   CalDef.h
 //  \brief
 // 
 //  \author czhou
 //  \date   14 janv. 2015 
 ***********************/
#ifndef CALDEF_H_
#define CALDEF_H_

#include "../TypeDef.h"
#include "../Grid/SPTree.h"
#include "../Grid/SPTreeNode.h"
#include "../Grid/Forest.h"

namespace Larus {

const LarusDef::size_type Idx_IDX = 0;
const LarusDef::size_type LS_IDX = 0;
const LarusDef::size_type VOF_IDX = 0;
const LarusDef::size_type VOFt_IDX = 0;
const LarusDef::size_type U_IDX = 2;
const LarusDef::size_type V_IDX = 3;

template<class NODE>
inline int getIDX(NODE* pn) {
	ASSERT(pn!=NULL_PTR);
	ASSERT(pn->data!=NULL_PTR);
	return int(pn->data->aCenterData[Idx_IDX]);
}
inline int getIDX(const Forest2D::iterator& iter) {
	ASSERT(iter->data!=NULL_PTR);
	return iter->data->aCenterData[Idx_IDX];
}
inline int getIDX(const Forest2D::const_iterator& iter) {
	ASSERT(iter->data!=NULL_PTR);
	return iter->data->aCenterData[Idx_IDX];
}
// get center data on array
template<class NODE>  //get center value
inline Float getcVal(const NODE* pn, const LarusDef::size_type& idx) {
	ASSERT(pn!=NULL_PTR);
	ASSERT(pn->data!=NULL_PTR);
	return pn->data->aCenterData[idx];
}

template<class NODE>
inline Float& refcVal(const NODE* pn, const LarusDef::size_type& idx) {
	ASSERT(pn!=NULL_PTR);
	ASSERT(pn->data!=NULL_PTR);
	return pn->data->aCenterData[idx];
}
inline Float& refcVal(Forest2D::iterator& iter, const LarusDef::size_type idx) {
	ASSERT(iter->data!=NULL_PTR);
	return iter->data->aCenterData[idx];
}

template<typename NODE>
void pfun_visit_is_leaf(NODE* pn, utPointer utp) {
	typedef NODE* pNode;
	if (condition_is_leaf(pn)) {
		arrayListT<utPointer>& arrutp = (*CAST(arrayListT<utPointer>*, utp));
		ListT<pNode>& ll = (*(CAST(ListT<pNode>*, arrutp[1])));
		ll.push_back(pn);
	}
}
template<typename NODE>
void pfun_visit_is_leaf_push_to_list(NODE* pn, utPointer utp) {
	typedef NODE* pNode;
	if (condition_is_leaf(pn)) {
		arrayListT<utPointer>& arrutp = (*CAST(arrayListT<utPointer>*, utp));
		ListT<pNode>& ll = (*(CAST(ListT<pNode>*, arrutp[0])));
		ll.push_back(pn);
	}
}

const short DIRvsNODEIDX[26][8] = {		//
		{ 1, 0, 0, 0, 1, 0, 0, 0 },     //
				{ 0, 1, 0, 0, 0, 1, 0, 0 },  //
				{ 0, 0, 1, 0, 0, 0, 1, 0 },  //
				{ 0, 0, 0, 1, 0, 0, 0, 1 },  //
				{ 1, 1, 0, 0, 1, 1, 0, 0 },  //
				{ 1, 0, 1, 0, 1, 0, 1, 0 },  //
				{ 0, 0, 1, 1, 0, 0, 1, 1 },  //
				{ 0, 1, 0, 1, 0, 1, 0, 1 },  //
				{ 0, 0, 0, 0, 1, 1, 1, 1 },  //
				{ 1, 1, 1, 1, 0, 0, 0, 0 },  //
				{ 0, 0, 0, 0, 0, 1, 0, 1 },  //
				{ 0, 1, 0, 1, 0, 0, 0, 0 },  //
				{ 0, 0, 0, 0, 1, 0, 1, 0 },  //
				{ 1, 0, 1, 0, 0, 0, 0, 0 },  //
				{ 0, 0, 1, 1, 0, 0, 0, 0 },  //
				{ 1, 1, 0, 0, 0, 0, 0, 0 },  //
				{ 0, 0, 0, 0, 0, 0, 1, 1 },  //
				{ 0, 0, 0, 0, 1, 1, 0, 0 },  //
				{ 0, 1, 0, 0, 0, 0, 0, 0 },  //
				{ 0, 0, 0, 1, 0, 0, 0, 0 },  //
				{ 1, 0, 0, 0, 0, 0, 0, 0 },  //
				{ 0, 0, 1, 0, 0, 0, 0, 0 },  //
				{ 0, 0, 0, 0, 0, 1, 0, 0 },  //
				{ 0, 0, 0, 0, 0, 0, 0, 1 },  //
				{ 0, 0, 0, 0, 1, 0, 0, 0 },  //
				{ 0, 0, 0, 0, 0, 0, 1, 0 },  //
		};
template<typename NODE>
void pfun_condition_is_on_direction(arrayList& arr, NODE* pn, utPointer utp) {
	arrayListT<utPointer>& arrutp = (*CAST(arrayListT<utPointer>*, utp));
	SPDirection& dir = (*(CAST(SPDirection*, arrutp[0])));
	for (short i = 0; i < NODE::NUM_CELLS; ++i) {
		arr[i] = (DIRvsNODEIDX[short(dir)][i] == 1) ? 1 : 0;
	}
}
/*
 *   get List pNode on direction
 */
template<typename NODE>
void getListpNode_direction(NODE* pn, SPDirection dir, ListT<NODE*>& ll) {
	arrayListT<utPointer> arrutp(2);
	arrutp[0] = &dir;
	arrutp[1] = &ll;
	pn->Traversal_conditional(pfun_condition_is_on_direction,
			pfun_visit_is_leaf, &arrutp);
}
/*
 *   get List pNode all the children
 */
template<typename NODE>
void getListpNode_children(NODE* pn, ListT<NODE*>& ll) {
	arrayListT<utPointer> arrutp(1);
	arrutp[0] = &ll;
	pn->Traversal(pfun_visit_is_leaf_push_to_list, &arrutp);
}

template<typename NODE>
void visit_average_value_from_leafs_1(NODE* pn, utPointer up) {
	if (condition_is_leaf(pn)) {
		arrayListT<utPointer> &arrp = (*CAST(arrayListT<utPointer>*, up));
		int& count = (*CAST(int*, arrp[0]));
		int& idx = (*CAST(int*, arrp[1]));
		Float& res = (*CAST(Float*, arrp[2]));
		if (pn->data != NULL_PTR) {
			res += pn->data->aCenterData[idx];
			count = count + 1;
		}
	}
}
template<class NODE>  //get average value from leafs
inline Float getAverageVal(NODE* pn, LarusDef::size_type idx) {
	ASSERT(pn!=NULL_PTR);
	if (condition_is_leaf(pn)) {
		return getcVal(pn, idx);
	} else {
		Float res = 0.0;
		int count = 0;
		arrayListT<utPointer> arrp(3);
		arrp[0] = &count;
		arrp[1] = &idx;
		arrp[2] = &res;
		pn->Traversal(visit_average_value_from_leafs_1, &arrp);
		res /= count;
	}

}

//work with Forest===============================
//===============================================

//work with tree=================================
//===============================================

}

#endif /* CALDEF_H_ */
