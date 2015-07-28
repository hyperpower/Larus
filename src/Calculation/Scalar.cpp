/*
 * Scalar.cpp
 *
 *  Created on: Feb 7, 2015
 *      Author: zhou
 */

#include "Scalar.h"
#include "../Utility/ArrayList.h"
#include "../IO/IO.h"
#include "../IO/IO_vtk.h"
#include "../IO/Gnuplot.h"
#include "../Algebra/Arithmetic.h"
#include "../Algebra/Interpolation.h"

namespace Larus {

template<typename NODE>
void visit_new_array_on_center_leaf(NODE* node, utPointer pt) {
	if (condition_is_leaf(node)) {
		LarusDef::size_type len = (*CAST(LarusDef::size_type*, pt));
		if (node->data == NULL_PTR) {
			node->data = new typename NODE::Data_type(len);
		} else {
			node->data->aCenterData.reconstruct(len);
		}
	}
}

void new_array_on_center_leaf(pQuadTree qtt, LarusDef::size_type len) { //2D
	if (len <= 0) {
		return;
	}
	qtt->Traversal(visit_new_array_on_center_leaf, &len);
}
void new_array_on_center_leaf(pOCTree oct, LarusDef::size_type len) {
	if (len <= 0) {
		return;
	}
	oct->Traversal(visit_new_array_on_center_leaf, &len);
}

template<typename NODE>
void visit_resize_array_on_center_leaf(NODE* node, utPointer pt) {
	if (condition_is_leaf(node)) {
		LarusDef::size_type len = (*CAST(LarusDef::size_type*, pt));
		if (node->data == NULL_PTR) {
			node->data = new typename NODE::Data_type(len);
		} else {
			node->data->aCenterData.resize(len);
		}
	}
}

void resize_array_on_center_leaf(pOCTree pt, LarusDef::size_type len) { //3D
	if (len <= 0) {
		return;
	}
	pt->Traversal(visit_resize_array_on_center_leaf, &len);
}
void resize_array_on_center_leaf(pQuadTree pt, LarusDef::size_type len) { //2D
	if (len <= 0) {
		return;
	}
	pt->Traversal(visit_resize_array_on_center_leaf, &len);
}
void resize_array_on_center_leaf(Forest2D& forest, LarusDef::size_type len) { //2D
	if (len <= 0) {
		return;
	}
	for (LarusDef::size_type i = 0; i < forest.iLen(); ++i) {
		for (LarusDef::size_type j = 0; j < forest.jLen(); ++j) {
			if (forest.get_attribution(i, j) == ATT_ENABLE) {
				pQuadTree pt = forest.getpTree(i, j);
				resize_array_on_center_leaf(pt, len);
			}
		}
	}
}
void resize_array_on_center_leaf(Forest3D& forest, LarusDef::size_type len) { //3D
	if (len <= 0) {
		return;
	}
	for (LarusDef::size_type i = 0; i < forest.iLen(); ++i) {
		for (LarusDef::size_type j = 0; j < forest.jLen(); ++j) {
			for (LarusDef::size_type k = 0;
					k < (Forest3D::Dim == 3 ? forest.kLen() : 1); ++k) {
				if (forest.get_attribution(i, j, k) == ATT_ENABLE) {
					pOCTree pt = forest.getpTree(i, j, k);
					resize_array_on_center_leaf(pt, len);
				}
			}
		}
	}
}

template<typename NODE>
void visit_plus_scalar_on_leaf(NODE* node, utPointer pt) {
	if (condition_is_leaf(node)) {
		if (node->data == NULL_PTR) {
			std::cerr << " >! Node Data == NULL \n";
			return;
		}
		arrayListT<utPointer>& arr = (*CAST(arrayListT<utPointer>*, pt));
		arrayList_st& arridx = (*CAST(arrayList_st*, arr[0]));
		arrayList& arrval = (*CAST(arrayList*, arr[1]));
		for (int i = 0; i < arridx.size(); ++i) {
			node->data->aCenterData[arridx[i]] += arrval[i];
		}
	}
}

void plus_scalar_on_leaf( // 2D Forest
		pQuadTree ptree,  //pQuadTree
		arrayList_st& arridx,  //data index
		arrayList& arrval   //data plus
		) {
	arrayListT<utPointer> arrutp(2);
	arrutp[0] = &arridx;
	arrutp[1] = &arrval;
	ptree->Traversal(visit_plus_scalar_on_leaf, &arrutp);
}

void plus_scalar_on_leaf( // 3D Forest
		pOCTree ptree,  //pQuadTree
		arrayList_st& arridx,  //data index
		arrayList& arrval   //data plus
		) {
	arrayListT<utPointer> arrutp(2);
	arrutp[0] = &arridx;
	arrutp[1] = &arrval;
	ptree->Traversal(visit_plus_scalar_on_leaf, &arrutp);
}

void plus_scalar_on_leaf( // 2D Forest
		Forest2D& forest,  //pQuadTree
		arrayList_st& arridx,  //data index
		arrayList& arrval   //data plus
		) {
	for (LarusDef::size_type i = 0; i < forest.iLen(); ++i) {
		for (LarusDef::size_type j = 0; j < forest.jLen(); ++j) {
			if (forest.get_attribution(i, j) == ATT_ENABLE) {
				pQuadTree pt = forest.getpTree(i, j);
				plus_scalar_on_leaf(pt, arridx, arrval);
			}
		}
	}
}

void plus_scalar_on_leaf( // 3D Forest
		Forest3D& forest,  //pQuadTree
		arrayList_st& arridx,  //data index
		arrayList& arrval   //data plus
		) {
	for (LarusDef::size_type i = 0; i < forest.iLen(); ++i) {
		for (LarusDef::size_type j = 0; j < forest.jLen(); ++j) {
			for (LarusDef::size_type k = 0;
					k < (Forest3D::Dim == 3 ? forest.kLen() : 1); ++k) {
				if (forest.get_attribution(i, j, k) == ATT_ENABLE) {
					pOCTree pt = forest.getpTree(i, j, k);
					plus_scalar_on_leaf(pt, arridx, arrval);
				}
			}
		}
	}
}

template<typename NODE>
void visit_set_scalar_on_leaf_by_function(NODE* node, utPointer pt) {
	if (condition_is_leaf(node)) {
		if (node->data == NULL_PTR) {
			std::cerr << " >! Node Data == NULL \n";
			return;
		}
		arrayListT<utPointer>& arr = (*CAST(arrayListT<utPointer>*, pt));
		LarusDef::size_type& idx = (*CAST(LarusDef::size_type*, arr[0]));
		scalar_pfun& pf = (*CAST(scalar_pfun*, arr[1]));
		if (NODE::DIM == 2) {
			node->data->aCenterData[idx] = pf(node->cell->getCPX(),
					node->cell->getCPY(), 0);
		} else {
			node->data->aCenterData[idx] = pf(node->cell->getCPX(),
					node->cell->getCPY(), node->cell->getCPZ());
		}
	}
}

template<typename NODE>
void visit_plus_scalar_on_leaf_by_function(NODE* node, utPointer pt) {
	if (condition_is_leaf(node)) {
		if (node->data == NULL_PTR) {
			std::cerr << " >! Node Data == NULL \n";
			return;
		}
		arrayListT<utPointer>& arr = (*CAST(arrayListT<utPointer>*, pt));
		LarusDef::size_type& idx = (*CAST(LarusDef::size_type*, arr[0]));
		scalar_pfun& pf = (*CAST(scalar_pfun*, arr[1]));
		if (NODE::DIM == 2) {
			node->data->aCenterData[idx] += pf(node->cell->getCPX(),
					node->cell->getCPY(), 0);
		} else {
			node->data->aCenterData[idx] += pf(node->cell->getCPX(),
					node->cell->getCPY(), node->cell->getCPZ());
		}
	}
}

void set_scalar_on_leaf_by_function( // 2D tree
		pQuadTree ptree,  //pQuadTree
		LarusDef::size_type idx,  //data index
		scalar_pfun pf      //data plus
		) {
	arrayListT<utPointer> arrutp(2);
	arrutp[0] = &idx;
	arrutp[1] = &pf;
	ptree->Traversal(visit_set_scalar_on_leaf_by_function, &arrutp);
}

void set_scalar_on_leaf_by_function( // 2D Forest
		Forest2D& forest,  //pQuadTree
		LarusDef::size_type idx,  //data index
		scalar_pfun pf      //data plus
		) {
	for (LarusDef::size_type i = 0; i < forest.size(); ++i) {
		if (forest.get_attribution_1d(i) == ATT_ENABLE) {
			pQuadTree pt = forest.getpTree_1d(i);
			set_scalar_on_leaf_by_function(pt, idx, pf);
		}
	}
}

void plus_scalar_on_leaf_by_function( // 2D tree
		pQuadTree ptree,  //pQuadTree
		LarusDef::size_type idx,  //data index
		scalar_pfun pf      //data plus
		) {
	arrayListT<utPointer> arrutp(2);
	arrutp[0] = &idx;
	arrutp[1] = &pf;
	ptree->Traversal(visit_plus_scalar_on_leaf_by_function, &arrutp);
}
void plus_scalar_on_leaf_by_function( // 2D Forest
		Forest2D& forest,  //pQuadTree
		LarusDef::size_type idx,  //data index
		scalar_pfun pf      //data plus
		) {
	for (LarusDef::size_type i = 0; i < forest.size(); ++i) {
		if (forest.get_attribution_1d(i) == ATT_ENABLE) {
			pQuadTree pt = forest.getpTree_1d(i);
			plus_scalar_on_leaf_by_function(pt, idx, pf);
		}
	}
}

template<typename NODE>
void visit_set_const_on_center_leaf(NODE* node, utPointer pt) {
	if (condition_is_leaf(node)) {
		arrayList a = (*CAST(arrayList*, pt));
		if (node->data == NULL_PTR) {
			std::cerr << " >! Node Data == NULL \n";
			return;
		}
		if (node->data->aCenterData.size() <= a[1]) {
			std::cerr << " >! aCenterData.size() = "
					<< node->data->aCenterData.size() << "\n";
			std::cerr << " >! NO place for value \n";
			return;
		}
		node->data->aCenterData[int(a[1])] = a[0];
	}
}

void set_const_on_center_leaf(pOCTree oct, Float f, LarusDef::size_type idx) {
	if (idx < 0) {
		return;
	}
	arrayList a(2);
	a[0] = f;
	a[1] = idx;
	oct->Traversal(visit_set_const_on_center_leaf, &a);
}
template<typename NODE>
void visit_set_const_on_center_leaf_array(NODE* node, utPointer pt) {
	if (condition_is_leaf(node)) {
		arrayListT<utPointer> &a = (*CAST(arrayListT<utPointer>*, pt));
		arrayList& arr = *(CAST(arrayList*, a[0]));
		LarusDef::size_type idx_bt = *(CAST(LarusDef::size_type*, a[1]));
		if (!node->data->aCenterData.check_idx(idx_bt)) {
			std::cerr << " >! aCenterData.size() = "
					<< node->data->aCenterData.size() << "\n";
			std::cerr << " >! NO place for value \n";
			return;
		}
		LarusDef::size_type newsize = arr.size() + idx_bt;
		if (newsize > node->data->aCenterData.size()) {
			node->data->aCenterData.resize(newsize);
		}
		for (LarusDef::size_type i = idx_bt; i < newsize; ++i) {
			node->data->aCenterData[i] = arr[i - idx_bt];
		}
	}
}
void set_const_on_center_leaf(pOCTree oct, arrayList& arr,
		LarusDef::size_type idx_bt) {
	if (arr.non_Empty()) {
		arrayListT<utPointer> a(2);
		a[0] = &arr;
		a[1] = &idx_bt;
		oct->Traversal(visit_set_const_on_center_leaf_array, &a);
	} else {
		return;
	}
}

template<typename NODE>
void visit_set_index_on_center_leaf(NODE* node, utPointer pt) {
	if (condition_is_leaf(node)) {
		arrayListT<utPointer> &arr = *(CAST(arrayListT<utPointer>*, pt));
		LarusDef::size_type& idx_idx = *(CAST(LarusDef::size_type*, arr[0]));
		LarusDef::size_type& idx = *(CAST(LarusDef::size_type*, arr[1]));
		if (!node->data->aCenterData.check_idx(idx_idx)) {
			std::cerr << " >! aCenterData.size() = "
					<< node->data->aCenterData.size();
			std::cerr << " idx = " << idx_idx << "\n";
			std::cerr << " >! NO place for value \n";
			return;
		}
		node->data->aCenterData[idx_idx] = idx;
		idx++;
	}
}

void set_index_on_center_leaf(Forest2D& forest, LarusDef::size_type idx_idx) {
	LarusDef::size_type idx = 0;
	arrayListT<utPointer> arr(2);
	arr[0] = &idx_idx;
	arr[1] = &idx;
	for (LarusDef::size_type i = 0; i < forest.iLen(); ++i) {
		for (LarusDef::size_type j = 0; j < forest.jLen(); ++j) {
			if (forest.get_attribution(i, j) == ATT_ENABLE)
				forest.getpTree(i, j)->Traversal(visit_set_index_on_center_leaf,
						&arr);
		}
	}
}

template<typename NODE>
void visit_set_value_function_on_center_leaf(NODE* node, utPointer pt) {
	if (condition_is_leaf(node)) {
		arrayListT<utPointer> &arr = *(CAST(arrayListT<utPointer>*, pt));
		LarusDef::size_type& idx_idx = *(CAST(LarusDef::size_type*, arr[0]));
		pFun_value* pfun = CAST(pFun_value*, arr[1]);

		if (!node->data->aCenterData.check_idx(idx_idx)) {
			std::cerr << " >! aCenterData.size() = "
					<< node->data->aCenterData.size();
			std::cerr << " idx = " << idx_idx << "\n";
			std::cerr << " >! NO place for value \n";
			return;
		}
		typename NODE::Cell_type::Point cp = node->cell->getCenterPoint();
		if (NODE::DIM == 2) {
			node->data->aCenterData[idx_idx] = (*pfun)(arr[2], cp[0], cp[1], 0);
		} else {
			node->data->aCenterData[idx_idx] = (*pfun)(arr[2], cp[0], cp[1],
					cp[2]);
		}
	}
}

void set_value_function_on_leaf( //
		Forest2D& forest,    //
		LarusDef::size_type idx, //
		pFun_value pfun,    //
		utPointer utp) {
	arrayListT<utPointer> arr(3);
	arr[0] = &idx;
	arr[1] = &pfun;
	arr[2] = utp;
	for (LarusDef::size_type i = 0; i < forest.iLen(); ++i) {
		for (LarusDef::size_type j = 0; j < forest.jLen(); ++j) {
			if (forest.get_attribution(i, j) == ATT_ENABLE)
				forest.getpTree(i, j)->Traversal(
						visit_set_value_function_on_center_leaf, &arr);
		}
	}

}
void set_index_on_center_leaf(Forest3D& forest, LarusDef::size_type idx_idx) {
	LarusDef::size_type idx = 0;
	arrayListT<utPointer> arr(2);
	arr[0] = &idx_idx;
	arr[1] = &idx;
	for (LarusDef::size_type i = 0; i < forest.iLen(); ++i) {
		for (LarusDef::size_type j = 0; j < forest.jLen(); ++j) {
			for (LarusDef::size_type k = 0; k < forest.kLen(); ++k) {
				if (forest.get_attribution(i, j, k) == ATT_ENABLE)
					forest.getpTree(i, j, k)->Traversal(
							visit_set_index_on_center_leaf, &arr);
			}
		}
	}
}

void set_value_on_center_leaf(pOCTree oct, OcTree::pFun_SPTree pfun,
		utPointer utp) {
	oct->Traversal(pfun, utp);
}

//get value======================================
//===3D==========================================
template<typename NODE>
void visit_getListpNode_leaf_center_data_in_range(NODE* node, utPointer pt) {
	if (condition_is_leaf(node)) {
		arrayListT<utPointer> arr = (*CAST(arrayListT<utPointer>*, pt));
		LarusDef::size_type idx = *(CAST(LarusDef::size_type*, arr[1]));
		Float min = *(CAST(Float*, arr[2]));
		Float max = *(CAST(Float*, arr[3]));
		TYPE_Range range = *(CAST(TYPE_Range*, arr[4]));
		if (isInRange(min, node->data->aCenterData[idx], max, range)) {
			ListT<NODE*>* list = CAST(ListT<NODE*>*, arr[0]);
			list->push_back(node);
		}
	}
}

void getListpNode_leaf_center_data_in_range( //
		ListT<pOCNode>& listnode, //as output
		pOCTree pt, // ptree
		LarusDef::size_type idx, //data index
		Float min, //min value
		Float max, //max value
		TYPE_Range range) {
	if (idx < 0) {
		return;
	}
	arrayListT<utPointer> a(5);
	a[0] = &listnode;
	a[1] = &idx;
	a[2] = &min;
	a[3] = &max;
	a[4] = &range;
	pt->Traversal(visit_getListpNode_leaf_center_data_in_range, &a);
}

void getListpNode_leaf_center_data_in_range( //
		ListT<pOCNode>& listnode, //as output
		Forest3D& forest, // ptree
		LarusDef::size_type idx, //data index
		Float min, //min value
		Float max, //max value
		TYPE_Range range) {
	arrayListT<utPointer> a(5);
	a[0] = &listnode;
	a[1] = &idx;
	a[2] = &min;
	a[3] = &max;
	a[4] = &range;
	for (LarusDef::size_type i = 0; i < forest.iLen(); ++i) {
		for (LarusDef::size_type j = 0; j < forest.jLen(); ++j) {
			for (LarusDef::size_type k = 0; k < forest.kLen(); ++k) {
				if (forest.get_attribution(i, j, k) == ATT_ENABLE)
					forest.getpTree(i, j, k)->Traversal(
							visit_getListpNode_leaf_center_data_in_range, &a);
			}
		}
	}
}

template<typename NODE>
void visit_getListpNode_leaf_on_line(NODE* node, utPointer pt) {
	if (condition_is_leaf(node)) {
		arrayListT<utPointer> arr = (*CAST(arrayListT<utPointer>*, pt));
		Float& loc = *(CAST(Float*, arr[1]));
		CSAxis& axi = *(CAST(CSAxis*, arr[2]));
		if (isInRange(node->cell->get(axi, eCPL_M), loc,
				node->cell->get(axi, eCPL_P), Range_co)) {
			ListT<NODE*>* list = CAST(ListT<NODE*>*, arr[0]);
			list->push_back(node);
		}
	}
}

void getListpNode_leaf_on_line( // 2D Forest
		ListT<pQTNode>& listnode,           //as output
		Forest2D& forest,                  // Forest
		Float loc, CSAxis axi) {
	arrayListT<utPointer> a(3);
	a[0] = &listnode;
	a[1] = &loc;
	a[2] = &axi;
	for (LarusDef::size_type i = 0; i < forest.iLen(); ++i) {
		for (LarusDef::size_type j = 0; j < forest.jLen(); ++j) {
			if (forest.get_attribution(i, j) == ATT_ENABLE)
				forest.getpTree(i, j)->Traversal(
						visit_getListpNode_leaf_on_line, &a);
		}
	}
}

void getListpNode_leaf_center_data_in_range( // 2D Forest
		ListT<pQTNode>& listnode, //as output
		Forest2D& forest, // Forest
		LarusDef::size_type idx, //data index
		Float min, //min value
		Float max, //max value
		TYPE_Range range) {
	arrayListT<utPointer> a(5);
	a[0] = &listnode;
	a[1] = &idx;
	a[2] = &min;
	a[3] = &max;
	a[4] = &range;
	for (LarusDef::size_type i = 0; i < forest.iLen(); ++i) {
		for (LarusDef::size_type j = 0; j < forest.jLen(); ++j) {
			if (forest.get_attribution(i, j) == ATT_ENABLE)
				forest.getpTree(i, j)->Traversal(
						visit_getListpNode_leaf_center_data_in_range, &a);
		}
	}
}

//interpolate====================================
//===============================================
template<typename NODE>
int T_neighbor_find_adj( //
		NODE* pn,    //node
		CSAxis axis, //
		NODE*& pm, NODE*& pp //two neighbor as output
		) {
	if (axis == ErrCSAxis || pn == NULL_PTR) {
		return -3;
	}
	if (axis == CSAxis_X) {
		pp = pn->getNeighborFast(SPD_IP);
		pm = pn->getNeighborFast(SPD_IM);
	} else if (axis == CSAxis_Y) {
		pp = pn->getNeighborFast(SPD_JP);
		pm = pn->getNeighborFast(SPD_JM);
	} else if (axis == CSAxis_Z) {
		pp = pn->getNeighborFast(SPD_KP);
		pm = pn->getNeighborFast(SPD_KM);
	}
	if (pp == NULL_PTR && pm == NULL_PTR) {
		return 0;
	}
	if (pp == NULL_PTR) {
		return -1;
	}
	if (pm == NULL_PTR) {
		return 1;
	}
	return 2;
}

template<typename NODE>
int T_find_neighbor_f1( //
		NODE* pc,    //node
		SPDirection direction, //direction
		NODE*& pn    //one neighbor as output
		) {
	if (direction == ErrSPDirection || pn == NULL_PTR) {
		return -3;
	}
	pn = pc->getNeighborFast(direction);
	if (pn == NULL_PTR) {
		return 0;
	}
	return 1;
}

template<typename NODE>
int T_neighbor_find_plane( //
		NODE* pn,    //node
		SPDirection direction, //
		NODE*& pm, NODE*& pp //two neighbor as output
		) {
//direction 0 to 4, 10 to 13  14 to 17
	if (direction == ErrSPDirection || pn == NULL_PTR) {
		return -3;
	}
	pp = pn->getNeighborFast(direction);
	pm = pn->getNeighborFast(oppositeDirection(direction));
	if (pp == NULL_PTR && pm == NULL_PTR) {
		return 0;
	}
	if (pp == NULL_PTR) {
		return -1;
	}
	if (pm == NULL_PTR) {
		return 1;
	}
	return 2;
}

template<typename NODE>
bool T_distance_check( //
		NODE* pn,        //node
		CSAxis axis,     //axix
		Float dis        //distance to center of pn
		) {
	switch (axis) {
	case ErrCSAxis:
		return false;
	case CSAxis_X: {
		if (abs(dis) > pn->cell->gethDx()) {
			return false;
		} else {
			return true;
		}
	}
	case CSAxis_Y: {
		if (abs(dis) > pn->cell->gethDy()) {
			return false;
		} else {
			return true;
		}
	}
	case CSAxis_Z: {
		if (NODE::DIM == 3) {
			if (abs(dis) > pn->cell->gethDz()) {
				return false;
			} else {
				return true;
			}
		} else {
			return false;
		}
	}
	default:
		return false;
	}
}
template<typename NODE>
void visit_averange_value_from_leaf(NODE* pn, utPointer up) {
	if (condition_is_leaf(pn)) {
		arrayListT<utPointer> &arrp = (*CAST(arrayListT<utPointer>*, up));
		int& count = (*CAST(int*, arrp[0]));
		arrayList_st& arridx = (*CAST(arrayList_st*, arrp[1]));
		arrayList& tmp_arrres = (*CAST(arrayList*, arrp[2]));
		if (pn->data != NULL_PTR) {
			for (int i = 0; i < arridx.size(); ++i) {
				tmp_arrres[i] += pn->data->aCenterData[arridx[i]];
			}
			count = count + 1;
		}
	}
}
template<typename NODE>
void T_cal_averange_value_from_leaf(        //
		NODE* pn,  //pNode
		arrayList_st& arridx,  //data index
		arrayList& tmp_arrres //tmp_centerdata
		) {
//improve for special case
	if (condition_is_leaf(pn)) {
		for (int i = 0; i < arridx.size(); ++i) {
			tmp_arrres[i] = pn->data->aCenterData[arridx[i]];
		}
		return;
	}
//========================
	tmp_arrres.assign(0.0);
	int count = 0;
	arrayListT<utPointer> arrp(3);
	arrp[0] = &count;
	arrp[1] = &arridx;
	arrp[2] = &tmp_arrres;
	pn->Traversal(visit_averange_value_from_leaf, &arrp);
	for (int i = 0; i < arridx.size(); ++i) {
		tmp_arrres[i] /= count;
	}
}

template<typename NODE>
void visit_averange_expression_from_leaf(NODE* pn, utPointer up) {
	if (condition_is_leaf(pn)) {
		arrayListT<utPointer> &arrp = (*CAST(arrayListT<utPointer>*, up));
		int& count = (*CAST(int*, arrp[0]));
		ListT<NODE*>& ln = (*CAST(ListT<NODE*>*, arrp[1]));
		if (pn->data != NULL_PTR) {
			ln.push_back(pn);
			count = count + 1;
		}
	}
}

template<typename NODE>
void T_cal_averange_expression_from_leaf(        //
		NODE* pn,  //pNode
		Expression& exp //Expression
		) {
//improve for special case
	if (condition_is_leaf(pn)) {
		ExpTerm et(int(pn->data->aCenterData[Idx_IDX]), 1.0);
		exp.Insert(et);
		return;
	}
//========================
	int count = 0;
	ListT<NODE*> listnode;
	arrayListT<utPointer> arrp(3);
	arrp[0] = &count;
	arrp[1] = &listnode;
	pn->Traversal(visit_averange_expression_from_leaf, &arrp);
	for (typename ListT<NODE*>::iterator iter = listnode.begin();
			iter != listnode.end(); iter++) {
		ExpTerm et(int((*iter)->data->aCenterData[Idx_IDX]), 1.0 / count);
		exp.Insert(et);
	}
}

const int _CHILD_LOOKUP_2D[2][2][2] =   //
		//
		{ { { 0, 1 }, { 2, 3 } }, //
				{ { 1, 3 }, { 0, 2 } } };

void _2_node_interplolate_case3( //
		pQTNode po,  //node
		pQTNode ps,  //node
		int dir,     //dir = 1 or -1
		CSAxis axis, //axix
		arrayList_st& arridx,  //data index
		Point2D& tmp_cpoint,        //point
		arrayList& tmp_arrres   //data res
		) {
	int ch1 = _CHILD_LOOKUP_2D[axis][(dir == 1) ? 0 : 1][0];
	int ch2 = _CHILD_LOOKUP_2D[axis][(dir == 1) ? 0 : 1][1];
	if (ps->hasChild(ch1) && ps->hasChild(ch2)) {   //Normal case
		arrayList tmp_cd1(tmp_arrres.size());
		arrayList tmp_cd2(tmp_arrres.size());
		pQTNode pc1 = ps->child[ch1];
		pQTNode pc2 = ps->child[ch2];
		T_cal_averange_value_from_leaf(pc1, arridx, tmp_cd1);
		T_cal_averange_value_from_leaf(pc2, arridx, tmp_cd2);
		tmp_cpoint = calMid(pc1->cell->getCenterPoint(),
				pc2->cell->getCenterPoint());
		tmp_arrres = (tmp_cd1 + tmp_cd2) * 0.5;
	} else {   //special case
		T_cal_averange_value_from_leaf(ps, arridx, tmp_arrres);
		tmp_cpoint = ps->cell->getCenterPoint();
	}
}

void _2_node_expression_interplolate_case3( //
		pQTNode po,  //node
		pQTNode ps,  //node
		int dir,     //dir = 1 or -1
		CSAxis axis, //axix
		Point2D& tmp_cpoint,        //point
		Expression& exp  //expression
		) {
	int ch1 = _CHILD_LOOKUP_2D[axis][(dir == 1) ? 0 : 1][0];
	int ch2 = _CHILD_LOOKUP_2D[axis][(dir == 1) ? 0 : 1][1];
	if (ps->hasChild(ch1) && ps->hasChild(ch2)) {   //Normal case
		Expression tmp_exp1, tmp_exp2;
		pQTNode pc1 = ps->child[ch1];
		pQTNode pc2 = ps->child[ch2];
		T_cal_averange_expression_from_leaf(pc1, tmp_exp1);
		T_cal_averange_expression_from_leaf(pc2, tmp_exp2);
		tmp_cpoint = calMid(pc1->cell->getCenterPoint(),
				pc2->cell->getCenterPoint());
		tmp_exp1.plus(tmp_exp2); //tmp_exp1+tmp_exp2;
		tmp_exp1.times(0.5);     //(tmp_exp1+tmp_exp2)*0.5;
		exp.plus(tmp_exp1);
	} else { //special case
		Expression tmp_exp;
		T_cal_averange_expression_from_leaf(ps, tmp_exp);
		tmp_cpoint = ps->cell->getCenterPoint();
		exp.plus(tmp_exp);
	}
}

int _2_node_relation( //
		pQTNode po,  //node
		pQTNode ps,  //node
		int dir,     //dir = 1 or -1
		CSAxis axis, //axix
		Point2D& point, //tmp_point
		arrayList_st& arridx,  //data index
		arrayList& arrres   //data res
		) {
	if (po->getLevel() == ps->getLevel() && ps->hasChild() == false) {
		//case 1  on the same level
		point.x = ps->cell->getCenterPoint()[CSAxis_X];
		point.y = ps->cell->getCenterPoint()[CSAxis_Y];
		for (int i = 0; i < arridx.size(); ++i) {
			arrres[i] = ps->data->aCenterData[arridx[i]];
		}
		return 1;
	}
	if (po->getLevel() > ps->getLevel()) {
		//case 2  neighber is coarse
		if (axis == CSAxis_X) {
			point.x = ps->cell->getCenterPoint()[CSAxis_X];
			point.y = po->cell->getCenterPoint()[CSAxis_Y];
			Float tmp_dis = point.y - ps->cell->getCenterPoint().y;
			_interpolate_on_axis( // 2D QuadTree Node
					ps,       //node
					CSAxis_Y, //axix
					tmp_dis,  //distance to center of pn
					arridx,   //data index
					arrres    //data res
					);
		} else { // axis== CSAxis_Y
			point.x = po->cell->getCenterPoint()[CSAxis_X];
			point.y = ps->cell->getCenterPoint()[CSAxis_Y];
			Float tmp_dis = point.x - ps->cell->getCenterPoint().x;
			_interpolate_on_axis( // 2D QuadTree Node
					ps,       //node
					CSAxis_X, //axix
					tmp_dis, //distance to center of pn
					arridx,  //data index
					arrres   //data res
					);
		}
		return 2;
	}
	if (po->getLevel() == ps->getLevel()) {
		//case 3  neighber is fine
		arrayList tmp_arrres(arrres.size());
		_2_node_interplolate_case3( //
				po,       //node
				ps,       //node
				dir,     //dir = 1 or -1
				axis,    //axix
				arridx,  //data index
				point,   //point
				arrres   //data res
				);
		return 3;
	}
	return -1;
}
int _2_node_relation_expression( //
		pQTNode po,     //node
		pQTNode ps,     //node
		int dir,        //dir = 1 or -1
		CSAxis axis,    //axix
		Point2D& point, //tmp_point
		Expression& exp //expression
		) {
	if (po->getLevel() == ps->getLevel() && ps->hasChild() == false) {
		//case 1
		point.x = ps->cell->getCenterPoint()[CSAxis_X];
		point.y = ps->cell->getCenterPoint()[CSAxis_Y];
		ExpTerm term(ps->data->aCenterData[Idx_IDX], 1.0);
		exp.Insert(term);
		return 1;
	}
	if (po->getLevel() > ps->getLevel()) {
		//case 2
		if (axis == CSAxis_X) {
			point.x = ps->cell->getCenterPoint()[CSAxis_X];
			point.y = po->cell->getCenterPoint()[CSAxis_Y];
			Float tmp_dis = point.y - ps->cell->getCenterPoint().y;
			interpolate_expression_on_axis( // 2D QuadTree Node
					ps,       //node
					CSAxis_Y, //axix
					tmp_dis,  //distance to center of pn
					exp);
		} else { // axis== CSAxis_Y
			point.x = po->cell->getCenterPoint()[CSAxis_X];
			point.y = ps->cell->getCenterPoint()[CSAxis_Y];
			Float tmp_dis = point.x - ps->cell->getCenterPoint().x;
			interpolate_expression_on_axis( // 2D QuadTree Node
					ps,       //node
					CSAxis_X, //axix
					tmp_dis, //distance to center of pn
					exp);
		}
		return 2;
	}
	if (po->getLevel() == ps->getLevel()) {
		//case 3
		_2_node_expression_interplolate_case3( //
				po,    //node
				ps,    //node
				dir,   //dir = 1 or -1
				axis,  //axix
				point,  //point
				exp);
		return 3;
	}
	return -1;
}

template<typename NODE>
void T_1_node_interplolate( //
		NODE* po,  //node
		arrayList_st& arridx,  //data index
		arrayList& arrres   //data res
		) {
	for (int i = 0; i < arridx.size(); ++i) {
		arrres[i] = po->data->aCenterData[arridx[i]];
	}
}

template<typename NODE>
void T_expression_1_node_interplolate( //
		NODE* po,  //node
		Expression& exp) {
	ExpTerm et(int(po->data->aCenterData[Idx_IDX]), 1.0);
	exp.Insert(et);
}

void _2_node_interplolate( //
		pQTNode po,  //node
		pQTNode ps,  //node
		int dir,     //dir = 1 or -1
		CSAxis axis, //axix
		Float dis,   //distance to center of pn
		arrayList_st& arridx,  //data index
		arrayList& arrres   //data res
		) {
	Point2D tmp_point;
	arrayList tmp_arrres(arridx.size());
	_2_node_relation( //
			po,    //node
			ps,    //node
			dir,   //dir = 1 or -1
			axis,  //axix
			tmp_point, //tmp_point
			arridx,  //data index
			tmp_arrres   //data res
			);
	for (int i = 0; i < arridx.size(); ++i) {
		Float x = po->cell->getCenterPoint()[axis] + dis;
		Float x1 = po->cell->getCenterPoint()[axis];
		Float y1 = po->data->aCenterData[arridx[i]];
		Float x2 = tmp_point[axis];
		Float y2 = tmp_arrres[i];
		arrres[i] = linear_interpolation(x, x1, y1, x2, y2);
	}
}

void _2_node_gradient( //
		pQTNode po,  //node
		pQTNode ps,  //node
		int dir,     //dir = 1 or -1
		CSAxis axis, //axix
		arrayList_st& arridx,  //data index
		arrayList& arrres   //data res
		) {
	Point2D tmp_point;
	arrayList tmp_arrres(arridx.size());
	_2_node_relation(po, ps, dir, axis, tmp_point, arridx, tmp_arrres);
	for (int i = 0; i < arridx.size(); ++i) {
		Float x1 = po->cell->getCenterPoint()[axis];
		Float y1 = po->data->aCenterData[arridx[i]];
		Float x2 = tmp_point[axis];
		Float y2 = tmp_arrres[i];
		arrres[i] = linear_gradient(x1, y1, x2, y2);
	}
}

void _expression_2_node_interplolate( //
		pQTNode po,  //node
		pQTNode ps,  //node
		int dir,     //dir = 1 or -1
		CSAxis axis, //axix
		Float dis,   //distance to center of pn
		Expression& exp) {
	Point2D tmp_point;
	Expression tmp_exp;
	_2_node_relation_expression( //
			po,    //node
			ps,    //node
			dir,   //dir = 1 or -1
			axis,  //axix
			tmp_point, //tmp_point
			tmp_exp);
	Float x = po->cell->getCenterPoint()[axis] + dis;
	Float x1 = po->cell->getCenterPoint()[axis];
//Float y1 = po->data->aCenterData[arridx[i]];
	Float x2 = tmp_point[axis];
//Float y2 = tmp_arrres[i];
	Expression y1;
	y1.Insert(ExpTerm(1.0, int(po->data->aCenterData[Idx_IDX])));
	exp = linear_interpolation_expression(x, x1, y1, x2, tmp_exp);
}

void _expression_2_node_gradient( //
		pQTNode po,  //node
		pQTNode ps,  //node
		int dir,     //dir = 1 or -1
		CSAxis axis, //axix
		Expression& exp) {
	Point2D tmp_point;
	Expression tmp_exp;
	_2_node_relation_expression( //
			po,    //node
			ps,    //node
			dir,   //dir = 1 or -1
			axis,  //axix
			tmp_point, //tmp_point
			tmp_exp);
//Float x = po->cell->getCenterPoint()[axis] + dis;
	Float x1 = po->cell->getCenterPoint()[axis];
	Float x2 = tmp_point[axis];
	Expression exp1;
	exp1.Insert(ExpTerm(1.0, int(po->data->aCenterData[Idx_IDX])));
	exp = linear_gradient_expression(x1, exp1, x2, tmp_exp);
}

void _3_node_interplolate( //
		pQTNode pm,  //node
		pQTNode po,  //node
		pQTNode pp,  //node
		CSAxis axis, //axix
		Float dis,   //distance to center of pn
		arrayList_st& arridx,  //data index
		arrayList& arrres   //data res
		) {
	Point2D point_m; //tmp_point
	Point2D point_p; //tmp_point
	arrayList arrres_m(arridx.size());   //data res
	arrayList arrres_p(arridx.size());   //data res
	_2_node_relation(po,  //node
			pm,  //node
			-1,     //dir = 1 or -1
			axis, //axix
			point_m, //tmp_point
			arridx,  //data index
			arrres_m   //data res
			);
	_2_node_relation(po,  //node
			pp,  //node
			1,    //dir = 1 or -1
			axis, //axix
			point_p, //tmp_point
			arridx,  //data index
			arrres_p   //data res
			);
	for (int i = 0; i < arridx.size(); ++i) {
		Float x = po->cell->getCenterPoint()[axis] + dis;
		Float x1 = po->cell->getCenterPoint()[axis];
		Float y1 = po->data->aCenterData[arridx[i]];
		Float x2 = point_m[axis];
		Float y2 = arrres_m[i];
		Float x3 = point_p[axis];
		Float y3 = arrres_p[i];
		arrres[i] = second_order_interpolation(x, x1, y1, x2, y2, x3, y3);
	}
}

void _3_node_gradient( //
		pQTNode pm,  //node
		pQTNode po,  //node
		pQTNode pp,  //node
		CSAxis axis, //axix
		Float dis,   //distance to center of pn
		arrayList_st& arridx,  //data index
		arrayList& arrres   //data res
		) {
	Point2D point_m; //tmp_point
	Point2D point_p; //tmp_point
	arrayList arrres_m(arridx.size());   //data res
	arrayList arrres_p(arridx.size());   //data res
	_2_node_relation(po, pm, -1, axis, point_m, arridx, arrres_m);
	_2_node_relation(po, pp, 1, axis, point_p, arridx, arrres_p);
	for (int i = 0; i < arridx.size(); ++i) {
		Float x = po->cell->getCenterPoint()[axis] + dis;
		Float x1 = po->cell->getCenterPoint()[axis];
		Float y1 = po->data->aCenterData[arridx[i]];
		Float x2 = point_m[axis];
		Float y2 = arrres_m[i];
		Float x3 = point_p[axis];
		Float y3 = arrres_p[i];
		arrres[i] = second_order_gradient(x, x1, y1, x2, y2, x3, y3);
	}
}

void _expression_3_node_interplolate( //
		pQTNode pm,  //node
		pQTNode po,  //node
		pQTNode pp,  //node
		CSAxis axis, //axix
		Float dis,   //distance to center of pn
		Expression& exp) {
	Point2D point_m; //tmp_point
	Point2D point_p; //tmp_point
	Expression exp_m;   //data res
	Expression exp_p;   //data res
	_2_node_relation_expression(po,  //node
			pm,  //node
			-1,     //dir = 1 or -1
			axis, //axix
			point_m, //tmp_point
			exp_m   //data res
			);
	_2_node_relation_expression(po,  //node
			pp,  //node
			1,    //dir = 1 or -1
			axis, //axix
			point_p, //tmp_point
			exp_p   //data res
			);
	Float x = po->cell->getCenterPoint()[axis] + dis;
	Float x1 = po->cell->getCenterPoint()[axis];
	Expression y1;
	y1.Insert(ExpTerm(po->data->aCenterData[Idx_IDX], 1.0));
	Float x2 = point_m[axis];
	Float x3 = point_p[axis];
	exp = second_order_interpolation_expression(x, x1, y1, x2, exp_m, x3,
			exp_p);
}

void _expression_3_node_gradient( //
		pQTNode pm,  //node
		pQTNode po,  //node
		pQTNode pp,  //node
		CSAxis axis, //axix
		Float dis,   //distance to center of pn
		Expression& exp) {
	Point2D point_m; //tmp_point
	Point2D point_p; //tmp_point
	Expression exp_m;   //data res
	Expression exp_p;   //data res
	_2_node_relation_expression(po,  //node
			pm,  //node
			-1,     //dir = 1 or -1
			axis, //axix
			point_m, //tmp_point
			exp_m   //data res
			);
	_2_node_relation_expression(po,  //node
			pp,  //node
			1,    //dir = 1 or -1
			axis, //axix
			point_p, //tmp_point
			exp_p   //data res
			);
//exp_p.show();
//exp_m.show();
	Float x = po->cell->getCenterPoint()[axis] + dis;
	Float x1 = po->cell->getCenterPoint()[axis];
	Expression y1;
	y1.Insert(ExpTerm(po->data->aCenterData[Idx_IDX], 1.0));
	Float x2 = point_m[axis];
	Float x3 = point_p[axis];
	exp = second_order_gradient_expression(x, x1, y1, x2, exp_m, x3, exp_p);
}

void _interpolate_on_axis( // 2D QuadTree Node
		pQTNode pn, //node
		CSAxis axis, //axix
		Float dis, //distance to center of pn
		arrayList_st& arridx,  //data index
		arrayList& arrres   //data res
		) {
	_IF_TRUE_RETRUN(pn==NULL);

	if (T_distance_check(pn, axis, dis)) {
		pQTNode pm = NULL_PTR;
		pQTNode pp = NULL_PTR;
		int numnei = T_neighbor_find_adj(pn, axis, pm, pp);
		switch (numnei) {
		case 0: {
			T_1_node_interplolate(pn, arridx, arrres);
			return; //unfinish
		}
		case -1: {
			_2_node_interplolate(pn,  //node
					pm,      //node
					numnei,     //dir = 1 or -1
					axis, //axix
					dis,   //distance to center of pn
					arridx,  //data index
					arrres);
			return;  //unfinish
		}
		case 1: {
			_2_node_interplolate(pn,  //node
					pp,  //node
					numnei,     //dir = 1 or -1
					axis, //axix
					dis,   //distance to center of pn
					arridx,  //data index
					arrres);
			return; //unfinish
		}
		case 2: {
			_3_node_interplolate( //
					pm,  //node
					pn,  //node
					pp,  //node
					axis, //axix
					dis,   //distance to center of pn
					arridx,  //data index
					arrres   //data res
					);
			return; //unfinish
		}
		default: {
			return; //unfinish
		}
		}
	} else {
		std::cerr << " >! Interploate distance is not in node \n";
		return;
	}
}

const SPDirection AXIS_TO_2_DIRECTION[3][2] = { //
												//
		{ SPD_IM, SPD_IP },		//
				{ SPD_JM, SPD_JP },		//
				{ SPD_KM, SPD_KP } //
		};

void _interpolate_1order_on_axis( // 2D QuadTree Node
		pQTNode pn, //node
		CSAxis axis, //axix
		Float dis, //distance to center of pn
		arrayList_st& arridx,  //data index
		arrayList& arrres   //data res
		) {
	_IF_TRUE_RETRUN(pn==NULL);
	_IF_TRUE_RETRUN(axis == ErrCSAxis);
	if (dis == 0) {
		T_1_node_interplolate(pn, arridx, arrres);
		return;
	}
//transfer to direction
	SPDirection spd = AXIS_TO_2_DIRECTION[axis][(dis < 0) ? 0 : 1];
	if (T_distance_check(pn, axis, dis)) {
		pQTNode ps = pn->getNeighborFast(spd);
		if (ps == NULL_PTR) {
			T_1_node_interplolate(pn, arridx, arrres);
			return; //unfinish
		} else {
			_2_node_interplolate(pn,  //node
					ps,  //node
					(dis < 0) ? -1 : 1, //dir = 1 or -1
					axis, //axix
					dis,   //distance to center of pn
					arridx,  //data index
					arrres);
		}
	} else {
		std::cerr << " >! Interploate distance is not in node \n";
		return;
	}
}

void _interpolate_gradient_on_axis( // 2D QuadTree Node
		pQTNode pn, //node
		CSAxis axis, //axix
		Float dis, //distance to center of pn
		arrayList_st& arridx,  //data index
		arrayList& arrres   //data res
		) {
	_IF_TRUE_RETRUN(pn == NULL_PTR);
	if (T_distance_check(pn, axis, dis)) {
		pQTNode pm = NULL_PTR;
		pQTNode pp = NULL_PTR;
		int numnei = T_neighbor_find_adj(pn, axis, pm, pp);
		switch (numnei) {
		case 0: {
			for (int i = 0; i < arridx.size(); ++i) {
				arrres[i] = 0.0;
			}
			return;
		}
		case -1: {
			_2_node_gradient(pn,  //node
					pm,  //node
					numnei,     //dir = 1 or -1
					axis, //axix
					arridx,  //data index
					arrres);
			return;  //unfinish
		}
		case 1: {
			_2_node_gradient(pn,  //node
					pp,  //node
					numnei,     //dir = 1 or -1
					axis, //axix
					arridx,  //data index
					arrres);
			return; //unfinish
		}
		case 2: {
			_3_node_gradient( //
					pm,  //node
					pn,  //node
					pp,  //node
					axis, //axix
					dis,   //distance to center of pn
					arridx,  //data index
					arrres   //data res
					);
			return; //unfinish
		}
		default: {
			return; //unfinish
		}
		}
	} else {
		std::cerr << " >! Interploate distance is not in node \n";
		return;
	}
}

void interpolate_expression_on_axis( // 2D QuadTree Node
		pQTNode pn, //node
		CSAxis axis, //axix
		Float dis, //distance to center of pn
		Expression& exp   //Expression
		) {
	_IF_TRUE_RETRUN(pn == NULL_PTR);
	if (T_distance_check(pn, axis, dis)) {
		pQTNode pm = NULL_PTR;
		pQTNode pp = NULL_PTR;
		int numnei = T_neighbor_find_adj(pn, axis, pm, pp);
		switch (numnei) {
		case 0: {
			T_expression_1_node_interplolate(pn, exp);
			return;
		}
		case -1: {
			_expression_2_node_interplolate(pn,  //node
					pm,  //node
					numnei,     //dir = 1 or -1
					axis, //axix
					dis,   //distance to center of pn
					exp);
			return;  //unfinish
		}
		case 1: {
			_expression_2_node_interplolate(pn,  //node
					pp,  //node
					numnei,     //dir = 1 or -1
					axis, //axix
					dis,   //distance to center of pn
					exp);
			return; //unfinish
		}
		case 2: {
			_expression_3_node_interplolate( //
					pm,  //node
					pn,  //node
					pp,  //node
					axis, //axix
					dis,   //distance to center of pn
					exp);
			return; //unfinish
		}
		default: {
			return; //unfinish
		}
		}
	} else {
		std::cerr << " >! Interploate distance is not in node \n";
		return;
	}
}

void interpolate_expression_gradient_on_axis( // 2D QuadTree Node
		pQTNode pn, //node
		CSAxis axis, //axix
		Float dis, //distance to center of pn
		Expression& exp   //Expression
		) {
	_IF_TRUE_RETRUN(pn == NULL_PTR);
	if (T_distance_check(pn, axis, dis)) {
		pQTNode pm = NULL_PTR;
		pQTNode pp = NULL_PTR;
		int numnei = T_neighbor_find_adj(pn, axis, pm, pp);
		switch (numnei) {
		case 0: {
			exp.clear();
			return;
		}
		case -1: {
			_expression_2_node_gradient(pn,  //node
					pm,  //node
					numnei,     //dir = 1 or -1
					axis, //axix
					exp);
			return;  //unfinish
		}
		case 1: {
			_expression_2_node_gradient(pn,  //node
					pp,  //node
					numnei,     //dir = 1 or -1
					axis, //axix
					exp);
			return; //unfinish
		}
		case 2: {
			_expression_3_node_gradient( //
					pm,  //node
					pn,  //node
					pp,  //node
					axis, //axix
					dis,   //distance to center of pn
					exp);
			return; //unfinish
		}
		default: {
			return; //unfinish
		}
		}
	} else {
		std::cerr << " >! Interploate distance is not in node \n";
		return;
	}
}

bool _trans_face_to_axis_and_dis( // 2D QuadTree Node
		pQTNode pn,        // Node
		SPDirection face,  // Face
		Float& dis,        // [out] Distance
		CSAxis& axis       // [out] Aixs
		) {
	const CSAxis FACE_TO_AIXS[4] = { CSAxis_X, CSAxis_Y, CSAxis_X, CSAxis_Y };
	axis = FACE_TO_AIXS[face - 4];
	if (4 <= face && face <= 7) {
		if (axis == CSAxis_X) {
			dis = (face == SPD_IP) ?
					(pn->cell->get(CSAxis_X, eCPL_P)
							- pn->cell->getCenterPoint()[CSAxis_X]) :
					(pn->cell->get(CSAxis_X, eCPL_M)
							- pn->cell->getCenterPoint()[CSAxis_X]);
		} else {
			dis = (face == SPD_JP) ?
					(pn->cell->get(CSAxis_Y, eCPL_P)
							- pn->cell->getCenterPoint()[CSAxis_Y]) :
					(pn->cell->get(CSAxis_Y, eCPL_M)
							- pn->cell->getCenterPoint()[CSAxis_Y]);
		}
		return true;
	} else {
		std::cerr << " >! Error face direction\n";
		return false;
	}
}

void interpolate_on_face( // 2D QuadTree Node
		pQTNode pn,           //node
		SPDirection face,  //face
		arrayList_st& arridx, //data index
		arrayList& arrres     //data res
		) {
//face 4 to 7
	Float dis;
	CSAxis axis;
	_IF_FALSE_RETRUN(_trans_face_to_axis_and_dis(pn, face, dis, axis));
	_interpolate_on_axis( // 2D QuadTree Node
			pn, //node
			axis, //axix
			dis, //distance to center of pn
			arridx,  //data index
			arrres   //data res
			);
}

void interpolate_1order_on_face( // 2D QuadTree Node
		pQTNode pn,                  //node
		SPDirection face,                //face
		arrayList_st& arridx,              //data index
		arrayList& arrres               //data res
		) {
	Float dis;
	CSAxis axis;
	_IF_FALSE_RETRUN(_trans_face_to_axis_and_dis(pn, face, dis, axis));
	_interpolate_1order_on_axis( // 2D QuadTree Node
			pn, //node
			axis, //axix
			dis, //distance to center of pn
			arridx,  //data index
			arrres   //data res
			);
}

void interpolate_1order_on_face( // 2D QuadTree face
		QTNodeFace& face,                //node
		arrayList_st& arridx,            //data index
		arrayList& arrres                //data res
		) {
	pQTNode pn = face.pnode;
	_IF_TRUE_RETRUN(pn == NULL_PTR);
	Float dis;
	CSAxis axis;
	_IF_FALSE_RETRUN(_trans_face_to_axis_and_dis(pn, face.direction, dis, axis));
	pQTNode ps = face.pneighbor;      //don't search for neighbor
// do not check distance
	if (ps == NULL_PTR) {
		T_1_node_interplolate(pn, arridx, arrres);
		return;
	} else {
		_2_node_interplolate(pn,  //node
				ps,  //node
				(dis < 0) ? -1 : 1, //dir = 1 or -1
				axis, //axix
				dis,   //distance to center of pn
				arridx,  //data index
				arrres);
	}
}

Float interpolate_1order_on_face( // 2D QuadTree Node
		QTNodeFace& face,                      //node
		LarusDef::size_type idx          //data index
		) {
	arrayList_st arridx(1);
	arridx[0] = idx;
	arrayList arrres(1);
	interpolate_1order_on_face(face, arridx, arrres);
	return arrres[0];
}

Float interpolate_1order_on_face( // 2D QuadTree Node
		pQTNode pn,                      //node
		SPDirection face,                //face
		LarusDef::size_type idx            //data index
		) {
	arrayList_st arridx(1);
	arridx[0] = idx;
	arrayList arrres(1);
	interpolate_1order_on_face(pn, face, arridx, arrres);
	return arrres[0];
}

void interpolate_1order_weight_on_face( // 2D QuadTree Node
		pQTNode pn,                      //node
		SPDirection face,                //face
		arrayList_st& arridx,            //data index
		arrayList& arrres                //data res
		) {
	Float dis;
	CSAxis axis;
	_IF_FALSE_RETRUN(_trans_face_to_axis_and_dis(pn, face, dis, axis));
	pQTNode ps = pn->getNeighborFast(face);
	if (ps == NULL_PTR) {
		T_1_node_interplolate(pn, arridx, arrres);
		return; //unfinish
	} else {
		Point2D tmp_point;
		arrayList tmp_arrres(arridx.size());
		_2_node_relation( //
				pn,    //node
				ps,    //node
				(dis < 0) ? -1 : 1, //dir = 1 or -1
				axis,  //axix
				tmp_point, //tmp_point
				arridx,  //data index
				tmp_arrres   //data res
				);
		for (int i = 0; i < arridx.size(); ++i) {
			Float x1 = pn->cell->getD(axis);  //the width of the cell
			Float y1 = pn->data->aCenterData[arridx[i]];
			Float x2 = ps->cell->getD(axis);  //the width of the cell
			Float y2 = tmp_arrres[i];
			arrres[i] = linear_weight_interpolation(x1, y1, x2, y2);
		}
	}
}

void interpolate_gradient_on_face( // 2D QuadTree Node
		pQTNode pn,           //node
		SPDirection face,  //face
		arrayList_st& arridx, //data index
		arrayList& arrres     //data res
		) {
//face 4 to 7
	Float dis;
	CSAxis axis;
	_IF_FALSE_RETRUN(_trans_face_to_axis_and_dis(pn, face, dis, axis));
	_interpolate_gradient_on_axis(pn, axis, dis, arridx, arrres);
}

void interpolate_gradient_at_center( // 2D QuadTree Node
		pQTNode pn,                         //node
		CSAxis axis, arrayList_st& arridx, //data index
		arrayList& arrres                   //data res
		) {
	_IF_TRUE_RETRUN(pn==NULL_PTR);
	arrayList arrres_m(arridx.size());
	arrayList arrres_p(arridx.size());
	interpolate_1order_weight_on_face(pn, AXIS_TO_2_DIRECTION[axis][0], arridx,
			arrres_m);
	interpolate_1order_weight_on_face(pn, AXIS_TO_2_DIRECTION[axis][1], arridx,
			arrres_p);
	for (int i = 0; i < arridx.size(); ++i) {
		Float x1 = pn->cell->get(axis, eCPL_M);  //the width of the cell
		Float y1 = arrres_m[i];
		Float x2 = pn->cell->get(axis, eCPL_P);  //the width of the cell 1!!!!!!
		Float y2 = arrres_p[i];
		arrres[i] = linear_gradient(x1, y1, x2, y2);
	}
}

void interpolate_expression_on_face( // 2D QuadTree Node
		pQTNode pn,                         //node
		SPDirection face,                //face
		Expression& exp   //Expression
		) {
	Float dis;
	CSAxis axis;
	_IF_FALSE_RETRUN(_trans_face_to_axis_and_dis(pn, face, dis, axis));
	interpolate_expression_on_axis( // 2D QuadTree Node
			pn, //node
			axis, //axix
			dis, //distance to center of pn
			exp);
}

void expression_interpolate_gradient_on_face( // 2D QuadTree Node
		pQTNode pn,                         //node
		SPDirection face,                //face
		Expression& exp   //Expression
		) {
	Float dis;
	CSAxis axis;
	_IF_FALSE_RETRUN(_trans_face_to_axis_and_dis(pn, face, dis, axis));
	interpolate_expression_gradient_on_axis( // 2D QuadTree Node
			pn, //node
			axis, //axix
			dis, //distance to center of pn
			exp);
}

int _2_node_relation_vertex( //
		pQTNode po,            //node
		pQTNode ps,            //node
		SPDirection vertex,    //vertex
		Point2D& point,        //tmp_point
		arrayList_st& arridx,  //data index
		arrayList& arrres      //data res
		) {
	if (po->getLevel() == ps->getLevel() && ps->hasChild() == false) {
		//case 1 the two node on the same level
		point = ps->cell->getCenterPoint();
		for (int i = 0; i < arridx.size(); ++i) {
			arrres[i] = ps->data->aCenterData[arridx[i]];
		}
		return 1;
	}
	if (po->getLevel() > ps->getLevel() && ps->hasChild() == false) { //fine to coarse
//case 2 the neighbor node is lower than the original one
		Point2D centerp = ps->cell->getCenterPoint();
		Point2D vertexp = ps->getPoint(vertex);
		point.x = vertexp.x + (vertexp.x - centerp.x);
		point.y = vertexp.y + (vertexp.y - centerp.y);
		//point is point-symmetrical point of vertexp
		for (int i = 0; i < arridx.size(); ++i) {
			arrres[i] = ps->data->aCenterData[arridx[i]];
		}
		return 0;
	}
	if (po->getLevel() == ps->getLevel()) {
		//Coarse to fine
		//case 3 the neighbor has children
		// method 1 ==================
		//arrayList tmp_arrres(arrres.size());
		//T_cal_averange_value_from_leaf(ps, arridx, tmp_arrres);
		//point = ps->cell->getCenterPoint();
		//case 3 method2 =============
		arrayList tmp_arrres(arrres.size());
		pQTNode psn = getChild_vertex(ps, vertex);
		for (int i = 0; i < arridx.size(); ++i) {
			tmp_arrres[i] = psn->data->aCenterData[arridx[i]];
		}
		point = psn->cell->getCenterPoint();
		return 2;
	}
	return -1;
}

void _2_node_interplolate_vertex( //
		pQTNode po,  //node
		pQTNode ps,  //node
		SPDirection vertex, //vertex
		arrayList_st& arridx,  //data index
		arrayList& arrres   //data res
		) {
	Point2D tmp_point;
	arrayList tmp_arrres(arridx.size());
	_2_node_relation_vertex( //
			po,        //node
			ps,        //node
			vertex, tmp_point, //tmp_point
			arridx,    //data index
			tmp_arrres //data res
			);
//interpolate between "po" and "tmp_point" with data
	Point2D p = po->getPoint(vertex);
	Point2D p1 = po->cell->getCenterPoint();
	for (int i = 0; i < arridx.size(); ++i) {
		arrres[i] = linear_interpolation(p, p1,
				po->data->aCenterData[arridx[i]], tmp_point, tmp_arrres[i]);
	}
}

void _3_node_interplolate_vertex( //
		pQTNode pm,  //node
		pQTNode po,  //node
		pQTNode pp,  //node
		SPDirection vertex, //vertex
		arrayList_st& arridx,  //data index
		arrayList& arrres   //data res
		) {
	Point2D point_m; //tmp_point
	Point2D point_p; //tmp_point
	arrayList arrres_m(arridx.size());   //data res
	arrayList arrres_p(arridx.size());   //data res
	_2_node_relation_vertex(   //
			po,  //node
			pm,  //node
			vertex, point_m, //tmp_point
			arridx,  //data index
			arrres_m   //data res
			);
	_2_node_relation_vertex(   //
			po,  //node
			pp,  //node
			vertex, point_p, //tmp_point
			arridx,  //data index
			arrres_p   //data res
			);
	Point2D p = po->getPoint(vertex);
	for (int i = 0; i < arridx.size(); ++i) {
		arrres[i] = second_order_interpolation(p, point_m, arrres_m[i],
				po->cell->getCenterPoint(), po->data->aCenterData[arridx[i]],
				point_p, arrres_p[i]);
	}
}

void _interpolate_on_vetex( // 2D QuadTree Node
		pQTNode pn, //node
		SPDirection vertex, //vertex
		arrayList_st& arridx, //data index
		arrayList& arrres   //data res
		) {
	if (pn == NULL_PTR) {
		return;
	}
	if (0 <= vertex && vertex <= 3) {
		pQTNode pm = NULL_PTR;
		pQTNode pp = NULL_PTR;
		int numnei = T_neighbor_find_plane(pn, vertex, pm, pp);
		switch (numnei) {
		case 0: {
			T_1_node_interplolate(pn, arridx, arrres);
			return; //unfinish
		}
		case -1: {
			_2_node_interplolate_vertex(pn,  //node
					pm,      //node
					vertex,  //vertex
					arridx,  //data index
					arrres   //data res
					);
			return;
		}
		case 1: {
			_2_node_interplolate_vertex(pn,  //node
					pp,      //node
					vertex,  //vertex
					arridx,  //data index
					arrres   //data res
					);
			return;
		}
		case 2: {
			_3_node_interplolate_vertex( //
					pm, //node
					pn, //node
					pp, //node
					vertex, //vertex
					arridx, //data index
					arrres //data res
					);
			return;
		}
		default: {
			return;
		}
		}
	} else {
		std::cerr << " >! Interploate on error Direction \n";
		return;
	}
}

void interpolate_1order_on_vertex( // 2D QuadTree Node
		pQTNode pn,                      //node
		SPDirection vertex,              //vertex
		arrayList_st& arridx,            //data index
		arrayList& arrres                //data res
		) {
	_IF_TRUE_RETRUN(pn==NULL);
	_IF_TRUE_RETRUN(vertex == ErrSPDirection);
	pQTNode ps = pn->getNeighborFast(vertex);
	if (ps == NULL_PTR) {
		T_1_node_interplolate(pn, arridx, arrres);
		return; //unfinish
	} else {
		_2_node_interplolate_vertex( //
				pn,      //node original
				ps,      //node neighbor
				vertex,  //vertex direction
				arridx,  //data index
				arrres);
	}
}

void interpolate_1order_on_vertex_averange_neighbor( // 2D QuadTree Node
		pQTNode pn,                      //node
		SPDirection vertex,              //vertex
		arrayList_st& arridx,            //data index
		arrayList& arrres                //data res
		) {
	_IF_TRUE_RETRUN(pn==NULL);
	_IF_TRUE_RETRUN(vertex == ErrSPDirection);
	ASSERT(isInRange(SPD_MP, vertex, SPD_PM, Range_cc));
	SPDirection d1, d2, d3;
//d3 is useless, vertex decompose to d1 d2
	Direction_Decompose(vertex, d1, d2, d3);
//find the neighbor in three direction
	pQTNode pc = pn->getNeighborFast(vertex);    //neighbor corner
	arrayList arrresoc(arridx.size());
	pc = (pc != NULL_PTR) ?
			getChild_vertex(pc, oppositeDirection(vertex)) : NULL_PTR;
	if (pc != NULL_PTR) {
		Float d1 = pn->cell->getDvertex();
		Float d2 = pc->cell->getDvertex();
		for (int i = 0; i < arridx.size(); ++i) {
			arrresoc[i] = linear_interpolation_center(d1,
					pn->data->aCenterData[arridx[i]], d2,
					pc->data->aCenterData[arridx[i]]);
		}
	} else { //pc==NULL
		for (int i = 0; i < arridx.size(); ++i) {
			arrresoc[i] = pn->data->aCenterData[arridx[i]];
		}
	}
	pQTNode px = pn->getNeighborFast(d1);        //neighbor on X
	pQTNode py = pn->getNeighborFast(d2);        //neighbor on Y
	arrayList arrresxy(arridx.size());
	int flag = 1;
	if (px == NULL_PTR && py == NULL_PTR) {
		flag = 0;
	} else {
		px = (px != NULL_PTR) ?
				getChild_vertex(px,
						Direction_Compose(oppositeDirection(d1), d2)) :
				NULL_PTR;
		py = (py != NULL_PTR) ?
				getChild_vertex(py,
						Direction_Compose(d1, oppositeDirection(d2))) :
				NULL_PTR;
		if (px == NULL_PTR && py == NULL_PTR) {
			flag = 0;
		} else if (px != NULL_PTR && py == NULL_PTR) {
			for (int i = 0; i < arridx.size(); ++i) {
				arrresxy[i] = px->data->aCenterData[arridx[i]];
			}
		} else if (px == NULL_PTR && py != NULL_PTR) {
			for (int i = 0; i < arridx.size(); ++i) {
				arrresxy[i] = py->data->aCenterData[arridx[i]];
			}
		} else {
			Float d1 = px->cell->getDvertex();
			Float d2 = py->cell->getDvertex();
			for (int i = 0; i < arridx.size(); ++i) {
				arrresxy[i] = linear_interpolation_center(d1,
						px->data->aCenterData[arridx[i]], d2,
						py->data->aCenterData[arridx[i]]);
			}
		}
	}
//get averange of xy and co
	if (flag == 0) {
		for (int i = 0; i < arridx.size(); ++i) {
			arrres[i] = arrresoc[i];
		}
	} else {
		for (int i = 0; i < arridx.size(); ++i) {
			arrres[i] = (arrresoc[i] + arrresxy[i]) * 0.5;
		}
	}
}

int _interpolate_node(       // 2D QuadTree Node
		pQTNode pn,       //node
		const Point2D& point,       //point
		arrayList_st& arridx,       //data index
		arrayList& arrres       //data res
		) {
	if (pn->cell->getCenterPoint() == point) { //point equal special case
		for (int i = 0; i < arridx.size(); ++i) {
			arrres[i] = pn->data->aCenterData[arridx[i]];
		}
		return 1;
	}
	if (pn->cell->isInOnCell(point)) {
		arrayList arrres_cor(arridx.size()), arrres_x(arridx.size()), arrres_y(
				arridx.size());
		SPDirection child = toDirection(pn->whichChild(point)); //child from 0 to 3
		Point2D point_c = pn->getPoint(child);
		_interpolate_on_vetex(pn, child, arridx, arrres_cor);
		//arrres_cor.show();
		SPDirection d1, d2, d3;
		Direction_Decompose(child, d1, d2, d3);
		//Point2D point_x=pn->getPoint(d1);
		//Point2D point_y=pn->getPoint(d2);
		interpolate_on_face(pn, d1, arridx, arrres_x);
		//arrres_x.show();
		interpolate_on_face(pn, d2, arridx, arrres_y);
		//arrres_y.show();
		for (int i = 0; i < arridx.size(); ++i) {
			arrres[i] = bilinear_interpolation(point,
					pn->cell->getCenterPoint(), point_c,
					pn->data->aCenterData[arridx[i]], arrres_x[i], arrres_y[i],
					arrres_cor[i]);
		}
		return 1;
	} else {
		std::cerr << " >! Interploate point is not in cell\n";
		return -1;
	}
}

int _interpolate_node_1order( // 2D QuadTree Node
		pQTNode pn, //node
		const Point2D& point, //point
		arrayList_st& arridx, //data index
		arrayList& arrres //data res
		) {
	if (pn->cell->getCenterPoint() == point) { //point equal special case
		for (int i = 0; i < arridx.size(); ++i) {
			arrres[i] = pn->data->aCenterData[arridx[i]];
		}
		return 1;
	}
	if (pn->cell->isInOnCell(point)) {
		arrayList arrres_cor(arridx.size()), arrres_x(arridx.size()), arrres_y(
				arridx.size());
		SPDirection child = toDirection(pn->whichChild(point)); //child from 0 to 3
		Point2D point_c = pn->getPoint(child);
		interpolate_1order_on_vertex(pn, child, arridx, arrres_cor);
		//arrres_cor.show();
		SPDirection d1, d2, d3;
		Direction_Decompose(child, d1, d2, d3);
		//Point2D point_x=pn->getPoint(d1);
		//Point2D point_y=pn->getPoint(d2);
		interpolate_1order_on_face(pn, d1, arridx, arrres_x);
		//arrres_x.show();
		interpolate_1order_on_face(pn, d2, arridx, arrres_y);
		//arrres_y.show();
		for (int i = 0; i < arridx.size(); ++i) {
			arrres[i] = bilinear_interpolation(point,
					pn->cell->getCenterPoint(), point_c,
					pn->data->aCenterData[arridx[i]], arrres_x[i], arrres_y[i],
					arrres_cor[i]);
		}
		return 1;
	} else {
		std::cerr << " >! Interploate point is not in cell\n";
		return -1;
	}
}

// this function according to
// Ji H, Lien F S, Yee E.
// An efficient second order accurate cut cell method
// for solving the variable coefficient Poisson equation
// with jump conditions on irregular domains[J].
// International journal for numerical methods in fluids, 2006, 52(7): 723-748.
int _interpolate_node_LS(  // 2D QuadTree Node Least Square
		pQTNode pn,            //pnode
		const Point2D& point,  //point
		arrayList_st& arridx,  //data index
		arrayList& arrres   //data res
		) {
	if (pn->cell->getCenterPoint() == point) { //point equal special case
		for (int i = 0; i < arridx.size(); ++i) {
			arrres[i] = pn->data->aCenterData[arridx[i]];
		}
		return 1;
	}
	if (pn->cell->isInOnCell(point)) {
		// gradient of P
		arrayList arrgpx(arridx.size());
		arrayList arrgpy(arridx.size());
		// get neigbor node
		// method1 stencil  get 4 neighbor of P
		//    n
		//  n p n
		//    n
		ListT<pQTNode> lneip;

		// distance
		return 1;
	} else {
		std::cerr << " >! Interploate point is not in cell\n";
		return -1;
	}
}

int interpolate( // 2D QuadTree
		pQuadTree pqt, //pQuadTree
		const Point2D& point, //point
		arrayList_st& arridx, //data index
		arrayList& arrres //data res
		) {
	pQTNode pn = pqt->Find(point);
	if (pn != NULL_PTR) {
		return _interpolate_node(pn, point, arridx, arrres);
	} else {
		std::cerr << " >! Interploate point is not in tree\n";
		return -1;
	}
}

int interpolate_1order(pQuadTree pqt, const Point2D& point,
		arrayList_st& arridx, arrayList& arrres) {
	pQTNode pn = pqt->Find(point);
	if (pn != NULL_PTR) {
		return _interpolate_node_1order(pn, point, arridx, arrres);
	} else {
		std::cerr << " >! Interploate point is not in tree\n";
		return -1;
	}
}

int interpolate( // 2D Forest
		Forest2D& forest, //pQuadTree
		const Point2D& point, //point
		arrayList_st& arridx, //data index
		arrayList& arrres //data res
		) {
	pQuadTree pt = forest.getpTree(point);
	if (pt != NULL_PTR) {
		return interpolate( // 2D QuadTree
				pt, //pQuadTree
				point, //point
				arridx, //data index
				arrres //data res
				);
	} else {
		return -1; //not found
	}
}

int interpolate_1order(Forest2D& forest, const Point2D& point,
		arrayList_st& arridx, arrayList& arrres) {
	pQuadTree pt = forest.getpTree(point);
	if (pt != NULL_PTR) {
		return interpolate_1order( // 2D QuadTree
				pt, point, arridx, arrres);
	} else {
		return -1; //not found
	}
}
void interpolate( // 2D Forest
		Forest2D& forest, //pQuadTree
		const ListT<Point2D>& lp, //point
		arrayList_st& arridx, //data index
		ListT<Pair<int, arrayList> >& arrres //data res
		) {
	for (ListT<Point2D>::const_iterator iter = lp.begin(); iter != lp.end();
			++iter) {
		int tmpint;
		arrayList tmpres(arridx.size());
		Pair<int, arrayList> p(tmpint, tmpres);
		tmpint = interpolate( // 2D Forest
				forest, //pQuadTree
				(*iter), //point
				arridx, //data index
				p.second //data res
				);
		arrres.push_back(p);
	}
}

void interpolate( // 2D Forest
		Forest2D& forest, //Forest
		const MatrixT<Point2D>& lp, //point
		arrayList_st& arridx, //data index
		MatrixT<Pair<int, arrayList> >& arrres //data res if int=-1 result is wrong
		) {
	assert(lp.size() > 0);
	for (MatrixT<Point2D>::size_type i = 0; i < lp.iLen(); ++i) {
		for (MatrixT<Point2D>::size_type j = 0; j < lp.jLen(); ++j) {
			int tmpint;
			arrayList tmpres(arridx.size());
			Pair<int, arrayList> p(tmpint, tmpres);
			tmpint = interpolate( // 2D Forest
					forest, //pQuadTree
					lp[i][j], //point
					arridx, //data index
					p.second //data res
					);
			arrres[i][j] = p;
		}
	}
}

//===============================================
template<typename NODE>
void visit_get_average_value_on_level(NODE* node, utPointer pt) {
	arrayListT<utPointer> &arr = (*CAST(arrayListT<utPointer>*, pt));
	int& level = (*CAST(int*, arr[0]));
	if (condition_at_level(node, level)) {
		arrayList_st& arridx = (*CAST(arrayList_st*, arr[1]));
		typedef ListT<Pair<typename NODE::Cell_type::Point, arrayList> > listPPA;
		listPPA& listres = (*CAST(listPPA*, arr[2]));
		Pair<typename NODE::Cell_type::Point, arrayList> pair(
				node->cell->getCenterPoint(), arrayList(arridx.size()));
		T_cal_averange_value_from_leaf(node, arridx, pair.second);
		listres.push_back(pair);
	}
}

void get_average_value_on_level( // 2D QuadTree
		pQuadTree tree, //Forest
		int level, //level
		arrayList_st& arridx, //data index
		ListT<Pair<Point2D, arrayList> >& listres //data res if int=-1 result is wrong
		) {
	arrayListT<utPointer> arr(3);
	arr[0] = &level;
	arr[1] = &arridx;
	arr[2] = &listres;
	tree->Traversal(visit_get_average_value_on_level, &arr);
}

void get_average_value_on_level( // 2D Forest
		Forest2D& forest, //Forest
		int level, //level
		arrayList_st& arridx, //data index
		ListT<Pair<Point2D, arrayList> >& listres //data res if int=-1 result is wrong
		) {
	listres.clear();
	arrayListT<utPointer> arr(3);
	arr[0] = &level;
	arr[1] = &arridx;
	arr[2] = &listres;
	for (LarusDef::size_type i = 0; i < forest.iLen(); ++i) {
		for (LarusDef::size_type j = 0; j < forest.jLen(); ++j) {
			if (forest.get_attribution(i, j) == ATT_ENABLE) {
				forest.getpTree(i, j)->Traversal(
						visit_get_average_value_on_level, &arr);
			}
		}
	}
}

void get_max_value(   //
		Forest2D& forest,   //level
		arrayList_st& arridx,   //data index
		arrayList& arrres   //data res
		) {
	Forest2D::iterator iter = forest.begin();
	for (LarusDef::size_type i = 0; i < arridx.size(); ++i) {
		arrres[i] = iter->data->aCenterData[arridx[i]];
	}
	for (; iter != forest.end(); iter++) {
		for (LarusDef::size_type i = 0; i < arridx.size(); ++i) {
			if (iter->data->aCenterData[arridx[i]] > arrres[i]) {
				arrres[i] = iter->data->aCenterData[arridx[i]];
			}
		}
	}
}

void get_min_value(   //
		Forest2D& forest,   //level
		arrayList_st& arridx,   //data index
		arrayList& arrres   //data res
		) {
	Forest2D::iterator iter = forest.begin();
	for (LarusDef::size_type i = 0; i < arridx.size(); ++i) {
		arrres[i] = iter->data->aCenterData[arridx[i]];
	}
	for (; iter != forest.end(); iter++) {
		for (LarusDef::size_type i = 0; i < arridx.size(); ++i) {
			if (iter->data->aCenterData[arridx[i]] < arrres[i]) {
				arrres[i] = iter->data->aCenterData[arridx[i]];
			}
		}
	}
}

Float get_max_value(   //
		Forest2D& forest,   //level
		LarusDef::size_type idx   //data index
		) {
	arrayList_st arridx(1);  //data index
	arrayList arrres(1);  //data res
	arridx[0] = idx;
	get_max_value(forest, arridx, arrres);
	return arrres[0];
}

Float get_min_value(   //
		Forest2D& forest,   //level
		LarusDef::size_type idx   //data index
		) {
	arrayList_st arridx(1);  //data index
	arrayList arrres(1);  //data res
	arridx[0] = idx;
	get_min_value(forest, arridx, arrres);
	return arrres[0];
}

//draw===========================================
//===============================================
template<typename TREE>
void _draw_vtk_cell_leaf(FILE*& data, TREE* pt, int numleaf) {
	fprintf(data, "POINTS %d float\n", numleaf * 8);

	utPointer utp = data;
	pt->Traversal(visit_draw_to_vtk, utp);

	fprintf(data, "\n");
	fprintf(data, "CELLS %d %d \n", numleaf, numleaf * 9);

	for (int i = 0; i < numleaf * 8; ++i) {
		if ((i % 8) == 0) {
			fprintf(data, "%d ", 8);
		}
		fprintf(data, "%d ", i);
		if ((i % 8) == 7) {
			fprintf(data, "\n");
		}
	}
	fprintf(data, "\n");
	fprintf(data, "CELL_TYPES %d\n", numleaf);
	for (int i = 0; i < numleaf; ++i) {
		fprintf(data, "%d \n", 11);
	}
	fprintf(data, "\n");
}

template<typename NODE>
void visit_draw_to_vtu_centerdata(NODE* node, utPointer pt) {
	arrayListT<utPointer> arr = (*CAST(arrayListT<utPointer>*, pt));
	std::ofstream& fs = (*CAST(std::ofstream*, arr[0]));
	LarusDef::size_type& idx = (*CAST(LarusDef::size_type*, arr[1]));
	if (condition_is_leaf(node)) {
		fs << node->data->aCenterData[idx] << " ";
	}
}

template<typename NODE>
void visit_draw_to_vtu(NODE* node, utPointer pt) {
	std::ofstream& fs = (*CAST(std::ofstream*, pt));
	if (condition_is_leaf(node)) {
		node->cell->output_vertex_in_vtk_order(fs);
	}
}

template<typename TREE>
void T_draw_vtu_point_leaf_scalars( //
		std::string filename, //filename
		TREE* ptree, //tree
		std::string dataname, //data name
		LarusDef::size_type idx //data index
		) {
	LarusDef::size_type dim = TREE::Node::DIM;
	LarusDef::size_type vertexes = TREE::Node::Cell_type::NUM_VERTEXES;
	std::ofstream fs;
	open_file(fs, filename, 1);
	int numleaf = ptree->count_leaf();
	vtu_unstructured_grid_file_head(fs);
	fs << "<Piece NumberOfPoints=\"" << numleaf * vertexes
			<< "\" NumberOfCells=\"" << numleaf << "\">";
	fs << "<Points>";
	fs
			<< "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";
	ptree->Traversal(visit_draw_to_vtu, &fs);
	fs << "</DataArray>";
	fs << "</Points>";
	fs << "<Cells>";
	fs << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
	for (int i = 0; i < numleaf * vertexes; ++i) {
		fs << i << " ";
	}
	fs << "</DataArray>";
	fs << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">";
	for (int i = 1; i < numleaf + 1; ++i) {
		fs << i * vertexes << " ";
	}
	fs << "</DataArray>";
	fs << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">";
	for (int i = 0; i < numleaf; ++i) {
		if (dim == 3) {
			fs << 11 << " ";
		} else if (dim == 2) {
			fs << 8 << " ";
		}
	}
	fs << "</DataArray>";
	fs << "</Cells>";
	fs << "<CellData Scalars=\"scalars\">";
	fs << "<DataArray type=\"Float32\" Name=\"" << dataname
			<< "\" format=\"ascii\">";
	arrayListT<utPointer> arr(2);
	arr[0] = &fs;
	arr[1] = &idx;
	ptree->Traversal(visit_draw_to_vtu_centerdata, &arr);
	fs << "</DataArray>";
	fs << "</CellData>";
	fs << "</Piece>";
	vtu_unstructured_grid_file_end(fs);
	fs.close();
}

void draw_vtu_point_leaf_scalars( //
		std::string filename, //filename
		pOCTree oct, //tree
		std::string dataname, //data name
		LarusDef::size_type idx //data index
		) {
	T_draw_vtu_point_leaf_scalars(filename,  //filename
			oct,  //tree
			dataname,  //data name
			idx);
}

void draw_vtu_point_leaf_scalars( //
		std::string filename, //filename
		pQuadTree pqt, //tree
		std::string dataname, //data name
		LarusDef::size_type idx //data index
		) {
	T_draw_vtu_point_leaf_scalars(filename,  //filename
			pqt,  //tree
			dataname,  //data name
			idx);
}

template<typename TREE>
void T_draw_vtu_point_leaf_scalars( //
		std::string filename, //filename
		TREE* ptree, //tree
		arrayListT<Pair<std::string, LarusDef::size_type> > arrd // data pair
		) {
	LarusDef::size_type dim = TREE::Node::DIM;
	LarusDef::size_type vertexes = TREE::Node::Cell_type::NUM_VERTEXES;
	std::ofstream fs;
	open_file(fs, filename, 1);
	int numleaf = ptree->count_leaf();
	vtu_unstructured_grid_file_head(fs);
	fs << "<Piece NumberOfPoints=\"" << numleaf * vertexes
			<< "\" NumberOfCells=\"" << numleaf << "\">";
	fs << "<Points>";
	fs
			<< "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";
	ptree->Traversal(visit_draw_to_vtu, &fs);
	fs << "</DataArray>";
	fs << "</Points>";
	fs << "<Cells>";
	fs << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
	for (int i = 0; i < numleaf * vertexes; ++i) {
		fs << i << " ";
	}
	fs << "</DataArray>";
	fs << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">";
	for (int i = 1; i < numleaf + 1; ++i) {
		fs << i * vertexes << " ";
	}
	fs << "</DataArray>";
	fs << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">";
	for (int i = 0; i < numleaf; ++i) {
		if (dim == 3) {
			fs << 11 << " ";
		} else if (dim == 2) {
			fs << 8 << " ";
		}
	}
	fs << "</DataArray>";
	fs << "</Cells>";
	fs << "<CellData Scalars=\"scalars\">";
	for (int i = 0; i < arrd.size(); ++i) {
		fs << "<DataArray type=\"Float32\" Name=\"" << arrd[i].first
				<< "\" format=\"ascii\">";
		arrayListT<utPointer> arr(2);
		arr[0] = &fs;
		arr[1] = &arrd[i].second;
		ptree->Traversal(visit_draw_to_vtu_centerdata, &arr);
		fs << "</DataArray>";
	}
	fs << "</CellData>";
	fs << "</Piece>";
	vtu_unstructured_grid_file_end(fs);
	fs.close();
}

void draw_vtu_point_leaf_scalars( //
		std::string filename, //filename
		pOCTree oct, //tree
		arrayListT<Pair<std::string, LarusDef::size_type> > arrd // data pair
		) {
	T_draw_vtu_point_leaf_scalars( //
			filename, //filename
			oct, //tree
			arrd // data pair
			);
}

void draw_vtu_point_leaf_scalars( //
		std::string filename, //filename
		pQuadTree qt, //tree
		arrayListT<Pair<std::string, LarusDef::size_type> > arrd // data pair
		) {
	T_draw_vtu_point_leaf_scalars( //
			filename, //filename
			qt, //tree
			arrd // data pair
			);
}

template<typename NODE>
void T_draw_vtu_listnode_scalars( //
		std::string filename, //filename
		ListT<NODE*>& list, //tree
		arrayListT<Pair<std::string, LarusDef::size_type> > arrd // data pair
		) {
	LarusDef::size_type dim = NODE::DIM;
	LarusDef::size_type vertexes = NODE::Cell_type::NUM_VERTEXES;
	std::ofstream fs;
	open_file(fs, filename, 1);
	int numnode = list.size();
	vtu_unstructured_grid_file_head(fs);
	fs << "<Piece NumberOfPoints=\"" << numnode * vertexes
			<< "\" NumberOfCells=\"" << numnode << "\">";
	fs << "<Points>";
	fs
			<< "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";
	for (typename ListT<NODE*>::iterator iter = list.begin();
			iter != list.end(); iter++) {
		(*iter)->cell->output_vertex_in_vtk_order(fs);
	}
	fs << "</DataArray>";
	fs << "</Points>";
	fs << "<Cells>";
	fs << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
	for (int i = 0; i < numnode * vertexes; ++i) {
		fs << i << " ";
	}
	fs << "</DataArray>";
	fs << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">";
	for (int i = 1; i < numnode + 1; ++i) {
		fs << i * vertexes << " ";
	}
	fs << "</DataArray>";
	fs << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">";
	for (int i = 0; i < numnode; ++i) {
		if (dim == 3) {
			fs << 11 << " ";
		} else if (dim == 2) {
			fs << 8 << " ";
		}
	}
	fs << "</DataArray>";
	fs << "</Cells>";
	fs << "<CellData Scalars=\"scalars\">";
	for (int i = 0; i < arrd.size(); ++i) {
		fs << "<DataArray type=\"Float32\" Name=\"" << arrd[i].first
				<< "\" format=\"ascii\">";
		for (typename ListT<NODE*>::iterator iter = list.begin();
				iter != list.end(); iter++) {
			fs << (*iter)->data->aCenterData[arrd[i].second] << " ";
		}
		fs << "</DataArray>";
	}
	fs << "</CellData>";
	fs << "</Piece>";
	vtu_unstructured_grid_file_end(fs);
	fs.close();
}

void draw_vtu_point_leaf_scalars( // 2D Forest
		std::string filename, //filename
		Forest2D& forest, //forest
		std::string dataname, //data name
		LarusDef::size_type idx //data index
		) {
	ListT<pQTNode> list;
	ListT<Dimension_2D::pNode> &reflist = list;  //make eclipse happy
	for (LarusDef::size_type i = 0; i < forest.iLen(); ++i) {
		for (LarusDef::size_type j = 0; j < forest.jLen(); ++j) {
			if (forest.get_attribution(i, j) == ATT_ENABLE)
				toList_leaf(forest.getpTree(i, j), reflist);
		}
	}
	arrayListT<Pair<std::string, LarusDef::size_type> > arrd(1);
	arrd[0] = Pair<std::string, LarusDef::size_type>(dataname, idx);
	T_draw_vtu_listnode_scalars(filename, list, arrd);
}

void draw_vtu_point_leaf_scalars( // 3D Forest
		std::string filename, //filename
		Forest2D& forest, //forest
		arrayListT<Pair<std::string, LarusDef::size_type> > arrd // data pair
		) {
	ListT<pQTNode> list;
	ListT<Dimension_2D::pNode> &reflist = list;  //make eclipse happy
	toListpNode(forest, reflist);
	T_draw_vtu_listnode_scalars(filename, list, arrd);
}

void _draw_gnuplot_as_vector( // 2D QuadTree
		FILE *data, //filename
		ListT<pQTNode>& listpnode, //tree
		LarusDef::size_type idx_x, //idx x
		LarusDef::size_type idx_y //icx y
		) {
	for (ListT<pQTNode>::iterator iter = listpnode.begin();
			iter != listpnode.end(); iter++) {
		if (!(*iter)->hasChild()) {  //leaf node
			Float u = (*iter)->data->aCenterData[idx_x];
			Float v = (*iter)->data->aCenterData[idx_y];
			fprintf(data, "%lf %lf ", (*iter)->cell->getCenterPoint().x,
					(*iter)->cell->getCenterPoint().y);
			fprintf(data, "%lf %lf \n", u, v);
		} else {   //non-leaf node
			arrayList_st arridx(2);
			arridx[0] = idx_x;
			arridx[1] = idx_y;
			arrayList arrres(2);
			T_cal_averange_value_from_leaf((*iter), arridx, arrres);
			fprintf(data, "%lf %lf ", (*iter)->cell->getCenterPoint().x,
					(*iter)->cell->getCenterPoint().y);
			fprintf(data, "%lf %lf \n", arrres[0], arrres[1]);
		}
	}
}

void draw_gnuplot_as_vector(     // 2D QuadTree
		std::string filename,     //filename
		ListT<pQTNode>& listpnode,     //tree
		int mode,     //mode
		LarusDef::size_type idx_x,     //idx x
		LarusDef::size_type idx_y     //icx y
		) {
	FILE *data = open_file(filename, mode);
	_draw_gnuplot_as_vector(data, listpnode, idx_x, idx_y);
	fclose(data);
}

void draw_gnuplot_as_vector( // 2D QuadTree
		std::string filename, //filename
		pQuadTree qt, //tree
		int mode, LarusDef::size_type idx_x, //idx x
		LarusDef::size_type idx_y //icx y
		) {
	FILE *data = open_file(filename, mode);
	ListT<pQTNode> listpnode;
	toList_leaf(qt, listpnode);
	_draw_gnuplot_as_vector(data, listpnode, idx_x, idx_y);
	fclose(data);
}

void draw_gnuplot_as_vector( // 2D QuadTree
		std::string filename, //filename
		Forest2D& forest, //tree
		LarusDef::size_type idx_x, //idx x
		LarusDef::size_type idx_y //icx y
		) {
	FILE *data = open_file(filename, 1);
	ListT<pQTNode> listpnode;
	ListT<Dimension_2D::pNode> &reflist = listpnode;  //make eclipse happy
	toListpNode(forest, reflist);
	_draw_gnuplot_as_vector(data, listpnode, idx_x, idx_y);
	fclose(data);
}

void _draw_gnuplot_as_contour( // 2D QuadTree
		FILE *data, //data
		ListT<pQTNode>& listpnode, //list node
		LarusDef::size_type idx //idx x
		) {
	for (ListT<pQTNode>::iterator iter = listpnode.begin();
			iter != listpnode.end(); iter++) {
		if ((*iter)->data != NULL_PTR) {
			fprintf(data, "%f %f ", (*iter)->cell->getCenterPoint().x,
					(*iter)->cell->getCenterPoint().y);
			fprintf(data, "%f %f ", (*iter)->cell->getMM().x,
					(*iter)->cell->getPP().x);
			fprintf(data, "%f %f ", (*iter)->cell->getMM().y,
					(*iter)->cell->getPP().y);
			fprintf(data, "%f\n", (*iter)->data->aCenterData[idx]);
		}
	}
}

//marching square ===============================
inline Float _ms_interpolate_1(Float d, Float vo, Float vx, Float iso) {
	return (iso - vo) / (vx - vo) * d;
}

int _ms_which_case( //
		const arrayList& av, //array value  size=4
		Float val, std::bitset<4>& ab) //array level  size=n
		{
	for (int i = 0; i < 4; i++) {
		ab[i] = av[i] > val ? 1 : 0;
		//value is above the isovalue = 1
		//value is blew  the isovalue = 0
	}
	return int(ab.to_ulong());
}

int _ms_a_level(		//
		const arrayListT<Point2D> aloc,		//the points
		const arrayList& av,		//array value  size=4
		const Float& cd,		//center data
		const Float& iso,		//iso
		Segment2D& s1,		//level case
		Segment2D& s2)		//level case
		{
	std::bitset<4> ab(0);
	int wca = _ms_which_case(av, iso, ab);
	if (wca == 0 || wca == 15) {     //no segment
		return 0;
	} else {                         //has segments
		arrayListT<Point2D> arrp(4);
		for (int i = 0; i < 4; ++i) {
			int ip = (i + 1) % 4;
			if (ab[i] != ab[ip]) {
				int d = (i == 0 || i == 2) ? 0 : 1;       //direction choose
				arrp[i][d] = _ms_interpolate_1(       //
						(aloc[ip][d] - aloc[i][d]), av[i], av[ip], iso)
						+ aloc[i][d];
				arrp[i][(d + 1) % 2] = aloc[i][(d + 1) % 2];
			}
		}
		if (wca == 5 || wca == 10) { //two segment case
			const int CSEG2[2][4] = { { 0, 3, 2, 1 },         //center point  =0
					{ 0, 1, 2, 3 }   //center point  =1
			};
			int idx1 = cd > iso ? 1 : 0;
			s1[0] = arrp[CSEG2[idx1][wca == 5 ? 0 : 1]];
			s1[1] = arrp[CSEG2[idx1][wca == 5 ? 1 : 0]];
			s2[0] = arrp[CSEG2[idx1][wca == 5 ? 2 : 3]];
			s2[1] = arrp[CSEG2[idx1][wca == 5 ? 3 : 2]];
			return 2;
		} else {                     //one segment case
			const int CSEG1[8][2] = { { 0, 0 },  //no
					{ 0, 3 }, { 1, 0 }, { 1, 3 }, { 2, 1 }, { 0, 0 },  //no
					{ 2, 0 }, { 2, 3 }, };
			int idx1 = (wca <= 7) ? wca : int((~ab).to_ulong());
			s1[0] = arrp[CSEG1[idx1][wca <= 7 ? 0 : 1]];
			s1[1] = arrp[CSEG1[idx1][wca <= 7 ? 1 : 0]];
			return 1;
		}
	}
}

void _ms_many_levels( //
		const arrayListT<Point2D>& aloc, //the points
		const arrayList& av, //array value  size=4
		const Float& cd, //center data
		const arrayList& aiso, //array iso
		ListT<Segment2D>& ls, //Segment2D
		ListT<Pair<int, int> >& lidx) //list index segment
		{
	Segment2D s1, s2;
	int nums = 0;
	for (int i = 0; i < aiso.size(); ++i) {
		nums = _ms_a_level(aloc, av, cd, aiso[i], s1, s2);
		if (nums == 1) {
			ls.push_back(s1);
			lidx.push_back(Pair<int, int>(nums, i));
		} else if (nums == 2) {
			ls.push_back(s1);
			ls.push_back(s2);
			lidx.push_back(Pair<int, int>(nums, i));
		}
		nums = 0;
	}
}

void _output_ms_segments( //
		FILE *data, //
		const ListT<Segment2D>& ls, //
		const ListT<Pair<int, int> >& lidx, //
		const arrayList& alevel) {
	_IF_TRUE_RETRUN(ls.empty());
	ListT<Segment2D>::const_iterator seg_iter = ls.begin();
	int nums = 0;
	for (ListT<Pair<int, int> >::const_iterator iter = lidx.begin();
			iter != lidx.end();) {
		nums++;
		fprintf(data, "%f %f ", seg_iter->PSX(), seg_iter->PSY());
		fprintf(data, "%d %f \n", iter->second, alevel[iter->second]);
		fprintf(data, "%f %f ", seg_iter->PEX(), seg_iter->PEY());
		fprintf(data, "%d %f \n\n", iter->second, alevel[iter->second]);
		if (iter->first == nums) {
			++iter;
			nums = 0;
		}
		++seg_iter;
	}
}

void _draw_gnuplot_as_contour_line( // 2D QuadTree
		FILE *data, //data
		pQTNode pnode, //list node
		LarusDef::size_type idx, //idx x
		const arrayList& alevel) {
	if (pnode->data != NULL_PTR) {
		arrayList aval(4), ares(1);
		arrayList_st aidx(1);
		aidx[0] = idx;
		const SPDirection index_dir[4] = { SPD_MM, SPD_PM, SPD_PP, SPD_MP };
		//  0       1       2       3
		arrayListT<Point2D> aloc(4);
		for (int i = 0; i < 4; ++i) {
			//get value on vertex =====================
			interpolate_1order_on_vertex_averange_neighbor(pnode, index_dir[i],
					aidx, ares);
			aval[i] = ares[0];
			aloc[i] = pnode->getPoint(index_dir[i]); //get point location
		}
		ListT<Segment2D> ls;                   //list Segment2D
		ListT<Pair<int, int> > lidx;           //list index segment
		_ms_many_levels(aloc, aval, pnode->data->aCenterData[idx], alevel, ls,
				lidx);
		_output_ms_segments(data, ls, lidx, alevel);
	}
}
//end marching square ---------------------------

void draw_gnuplot_as_contour_line(           // 2D QuadTree
		std::string filename,           //filename
		Forest2D& forest,           //tree
		int mode,           //mode
		LarusDef::size_type idx,           //idx x
		arrayList alevel) {
	FILE *data = open_file(filename, mode);
	for (Forest2D::iterator iter = forest.begin(); iter != forest.end();
			iter++) {
		_draw_gnuplot_as_contour_line(data, &(*iter), idx, alevel);
	}
	fclose(data);
}

void draw_gnuplot_as_contour_line(  // 2D QuadTree
		std::string filename,  //filename
		pQuadTree qt,  //tree
		int mode,  //mode
		LarusDef::size_type idx,  //idx x
		arrayList alevel  //array level
		) {
	FILE *data = open_file(filename, mode);
	for (QuadTree::iterator iter = qt->begin(); iter != qt->end(); iter++) {
		_draw_gnuplot_as_contour_line(data, &(*iter), idx, alevel);
	}
	fclose(data);
}

void draw_gnuplot_as_contour_line(      // 2D QuadTree
		std::string filename,      //filename
		ListT<pQTNode> lnode,      //tree
		int mode,      //mode
		LarusDef::size_type idx,      //idx x
		arrayList alevel      //array level
		) {
	FILE *data = open_file(filename, mode);
	for (ListT<pQTNode>::iterator iter = lnode.begin(); iter != lnode.end();
			iter++) {
		_draw_gnuplot_as_contour_line(data, (*iter), idx, alevel);
	}
	fclose(data);
}

void _draw_gnuplot_as_contour( // 2D QuadTree
		FILE *data, //data
		pQTNode pnode, //list node
		LarusDef::size_type idx //idx x
		) {
	if (pnode->data != NULL_PTR) {
		fprintf(data, "%f %f ", pnode->cell->getCenterPoint().x,
				pnode->cell->getCenterPoint().y);
		fprintf(data, "%f %f ", pnode->cell->getMM().x, pnode->cell->getPP().x);
		fprintf(data, "%f %f ", pnode->cell->getMM().y, pnode->cell->getPP().y);
		fprintf(data, "%f\n", pnode->data->aCenterData[idx]);
	}
}

void draw_gnuplot_as_contour( // 2D QuadTree
		std::string filename, //filename
		pQuadTree qt, //tree
		int mode, //mode
		LarusDef::size_type idx) //idx x
		{
	FILE *data = open_file(filename, mode);
	for (QuadTree::iterator iter = qt->begin(); iter != qt->end(); iter++) {
		_draw_gnuplot_as_contour(data, &(*iter), idx);
	}
	fclose(data);
}

void draw_gnuplot_as_contour(  // 2D QuadTree
		std::string filename,  //filename
		Forest2D& forest,  //tree
		int mode,  //mode
		LarusDef::size_type idx  //idx x
		) {
	FILE *data = open_file(filename, mode);
	for (Forest2D::iterator iter = forest.begin(); iter != forest.end();
			iter++) {
		_draw_gnuplot_as_contour(data, &(*iter), idx);
	}
	fclose(data);
}
void show_gnuplot_as_contour( // 2D QuadTree
		Forest2D& forest,       //tree
		LarusDef::size_type idx //idx x
		) {
	typedef Float vt;
	ListT<vt> lxc, lyc, lxm, lxp, lym, lyp, lval;
	for (Forest2D::const_iterator iter = forest.begin(); iter != forest.end();
			iter++) {
		const Forest2D::Node* pnode = iter.get_pointer();
		lxc.push_back(pnode->cell->getCenterPoint().x);
		lyc.push_back(pnode->cell->getCenterPoint().y);
		lxm.push_back(pnode->cell->getMM().x);
		lxp.push_back(pnode->cell->getPP().x);
		lym.push_back(pnode->cell->getMM().y);
		lyp.push_back(pnode->cell->getPP().y);
		lval.push_back(pnode->data->aCenterData[idx]);
	}
	ListT<vt>::const_iterator iter = lval.begin();
	vt max = (*iter);
	vt min = (*iter);
	for (++iter; iter != lval.end(); ++iter) {
		if ((*iter) > max) {
			max = (*iter);
		}
		if ((*iter) < min) {
			min = (*iter);
		}
	}
	Gnuplot gp("boxes");
	std::string cmdstr = "with boxxy title \"\" fs solid palette";
	gp.set_palette_blue_red();
	if (max == min) {
		gp.set_cbrange(ABS(max), -ABS(min));
	} else {
		gp.set_cbrange(min, max);
	}

	std::ostringstream ss;
	gp.set_xlabel(ss.str());
	gp.plot_7(lxc, lyc, lxm, lxp, lym, lyp, lval, cmdstr);
}
}  //end of namespace
