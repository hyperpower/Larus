/************************
 //  \file   Octree.h
 //  \brief
 // 
 //  \author zhou
 //  \date   3 sept. 2014 
 ***********************/
#ifndef OCTREE_H_
#define OCTREE_H_

#include "../TypeDef.h"
#include "../Algebra/Arithmetic.h"
#include "../Utility/Array.h"
#include "../Utility/List.h"
#include "../IO/IO.h"
#include "../IO/IO_vtk.h"
#include "Cell.h"
#include "SPTreeNode.h"

#include <string>

namespace Larus {

template<class NODE, LarusDef::size_type Dim>
class SPTree {
public:
	//def========================================
	typedef SPTree<NODE, Dim> Tree;
	typedef SPTree<NODE, Dim>* pTree;
	typedef NODE* pNode;
	typedef NODE Node;
	typedef const NODE* const_pNode;
	typedef const NODE const_Node;
	typedef void (*pFun_SPTree)(pNode, utPointer);
	typedef void (*pFun_SPTree_const)(const_pNode, const_utPointer);
	typedef void (*pFun_SPTree_condition)(arrayList&, pNode, utPointer);
	typedef int (*pFun_SPTree_condition_re)(arrayList&, pNode, utPointer);
	typedef LarusDef::size_type size_t;

	typedef _SPNode_iterator<typename NODE::Cell_type, typename NODE::Data_type,
			Dim, NODE&, NODE*> iterator;
	typedef _SPNode_iterator<typename NODE::Cell_type, typename NODE::Data_type,
			Dim, const NODE&, const NODE*> const_iterator;
protected:
	pNode _root;
	size_t _maxlevel;
	pTree neighbor[Dim + Dim];

	void _delete(pNode Current);

	int _creatchildren_partial(pNode, const arrayList&);  //test

	//level_0====================================
	void _traversal_level(size_t, pNode, pFun_SPTree, utPointer);
	void _traversal(pNode, pFun_SPTree, utPointer);
	void _traversal_b(pNode, pFun_SPTree, utPointer);
	void _conditional_traversal(pNode, pFun_SPTree_condition, pFun_SPTree,
			utPointer);
	void _refine_leaf(pNode, pFun_SPTree_condition_re, utPointer, pFun_SPTree,
			utPointer); //test
	//end level 0================================

	pNode _father(pNode);
	pNode _child(pNode, int);
	const_pNode _father(pNode) const;
	const_pNode _child(pNode, int) const;
	//pNode _find(pNode, const Point2D &);

public:
	SPTree(const typename NODE::Cell_type& c, LarusDef::size_type = 3);
	~SPTree();

	inline size_t getMaxLevel() const {
		return _maxlevel;
	}

	//Traversal==================================
	void Traversal_partial(pNode, pFun_SPTree, utPointer);
	void Traversal_level_RtL(pFun_SPTree, utPointer);  //Root to Leaf
	void Traversal_level_LtR(pFun_SPTree, utPointer);  //Leaf to Root
	void Traversal(pFun_SPTree, utPointer);
	void Traversal_b(pFun_SPTree, utPointer);
	void Conditional_Traversal(pFun_SPTree_condition, pFun_SPTree, utPointer);
	void Traversal_leaf_to_list(ListT<pNode>&);
	//Refine=======================================
	void Refine(pFun_SPTree_condition_re, utPointer, pFun_SPTree, utPointer);
	//Creat=======================================
	int Creatchildren(pNode Current, SPNodeIdx idx_child);
	int Creatchildren_full(pNode);
	int Creatchildren_partial(pNode, bool, bool, bool, bool, bool = false,
			bool = false, bool = false, bool = false);
	void CreatFullTree();
	void CreatFullTree_to_Level(size_t);
	//Connect ====================================
	void ConnectNodes();
	//============================================
	bool isEmpty() const {
		return _root == NULL ? true : false;
	}
	bool isNewTree() const {
		return _root->hasChild() ? false : true;
	}
	bool isInTree(const typename NODE::Cell_type::Point&) const;
	bool isInOnTree(const typename NODE::Cell_type::Point&) const;
	//count=======================================
	int count_all();
	int count_level(int);
	int count_leaf();
	int count_leaf_at_level(int);

	void show_info();
	//Other ======================================
	const typename NODE::Cell_type* getpRootCell() const {
		return _root->cell;
	}
	pNode getpRootNode() {
		return _root;
	}
	const pNode getpRootNode() const {
		return _root;
	}
	pNode getpNode(const typename NODE::Cell_type::Point&);

	//Neighbor tree ==============================
	void setNeighborTree( //
			pTree xm, pTree xp, //x
			pTree ym, pTree yp, //y
			pTree = NULL_PTR, pTree = NULL_PTR); //z
	pTree getNeighborpTree(SPDirection);
	const pTree getNeighborpTree(SPDirection) const;
	//iterator ====================================
public:
	iterator begin();
	const_iterator begin() const;
	iterator end();
	const_iterator end() const;
protected:
	pTree _walk_path(SPDirection d1, SPDirection d2);
	pTree _walk_path(SPDirection d1, SPDirection d2, SPDirection d3);
	pTree _getpNeighborTree_path(SPDirection, SPDirection);
	pTree _getpNeighborTree_path(SPDirection, SPDirection, SPDirection);
public:
	pTree getpNeighborTree(SPDirection);

	//Neighbor Finding============================
	pNode getAdjNeighbor(pNode Current, SPDirection d);
	pNode getCorNeighbor(pNode Current, SPDirection d);
	pNode getNeighbor_Adj_Cor(pNode Current, SPDirection d);
	pNode getCorNeighbor_XYZ(pNode Current, SPDirection d);
	pNode getNeighbor(pNode Current, SPDirection d);

	//IO==========================================
protected:
	void _draw_to_vtk_leaf(FILE*& data, int numleaf);
public:
	void draw_to_vtk_leaf(std::string filename);
	void draw_to_vtk(FILE*&);
	void draw_gnuplot(std::string filename, int mode);  //for 2D
};

//function out of class===========================
template<class NODE>
bool is_height_1(NODE* pn) {
	if (pn != NULL_PTR) {
		if (!pn->hasChild()) {
			return false;
		}
		for (int i = 0; i < NODE::NUM_CELLS; i++) {
			if (pn->child[i] != NULL_PTR) {
				if (pn->child[i]->hasChild()) {
					return false;
				}
			}
		}
		return true;
	} else {
		return false;
	}
}
template<class NODE>
bool condition_is_root(NODE* pnode) {
	if (pnode->getLevel() == 0) {
		return true;
	} else {
		return false;
	}
}
template<class NODE>
bool condition_is_leaf(NODE* pnode) {
	if (!pnode->hasChild()) {
		return true;
	} else {
		return false;
	}
}
template<class NODE>
bool condition_at_level(NODE* pn, int level) {
	if (pn->getLevel() == level) {
		return true;
	} else {
		return false;
	}
}
template<class NODE>
bool condition_is_node(NODE* pnode) {
	if (pnode != NULL_PTR) {
		return true;
	} else {
		return false;
	}
}
template<class NODE>
bool condition_is_le_level(NODE* pnode, LarusDef::size_type level) {
	if (pnode->getLevel() <= level) {
		return true;
	} else {
		return false;
	}
}
template<class NODE>
bool condition_is_eq_level(NODE* pnode, LarusDef::size_type level) {
	if (pnode->getLevel() == level) {
		return true;
	} else {
		return false;
	}
}
//visit function ================================
template<class NODE>
void visit_empty(NODE* pn, utPointer p) {
}

template<class NODE>
void visit_leaf_to_list(NODE* pn, utPointer p) {
	if (condition_is_leaf(pn)) {
		ListT<NODE*>* listp = CAST(ListT<NODE*>*, p);
		listp->push_back(pn);
	}
}
template<class NODE>
void visit_conunt_level(NODE* pn, utPointer p) {
	arrayListT<size_t*>* arg = CAST(arrayListT<size_t*>*, p);
	size_t& level = *(*arg)[0];
	size_t& res = *(*arg)[1];
	if (condition_at_level(pn, level)) {
		res++;
	}
}
template<class NODE>
void visit_conunt_all(NODE* pn, utPointer p) {
	size_t& res = *(CAST(size_t*, p));
	res++;
}
template<class NODE>
void visit_conunt_leaf(NODE* pn, utPointer p) {
	size_t& res = *(CAST(size_t*, p));
	if (condition_is_leaf(pn)) {
		res++;
	}
}
template<class NODE>
void visit_conunt_leaf_at_level(NODE* pn, utPointer p) {
	arrayListT<size_t*>* arg = CAST(arrayListT<size_t*>*, p);
	size_t& level = *(*arg)[0];
	size_t& res = *(*arg)[1];
	if (condition_at_level(pn, level) && condition_is_leaf(pn)) {
		res++;
	}
}
template<class NODE>
void visit_CreatFullTree(NODE* pn, utPointer p) {
	typedef SPTree<NODE, NODE::DIM>* ptree_t;
	ptree_t pt = CAST(ptree_t, p);
	if (condition_is_node(pn)) {
		pt->Creatchildren_full(pn);
	}
}
template<class NODE>
void visit_CreatFullTree_to_Level(NODE* pn, utPointer p) {
	typedef SPTree<NODE, NODE::DIM>* ptree_t;
	arrayListT<utPointer>& arrp = (*CAST(arrayListT<utPointer>*, p));
	ptree_t pt = CAST(ptree_t, arrp[0]);
	int& level = (*CAST(int*, arrp[1]));
	if (condition_is_node(pn) && pn->getLevel()<level) {
		pt->Creatchildren_full(pn);
	}
}
template<class NODE>
void visit_to_list_leaf(NODE* pn, utPointer p) {
	if (condition_is_leaf(pn)) {
		ListT<NODE*>* list = CAST(ListT<NODE*>*, p);
		list->push_back(pn);
	}
}
template<class NODE>
void visit_to_list_leaf_slice(NODE* pn, utPointer p) {
	if (condition_is_leaf(pn)) {
		arrayListT<utPointer>& arr = (*CAST(arrayListT<utPointer>*, p));
		CSAxis& axis = (*CAST(CSAxis*, arr[0]));
		Float& loc = (*CAST(Float*, arr[1]));
		ListT<NODE*>* list = CAST(ListT<NODE*>*, arr[2]);
		if (isInRange(pn->cell->get(axis, eCPL_M), loc,
				pn->cell->get(axis, eCPL_P), Range_oc)) {
			list->push_back(pn);
		}
	}
}
template<class NODE>
void visit_to_list_leaf_slice_level(NODE* pn, utPointer p) {
	arrayListT<utPointer>& arr = (*CAST(arrayListT<utPointer>*, p));
	int& level = (*CAST(int*, arr[3]));
	if (condition_is_eq_level(pn, level)) {
		CSAxis& axis = (*CAST(CSAxis*, arr[0]));
		Float& loc = (*CAST(Float*, arr[1]));
		ListT<NODE*>* list = CAST(ListT<NODE*>*, arr[2]);
		if (isInRange(pn->cell->get(axis, eCPL_M), loc,
				pn->cell->get(axis, eCPL_P), Range_oc)) {
			list->push_back(pn);
		}
	}
}

template<class NODE>
void visit_find(NODE* qnode, utPointer p) {
	if (condition_is_leaf(qnode)) {
		CAST(arrayListT<utPointer>*,p)->at(1) = qnode;
	}
}
template<class NODE>
void visit_draw_to_vtk(NODE* node, utPointer p) {
	FILE *data = CAST(FILE*, p);
	if (condition_is_leaf(node)) {
		node->cell->output_vertex_in_vtk_order(data);
	}
}
template<class NODE>
void condition_point_at_which_child(arrayList& arrt, NODE* pnode, utPointer p) {
	if (pnode != NULL && p != NULL) {
		utPointer ap = CAST(arrayListT<utPointer>*,p)->at(0);
		arrt[pnode->whichChild((*CAST(typename NODE::Cell_type::Point*, ap)))] =
				1;
	}
}
//================================================
//Class member function ==========================
template<class NODE, LarusDef::size_type Dim>
void SPTree<NODE, Dim>::_delete(NODE* Current) {
	if (Current == NULL_PTR) {
		return;
	}
	if (Current->hasChild()) {
		for (int i = 0; i < NODE::NUM_CELLS; i++) {
			pNode c = Current->child[i];
			if (c != NULL_PTR) {
				_delete(c);
			}
		}
	}
	if (is_height_1(Current)) {
		//Current->sethasChild(false);
		for (int i = 0; i < NODE::NUM_CELLS; i++) {
			if (Current->child[i] != NULL_PTR) {
				delete Current->child[i];
				Current->child[i] = NULL_PTR;
			}
		}
	}
}

template<class NODE, LarusDef::size_type Dim>
SPTree<NODE, Dim>::SPTree(const typename NODE::Cell_type& c, int maxlevel) {
	_maxlevel = maxlevel;
	_root = new NODE(NULL, SPT_normal, 0, -1, c);
}
template<class NODE, LarusDef::size_type Dim>
SPTree<NODE, Dim>::~SPTree() {
	_delete(_root);
}
template<class NODE, LarusDef::size_type Dim>
void SPTree<NODE, Dim>::setNeighborTree(
		//
		typename SPTree<NODE, Dim>::pTree xm,
		typename SPTree<NODE, Dim>::pTree xp, //x
		typename SPTree<NODE, Dim>::pTree ym,
		typename SPTree<NODE, Dim>::pTree yp, //y
		typename SPTree<NODE, Dim>::pTree zm,
		typename SPTree<NODE, Dim>::pTree zp) { //z
	//       yp 1
	//      ______
	//     |      |
	//xm 0 |      | xp 2
	//     |______|
	//       ym 3
	neighbor[0] = xm;
	neighbor[1] = yp;
	neighbor[2] = xp;
	neighbor[3] = ym;
	if (Dim == 3) {
		neighbor[4] = zp;
		neighbor[5] = zm;
	}
}
template<class NODE, LarusDef::size_type Dim>
typename SPTree<NODE, Dim>::pTree SPTree<NODE, Dim>::getNeighborpTree(
		SPDirection dir) {
	ASSERT(isFaceDirection(dir, Dim));
	return neighbor[int(dir) - 4];
}
template<class NODE, LarusDef::size_type Dim>
const typename SPTree<NODE, Dim>::pTree SPTree<NODE, Dim>::getNeighborpTree(
		SPDirection dir) const {
	ASSERT(isFaceDirection(dir, Dim));
	return neighbor[int(dir) - 4];
}

template<class NODE, LarusDef::size_type Dim>
typename SPTree<NODE, Dim>::pTree SPTree<NODE, Dim>::_walk_path(SPDirection d1,
		SPDirection d2) {
	pTree pstep1 = this->getpNeighborTree(d1);
	if (pstep1 == NULL_PTR) {
		return NULL_PTR;
	} else {
		return pstep1->getpNeighborTree(d2);
	}
}
template<class NODE, LarusDef::size_type Dim>
typename SPTree<NODE, Dim>::pTree SPTree<NODE, Dim>::_walk_path(SPDirection d1,
		SPDirection d2, SPDirection d3) {
	pTree pstep = this->getpNeighborTree(d1);
	if (pstep == NULL_PTR) {
		return NULL_PTR;
	} else {
		pstep = pstep->getpNeighborTree(d2);
		if (pstep == NULL_PTR) {
			return NULL_PTR;
		} else {
			return pstep->getpNeighborTree(d3);
		}
	}
}

template<class NODE, LarusDef::size_type Dim>
typename SPTree<NODE, Dim>::pTree SPTree<NODE, Dim>::_getpNeighborTree_path(
		SPDirection d1, SPDirection d2) {
	pTree pt = this->_walk_path(d1, d2);
	if (pt != NULL_PTR) {
		return pt;
	} else {
		return this->_walk_path(d2, d1);
	}
}
template<class NODE, LarusDef::size_type Dim>
typename SPTree<NODE, Dim>::pTree SPTree<NODE, Dim>::_getpNeighborTree_path(
		SPDirection d1, SPDirection d2, SPDirection d3) {
	pTree pt = this->_walk_path(d1, d2, d3);
	if (pt != NULL_PTR) {
		return pt;
	}
	pt = this->_walk_path(d1, d3, d2);
	if (pt != NULL_PTR) {
		return pt;
	}
	pt = _walk_path(d2, d1, d3);
	if (pt != NULL_PTR) {
		return pt;
	}
	pt = _walk_path(d2, d3, d1);
	if (pt != NULL_PTR) {
		return pt;
	}
	pt = _walk_path(d3, d1, d2);
	if (pt != NULL_PTR) {
		return pt;
	}
	pt = _walk_path(d3, d2, d1);
	if (pt != NULL_PTR) {
		return pt;
	}
	return NULL_PTR;
}
template<class NODE, LarusDef::size_type Dim>
typename SPTree<NODE, Dim>::pTree SPTree<NODE, Dim>::getpNeighborTree(
		SPDirection d) {
	if (d == ErrSPDirection) {
		return NULL_PTR;
	}
	if (Dim == 2) {
		if (d > 7) {
			return NULL_PTR;
		}
	}
	if (d >= 4 && d <= 9) {
		return neighbor[d - 4];
	}
	SPDirection d1 = ErrSPDirection;
	SPDirection d2 = ErrSPDirection;
	SPDirection d3 = ErrSPDirection;
	int numd = Direction_Decompose(d, d1, d2, d3);
	if (numd == 2) {
		return _getpNeighborTree_path(d1, d2);
	} else if (numd == 3) {
		return _getpNeighborTree_path(d1, d2, d3);
	}
	return NULL_PTR;
}
//level_0====================================
template<class NODE, LarusDef::size_type Dim>
void SPTree<NODE, Dim>::_traversal_level(SPTree<NODE, Dim>::size_t le, pNode pn,
		pFun_SPTree visit, utPointer p) {
	if (pn == NULL_PTR) {
		return;
	} else {
		if (pn->getLevel() == le) {
			(*visit)(pn, p);
		} else {
			if (pn->hasChild()) {
				for (int i = 0; i < NODE::NUM_CELLS; i++) {
					pNode c = pn->child[i];
					if (c != NULL_PTR) {
						_traversal_level(le, c, visit, p);
					}
				}
			}
		}
	}
}

template<class NODE, LarusDef::size_type Dim>
void SPTree<NODE, Dim>::_traversal(SPTree<NODE, Dim>::pNode pn,
		SPTree<NODE, Dim>::pFun_SPTree visit, utPointer p) {
	if (pn == NULL_PTR) {
		return;
	} else {
		(*visit)(pn, p);
		if (pn->hasChild()) {
			for (int i = 0; i < NODE::NUM_CELLS; i++) {
				pNode c = pn->child[i];
				if (c != NULL_PTR) {
					_traversal(c, visit, p);
				}
			}
		}
	}
}

template<class NODE, LarusDef::size_type Dim>
void SPTree<NODE, Dim>::_traversal_b(SPTree<NODE, Dim>::pNode pn,
		SPTree<NODE, Dim>::pFun_SPTree visit, utPointer p) {
	if (pn == NULL_PTR) {
		return;
	}
	if (pn->hasChild()) {
		for (int i = 0; i < NODE::NUM_CELLS; i++) {
			pNode c = pn->child[i];
			if (c != NULL_PTR) {
				_traversal_b(c, visit, p);
			}
		}
	}
	(*visit)(pn, p);
}

template<class NODE, LarusDef::size_type Dim>
void SPTree<NODE, Dim>::_conditional_traversal(pNode pn,
		pFun_SPTree_condition t_condition, pFun_SPTree visit, utPointer p) {
	if (pn == NULL_PTR) {
		return;
	} else {
		(*visit)(pn, p);
		if (pn->hasChild()) {
			arrayList avt(NODE::NUM_CELLS);
			t_condition(avt, pn, p);
			for (int i = 0; i < NODE::NUM_CELLS; i++) {
				pNode c = pn->child[i];
				if (c != NULL_PTR && avt[i] == 1) {
					_conditional_traversal(c, t_condition, visit, p);
				}
			}
		}
	}
}

template<class NODE, LarusDef::size_type Dim>
void SPTree<NODE, Dim>::_refine_leaf( //
		pNode pn,  //Current node
		pFun_SPTree_condition_re t_condition, // condition function with return value
		utPointer utp_condition,  //untype pointer
		pFun_SPTree visit, //visit function
		utPointer utp_visit) {  //untype pointer
	if (pn == NULL_PTR || !condition_is_leaf(pn)) {
		return;
	} else {
		arrayList avt(NODE::NUM_CELLS);
		int ir = (*t_condition)(avt, pn, utp_condition);
		if (ir > 0) {
			_creatchildren_partial(pn, avt);
			for (int i = 0; i < NODE::NUM_CELLS; i++) {
				pNode c = pn->child[i];
				if (c != NULL_PTR && avt[i] == 1) {
					visit(c, utp_visit);
				}
			}
			for (int i = 0; i < NODE::NUM_CELLS; i++) {
				pNode c = pn->child[i];
				if (c != NULL_PTR && avt[i] == 1) {
					_refine_leaf(c, t_condition, utp_condition, visit,
							utp_visit);
				}
			}
		}
	}
}

//end level 0================================
template<class NODE, LarusDef::size_type Dim>
NODE* SPTree<NODE, Dim>::_father(pNode Current) {
	return Current->father;
}
template<class NODE, LarusDef::size_type Dim>
NODE* SPTree<NODE, Dim>::_child(pNode Current, int idx) {
	ASSERT(idx >= 0 && idx < NODE::NUM_CELLS);
	return Current->child[idx];
}
template<class NODE, LarusDef::size_type Dim>
const NODE* SPTree<NODE, Dim>::_father(NODE* Current) const {
	return Current->father;
}
template<class NODE, LarusDef::size_type Dim>
const NODE* SPTree<NODE, Dim>::_child(NODE* Current, int idx) const {
	ASSERT(idx >= 0 && idx < NODE::NUM_CELLS);
	return Current->child[idx];
}

//Traversal==================================
template<class NODE, LarusDef::size_type Dim>
void SPTree<NODE, Dim>::Traversal_partial(SPTree<NODE, Dim>::pNode Current,
		SPTree<NODE, Dim>::pFun_SPTree visit, utPointer p) {
	_traversal(Current, visit, p);
}
template<class NODE, LarusDef::size_type Dim>
void SPTree<NODE, Dim>::Traversal_level_RtL(
		SPTree<NODE, Dim>::pFun_SPTree visit, utPointer p) {
	for (size_t i = 0; i <= _maxlevel; i++) {
		_traversal_level(i, _root, visit, p);
	}
}
template<class NODE, LarusDef::size_type Dim>
void SPTree<NODE, Dim>::Traversal_level_LtR(
		SPTree<NODE, Dim>::pFun_SPTree visit, utPointer p) {
	for (int i = _maxlevel; i >= 0; i--) {
		_traversal_level(i, _root, visit, p);
	}
}
template<class NODE, LarusDef::size_type Dim>
void SPTree<NODE, Dim>::Traversal(SPTree<NODE, Dim>::pFun_SPTree visit,
		utPointer p) {
	ASSERT(_root!=NULL);
	_traversal(_root, visit, p);
}
template<class NODE, LarusDef::size_type Dim>
void SPTree<NODE, Dim>::Traversal_b(SPTree<NODE, Dim>::pFun_SPTree visit,
		utPointer p) {
	ASSERT(_root!=NULL);
	_traversal_b(_root, visit, p);
}
template<class NODE, LarusDef::size_type Dim>
void SPTree<NODE, Dim>::Conditional_Traversal(
		SPTree<NODE, Dim>::pFun_SPTree_condition t_condition,
		SPTree<NODE, Dim>::pFun_SPTree visit,
		utPointer p) {
	ASSERT(_root!=NULL);
	_conditional_traversal(_root, t_condition, visit, p);
}
template<class NODE, LarusDef::size_type Dim>
void SPTree<NODE, Dim>::Traversal_leaf_to_list(
		ListT<SPTree<NODE, Dim>::pNode>& listp) {
	ASSERT(_root!=NULL_PTR);
	_traversal(_root, visit_leaf_to_list, &listp);
}
//============================================
template<class NODE, LarusDef::size_type Dim>
int SPTree<NODE, Dim>::count_all() {
	ASSERT(_root!=NULL_PTR);
	size_t res = 0;
	_traversal(_root, visit_conunt_all, &res);
	return res;
}

template<class NODE, LarusDef::size_type Dim>
int SPTree<NODE, Dim>::count_level(SPTree<NODE, Dim>::size_t l) {
	ASSERT(_root!=NULL_PTR);
	size_t res = 0;
	arrayListT<size_t*> vp(2);
	vp[0] = &l;
	vp[1] = &res;
	_traversal(_root, visit_conunt_level, &vp);
	return res;
}
template<class NODE, LarusDef::size_type Dim>
int SPTree<NODE, Dim>::count_leaf() {
	ASSERT(_root!=NULL_PTR);
	size_t res = 0;
	_traversal(_root, visit_conunt_leaf, &res);
	return res;
}
template<class NODE, LarusDef::size_type Dim>
int SPTree<NODE, Dim>::count_leaf_at_level(size_t l) {
	ASSERT(_root!=NULL_PTR);
	size_t res = 0;
	arrayListT<size_t*> vp(2);
	vp[0] = &l;
	vp[1] = &res;
	_traversal(_root, visit_conunt_leaf_at_level, &vp);
	return res;
}
const CellPointLocation CHILD_PL[8][3] = {
//
		{ eCPL_M, eCPL_C, eCPL_M },  //
		{ eCPL_M, eCPL_M, eCPL_M },  //
		{ eCPL_C, eCPL_C, eCPL_M },  //
		{ eCPL_C, eCPL_M, eCPL_M },  //
		{ eCPL_M, eCPL_C, eCPL_C },  //
		{ eCPL_M, eCPL_M, eCPL_C },  //
		{ eCPL_C, eCPL_C, eCPL_C },  //
		{ eCPL_C, eCPL_M, eCPL_C },  //
		};

template<class NODE, LarusDef::size_type Dim>
void SPTree<NODE, Dim>::Refine(
		typename SPTree<NODE, Dim>::pFun_SPTree_condition_re fun_condition,
		utPointer utp_condition,
		typename SPTree<NODE, Dim>::pFun_SPTree fun_visit,
		utPointer utp_visit) {
	ASSERT(_root!=NULL);
	_refine_leaf(_root, fun_condition, utp_condition, fun_visit, utp_visit);
}

template<class NODE, LarusDef::size_type Dim>
int SPTree<NODE, Dim>::Creatchildren(typename SPTree<NODE, Dim>::pNode Current,
		SPNodeIdx idx_child) {
	if (idx_child == ErrSPIdx) {
		return 0;
	}
	int le = Current->getLevel();
	if ((Current->hasChild(idx_child) == false) && le < _maxlevel) {
		int ltmp = le + 1;
		Current->child[idx_child] = new NODE(Current, SPT_normal, ltmp,
				idx_child,
				typename NODE::Cell_type(
						Current->cell->getPoint(CHILD_PL[idx_child][0],
								CHILD_PL[idx_child][1], CHILD_PL[idx_child][2]),
						Current->cell->getPoint(
								CellPointLocation(CHILD_PL[idx_child][0] + 1),
								CellPointLocation(CHILD_PL[idx_child][1] + 1),
								CellPointLocation(
										CHILD_PL[idx_child][2] + 1))));
		Current->child[idx_child]->father = Current;
		return 1;
	}
	return 0;
}

template<class NODE, LarusDef::size_type Dim>
int SPTree<NODE, Dim>::Creatchildren_full(
		typename SPTree<NODE, Dim>::pNode Current) {
	int le = Current->getLevel();
	if (Current->hasChild() == false && le < _maxlevel) {
		int ltmp = le + 1;
		Current->child[0] = new NODE(Current, SPT_normal, ltmp, 0,
				typename NODE::Cell_type(
						Current->cell->getPoint(eCPL_M, eCPL_C, eCPL_M),
						Current->cell->getPoint(eCPL_C, eCPL_P, eCPL_C)));
		Current->child[0]->father = Current;
		Current->child[1] = new NODE(Current, SPT_normal, ltmp, 1,
				typename NODE::Cell_type(
						Current->cell->getPoint(eCPL_M, eCPL_M, eCPL_M),
						Current->cell->getPoint(eCPL_C, eCPL_C, eCPL_C)));
		Current->child[1]->father = Current;
		Current->child[2] = new NODE(Current, SPT_normal, ltmp, 2,
				typename NODE::Cell_type(
						Current->cell->getPoint(eCPL_C, eCPL_C, eCPL_M),
						Current->cell->getPoint(eCPL_P, eCPL_P, eCPL_C)));
		Current->child[2]->father = Current;
		Current->child[3] = new NODE(Current, SPT_normal, ltmp, 3,
				typename NODE::Cell_type(
						Current->cell->getPoint(eCPL_C, eCPL_M, eCPL_M),
						Current->cell->getPoint(eCPL_P, eCPL_C, eCPL_C)));
		Current->child[3]->father = Current;
		if (Dim == 3) {
			Current->child[4] = new NODE(Current, SPT_normal, ltmp, 4,
					typename NODE::Cell_type(
							Current->cell->getPoint(eCPL_M, eCPL_C, eCPL_C),
							Current->cell->getPoint(eCPL_C, eCPL_P, eCPL_P)));
			Current->child[4]->father = Current;
			Current->child[5] = new NODE(Current, SPT_normal, ltmp, 5,
					typename NODE::Cell_type(
							Current->cell->getPoint(eCPL_M, eCPL_M, eCPL_C),
							Current->cell->getPoint(eCPL_C, eCPL_C, eCPL_P)));
			Current->child[5]->father = Current;
			Current->child[6] = new NODE(Current, SPT_normal, ltmp, 6,
					typename NODE::Cell_type(
							Current->cell->getPoint(eCPL_C, eCPL_C, eCPL_C),
							Current->cell->getPoint(eCPL_P, eCPL_P, eCPL_P)));
			Current->child[6]->father = Current;
			Current->child[7] = new NODE(Current, SPT_normal, ltmp, 7,
					typename NODE::Cell_type(
							Current->cell->getPoint(eCPL_C, eCPL_M, eCPL_C),
							Current->cell->getPoint(eCPL_P, eCPL_C, eCPL_P)));
			Current->child[7]->father = Current;
		}
		return 1;
	}
	return 0;
}
template<class NODE, LarusDef::size_type Dim>
int SPTree<NODE, Dim>::Creatchildren_partial(
		typename SPTree<NODE, Dim>::pNode Current, bool b0, bool b1, bool b2,
		bool b3, bool b4, bool b5, bool b6, bool b7) {
	if (Current == NULL_PTR) {
		return 0;
	}
	if ((!b0) && (!b1) && (!b2) && (!b3) && (!b4) && (!b5) && (!b6) && (!b7)) {
		return 0;
	}
	if (b0) {
		Creatchildren(Current, MPM);
	}
	if (b1) {
		Creatchildren(Current, MMM);
	}
	if (b2) {
		Creatchildren(Current, PPM);
	}
	if (b3) {
		Creatchildren(Current, PMM);
	}
	if (Dim == 3) {
		if (b4) {
			Creatchildren(Current, MPP);
		}
		if (b5) {
			Creatchildren(Current, MMP);
		}
		if (b6) {
			Creatchildren(Current, PPP);
		}
		if (b7) {
			Creatchildren(Current, PMP);
		}
	}
	return 1;
}
template<class NODE, LarusDef::size_type Dim>
int SPTree<NODE, Dim>::_creatchildren_partial( //
		typename SPTree<NODE, Dim>::pNode pnode, //
		const arrayList& arr) {  //
	if (Dim == 2) {
		return Creatchildren_partial(pnode, (int(arr[0]) == 1),
				(int(arr[1]) == 1), (int(arr[2]) == 1), (int(arr[3]) == 1));
	} else { //dim ==3
		return Creatchildren_partial(pnode, (int(arr[0]) == 1),
				(int(arr[1]) == 1), (int(arr[2]) == 1), (int(arr[3]) == 1),
				(int(arr[4]) == 1), (int(arr[5]) == 1), (int(arr[6]) == 1),
				(int(arr[7]) == 1));
	}
}

template<class NODE, LarusDef::size_type Dim>
void SPTree<NODE, Dim>::CreatFullTree() {
	assert(_root!=NULL);
	_traversal(_root, visit_CreatFullTree, this);
}

template<class NODE, LarusDef::size_type Dim>
void SPTree<NODE, Dim>::CreatFullTree_to_Level(int level) {
	assert(_root!=NULL);
	assert(level >= 1);
	arrayListT<utPointer> arrp(2);
	arrp[0] = this;
	arrp[1] = &level;
	_traversal(_root, visit_CreatFullTree_to_Level, &arrp);
}

const SPDirection NODEIDX_NEIGBOR[8][7] = {
//
		{ SPD_IM, SPD_JP, SPD_KM, SPD_MP_XY, SPD_PM_YZ, SPD_MM_ZX }, //
		{ SPD_IM, SPD_JM, SPD_KM, SPD_MM_XY, SPD_MM_YZ, SPD_MM_ZX }, //
		{ SPD_IP, SPD_JP, SPD_KM, SPD_PP_XY, SPD_PM_YZ, SPD_MP_ZX }, //
		{ SPD_IP, SPD_JM, SPD_KM, SPD_PM_XY, SPD_MM_YZ, SPD_MP_ZX }, //
		{ SPD_IM, SPD_JP, SPD_KP, SPD_MP_XY, SPD_PP_YZ, SPD_PM_ZX }, //
		{ SPD_IM, SPD_JM, SPD_KP, SPD_MM_XY, SPD_MP_YZ, SPD_PM_ZX }, //
		{ SPD_IP, SPD_JP, SPD_KP, SPD_PP_XY, SPD_PP_YZ, SPD_PP_ZX }, //
		{ SPD_IP, SPD_JM, SPD_KP, SPD_PM_XY, SPD_MP_YZ, SPD_PP_ZX }, //
		};
template<class NODE>
void visit_ConnectNodes(NODE* pn, utPointer p) {
	if (!condition_is_root(pn)) {
		typedef SPTree<NODE, NODE::DIM>* pTree;
		typedef NODE* pNode;
		pTree ptree = CAST(pTree, p);
		int idx = pn->getIdx();
		pNode px = NULL_PTR, py = NULL_PTR, pz = NULL_PTR, pxy = NULL_PTR, pyz =
		NULL_PTR, pzx = NULL_PTR, pxyz = NULL_PTR;
		if (NODE::DIM == 2) {
			px = ptree->getNeighbor(pn, NODEIDX_NEIGBOR[idx][0]);
			py = ptree->getNeighbor(pn, NODEIDX_NEIGBOR[idx][1]);
			pxy = ptree->getNeighbor(pn, NODEIDX_NEIGBOR[idx][3]);
			pn->setNeighbors(px, py, pz, pxy, pyz, pzx, pxyz);
		} else {
			px = ptree->getNeighbor(pn, NODEIDX_NEIGBOR[idx][0]);
			py = ptree->getNeighbor(pn, NODEIDX_NEIGBOR[idx][1]);
			pz = ptree->getNeighbor(pn, NODEIDX_NEIGBOR[idx][2]);
			pxy = ptree->getNeighbor(pn, NODEIDX_NEIGBOR[idx][3]);
			pyz = ptree->getNeighbor(pn, NODEIDX_NEIGBOR[idx][4]);
			pzx = ptree->getNeighbor(pn, NODEIDX_NEIGBOR[idx][5]);
			pxyz = ptree->getNeighbor(pn, NODEIDX_NEIGBOR[idx][6]);
			pn->setNeighbors(px, py, pz, pxy, pyz, pzx, pxyz);
		}
	}
}

template<class NODE, LarusDef::size_type Dim>
void SPTree<NODE, Dim>::ConnectNodes() {
	assert(_root!=NULL);
	_traversal(_root, visit_ConnectNodes, this);
}

template<class NODE, LarusDef::size_type Dim>
void SPTree<NODE, Dim>::show_info() {
	std::cout << "SPTree ------\n";
	std::cout << "Dimension = " << Dim << "D\n";
	std::cout.precision(4);
	std::cout << "Root: (" << _root->cell->getPoint(eCPL_M, eCPL_M, eCPL_M).x
			<< ", " << _root->cell->getPoint(eCPL_M, eCPL_M, eCPL_M).y;
	if (Dim == 3) {
		std::cout << ", " << _root->cell->getPoint(eCPL_M, eCPL_M, eCPL_M)[2]
				<< ")\n";
	} else {
		std::cout << ")\n";
	}
	std::cout << "      (" << _root->cell->getPoint(eCPL_P, eCPL_P, eCPL_P).x
			<< ", " << _root->cell->getPoint(eCPL_P, eCPL_P, eCPL_P).y;
	if (Dim == 3) {
		std::cout << ", " << _root->cell->getPoint(eCPL_P, eCPL_P, eCPL_P)[2]
				<< ")\n";
	} else {
		std::cout << ")\n";
	}
	std::cout << "max_level set :" << getMaxLevel() << std::endl;
//std::cout<<"max_level real:"<<getMaxLevel()<<std::endl; //===
	int totalnode = 0;
	int totalleaf = 0;
	std::cout << "level  Num  leaf  ratio%\n";
	for (int i = 0; i <= _maxlevel; i++) {
		int nn = count_level(i);
		std::cout.flags(std::ios::right);
		std::cout.width(4);
		std::cout << i;
		std::cout.width(6);
		std::cout << nn;
		std::cout.width(6);
		int nl = count_leaf_at_level(i);
		std::cout << nl;
		std::cout.width(7);
		std::cout.precision(1);
		std::cout.setf(std::ios::fixed, std::ios::floatfield);
		std::cout << Float(nl) / Float(nn) * 100 << std::endl;
		totalnode += nn;
		totalleaf += nl;
	}
	std::cout.flags(std::ios::right);
	std::cout.width(10);
	std::cout << totalnode;
	std::cout.width(6);
	std::cout << totalleaf;
	std::cout.width(7);
	std::cout.precision(1);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout << Float(totalleaf) / Float(totalnode) * 100 << std::endl;
}
template<class NODE, LarusDef::size_type Dim>
typename SPTree<NODE, Dim>::pNode SPTree<NODE, Dim>::getpNode(
		const typename NODE::Cell_type::Point& p) {
	ASSERT(this->_root != NULL_PTR);
	if (!isInOnTree(p)) {
		return NULL_PTR;
	}
	typename NODE::Cell_type::Point inp = p;
	arrayListT<utPointer> vp(2);
	vp[0] = &inp;
	vp[1] = NULL_PTR;
	_conditional_traversal(_root, condition_point_at_which_child, visit_find,
			&vp);
	if (vp[1] == NULL_PTR) {
		return NULL_PTR;
	} else {
		return CAST(pNode, vp[1]);
	}
}
template<class NODE, LarusDef::size_type Dim>
bool SPTree<NODE, Dim>::isInTree(
		const typename NODE::Cell_type::Point& p) const {
	if (isEmpty()) {
		return false;
	} else {
		return _root->cell->isInCell(p);
	}
}

template<class NODE, LarusDef::size_type Dim>
bool SPTree<NODE, Dim>::isInOnTree(
		const typename NODE::Cell_type::Point& p) const {
	if (isEmpty()) {
		return false;
	} else {
		return _root->cell->isInOnCell(p);
	}
}

//===============================================
template<class NODE, LarusDef::size_type Dim>
typename SPTree<NODE, Dim>::pNode SPTree<NODE, Dim>::getAdjNeighbor(
		SPTree<NODE, Dim>::pNode Current, SPDirection d) {
// d 4 to 9
	pNode ca = NULL_PTR;    //common ancestor
	if (Current->father != NULL_PTR
			&& Current->isAdjacent(toSPNodeBoundary(d))) {
		ca = getAdjNeighbor(Current->father, d);
	} else {
		ca = Current->father;
	}
	pNode pt = NULL_PTR;
	if (ca != NULL_PTR && ca->hasChild()) {
		pt = _child(ca, Current->reflectIdx(d));
	} else if (ca == NULL_PTR) {
		pTree ptree = this->getpNeighborTree(d);
		if (ptree != NULL_PTR) {
			pt = ptree->getpRootNode();
		}
	} else {
		pt = ca;
	}
	return pt;
}
template<class NODE, LarusDef::size_type Dim>
typename SPTree<NODE, Dim>::pNode SPTree<NODE, Dim>::getCorNeighbor(
		SPTree<NODE, Dim>::pNode Current, SPDirection d) {
	pNode ca = NULL_PTR;    //common ancestor
	int flag = 0;
	if (Current->father != NULL_PTR && !dia_Sibling(Current->getEIdx(), d)) {
//Find a common ancestor
		if (out_cor(Current->getEIdx(), d)) {
			ca = getCorNeighbor(Current->father, d);
		} else {
			ca = getAdjNeighbor(Current->father,
					commonDirection(Current->getEIdx(), d));
		}
	} else {
		flag = 1;
		ca = Current->father;
	}
//Follow opposite path to locate the neighbor
	pNode pt = NULL_PTR;
	if (ca != NULL_PTR && ca->hasChild()) {
		pt = _child(ca, diagonalIdx(Current->getEIdx(), d));
	} else if (ca == NULL_PTR && flag == 1) {
		pTree ptree = this->getpNeighborTree(d);
		if (ptree != NULL_PTR) {
			pt = ptree->getpRootNode();
		}
	} else {
		pt = ca;
	}
	return pt;
}
template<class NODE, LarusDef::size_type Dim>
typename SPTree<NODE, Dim>::pNode SPTree<NODE, Dim>::getNeighbor_Adj_Cor(
		SPTree<NODE, Dim>::pNode Current, //
		SPDirection d) { //
// AdjNeighbor
	if (d >= 4 && d <= 9) {
		return getAdjNeighbor(Current, d);
	}
	if (d != ErrSPDirection && d <= 18) {
		return getCorNeighbor(Current, d);
	}
	return NULL_PTR;
}

template<class NODE, LarusDef::size_type Dim>
typename SPTree<NODE, Dim>::pNode SPTree<NODE, Dim>::getCorNeighbor_XYZ(
		SPTree<NODE, Dim>::pNode Current, SPDirection d) {
// d  17+1 to 17+8
	pNode ca = NULL_PTR;    //common ancestor
	if (Current->father != NULL_PTR && Current->isAdj_Vertex(d)) {
		ca = getCorNeighbor_XYZ(Current->father, d);
	} else {
		SPDirection newd = transfer_3Ddir_to_2Ddir(Current->getEIdx(), d);
		if (newd != ErrSPDirection) {
			ca = getNeighbor_Adj_Cor(Current->father, newd);
		} else {
			ca = Current->father;
		}
	}
	pNode pt = NULL_PTR;
	if (ca != NULL_PTR && ca->hasChild()) {
		pt = _child(ca, Current->reflectIdx_Vertex());
	} else if (ca == NULL_PTR) {
		pTree ptree = this->getpNeighborTree(d);
		if (ptree != NULL_PTR) {
			pt = ptree->getpRootNode();
		}
	} else {
		pt = ca;
	}
	return pt;
}
template<class NODE, LarusDef::size_type Dim>
typename SPTree<NODE, Dim>::pNode SPTree<NODE, Dim>::getNeighbor(
		SPTree<NODE, Dim>::pNode Current, SPDirection d) {
// AdjNeighbor
	if (d != ErrSPDirection && d <= 18) {
		return getNeighbor_Adj_Cor(Current, d);
	} else if (d >= 18 && d <= 25) {
		return getCorNeighbor_XYZ(Current, d);
	}
	return NULL_PTR;
}
//IO=============================================
//===============================================
template<class NODE, LarusDef::size_type Dim>
void SPTree<NODE, Dim>::_draw_to_vtk_leaf(FILE*& data, int numleaf) {
	fprintf(data, "POINTS %d float\n", numleaf * 8);

	utPointer utp = data;
	_traversal(this->_root, visit_draw_to_vtk, utp);

	fprintf(data, "\n");
	fprintf(data, "CELLS %d %d \n", numleaf, numleaf * 9);

	for (int i = 0; i < numleaf * 8; i++) {
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
	for (int i = 0; i < numleaf; i++) {
		fprintf(data, "%d \n", 11);
	}
	fprintf(data, "\n");
}
template<class NODE, LarusDef::size_type Dim>
void SPTree<NODE, Dim>::draw_to_vtk_leaf(std::string filename) {
	arrayList arr;
	int numleaf = count_leaf();
	FILE *data = open_file(filename, 1);
	vtk_unstructured_grid_head(data, "OcTree leafs");
	_draw_to_vtk_leaf(data, numleaf);
	fclose(data);
}
template<class NODE, LarusDef::size_type Dim>
void SPTree<NODE, Dim>::draw_to_vtk(FILE*& data) {
	utPointer utp = data;
	_traversal(this->_root, visit_draw_to_vtk, utp);
}

void visit_draw_gnuplot(pQTNode pn, utPointer p);
//void visit_draw_gnuplot(pOCNode pn, utPointer p); //don't use

template<class NODE, LarusDef::size_type Dim>
void SPTree<NODE, Dim>::draw_gnuplot(std::string filename, int mode) { //2D
	ASSERT(Dim == 2);
	FILE *data = open_file(filename, mode);
	//output root
	draw_boundary(_root, data);
	//Traversal_partial(_root, data);
	_traversal(_root, visit_draw_gnuplot, data);
	fclose(data);
}

template<class NODE, LarusDef::size_type Dim>
typename SPTree<NODE, Dim>::iterator SPTree<NODE, Dim>::begin() {
	pNode c = getFirstChild(_root);
	if (c == NULL_PTR) {
		return _root;
	} else {
		pNode resc = c;
		while (c != NULL_PTR) {
			resc = c;
			c = getFirstChild(c);
		}
		return resc;
	}
}
template<class NODE, LarusDef::size_type Dim>
typename SPTree<NODE, Dim>::const_iterator SPTree<NODE, Dim>::begin() const {
	pNode c = getFirstChild(_root);
	if (c == NULL_PTR) {
		return _root;
	} else {
		pNode resc = c;
		while (c != NULL_PTR) {
			resc = c;
			c = getFirstChild(c);
		}
		return resc;
	}
}
template<class NODE, LarusDef::size_type Dim>
typename SPTree<NODE, Dim>::iterator SPTree<NODE, Dim>::end() {
	return _root;
}

template<class NODE, LarusDef::size_type Dim>
typename SPTree<NODE, Dim>::const_iterator SPTree<NODE, Dim>::end() const {
	return _root;
}

//===============================================
//QuadTree=======================================
//===============================================
typedef SPTree<QTNode, 2> QuadTree;
typedef QuadTree* pQuadTree;

//===============================================
//OCTree=========================================
//===============================================
typedef SPTree<OCNode, 3> OcTree;
typedef OcTree* pOCTree;

//Function out of class==========================
template<class NODE, LarusDef::size_type Dim>
void toList_leaf(SPTree<NODE, Dim>* ptree,  //
		ListT<typename SPTree<NODE, Dim>::pNode>& list) {
	_IF_TRUE_RETRUN(ptree==NULL_PTR || ptree->isEmpty());
	ptree->Traversal(visit_to_list_leaf, &list);
}

template<class NODE, LarusDef::size_type Dim>
void toList_leaf(SPTree<NODE, Dim>* ptree,
		ListT<typename SPTree<NODE, Dim>::pNode>& list, CSAxis axis, //
		Float loc) {  //slice
	_IF_TRUE_RETRUN(ptree==NULL_PTR || ptree->isEmpty());
	arrayListT<utPointer> arr(3);
	arr[0] = &axis;
	arr[1] = &loc;
	arr[2] = &list;
	ptree->Traversal(visit_to_list_leaf_slice, &arr);
}

template<class NODE, LarusDef::size_type Dim>
void toList_leaf(SPTree<NODE, Dim>* ptree,
		ListT<typename SPTree<NODE, Dim>::pNode>& list, CSAxis axis, Float loc,
		int level) {
//slice
	_IF_TRUE_RETRUN(ptree==NULL_PTR || ptree->isEmpty());
	arrayListT<utPointer> arr(4);
	arr[0] = &axis;
	arr[1] = &loc;
	arr[2] = &list;
	arr[3] = &level;
	ptree->Traversal(visit_to_list_leaf_slice_level, &arr);
}

}

#endif /* OCTREE_H_ */
