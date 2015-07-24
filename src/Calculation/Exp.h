/*
 * Exp.h
 *
 *  Created on: May 15, 2015
 *      Author: zhou
 */

#ifndef CALCULATION_EXP_H_
#define CALCULATION_EXP_H_

#include "../TypeDef.h"
#include "CalDef.h"
#include "../Grid/SPTree.h"
#include "../Grid/SPTreeNode.h"
#include "../Grid/Forest.h"
#include "../Algebra/Arithmetic.h"
//#include "../Algebra/Expression.h"
#include "Boundary.h"

namespace Larus {

template<typename VALUE, typename NODE>
class ExpTNode {
public:
	//typedef======
	static const int IDX_CONST = -100000;
	//Data=========
	int idx;
	NODE* pnode;
	VALUE val;
	//Function=====
	ExpTNode();
	//ExpTNode(int, const VALUE&);
	ExpTNode(int i, NODE* pn, const VALUE& v) { //
		pnode = pn;
		idx = i;
		val = v;
	}
	ExpTNode(VALUE);
	ExpTNode(const ExpTNode<VALUE, NODE>&);
	ExpTNode<VALUE, NODE>& operator=(const ExpTNode<VALUE, NODE> &);
	~ExpTNode() {
	}

	void reconstract(int, const VALUE&, NODE*);

	bool operator==(const ExpTNode<VALUE, NODE> &) const;
	bool operator!=(const ExpTNode<VALUE, NODE> &) const;
	bool operator<(const ExpTNode<VALUE, NODE> &) const;
	bool operator>(const ExpTNode<VALUE, NODE> &) const;
	bool isConstNode() const;
	ExpTNode<VALUE, NODE> operator+(const ExpTNode<VALUE, NODE>&);
	ExpTNode<VALUE, NODE> operator-();
	ExpTNode<VALUE, NODE> operator-(const ExpTNode<VALUE, NODE>&);
	void show() const;
};

template<typename VALUE, typename NODE>
ExpTNode<VALUE, NODE>::ExpTNode() {
	idx = IDX_CONST;
	pnode = NULL_PTR;
	val = 0.0;
}
//template<typename VALUE, typename NODE>
//ExpTNode<VALUE, NODE>::ExpTNode(int i, const VALUE& v) {
//	pnode = NULL_PTR;
//	idx = i;
//	val = v;
//}
template<typename VALUE, typename NODE>
ExpTNode<VALUE, NODE>::ExpTNode(VALUE v) {
	idx = IDX_CONST;
	pnode = NULL_PTR;
	val = v;
}
template<typename VALUE, typename NODE>
ExpTNode<VALUE, NODE>::ExpTNode(const ExpTNode<VALUE, NODE>& en) {
	idx = en.idx;
	val = en.val;
	pnode = en.pnode;
}
template<typename VALUE, typename NODE>
ExpTNode<VALUE, NODE>& ExpTNode<VALUE, NODE>::operator=(
		const ExpTNode<VALUE, NODE>& en) {
	idx = en.idx;
	val = en.val;
	pnode = en.pnode;
	return *this;
}
template<typename VALUE, typename NODE>
void ExpTNode<VALUE, NODE>::reconstract(int i, const VALUE& v, NODE* pn) {
	idx = i;
	val = v;
	pnode = pn;
}
template<typename VALUE, typename NODE>
bool ExpTNode<VALUE, NODE>::operator==(const ExpTNode<VALUE, NODE> & en) const {
	return en.idx == this->idx ? true : false;
}
template<typename VALUE, typename NODE>
bool ExpTNode<VALUE, NODE>::operator!=(const ExpTNode<VALUE, NODE> & en) const {
	return en.idx != this->idx ? true : false;
}
template<typename VALUE, typename NODE>
bool ExpTNode<VALUE, NODE>::operator<(const ExpTNode<VALUE, NODE> & en) const {
	return this->idx < en.idx ? true : false;
}
template<typename VALUE, typename NODE>
bool ExpTNode<VALUE, NODE>::operator>(const ExpTNode<VALUE, NODE> & en) const {
	return this->idx > en.idx ? true : false;
}
template<typename VALUE, typename NODE>
bool ExpTNode<VALUE, NODE>::isConstNode() const {
	return this->idx == IDX_CONST ? true : false;
}
template<typename VALUE, typename NODE>
ExpTNode<VALUE, NODE> ExpTNode<VALUE, NODE>::operator+(
		const ExpTNode<VALUE, NODE> & en) {
	ASSERT(this->idx == en.idx);
	this->val += en.val;
	return *this;
}
template<typename VALUE, typename NODE>
ExpTNode<VALUE, NODE> ExpTNode<VALUE, NODE>::operator-() {
	this->val = -this->val;
	return *this;
}

template<typename VALUE, typename NODE>
ExpTNode<VALUE, NODE> ExpTNode<VALUE, NODE>::operator-(
		const ExpTNode<VALUE, NODE> & en) {
	ASSERT(this->idx == en.idx);
	this->val -= en.val;
	return *this;
}
template<typename VALUE, typename NODE>
void ExpTNode<VALUE, NODE>::show() const {
	std::cout.flags(std::ios::right);
	std::cout.width(9);
	if (idx == IDX_CONST) {
		std::cout << "CONST" << "  ";
	} else {
		std::cout << idx << "  ";
	}
	std::cout.precision(5);
	std::cout.flags(std::ios::right);
	std::cout.width(12);
	std::cout << val << "\n";
}

//end of ExpTNode class =========================
//===============================================
//===============================================
template<typename TERM>
class Exp: public BSTree<TERM> {
public:
	typedef TERM ExpTerm;
	typedef typename BSTree<ExpTerm>::reference reference;
	typedef typename BSTree<ExpTerm>::const_reference const_reference;
	typedef typename BSTree<ExpTerm>::pointer pointer;
	typedef typename BSTree<ExpTerm>::const_pointer const_pointer;
	typedef typename BSTree<ExpTerm>::iterator iterator;
	typedef typename BSTree<ExpTerm>::const_iterator const_iterator;

	typedef typename BSTree<ExpTerm>::difference_type difference_type;
	typedef typename BSTree<ExpTerm>::size_type size_type;

protected:
	void _insert(BinaryTreeNode<ExpTerm>*& Current, const ExpTerm &data) {
		if (Current == NULL) {
			Current = new BinaryTreeNode<ExpTerm>(data);
			assert(Current != NULL);
		} else if (data == Current->m_value) {           //operation ==
			Current->m_value = Current->m_value + data;  //operation +
		} else if (data < Current->m_value) {            //operation <
			_insert(Current->lchild, data);
			Current->lchild->father = Current;
		} else {
			_insert(Current->rchild, data);
			Current->rchild->father = Current;
		}
	}
public:
	Exp() :
			BSTree<ExpTerm>() {
	}
	//Expression operator=(const Expression &);
	void Insert(const ExpTerm &x) {    //override
		if (x.val == 0.0) {
			return;
		}
		if (this->_root->lchild == NULL_PTR) {
			_insert(this->_root->lchild, x);
			this->_root->lchild->father = this->_root;
		} else {
			_insert(this->_root->lchild, x);
		}
	}
	//bool operator==(const Expression&) const;
	Exp operator+(const Exp &exp) {
		Exp res(*this);
		for (const_iterator iter = exp.begin(); iter != exp.end(); iter++) {
			res.Insert((*iter));
		}
		return res;
	}
	Exp operator-(const Exp &exp) {
		Exp res(*this);
		for (const_iterator iter = exp.begin(); iter != exp.end(); iter++) {
			iterator itertmp = res.Find((*iter));
			if (itertmp.isExist()) {
				(*itertmp) = (*itertmp) - (*iter);
			} else {
				ExpTerm tmp((*iter));
				res.Insert(-tmp);
			}
		}
		res.trim_zero();
		return res;
	}
	Exp operator*(const Float &a) {
		if (a == 0.0) {
			this->clear();
			return *this;
		}
		for (iterator iter = this->begin(); iter != this->end(); iter++) {
			iter->val = iter->val * a;
		}
		return *this;
	}
	Exp operator/(const Float &a) {
		for (iterator iter = this->begin(); iter != this->end(); iter++) {
			iter->val = iter->val / a;
		}
		return *this;
	}

	void plus(const Exp &exp) {
		if (exp.empty()) {
			return;
		}
		for (const_iterator iter = exp.begin(); iter != exp.end(); iter++) {
			this->Insert((*iter));
		}
	}
	void minus(const Exp &exp) {
		for (const_iterator iter = exp.begin(); iter != exp.end(); iter++) {
			iterator itertmp = this->Find((*iter));
			if (itertmp.isExist()) {
				(*itertmp) = (*itertmp) - (*iter);
			} else {
				ExpTerm tmp((*iter));
				this->Insert(-tmp);
			}
		}
		this->trim_zero();
	}
	void times(const Float &a) {
		if (a == 0.0) {
			this->clear();
			return;
		}
		for (iterator iter = this->begin(); iter != this->end(); iter++) {
			iter->val = iter->val * a;
		}
	}
	void divide(const Float &a) {
		for (iterator iter = this->begin(); iter != this->end(); iter++) {
			iter->val = iter->val / a;
		}
	}

	void trim_zero() { //can be improved
		if (this->_root->lchild == NULL_PTR) {
			return;
		}
		int num = this->size();
		int numz = count_zero();
		if (numz == 0) {
			return;
		}
		BinaryTreeNode<ExpTerm>* pnode = this->MinValNode();
		if (pnode != NULL_PTR) {
			BinaryTreeNode<ExpTerm>* cur = pnode;
			do {
				if (cur->m_value.val == 0.0) {
					this->_delete(cur);
					continue;
				}
				cur = this->NextNode(cur);
			} while (cur != this->_root);
		}
		while (this->size() != num - numz) {
			trim_zero();
		}
	}
	int count_zero() const {
		int res = 0;
		for (const_iterator iter = this->begin(); iter != this->end(); iter++) {
			if (iter->val == 0.0) {
				res += 1;
			}
		}
		return res;
	}

	void show() const {
		for (const_iterator iter = this->begin(); iter != this->end(); iter++) {
			iter->show();
		}
	}
};

//===============================================
template<typename TERM>
Exp<TERM> operator-(const Exp<TERM> & exp1, const Exp<TERM> & exp2) {
	Exp<TERM> res(exp1);
	for (typename Exp<TERM>::const_iterator iter = exp2.begin();
			iter != exp2.end(); iter++) {
		typename Exp<TERM>::iterator itertmp = res.Find((*iter));
		if (itertmp.isExist()) {
			(*itertmp) = (*itertmp) - (*iter);
		} else {
			TERM tmp((*iter));
			res.Insert(-tmp);
		}
	}
	res.trim_zero();
	return res;
}
template<typename TERM>
Exp<TERM> operator+(const Exp<TERM> & exp1, const Exp<TERM> & exp2) {
	Exp<TERM> res(exp1);
	for (typename Exp<TERM>::const_iterator iter = exp2.begin();
			iter != exp2.end(); iter++) {
		res.Insert((*iter));
	}
	return res;
}
template<typename TERM>
Exp<TERM> operator*(const Float& a, const Exp<TERM>& exp) {
	if (a == 0.0) {
		return Exp<TERM>();
	}
	Exp<TERM> res(exp);
	for (typename Exp<TERM>::iterator iter = res.begin(); iter != res.end();
			iter++) {
		iter->val = iter->val * a;
	}
	return res;
}

typedef ExpTNode<Float, QTNode> ExpTerm2D;
typedef ExpTNode<Float, OCNode> ExpTerm3D;

typedef Exp<ExpTerm2D> Exp2D;
typedef Exp<ExpTerm3D> Exp3D;

inline Exp2D linear_interpolation( //
		const Float& x,  //
		const Float& x0, const Exp2D& y0, //
		const Float& x1, const Exp2D& y1 //
		) {
	Exp2D res = y1 - y0;
	res.times((x - x0) / (x1 - x0));
	res.plus(y0);
	return res;
}
template<typename TERM>
Exp<TERM> second_order_interpolation( //
		const Float& x,  //
		const Float& x1, const Exp<TERM>& y1, //
		const Float& x2, const Exp<TERM>& y2, //
		const Float& x3, const Exp<TERM>& y3  //
		) {
	Exp<TERM> res(y1);
	res.times((x - x2) * (x - x3) / (x1 - x2) / (x1 - x3));
	Exp<TERM> res2(y2);
	res2.times((x - x1) * (x - x3) / (x2 - x1) / (x2 - x3));
	Exp<TERM> res3(y3);
	res3.times((x - x1) * (x - x2) / (x3 - x1) / (x3 - x2));
	res.plus(res2);
	res.plus(res3);
	return res;
}
template<typename TERM>
Exp<TERM> linear_gradient( //
		const Float& x0, const Exp<TERM>& y0, //
		const Float& x1, const Exp<TERM>& y1 //
		) {
	Exp<TERM> res = y1 - y0;
	res.times(1.0 / (x1 - x0));
	return res;
}
template<typename TERM>
Exp<TERM> second_order_gradient( //
		const Float& x,  //
		const Float& x1, const Exp<TERM>& y1, //
		const Float& x2, const Exp<TERM>& y2, //
		const Float& x3, const Exp<TERM>& y3  //
		) {
	Exp<TERM> res(y1);
	res.times((2.0 * x - x2 - x3) / (x1 - x2) / (x1 - x3));
	Exp<TERM> res2(y2);
	res2.times((2.0 * x - x1 - x3) / (x2 - x1) / (x2 - x3));
	Exp<TERM> res3(y3);
	res3.times((2.0 * x - x1 - x2) / (x3 - x1) / (x3 - x2));
	res.plus(res2);
	res.plus(res3);
	return res;
}

// SPDirection
//                         4        5         6       7
const int PoN_table[] = { -1, 1, -1, 1 };  // positive or nagative
const CSAxis D_table[] = { CSAxis_X, CSAxis_Y, CSAxis_X, CSAxis_Y };
const CSAxis P_table[] = { CSAxis_Y, CSAxis_X, CSAxis_Y, CSAxis_X };

void interpolate_expression_on_axis_2f( // 2D QuadTree Node
		pQTNode pc,  //node
		CSAxis axis, //axix
		Float dis,   //distance to center of pn with sign
		Exp2D& exp   //Expression
		);
void interpolate_expression_on_axis_2f( // 2D QuadTree Node
		pQTNode pc,       //node
		CSAxis axis,      //axix
		Float dis,        //distance to center of pn with sign
		Exp2D& exp,        //Expression
		BCManager<Dimension_2D>& bcm //
		);
void interpolate_expression_on_axis_2b( // 2D QuadTree Node
		pQTNode pc,  //node
		CSAxis axis, //axix
		Float dis,   //distance to center of pn with sign
		Exp2D& exp   //Expression
		);
void interpolate_expression_on_axis_3fb( // 2D QuadTree Node
		pQTNode pc,  //node
		CSAxis axis, //axix
		Float dis,   //distance to center of pn with sign
		Exp2D& exp   //Expression
		);
void interpolate_expression_on_axis_3fb( // 2D QuadTree Node
		pQTNode pc,  //node
		CSAxis axis, //axix
		Float dis,   //distance to center of pn with sign
		Exp2D& exp,   //Expression
		BCManager<Dimension_2D>& bcm);
// =======================
// =======================
template<typename NODE>
int T_find_neighbor_2( //
		NODE* pc,    //node
		SPDirection direction, //direction
		NODE*& pn    //one neighbor as output
		) {
	if (direction == ErrSPDirection || pc == NULL_PTR) {
		return -3;
	}
	pn = pc->getNeighborFast(direction);
	if (pn == NULL_PTR) {
		return 0;
	}
	return 1;
}

inline int find_neighbor_2( //
		pQTNode pc,    //node
		SPDirection direction, //direction
		pQTNode& pn,    //one neighbor as output
		BCManager<Dimension_2D>& bcm) {
	int flag = T_find_neighbor_2(pc, direction, pn);
	if (flag == 0) {
		pn = bcm.find_ghost(getIDX(pc), direction);
	}
	if (pn == NULL_PTR) {
		return 0;
	}
	return 1;
}

template<typename NODE>
int T_find_neighbor_3ff( //
		NODE* pc,    //node
		SPDirection df, //direction forward
		NODE*& pn,   //one neighbor as output
		NODE*& pnn) {
	if (df == ErrSPDirection || pc == NULL_PTR) {
		return -3;
	}
	pn = pc->getNeighborFast(df);
	if (pn == NULL_PTR) {
		pnn = NULL_PTR;
		return 0;
	}
	pnn = pn->getNeighborFast(df);
	if (pnn == NULL_PTR) {
		return 1;
	}
	return 2;
}
template<typename NODE>
int T_find_neighbor_3fb( //
		NODE* pc,    //node
		SPDirection df, //direction forward
		NODE*& pf,   //one neighbor as output
		NODE*& pb) {
	if (df == ErrSPDirection || pc == NULL_PTR) {
		return -3;
	}
	pf = pc->getNeighborFast(df);
	SPDirection opd = oppositeDirection(df);
	pb = pc->getNeighborFast(opd);
	if (pf == NULL_PTR && pb == NULL_PTR) {
		return 0;
	}
	if (pf == NULL_PTR) {
		return 1;
	}
	if (pb == NULL_PTR) {
		return -1;
	}
	return 3;
}
template<typename NODE>
int T_find_neighbor_3bb( //
		NODE* pc,    //node
		SPDirection df, //direction forward
		NODE*& pb,   //one neighbor as output
		NODE*& pbb) {
	if (df == ErrSPDirection || pc == NULL_PTR) {
		return -3;
	}
	SPDirection opd = oppositeDirection(df);
	return T_find_neighbor_3ff(pc, opd, pb, pbb);
}

template<typename NODE>
bool T_dis_check( //
		NODE* pn,        //node
		CSAxis axis,     //axix
		Float dis        //distance to center of pn
		) {
	switch (axis) {
	case ErrCSAxis:
		return false;
	case CSAxis_X: {
		if (ABS(dis) > pn->cell->gethDx()) {
			return false;
		} else {
			return true;
		}
	}
	case CSAxis_Y: {
		if (ABS(dis) > pn->cell->gethDy()) {
			return false;
		} else {
			return true;
		}
	}
	case CSAxis_Z: {
		if (NODE::DIM == 3) {
			if (ABS(dis) > pn->cell->gethDz()) {
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
void T_exp_1node( //
		NODE* po,  //node
		Exp<ExpTNode<Float, NODE> >& exp) {
	ExpTNode<Float, NODE> et(int(po->data->aCenterData[Idx_IDX]), po, 1.0);
	exp.Insert(et);
}
template<typename NODE>
void visit_averange_exp_from_leaf(NODE* pn, utPointer up) {
	if (condition_is_leaf(pn)) {
		arrayListT<utPointer>& arrp = (*CAST(arrayListT<utPointer>*, up));
		int& count = (*CAST(int*, arrp[0]));
		ListT<NODE*>& ln = (*CAST(ListT<NODE*>*, arrp[1]));
		if (pn->data != NULL_PTR) {
			ln.push_back(pn);
			count = count + 1;
		}
	}
}
template<typename NODE>
void T_averange_exp_from_leaf(        //
		NODE* pn,  //pNode
		Exp<ExpTNode<Float, NODE> >& exp //Expression
		) {
//improve for special case
	if (condition_is_leaf(pn)) {
		ExpTNode<Float, NODE> et(int(pn->data->aCenterData[Idx_IDX]), pn, 1.0);
		exp.Insert(et);
		return;
	}
//========================
	int count = 0;
	ListT<NODE*> listnode;
	arrayListT<utPointer> arrp(3);
	arrp[0] = &count;
	arrp[1] = &listnode;
	pn->Traversal(visit_averange_exp_from_leaf, &arrp);
	for (typename ListT<NODE*>::iterator iter = listnode.begin();
			iter != listnode.end(); iter++) {
		ExpTNode<Float, NODE> et(int((*iter)->data->aCenterData[Idx_IDX]),
				(*iter), 1.0 / count);
		exp.Insert(et);
	}
}

inline int _exp_2node_relation( //
		pQTNode po,     //node
		pQTNode ps,     //node
		SPDirection dir,     //direction
		Point2D& point, //tmp_point
		Exp2D& exp //expression
		) {
	//dir 4 5 6 7
	ASSERT(int(dir) >= 4 && int(dir) <= 7);
	int pon = PoN_table[int(dir) - 4];
	CSAxis ad = D_table[int(dir) - 4];
	CSAxis ap = P_table[int(dir) - 4];
	//case 1 two node are aqual
	if (po->getLevel() == ps->getLevel() && ps->hasChild() == false) {
		point.x = ps->cell->getCenterPoint(CSAxis_X);
		point.y = ps->cell->getCenterPoint(CSAxis_Y);
		ExpTerm2D term(ps->data->aCenterData[Idx_IDX], ps, 1.0);
		exp.Insert(term);
		return 1;
	}
	//case 2  fine to Coarse
	if (po->getLevel() > ps->getLevel()) {
		if (ad == CSAxis_X) {
			point.x = ps->cell->getCenterPoint(CSAxis_X);
			point.y = po->cell->getCenterPoint(CSAxis_Y);
		} else {   // axis== CSAxis_Y
			point.x = po->cell->getCenterPoint(CSAxis_X);
			point.y = ps->cell->getCenterPoint(CSAxis_Y);
		}
		Float tmp_dis = point[ap] - ps->cell->getCenterPoint()[ap];
		Exp2D tmp_exp;
		interpolate_expression_on_axis_2f( // 2D QuadTree Node
				ps,       //node
				ap,       //axix
				tmp_dis,  //distance to center of pn
				tmp_exp);
		exp.plus(tmp_exp);
		return 2;
	}
	//case 3  Coarse to fine
	if (po->getLevel() == ps->getLevel()) {
		const int _CHILD_LOOKUP_2D[2][2][2] =   //
				//  axis choose  x or y
				//  direction choose negative or positive
				//  choose the two child 0 and 1
				{ { { 0, 1 }, { 2, 3 } }, { { 1, 3 }, { 0, 2 } } }; //
		int ch1 = _CHILD_LOOKUP_2D[int(ad)][(-pon == 1) ? 0 : 1][0];
		int ch2 = _CHILD_LOOKUP_2D[int(ad)][(-pon == 1) ? 0 : 1][1];
		if (ps->hasChild(ch1) && ps->hasChild(ch2)) {
			//Normal case
			//There two children
			Exp2D tmp_exp1, tmp_exp2;
			pQTNode pc1 = ps->child[ch1];
			pQTNode pc2 = ps->child[ch2];
			T_averange_exp_from_leaf(pc1, tmp_exp1);
			T_averange_exp_from_leaf(pc2, tmp_exp2);
			point = calMid(pc1->cell->getCenterPoint(),
					pc2->cell->getCenterPoint());
			tmp_exp1.plus(tmp_exp2); //  tmp_exp1+tmp_exp2;
			tmp_exp1.times(0.5);     // (tmp_exp1+tmp_exp2)*0.5;
			exp.plus(tmp_exp1);
		} else {
			//special case the children is not complete
			Exp2D tmp_exp;
			T_averange_exp_from_leaf(ps, tmp_exp);
			point = ps->cell->getCenterPoint();
			exp.plus(tmp_exp);
		}
		return 3;
	}
	return -1;
}

inline void _exp_interplolate_3node( //
		pQTNode pm,  //node
		pQTNode po,  //node
		pQTNode pp,  //node
		CSAxis axis, //axix
		Float dis,   //distance to center of pn
		Exp2D& exp) {
	Point2D point_m; //tmp_point
	Point2D point_p; //tmp_point
	Exp2D exp_m;   //data res
	Exp2D exp_p;   //data res
	SPDirection dp =
			(axis == CSAxis_X) ?
					(dis < 0 ? SPD_IM : SPD_IP) : (dis < 0 ? SPD_JM : SPD_JP);
	SPDirection dm = oppositeDirection(dp);

	_exp_2node_relation(po,  //node
			pp,              //node p
			dp,              //dirction
			point_p,         //tmp_point
			exp_p            //data res
			);
	_exp_2node_relation(po,  //node
			pm,              //node
			dm,              //dirction
			point_m,         //tmp_point
			exp_m            //data res
			);
	Float x = po->cell->getCenterPoint()[axis] + dis;
	Float x1 = po->cell->getCenterPoint()[axis];
	Exp2D y1;
	y1.Insert(ExpTerm2D(po->data->aCenterData[Idx_IDX], po, 1.0));
	Float x2 = point_m[axis];
	Float x3 = point_p[axis];
	exp = second_order_interpolation(x, x1, y1, x2, exp_m, x3, exp_p);
}

inline void _exp_interplolate_2node( //
		pQTNode po,  //node
		pQTNode ps,  //node
		CSAxis axis, //aix
		Float dis,   //distance to center of pn
		Exp2D& exp) {
	Point2D tmp_point;
	Exp2D tmp_exp;
	// axis and dis to direction
	SPDirection dir =
			(axis == CSAxis_X) ?
					(dis < 0 ? SPD_IM : SPD_IP) : (dis < 0 ? SPD_JM : SPD_JP);
	_exp_2node_relation( //
			po,    //node
			ps,    //node
			dir,   //direction
			tmp_point, //tmp_point
			tmp_exp);
	Float x = po->cell->getCenterPoint(axis) + dis;
	Float x1 = po->cell->getCenterPoint(axis);
	Float x2 = tmp_point[axis];
	Exp2D y1;
	//po must be leaf node
	y1.Insert(ExpTerm2D(int(po->data->aCenterData[Idx_IDX]), po, 1.0));
	exp = linear_interpolation(x, x1, y1, x2, tmp_exp);
}

}

#endif /* CALCULATION_EXP_H_ */
