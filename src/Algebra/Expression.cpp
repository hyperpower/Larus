/************************
 //  \file   Expression.cpp
 //  \brief
 // 
 //  \author czhou
 //  \date   4 févr. 2015 
 ***********************/

#include "Expression.h"

namespace Larus {
void Expression::_insert(BinaryTreeNode<ExpTerm> *&Current,
		const ExpTerm &data) {
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

void Expression::Insert(const ExpTerm &x) {
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

Expression Expression::operator+(const Expression & exp) {
	Expression res(*this);
	for (const_iterator iter = exp.begin(); iter != exp.end(); iter++) {
		res.Insert((*iter));
	}
	return res;
}

Expression Expression::operator-(const Expression &exp) {
	Expression res(*this);
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

Expression Expression::operator*(const Float & a) {
	if (a == 0.0) {
		this->clear();
		return *this;
	}
	for (iterator iter = this->begin(); iter != this->end(); iter++) {
		iter->val = iter->val * a;
	}
	return *this;
}

Expression Expression::operator/(const Float & a) {
	for (iterator iter = this->begin(); iter != this->end(); iter++) {
		iter->val = iter->val / a;
	}
	return *this;
}

void Expression::show() const {
	if(this->empty()){
		return;
	}
	for (const_iterator iter = this->begin(); iter != this->end(); iter++) {
		iter->show();
	}
}

int Expression::count_zero() const {
	int res = 0;
	for (const_iterator iter = this->begin(); iter != this->end(); iter++) {
		if (iter->val == 0.0) {
			res += 1;
		}
	}
	return res;
}

void Expression::trim_zero() {
	if (this->_root->lchild == NULL_PTR) {
		return;
	}
	int num = size();
	int numz = count_zero();
	BinaryTreeNode<ExpTerm>* pnode = MinValNode();
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
void Expression::plus(const Expression & exp) {
	if(exp.empty()){
		return;
	}
	for (const_iterator iter = exp.begin(); iter != exp.end(); iter++) {
		this->Insert((*iter));
	}
}
void Expression::minus(const Expression & exp) {
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
void Expression::times(const Float & a) {
	if (a == 0.0) {
		this->clear();
		return;
	}
	for (iterator iter = this->begin(); iter != this->end(); iter++) {
		iter->val = iter->val * a;
	}
}
void Expression::divide(const Float & a) {
	for (iterator iter = this->begin(); iter != this->end(); iter++) {
		iter->val = iter->val / a;
	}
}

//===============================================
Expression operator-(const Expression & exp1, const Expression & exp2) {
	Expression res(exp1);
	for (Expression::const_iterator iter = exp2.begin(); iter != exp2.end();
			iter++) {
		Expression::iterator itertmp = res.Find((*iter));
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
Expression operator+(const Expression & exp1, const Expression & exp2) {
	Expression res(exp1);
	for (Expression::const_iterator iter = exp2.begin(); iter != exp2.end(); iter++) {
		res.Insert((*iter));
	}
	return res;
}
Expression operator*(const Float& a, const Expression& exp) {
	if (a == 0.0) {
		return Expression();
	}
	Expression res(exp);
	for (Expression::iterator iter = res.begin(); iter != res.end(); iter++) {
		iter->val = iter->val * a;
	}
	return res;
}

}
