/*
 * Expression.h
 *
 *  Created on: Feb 2, 2015
 *      Author: zhou
 */

#ifndef _EXPRESSION_H_
#define _EXPRESSION_H_

#include "../TypeDef.h"
#include "../Utility/BinaryTree.h"
#include "../Utility/BSTree.h"
#include <stdint.h>

namespace Larus {

template<typename VALUE>
class ExpNode {
public:
	//typedef======
	static const int IDX_CONST = -100000;
	//Data=========
	int idx;
	VALUE val;
	//Function=====
	ExpNode();
	ExpNode(int, const VALUE&);
	ExpNode(VALUE);
	ExpNode(const ExpNode<VALUE>&);
	ExpNode<VALUE>& operator=(const ExpNode<VALUE> &);
	~ExpNode() {
	}

	void reconstract(int, const VALUE&);

	bool operator==(const ExpNode<VALUE> &) const;
	bool operator!=(const ExpNode<VALUE> &) const;
	bool operator<(const ExpNode<VALUE> &) const;
	bool operator>(const ExpNode<VALUE> &) const;
	bool isConstNode() const;
	ExpNode<VALUE> operator+(const ExpNode<VALUE>&);
	ExpNode<VALUE> operator-();
	ExpNode<VALUE> operator-(const ExpNode<VALUE>&);
	void show() const;
};

template<typename VALUE>
ExpNode<VALUE>::ExpNode() {
	idx = IDX_CONST;
	val = 0.0;
}
template<typename VALUE>
ExpNode<VALUE>::ExpNode(int i, const VALUE& v) {
	idx = i;
	val = v;
}
template<typename VALUE>
ExpNode<VALUE>::ExpNode(VALUE v) {
	idx = IDX_CONST;
	val = v;
}
template<typename VALUE>
ExpNode<VALUE>::ExpNode(const ExpNode<VALUE>& en) {
	idx = en.idx;
	val = en.val;
}
template<typename VALUE>
ExpNode<VALUE>& ExpNode<VALUE>::operator=(const ExpNode<VALUE>& en) {
	idx = en.idx;
	val = en.val;
	return *this;
}
template<typename VALUE>
void ExpNode<VALUE>::reconstract(int i, const VALUE& v) {
	idx = i;
	val = v;
}
template<typename VALUE>
bool ExpNode<VALUE>::operator==(const ExpNode<VALUE> & en) const {
	return en.idx == this->idx ? true : false;
}
template<typename VALUE>
bool ExpNode<VALUE>::operator!=(const ExpNode<VALUE> & en) const {
	return en.idx != this->idx ? true : false;
}
template<typename VALUE>
bool ExpNode<VALUE>::operator<(const ExpNode<VALUE> & en) const {
	return this->idx < en.idx ? true : false;
}
template<typename VALUE>
bool ExpNode<VALUE>::operator>(const ExpNode<VALUE> & en) const {
	return this->idx > en.idx ? true : false;
}
template<typename VALUE>
bool ExpNode<VALUE>::isConstNode() const{
	return this->idx == IDX_CONST?true:false;
}
template<typename VALUE>
ExpNode<VALUE> ExpNode<VALUE>::operator+(const ExpNode<VALUE> & en) {
	ASSERT(this->idx == en.idx);
	this->val += en.val;
	return *this;
}
template<typename VALUE>
ExpNode<VALUE> ExpNode<VALUE>::operator-() {
	this->val = -this->val;
	return *this;
}

template<typename VALUE>
ExpNode<VALUE> ExpNode<VALUE>::operator-(const ExpNode<VALUE> & en) {
	ASSERT(this->idx == en.idx);
	this->val -= en.val;
	return *this;
}
template<typename VALUE>
void ExpNode<VALUE>::show() const{
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

//end of ExpNode class ==========================
//===============================================
//===============================================

typedef ExpNode<Float> ExpTerm;

class Expression: public BSTree<ExpTerm> {
public:
	typedef typename BSTree<ExpTerm>::reference reference;
	typedef typename BSTree<ExpTerm>::const_reference const_reference;
	typedef typename BSTree<ExpTerm>::pointer pointer;
	typedef typename BSTree<ExpTerm>::const_pointer const_pointer;
	typedef typename BSTree<ExpTerm>::iterator iterator;
	typedef typename BSTree<ExpTerm>::const_iterator const_iterator;

	typedef typename BSTree<ExpTerm>::difference_type difference_type;
	typedef typename BSTree<ExpTerm>::size_type size_type;

protected:
	void _insert(BinaryTreeNode<ExpTerm>*& Current, const ExpTerm &data);
public:
	Expression() :
			BSTree<ExpTerm>() {
	}
	//Expression operator=(const Expression &);
	void Insert(const ExpTerm &x);    //override

	//bool operator==(const Expression&) const;
	Expression operator+(const Expression &);
	Expression operator-(const Expression &);
	Expression operator*(const Float &);
	Expression operator/(const Float &);

	void plus(const Expression &);
	void minus(const Expression &);
	void times(const Float &);
    void divide(const Float &);

	void trim_zero();  //can be improved
	int  count_zero() const;

	void show() const;
};
Expression operator-(const Expression &, const Expression &);
Expression operator+(const Expression &, const Expression &);
Expression operator*(const Float&, const Expression&);

}

#endif /* ALGEBRA_EXPRESSION_H_ */
