/************************
 //  \file   ArrayT.h
 //  \brief
 // 
 //  \author zhou
 //  \date   23 janv. 2014 
 ***********************/
#ifndef _ARRAYLIST_H_
#define _ARRAYLIST_H_

#pragma once

#include "../TypeDef.h"

#include <stddef.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>

namespace Larus
{

template<typename T>
inline void _add_equal(
//+=
		LarusDef::size_type n, const T* src, T* dst)
{
	for (LarusDef::size_type i = 0; i < n; i++) {
		dst[i] += src[i];
	}
}

template<typename T>
inline int _copy_(LarusDef::size_type n, const T * a, T* b)
{
	//make sure: the length of a and b is equal
	//           and n>0
	typedef LarusDef::size_type size_type;
	size_type LN = 7;
	size_type m = (n - 1) % LN;
	for (size_type i = 0; i <= m; ++i) {
		b[i] = a[i];
	}
	if (n <= LN) {
		return 1;
	}
	size_type mp1 = m + 1;
	for (size_type i = mp1; i < n; i += LN) {
		b[i] = a[i];
		b[i + 1] = a[i + 1];
		b[i + 2] = a[i + 2];
		b[i + 3] = a[i + 3];
		b[i + 4] = a[i + 4];
		b[i + 5] = a[i + 5];
		b[i + 6] = a[i + 6];
	}
	return 2;
}

/*! \class Array
 \brief A class designed for array.
 It's the basic data structure class.
 */
template<typename TYPE>
class arrayListT: public ObjectBase
{
protected:
	LarusDef::size_type m_Len;
	TYPE *m_p;    //!< The heap for the real data.

public:
	// type definitions
	typedef TYPE value_type;
	typedef TYPE* iterator;
	typedef const TYPE* const_iterator;
	typedef TYPE& reference;
	typedef const TYPE& const_reference;
	typedef LarusDef::size_type size_type;
	typedef LarusDef::size_type difference_type;

	//constructor==================================
	arrayListT();
	arrayListT(const arrayListT<TYPE>& a);
	arrayListT(size_type Len);
	void reconstruct(size_type Len);
	arrayListT(size_type Len, const TYPE& nd);
	arrayListT(TYPE *nd, size_type Len);
	//=============================================
	~arrayListT();
	//=============================================
	TYPE* getPointer()
	{
		return m_p;
	}
	const TYPE* getPointer() const
	{
		return m_p;
	}
	//=============================================
	const_reference operator[](size_type index) const;  //overload []
	reference operator[](size_type index);
	const_reference operator()(size_type index) const;
	reference operator()(size_type index);
	const_reference at(size_type index) const;
	reference at(size_type index);
	//operator=====================================
	arrayListT<TYPE>& operator=(const arrayListT<TYPE> &a);
	// iterator support============================
	iterator begin()
	{
		return m_p;
	}
	const_iterator begin() const
	{
		return m_p;
	}
	iterator end()
	{
		return m_p + m_Len;
	}
	const_iterator end() const
	{
		return m_p + m_Len;
	}

	size_type size() const;  //!< get length of the array
	TYPE get(size_type i) const;
	void set(size_type i, const TYPE& value);
	void assign(const TYPE& nd);

	// front() and back()
	reference front();
	const_reference front() const;
	reference back();
	const_reference back() const;

	void swap(size_type i1, size_type i2);
	void reverse();
	void push_back(const TYPE& nd);
	void pop_back();
	void erase(size_type i);
	void resize(size_type new_size);

	bool empty() const;
	bool non_Empty() const;
	inline bool check_idx(size_type) const;

	bool has(const TYPE& nd) const;       //overload ==
	size_type find(const TYPE& nd) const; //overload ==
	size_type count_equal(const TYPE& nd) const; //overload ==
};

template<typename TYPE>
arrayListT<TYPE>::arrayListT()
{
	m_Len = 0;
	m_p = NULL_PTR;
}

template<typename TYPE>
arrayListT<TYPE>::arrayListT(const arrayListT<TYPE>& a)
{
	m_Len = a.size();
	m_p = new TYPE[m_Len];
	//unrolled loop
	_copy_(m_Len, a.getPointer(), m_p);
}

template<typename TYPE>
arrayListT<TYPE>::arrayListT(size_type Len)
{
	m_Len = Len;
	m_p = new TYPE[m_Len];
}

template<typename TYPE>
void arrayListT<TYPE>::reconstruct(size_type Len)
{
	assert(Len > 0);
	if (NULL != m_p)
		delete[] m_p;
	m_Len = Len;
	m_p = new TYPE[m_Len];
}

template<typename TYPE>
arrayListT<TYPE>::arrayListT(size_type Len, const TYPE& nd)
{
	m_Len = Len;
	m_p = new TYPE[Len];
	for (size_type i = 0; i < m_Len; i++) {
		m_p[i] = nd;
	}
}

template<typename TYPE>
arrayListT<TYPE>::arrayListT(TYPE *nd, size_type Len)
{
	m_Len = Len;
	m_p = new TYPE[Len];
	for (size_type i = 0; i < m_Len; i++) {
		m_p[i] = nd[i];
	}
}

template<typename TYPE>
arrayListT<TYPE>::~arrayListT()
{
	delete[] m_p;
}

template<typename TYPE>
const TYPE& arrayListT<TYPE>::operator[](size_type index) const
{
	assert(index >= 0 && index < m_Len);
	return m_p[index];
}

template<typename TYPE>
TYPE& arrayListT<TYPE>::operator[](size_type index)
{
	assert(index >= 0 && index < m_Len);
	return m_p[index];
}

template<typename TYPE>
const TYPE& arrayListT<TYPE>::operator()(size_type index) const
{
	assert(index >= 0 && index < m_Len);
	return m_p[index];
}

template<typename TYPE>
TYPE& arrayListT<TYPE>::operator()(size_type index)
{
	assert(index >= 0 && index < m_Len);
	return m_p[index];
}
template<typename TYPE>
const TYPE& arrayListT<TYPE>::at(size_type index) const
{
	assert(index >= 0 && index < m_Len);
	return m_p[index];
}
template<typename TYPE>
TYPE& arrayListT<TYPE>::at(size_type index)
{
	assert(index >= 0 && index < m_Len);
	return m_p[index];
}

template<typename TYPE>
arrayListT<TYPE>& arrayListT<TYPE>::operator=(const arrayListT<TYPE> &a)
{
	if (this == &a) {
		return *this;
	}
	if (m_Len == a.size()) {
		//unrolled loop
		_copy_(m_Len, a.getPointer(), m_p);
	} else {
		delete[] this->m_p;
		m_Len = a.size();
		m_p = new TYPE[m_Len];
		_copy_(m_Len, a.getPointer(), m_p);
	}
	return *this;
}

template<typename TYPE>
void arrayListT<TYPE>::set(size_type i, const TYPE& value)
{
	assert(i >= 0 && i < m_Len);
	m_p[i] = value;
}

template<typename TYPE>
void arrayListT<TYPE>::assign(const TYPE& nd)
{
	for (size_type i = 0; i < m_Len; i++) {
		m_p[i] = nd;
	}
}
template<typename TYPE>
LarusDef::size_type arrayListT<TYPE>::size() const
{
	return m_Len;
}
template<typename TYPE>
TYPE& arrayListT<TYPE>::front()
{
	ASSERT(m_Len > 0);
	return m_p[0];
}
template<typename TYPE>
const TYPE& arrayListT<TYPE>::front() const
{
	ASSERT(m_Len > 0);
	return m_p[0];
}
template<typename TYPE>
TYPE& arrayListT<TYPE>::back()
{
	ASSERT(m_Len > 0);
	return m_p[m_Len - 1];
}
template<typename TYPE>
const TYPE& arrayListT<TYPE>::back() const
{
	ASSERT(m_Len > 0);
	return m_p[m_Len - 1];
}

template<typename TYPE>
TYPE arrayListT<TYPE>::get(size_type i) const
{
	assert(i >= 0 && i < m_Len);
	return m_p[i];
}

template<typename TYPE>
void arrayListT<TYPE>::swap(size_type i1, size_type i2)
{
	assert(i1 >= 0 && i1 < m_Len && i2 >= 0 && i2 < m_Len);
	if (i1 == i2) {
		return;
	}
	TYPE tmp;
	tmp = m_p[i1];
	m_p[i1] = m_p[i2];
	m_p[i2] = tmp;
}

template<typename TYPE>
void arrayListT<TYPE>::reverse()
{
	if (empty()) {
		return;
	}
	TYPE tmp;
	for (int i1 = 0, i2 = size() - 1; i1 < i2; i1++, i2--) {
		tmp = m_p[i1];
		m_p[i1] = m_p[i2];
		m_p[i2] = tmp;
	}
}

template<typename TYPE>
void arrayListT<TYPE>::push_back(const TYPE& nd)
{
	if (m_p == NULL_PTR) {
		m_Len = 1;
		m_p = new TYPE[1];
		m_p[0] = nd;
	} else {
		TYPE *tmp = new TYPE[m_Len];
		for (size_type i = 0; i < m_Len; i++) {
			tmp[i] = m_p[i];
		}
		delete[] m_p;
		m_Len += 1;
		m_p = new TYPE[m_Len];
		for (size_type i = 0; i < m_Len - 1; i++) {
			m_p[i] = tmp[i];
		}
		m_p[m_Len - 1] = nd;
		delete[] tmp;
	}
}

template<typename TYPE>
void arrayListT<TYPE>::pop_back()
{
	if (m_p == NULL) {
		return;
	} else if (m_Len == 1) {
		m_Len = 0;
		delete[] m_p;
		m_p = NULL;
	} else {
		TYPE *tmp = new TYPE[m_Len];
		for (size_type i = 0; i < m_Len; i++) {
			tmp[i] = m_p[i];
		}
		delete[] m_p;
		m_Len -= 1;
		m_p = new TYPE[m_Len];
		for (size_type i = 0; i < m_Len; i++) {
			m_p[i] = tmp[i];
		}
		delete[] tmp;
	}
}

template<typename TYPE>
void arrayListT<TYPE>::erase(size_type idx)
{
	assert(idx >= 0 && idx < m_Len);
	if (m_Len == 1) {
		m_Len = 0;
		delete[] m_p;
		m_p = NULL;
	} else {
		TYPE *tmp = new TYPE[m_Len];
		for (size_type i = 0; i < m_Len; i++) {
			if (i != idx)
				tmp[i] = m_p[i];
		}
		delete[] m_p;
		m_Len--;
		m_p = new TYPE[m_Len];
		for (size_type i = 0; i < idx; i++) {
			m_p[i] = tmp[i];
		}
		for (size_type i = idx; i < m_Len; i++) {
			m_p[i] = tmp[i + 1];
		}
		delete[] tmp;
	}
}

template<typename TYPE>
void arrayListT<TYPE>::resize(size_type new_size)
{
	assert(new_size >= 0);
	if (empty()) {
		reconstruct(new_size);
	} else {
		TYPE *tmp = new TYPE[new_size];
		for (size_type i = 0; i < m_Len && i < new_size; i++) {
			tmp[i] = m_p[i];
		}
		delete[] m_p;
		m_Len = new_size;
		m_p = tmp;
	}
}

template<typename TYPE>
bool arrayListT<TYPE>::has(const TYPE& nd) const
{
	if (m_p == NULL) {
		return false;
	} else {
		for (size_type i = 0; i < m_Len; i++) {
			if (nd == m_p[i]) {
				return true;
			}
		}
		return false;
	}
}

template<typename TYPE>
LarusDef::size_type arrayListT<TYPE>::find(const TYPE& nd) const
{
	if (m_p == NULL) {
		return -1;
	} else {
		for (size_type i = 0; i < m_Len; i++) {
			if (nd == m_p[i]) {
				return i;
			}
		}
		return -1;
	}
}
template<typename TYPE>
LarusDef::size_type arrayListT<TYPE>::count_equal(const TYPE& nd) const
{ //overload ==
	size_type res = 0;
	if (m_p == NULL) {
		return 0;
	} else {
		for (size_type i = 0; i < m_Len; i++) {
			if (nd == m_p[i]) {
				res++;
			}
		}
	}
	return res;
}
template<typename TYPE>
bool arrayListT<TYPE>::empty() const
{
	if (m_p == NULL && m_Len == 0) {
		return true;
	} else {
		return false;
	}
}

template<typename TYPE>
bool arrayListT<TYPE>::non_Empty() const
{
	return !empty();
}
template<typename TYPE>
inline bool arrayListT<TYPE>::check_idx(arrayListT<TYPE>::size_type idx) const
{
	if (idx >= 0 && idx < m_Len) {
		return true;
	} else {
		return false;
	}

}
//end of Class arrayListT
//=========================================================
//=========================================================
//Function out of class====================================

//=========================================================

template<typename TYPE>
class arrayListV: public arrayListT<TYPE>
{
public:
	// type definitions
	typedef TYPE value_type;
	typedef TYPE* pointer;
	typedef const TYPE* const_pointer;
	typedef TYPE& reference;
	typedef const TYPE& const_reference;
	typedef LarusDef::size_type size_type;
	typedef LarusDef::size_type difference_type;
	//constructor==================================
	arrayListV();
	arrayListV(size_type Len);
	arrayListV(size_type Len, const TYPE& nd);
	//arrayListV(TYPE *nd, size_type Len);
	//opertator====================================
	arrayListV<TYPE> operator+(const arrayListV<TYPE> &a);
	arrayListV<TYPE> operator+(const TYPE &a);
	arrayListV<TYPE> operator-(const arrayListV<TYPE> &a);
	arrayListV<TYPE> operator-(const TYPE &a);
	arrayListV<TYPE> operator*(const arrayListV<TYPE> &a);
	arrayListV<TYPE> operator*(const TYPE &a);
	arrayListV<TYPE> operator/(const arrayListV<TYPE> &a);
	arrayListV<TYPE> operator/(const TYPE &a);

	//other functions==============================
	TYPE sum() const;
	TYPE findMin() const;
	TYPE findMax() const;
	size_type findMinIdx() const;
	size_type findMaxIdx() const;
	//fill ----------------------------------------
	void assign_forward(const TYPE& be, const TYPE& step);
	void assign_backward(const TYPE& be, const TYPE& step);
	//count ---------------------------------------
	size_type countEq(const TYPE &a) const;

	void show() const;
};
template<typename TYPE>
arrayListV<TYPE>& operator*=(arrayListV<TYPE> &x, const TYPE &a);
template<typename TYPE>
arrayListV<TYPE>& operator+=(arrayListV<TYPE> &x, const arrayListV<TYPE> &y);
template<typename TYPE>
arrayListV<TYPE>& operator-=(arrayListV<TYPE> &x, const arrayListV<TYPE> &y);
template<typename TYPE>
arrayListV<TYPE> operator*(const TYPE &a, const arrayListV<TYPE> &x);
template<typename TYPE>
arrayListV<TYPE> operator*(const arrayListV<TYPE> &x, const TYPE &a);
template<typename TYPE>
arrayListV<TYPE> operator-(const arrayListV<TYPE> &a,
		const arrayListV<TYPE> &b);
//=================================================

template<typename TYPE>
arrayListV<TYPE>::arrayListV() :
		arrayListT<TYPE>()
{
}
template<typename TYPE>
arrayListV<TYPE>::arrayListV(size_type Len) :
		arrayListT<TYPE>(Len)
{
	this->assign(0);
}
template<typename TYPE>
arrayListV<TYPE>::arrayListV(size_type Len, const TYPE& nd) :
		arrayListT<TYPE>(Len, nd)
{
}
template<typename TYPE>
arrayListV<TYPE> arrayListV<TYPE>::operator+(const arrayListV<TYPE> &a)
{
	assert(a.size() == this->size());
	arrayListV<TYPE> res(this->m_Len);
	_copy_(a.size(), this->m_p, res.m_p);
	_add_equal(a.size(), a.m_p, res.m_p);
	return res;
}
template<typename TYPE>
arrayListV<TYPE> arrayListV<TYPE>::operator-(const arrayListV<TYPE> &a)
{
	assert(a.size() == this->size());
	arrayListV<TYPE> tmp(this->m_Len);
	for (size_type i = 0; i < this->m_Len; i++) {
		tmp[i] = this->m_p[i] - a[i];
	}
	return tmp;
}
template<typename TYPE>
arrayListV<TYPE> arrayListV<TYPE>::operator*(const arrayListV<TYPE> &a)
{
	assert(a.size() == this->size());
	arrayListV<TYPE> tmp(this->m_Len);
	for (size_type i = 0; i < this->m_Len; i++) {
		tmp[i] = this->m_p[i] * a[i];
	}
	return tmp;
}
template<typename TYPE>
arrayListV<TYPE> arrayListV<TYPE>::operator/(const arrayListV<TYPE> &a)
{
	assert(a.size() == this->size());
	arrayListV<TYPE> tmp(this->m_Len);
	for (size_type i = 0; i < this->m_Len; i++) {
		tmp[i] = this->m_p[i] / a[i];
	}
	return tmp;
}
template<typename TYPE>
arrayListV<TYPE> arrayListV<TYPE>::operator+(const TYPE &a)
{
	arrayListV<TYPE> sum(this->m_Len);
	for (size_type i = 0; i < this->m_Len; i++) {
		sum[i] = this->m_p[i] + a;
	}
	return sum;
}
template<typename TYPE>
arrayListV<TYPE> arrayListV<TYPE>::operator-(const TYPE &a)
{
	arrayListV<TYPE> sum(this->m_Len);
	for (size_type i = 0; i < this->m_Len; i++) {
		sum[i] = this->m_p[i] - a;
	}
	return sum;
}
template<typename TYPE>
arrayListV<TYPE> arrayListV<TYPE>::operator*(const TYPE &a)
{
	arrayListV<TYPE> sum(this->m_Len);
	for (size_type i = 0; i < this->m_Len; i++) {
		sum[i] = this->m_p[i] * a;
	}
	return sum;
}
template<typename TYPE>
arrayListV<TYPE> arrayListV<TYPE>::operator/(const TYPE &a)
{
	arrayListV<TYPE> sum(this->m_Len);
	for (size_type i = 0; i < this->m_Len; i++) {
		sum[i] = this->m_p[i] / a;
	}
	return sum;
}
template<typename TYPE>
TYPE arrayListV<TYPE>::sum() const
{
	TYPE sum = 0;
	for (size_type i = 0; i < this->m_Len; i++) {
		sum += this->m_p[i];
	}
	return sum;
}
template<typename TYPE>
TYPE arrayListV<TYPE>::findMin() const
{
	assert(this->m_p!=NULL);
	TYPE min = this->m_p[0];
	for (size_type i = 0; i < this->m_Len; i++) {
		if (this->m_p[i] < min) {
			min = this->m_p[i];
		}
	}
	return min;
}
template<typename TYPE>
TYPE arrayListV<TYPE>::findMax() const
{
	assert(this->m_p!=NULL);
	TYPE max = this->m_p[0];
	for (size_type i = 0; i < this->m_Len; i++) {
		if (this->m_p[i] > max) {
			max = this->m_p[i];
		}
	}
	return max;
}
template<typename TYPE>
LarusDef::size_type arrayListV<TYPE>::findMinIdx() const
{
	assert(this->m_p!=NULL);
	TYPE min = this->m_p[0];
	size_type idx = 0;
	for (size_type i = 0; i < this->m_Len; i++) {
		if (this->m_p[i] < min) {
			min = this->m_p[i];
			idx = i;
		}
	}
	return idx;
}
template<typename TYPE>
LarusDef::size_type arrayListV<TYPE>::findMaxIdx() const
{
	assert(this->m_p!=NULL);
	TYPE max = this->m_p[0];
	size_type idx = 0;
	for (size_type i = 0; i < this->m_Len; i++) {
		if (this->m_p[i] > max) {
			max = this->m_p[i];
			idx = i;
		}
	}
	return idx;
}
template<typename TYPE>
void arrayListV<TYPE>::assign_forward(const TYPE& be, const TYPE& step)
{
	for (size_type i = 0; i < this->m_Len; i++) {
		this->m_p[i] = be + i * step;
	}
}
template<typename TYPE>
void arrayListV<TYPE>::assign_backward(const TYPE& be, const TYPE& step)
{
	for (size_type i = this->m_Len - 1; i >= 0; i++) {
		this->m_p[i] = be + (this->m_Len - 1 - i) * step;
	}
}

template<typename TYPE>
LarusDef::size_type arrayListV<TYPE>::countEq(const TYPE &a) const
{
	if (this->m_p != NULL) {
		return 0;
	}
	size_type count = 0;
	for (size_type i = 0; i < this->m_Len; i++) {
		if (this->m_p[i] == a) {
			count++;
		}
	}
	return count;
}
template<typename TYPE>
void arrayListV<TYPE>::show() const
{
	std::cout << "Size = " << this->m_Len << "\n";
	for (int i = 0; i < this->m_Len; i++) {
		std::cout.width(5);
		std::cout << i <<" ";
		std::cout.width(12);
		std::cout.precision(8);
		std::cout << this->m_p[i] << "\n";
	}
}
//=========================================================
template<typename TYPE>
arrayListV<TYPE>& operator*=(arrayListV<TYPE> &x, const TYPE &a)
{
	for (idx_t i = 0; i < x.size(); ++i) {
		x[i] *= a;
	}
	return x;
}
template<typename TYPE>
arrayListV<TYPE>& operator+=(arrayListV<TYPE> &x, const arrayListV<TYPE> &y)
{
	assert(x.size() == y.size());
	for (idx_t i = 0; i < x.size(); i++) {
		x[i] += y[i];
	}
	return x;
}
template<typename TYPE>
arrayListV<TYPE>& operator-=(arrayListV<TYPE> &x, const arrayListV<TYPE> &y)
{
	assert(x.size() == y.size());
	for (idx_t i = 0; i < x.size(); i++) {
		x[i] -= y[i];
	}
	return x;
}
template<typename TYPE>
arrayListV<TYPE> operator*(const TYPE& a, const arrayListV<TYPE>& x)
{
	arrayListV<TYPE> res(x);
	res *= a;
	return res;
}
template<typename TYPE>
arrayListV<TYPE> operator*(const arrayListV<TYPE>& x, const TYPE&a)
{
	arrayListV<TYPE> res(x);
	res *= a;
	return res;
}

template<typename TYPE>
arrayListV<TYPE> operator-(const arrayListV<TYPE> &a, const arrayListV<TYPE> &b)
{
	arrayListV<TYPE> res(a);
	res -= b;
	return res;
}

typedef arrayListV<Float> arrayList;
typedef arrayListV<LarusDef::size_type> arrayList_st;
typedef arrayListV<int> arrayList_int;
typedef arrayListT<bool> arrayList_bool;
typedef arrayListT<unsigned long int> arrayList_ul;

}

#endif /* ARRAYT_H_ */
