/************************
 //  \file   MatrixT.h
 //  \brief
 // 
 //  \author zhou
 //  \date   23 janv. 2014 
 ***********************/
#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
#include <assert.h>
#include "../TypeDef.h"
#include "../Utility/ArrayList.h"

namespace Larus {

template<typename TYPE>
class MatrixT {
public:
	// type definitions===================
	typedef TYPE value_type;
	typedef TYPE* pointer;
	typedef const TYPE* const_pointer;
	typedef TYPE& reference;
	typedef const TYPE& const_reference;
	typedef LarusDef::size_type size_type;
	typedef LarusDef::size_type difference_type;
protected:
	LarusDef::size_type m_iLen;
	LarusDef::size_type m_jLen;
	arrayListT<TYPE> *m_mp;
public:
	//constructor==========================
	MatrixT();
	MatrixT(const MatrixT<TYPE>& a);
	MatrixT(size_type iLen, size_type jLen);
	MatrixT(size_type iLen, size_type jLen, TYPE **value);
	void reconstruct(size_type iLen, size_type jLen);
	//=============================================
	MatrixT<TYPE>& operator=(const MatrixT<TYPE> &a);
	//=============================================
	~MatrixT();
	//Capacity=====================================
	size_type size() const;
	size_type iLen() const;
	size_type jLen() const;
	bool empty() const;
	//Element access===============================
	arrayListT<TYPE>& operator[](size_type index);
	const arrayListT<TYPE>& operator[](size_type index) const;
	reference operator()(size_type i, size_type j);
	const_reference operator()(size_type i, size_type j) const;
	reference at(size_type i, size_type j);
	const_reference at(size_type i, size_type j) const;
	TYPE get(size_type i, size_type j);
	TYPE* getpValue(size_type i, size_type j);
	void set(size_type i, size_type j, const TYPE& value);
	void set_row(size_type i, const TYPE& value);
	void set_col(size_type i, const TYPE& value);
	void assign(const TYPE& value);
	//=============================================

	void swapValue(size_type i1, size_type j1, size_type i2, size_type j2);

	size_type countNumEq(TYPE value); //overload ==

	inline bool testIdxI(size_type);
	inline bool testIdxJ(size_type);

	//Low efficient function for small matrix
	void appendRow(const arrayListT<TYPE> &);
	void appendCol(const arrayListT<TYPE> &);
	void insertRowAfter(size_type, const arrayListT<TYPE>&);
	void insertColAfter(size_type, const arrayListT<TYPE>&);
	void deleteRow(size_type);
	void deleteCol(size_type);
};

template<typename TYPE>
inline bool MatrixT<TYPE>::testIdxI(LarusDef::size_type i) {
	return (i >= 0 && i < m_iLen) ? true : false;
}

template<typename TYPE>
inline bool MatrixT<TYPE>::testIdxJ(LarusDef::size_type j) {
	return (j >= 0 && j < m_jLen) ? true : false;
}

template<typename TYPE>
MatrixT<TYPE>::MatrixT() {
	m_iLen = 0;
	m_jLen = 0;
	m_mp = NULL_PTR;
}

template<typename TYPE>
MatrixT<TYPE>::MatrixT(const MatrixT<TYPE>& a) {
	m_iLen = a.iLen();
	m_jLen = a.jLen();
	m_mp = new arrayListT<TYPE> [m_iLen];
	for (size_type i = 0; i < m_iLen; i++) {
		m_mp[i].reconstruct(m_jLen);
	}
	for (size_type i = 0; i < m_iLen; i++) {
		for (size_type j = 0; j < m_jLen; j++) {
			m_mp[i][j] = a[i][j];
		}
	}
}

template<typename TYPE>
MatrixT<TYPE>::MatrixT(LarusDef::size_type iLen, LarusDef::size_type jLen) {
	m_iLen = iLen;
	m_jLen = jLen;
	//m_total = m_iLen * m_jLen;
	m_mp = new arrayListT<TYPE> [m_iLen];
	for (size_type i = 0; i < m_iLen; i++) {
		m_mp[i].reconstruct(m_jLen);
	}
}

template<typename TYPE>
MatrixT<TYPE>::MatrixT(LarusDef::size_type iLen, LarusDef::size_type jLen,
		TYPE **value) {
	m_iLen = iLen;
	m_jLen = jLen;
	//m_total = m_iLen * m_jLen;
	m_mp = new arrayListT<TYPE> [m_iLen];
	for (size_type i = 0; i < m_iLen; i++) {
		m_mp[i].reconstruct(m_jLen);
	}
	for (size_type i = 0; i < m_iLen; i++) {
		for (size_type j = 0; j < m_jLen; j++) {
			m_mp[i][j] = value[i][j];
		}
	}
}
template<typename TYPE>
void MatrixT<TYPE>::reconstruct(LarusDef::size_type iLen,
		LarusDef::size_type jLen) {
	m_iLen = iLen;
	m_jLen = jLen;
	//m_total = m_iLen * m_jLen;
	if (m_mp != NULL_PTR) {
		delete[] m_mp;
	}
	m_mp = new arrayListT<TYPE> [m_iLen];
	for (size_type i = 0; i < m_iLen; i++) {
		m_mp[i].reconstruct(m_jLen);
	}
}
template<typename TYPE>
MatrixT<TYPE>::~MatrixT() {
	delete[] m_mp;
}
template<typename TYPE>
const arrayListT<TYPE>& MatrixT<TYPE>::operator[](
		LarusDef::size_type index) const {
	assert(index >= 0 && index < m_iLen);
	return m_mp[index];
}
template<typename TYPE>
arrayListT<TYPE>& MatrixT<TYPE>::operator[](LarusDef::size_type index) {
	assert(index >= 0 && index < m_iLen);
	return m_mp[index];
}
template<typename TYPE>
bool MatrixT<TYPE>::empty() const {
	return m_mp == NULL_PTR;
}

template<typename TYPE>
MatrixT<TYPE>& MatrixT<TYPE>::operator=(const MatrixT<TYPE> &a) {
	if (this == &a) {
		return *this;
	}
	if (m_iLen == a.iLen() && m_jLen == a.jLen()) {
		for (size_type i = 0; i < m_iLen; i++) {
			for (size_type j = 0; j < m_jLen; j++) {
				this->m_mp[i][j] = a[i][j];
			}
		}
	} else {
		delete[] this->m_mp;
		m_iLen = a.iLen();
		m_jLen = a.jLen();
		//m_total = m_iLen * m_jLen;
		m_mp = new arrayListT<TYPE> [m_iLen];
		for (size_type i = 0; i < m_iLen; i++) {
			m_mp[i].reconstruct(m_jLen);
		}
		for (size_type i = 0; i < m_iLen; i++) {
			for (size_type j = 0; j < m_jLen; j++) {
				this->m_mp[i][j] = a[i][j];
			}
		}
	}
	return *this;
}

template<typename TYPE>
LarusDef::size_type MatrixT<TYPE>::size() const {
	return m_iLen * m_jLen;
}
template<typename TYPE>
LarusDef::size_type MatrixT<TYPE>::iLen() const {
	return m_iLen;
}
template<typename TYPE>
LarusDef::size_type MatrixT<TYPE>::jLen() const {
	return m_jLen;
}
template<typename TYPE>
TYPE MatrixT<TYPE>::get(LarusDef::size_type i, LarusDef::size_type j) {
	assert(i >= 0 && i < m_iLen);
	assert(j >= 0 && j < m_jLen);
	return m_mp[i][j];
}
template<typename TYPE>
typename MatrixT<TYPE>::reference MatrixT<TYPE>::operator()(size_type i,
		size_type j) {
	assert(i >= 0 && i < m_iLen);
	assert(j >= 0 && j < m_jLen);
	return m_mp[i][j];
}
template<typename TYPE>
typename MatrixT<TYPE>::const_reference MatrixT<TYPE>::operator()(size_type i,
		size_type j) const {
	assert(i >= 0 && i < m_iLen);
	assert(j >= 0 && j < m_jLen);
	return m_mp[i][j];
}

template<typename TYPE>
TYPE& MatrixT<TYPE>::at(LarusDef::size_type i, LarusDef::size_type j) {
	assert(i >= 0 && i < m_iLen);
	assert(j >= 0 && j < m_jLen);
	return m_mp[i][j];
}

template<typename TYPE>
const TYPE& MatrixT<TYPE>::at(LarusDef::size_type i,
		LarusDef::size_type j) const {
	assert(i >= 0 && i < m_iLen);
	assert(j >= 0 && j < m_jLen);
	return m_mp[i][j];
}

template<typename TYPE>
TYPE* MatrixT<TYPE>::getpValue(LarusDef::size_type i, LarusDef::size_type j) {
	assert(i >= 0 && i < m_iLen);
	assert(j >= 0 && j < m_jLen);
	return &m_mp[i][j];
}

template<typename TYPE>
void MatrixT<TYPE>::set(LarusDef::size_type i, LarusDef::size_type j,
		const TYPE& value) {
	assert(i >= 0 && i < m_iLen);
	assert(j >= 0 && j < m_jLen);
	m_mp[i][j] = value;
}

template<typename TYPE>
void MatrixT<TYPE>::set_row(LarusDef::size_type i, const TYPE& value) {
	assert(i >= 0 && i < m_iLen);
	for (size_type j = 0; j < m_jLen; j++) {
		m_mp[i][j] = value;
	}
}

template<typename TYPE>
void MatrixT<TYPE>::set_col(LarusDef::size_type j, const TYPE& value) {
	assert(j >= 0 && j < m_jLen);
	for (size_type i = 0; i < m_iLen; i++) {
		m_mp[i][j] = value;
	}
}

template<typename TYPE>
void MatrixT<TYPE>::assign(const TYPE& value) {
	for (size_type i = 0; i < m_iLen; i++) {
		for (size_type j = 0; j < m_jLen; j++) {
			m_mp[i][j] = value;
		}
	}
}

template<typename TYPE>
void MatrixT<TYPE>::swapValue(LarusDef::size_type i1, LarusDef::size_type j1,
		LarusDef::size_type i2, LarusDef::size_type j2) {
	assert(i1 >= 0 && i1 < m_iLen && i2 >= 0 && i2 < m_iLen);
	assert(j1 >= 0 && j1 < m_jLen && j2 >= 0 && j2 < m_jLen);
	if (i1 == i2 && j1 == j2) {
		return;
	}
	TYPE tmp;
	tmp = m_mp[i1][j1];
	m_mp[i1][j1] = m_mp[i2][j2];
	m_mp[i2][j2] = tmp;
}

template<typename TYPE>
LarusDef::size_type MatrixT<TYPE>::countNumEq(TYPE value) {
	size_type num = 0;
	for (size_type i = 0; i < m_iLen; i++) {
		for (size_type j = 0; j < m_jLen; j++) {
			if (m_mp[i][j] == value) {
				num++;
			}
		}
	}
	return num;
}
template<typename TYPE>
void MatrixT<TYPE>::appendRow(const arrayListT<TYPE> &a) {
	assert(a.size()==m_jLen||m_mp==NULL);
	if (m_mp == NULL) {
		m_iLen = 1;
		m_jLen = a.size();
		m_mp = new arrayListT<TYPE> [m_iLen];
		for (size_type i = 0; i < m_iLen; i++) {
			m_mp[i].reconstruct(m_jLen);
		}
		for (size_type j = 0; j < m_jLen; j++) {
			m_mp[0][j] = a[j];
		}
	} else {
		m_iLen += 1;
		arrayListT<TYPE>* tmp = new arrayListT<TYPE> [m_iLen];
		for (size_type i = 0; i < m_iLen; i++) {
			tmp[i].reconstruct(m_jLen);
		}
		for (size_type i = 0; i < m_iLen - 1; i++) {
			for (size_type j = 0; j < m_jLen; j++) {
				tmp[i][j] = m_mp[i][j];
			}
		}
		for (size_type j = 0; j < m_jLen; j++) {
			tmp[m_iLen - 1][j] = a[j];
		}
		arrayListT<TYPE>* tmp2 = m_mp;
		m_mp = tmp;
		tmp = tmp2;
		delete[] tmp;
	}
}
template<typename TYPE>
void MatrixT<TYPE>::appendCol(const arrayListT<TYPE> &a) {
	assert(a.size()==m_iLen||m_mp==NULL);
	if (m_mp == NULL) {
		m_iLen = a.size();
		m_jLen = 1;
		m_mp = new arrayListT<TYPE> [m_iLen];
		for (size_type i = 0; i < m_iLen; i++) {
			m_mp[i].reconstruct(1);
		}
		for (size_type i = 0; i < m_iLen; i++) {
			m_mp[i][0] = a[i];
		}
	} else {
		m_jLen += 1;
		arrayListT<TYPE>* tmp = new arrayListT<TYPE> [m_iLen];
		for (size_type i = 0; i < m_iLen; i++) {
			tmp[i].reconstruct(m_jLen);
		}
		for (size_type i = 0; i < m_iLen; i++) {
			for (size_type j = 0; j < m_jLen - 1; j++) {
				tmp[i][j] = m_mp[i][j];
			}
		}
		for (size_type i = 0; i < m_iLen; i++) {
			tmp[i][m_jLen - 1] = a[i];
		}
		arrayListT<TYPE>* tmp2 = m_mp;
		m_mp = tmp;
		tmp = tmp2;
		delete[] tmp;
	}
}
template<typename TYPE>
void MatrixT<TYPE>::insertRowAfter(MatrixT<TYPE>::size_type ii,
		const arrayListT<TYPE>& a) {
	assert(ii >= 0 && ii < m_iLen);
	assert(a.Len() == m_jLen);
	m_iLen += 1;
	arrayListT<TYPE>* tmp = new arrayListT<TYPE> [m_iLen];
	for (size_type i = 0; i < m_iLen; i++) {
		tmp[i].reconstruct(m_jLen);
	}
	for (size_type i = 0; i <= ii; i++) {
		for (size_type j = 0; j < m_jLen; j++) {
			tmp[i][j] = m_mp[i][j];
		}
	}
	for (size_type j = 0; j < m_jLen; j++) {
		tmp[ii + 1][j] = a[j];
	}
	for (size_type i = ii + 2; i < m_iLen; i++) {
		for (size_type j = 0; j < m_jLen; j++) {
			tmp[i][j] = m_mp[i - 1][j];
		}
	}
	arrayListT<TYPE>* tmp2 = m_mp;
	m_mp = tmp;
	tmp = tmp2;
	delete[] tmp;
}
template<typename TYPE>
void MatrixT<TYPE>::insertColAfter(MatrixT<TYPE>::size_type jj,
		const arrayListT<TYPE>& a) {
	assert(jj >= 0 && jj < m_jLen);
	assert(a.Len() == m_iLen);
	m_jLen += 1;
	arrayListT<TYPE>* tmp = new arrayListT<TYPE> [m_iLen];
	for (size_type i = 0; i < m_iLen; i++) {
		tmp[i].reconstruct(m_jLen);
	}
	for (size_type i = 0; i < m_iLen; i++) {
		for (size_type j = 0; j <= jj; j++) {
			tmp[i][j] = m_mp[i][j];
		}
	}
	for (size_type i = 0; i < m_iLen; i++) {
		tmp[i][jj + 1] = a[i];
	}
	for (size_type i = 0; i < m_iLen; i++) {
		for (size_type j = jj + 2; j < m_jLen; j++) {
			tmp[i][j] = m_mp[i][j - 1];
		}
	}
	arrayListT<TYPE>* tmp2 = m_mp;
	m_mp = tmp;
	tmp = tmp2;
	delete[] tmp;
}
template<typename TYPE>
void MatrixT<TYPE>::deleteRow(MatrixT<TYPE>::size_type ii) {
	assert(ii >= 0 && ii < m_iLen);
	m_iLen -= 1;
	if (m_iLen == 0) {
		m_iLen = 0;
		m_jLen = 0;
		delete[] m_mp;
		m_mp = NULL;
	} else {
		arrayListT<TYPE>* tmp = new arrayListT<TYPE> [m_iLen];
		for (size_type i = 0; i < m_iLen; i++) {
			tmp[i].reconstruct(m_jLen);
		}
		for (size_type i = 0; i < ii; i++) {
			for (size_type j = 0; j < m_jLen; j++) {
				tmp[i][j] = m_mp[i][j];
			}
		}
		for (size_type i = ii + 1; i < m_iLen + 1; i++) {
			for (size_type j = 0; j < m_jLen; j++) {
				tmp[i - 1][j] = m_mp[i][j];
			}
		}
		arrayListT<TYPE>* tmp2 = m_mp;
		m_mp = tmp;
		tmp = tmp2;
		delete[] tmp;
	}
}
template<typename TYPE>
void MatrixT<TYPE>::deleteCol(MatrixT<TYPE>::size_type jj) {
	assert(jj >= 0 && jj < m_iLen);
	m_jLen -= 1;
	if (m_jLen == 0) {
		m_iLen = 0;
		m_jLen = 0;
		delete[] m_mp;
		m_mp = NULL;
	} else {
		arrayListT<TYPE>* tmp = new arrayListT<TYPE> [m_iLen];
		for (size_type i = 0; i < m_iLen; i++) {
			tmp[i].reconstruct(m_jLen);
		}
		for (size_type i = 0; i < m_iLen; i++) {
			for (size_type j = 0; j < jj; j++) {
				tmp[i][j] = m_mp[i][j];
			}
		}
		for (size_type i = 0; i < m_iLen; i++) {
			for (size_type j = jj + 1; j < m_jLen + 1; j++) {
				tmp[i][j - 1] = m_mp[i][j];
			}
		}
		arrayListT<TYPE>* tmp2 = m_mp;
		m_mp = tmp;
		tmp = tmp2;
		delete[] tmp;
	}
}

//end of class MatrixT===========================
//===============================================
//===============================================

template<typename TYPE>
class MatrixV: public MatrixT<TYPE> {
public:
	// type definitions===================
	typedef TYPE value_type;
	typedef TYPE* pointer;
	typedef const TYPE* const_pointer;
	typedef TYPE& reference;
	typedef const TYPE& const_reference;
	typedef LarusDef::size_type size_type;
	typedef LarusDef::size_type difference_type;
	//constructor==========================
	MatrixV();
	MatrixV(size_type iLen, size_type jLen);
	MatrixV(size_type iLen, size_type jLen, TYPE **value);
	//void reconstruct(size_type iLen, size_type jLen);
	//~MatrixV();
	//=============================================
	MatrixV<TYPE> operator+(const MatrixV<TYPE> &a);
	MatrixV<TYPE> operator-(const MatrixV<TYPE> &a);
	MatrixV<TYPE> operator*(const MatrixV<TYPE> &a);
	//show ========================================
	void show() const;

};

template<typename TYPE>
MatrixV<TYPE>::MatrixV() :
		MatrixT<TYPE>() {
}

template<typename TYPE>
MatrixV<TYPE>::MatrixV(LarusDef::size_type iLen, LarusDef::size_type jLen) :
		MatrixT<TYPE>(iLen, jLen) {
	this->assign(0);
}

template<typename TYPE>
MatrixV<TYPE>::MatrixV(LarusDef::size_type iLen, LarusDef::size_type jLen,
		TYPE **value) :
		MatrixT<TYPE>(iLen, jLen) {
	for (size_type i = 0; i < this->m_iLen; i++) {
		for (size_type j = 0; j < this->m_jLen; j++) {
			this->m_mp[i][j] = value[i][j];
		}
	}
}

template<typename TYPE>
MatrixV<TYPE> MatrixV<TYPE>::operator+(const MatrixV<TYPE> &a) {
	assert(a.iLen() == this->iLen());
	assert(a.jLen() == this->jLen());
	MatrixT<TYPE> sum(this->m_iLen, this->m_jLen);
	for (size_type i = 0; i < this->m_iLen; i++) {
		for (size_type j = 0; j < this->m_jLen; j++) {
			sum[i][j] = this->m_mp[i][j] + a[i][j];
		}
	}
	return sum;
}
template<typename TYPE>
MatrixV<TYPE> MatrixV<TYPE>::operator-(const MatrixV<TYPE> &a) {
	assert(a.iLen() == this->iLen());
	assert(a.jLen() == this->jLen());
	MatrixT<TYPE> sum(this->m_iLen, this->m_jLen);
	for (size_type i = 0; i < this->m_iLen; i++) {
		for (size_type j = 0; j < this->m_jLen; j++) {
			sum[i][j] = this->m_mp[i][j] - a[i][j];
		}
	}
	return sum;
}
template<typename TYPE>
MatrixV<TYPE> MatrixV<TYPE>::operator*(const MatrixV<TYPE> &a) {
	assert(a.iLen() == this->jLen());
	size_type nrow = this->m_iLen;
	size_type ncol = a.jLen();
	MatrixT<TYPE> res(nrow, ncol);
	for (size_type i = 0; i < nrow; i++) {
		for (size_type j = 0; j < ncol; j++) {
			res[i][j] = 0;
			for (size_type k = 0; k < this->m_jLen; k++) {
				res[i][j] += this->m_mp[i][k] * a[k][j];
			}
		}
	}
	return res;
}
template<typename TYPE>
void MatrixV<TYPE>::show() const {
	std::cout << "> Matrix " << this->m_iLen << " x " << this->m_jLen << "\n";
	std::cout << "> ";
	for (int i = 0; i < this->m_iLen; i++) {
		for (int j = 0; j < this->m_jLen; j++) {
			std::cout << std::scientific << this->m_mp[i][j] << "  ";
		}
		std::cout << std::endl;
		std::cout << "> ";
	}
	std::cout << "< ----------\n";
}

//===============================================
//Function outside of the class==================
//===============================================
template<typename TYPE>
arrayListV<TYPE> operator*(const MatrixV<TYPE> &m, const arrayListV<TYPE> &a) {
	assert(m.jLen() == a.size());
	arrayList res(a.Len());
	for (int i = 0; i < m.iLen(); i++) {
		for (int j = 0; j < m.jLen(); j++) {
			res[i] += m[i][j] * a[j];
		}
	}
	return res;
}

typedef MatrixV<Float> Matrix;

template<typename TYPE, LarusDef::size_type DIM1, LarusDef::size_type DIM2>
class MatrixS {
public:
	// type definitions===================
	typedef MatrixS<TYPE, DIM1, DIM2> self;
	typedef TYPE value_type;
	typedef TYPE* pointer;
	typedef const TYPE* const_pointer;
	typedef TYPE& reference;
	typedef const TYPE& const_reference;
	typedef LarusDef::size_type size_type;
	typedef LarusDef::size_type difference_type;

	TYPE elems[DIM1 > 0 ? DIM1 : 1][DIM2 > 0 ? DIM2 : 1];

public:
	//constructor==========================
	MatrixS() {
	}
	MatrixS(const TYPE& a) {
		this->assign(a);
	}
	MatrixS(const MatrixS<TYPE, DIM1, DIM2>& a) {
		typedef size_type st;
		for (st i = 0; i < DIM1; i++) {
			for (st j = 0; j < DIM2; j++) {
				elems[i][j] = a[i][j];
			}
		}
	}
	//=============================================
	self& operator=(const self &a) {
		if (this == &a) {
			return *this;
		} else {
			typedef size_type st;
			for (st i = 0; i < DIM1; i++) {
				for (st j = 0; j < DIM2; j++) {
					this->elems[i][j] = a[i][j];
				}
			}
		}
		return *this;
	}
	//=============================================
	~MatrixS() {
	}
	//Capacity=====================================
	static inline size_type size() {
		return DIM1 * DIM2;
	}
	static inline size_type iLen() {
		return DIM1;
	}
	static inline size_type jLen() {
		return DIM2;
	}
	//Element access===============================
	reference operator()(size_type i, size_type j) {
		ASSERT_MSG(i < DIM1 && j < DIM2, "out of range");
		return elems[i][j];
	}
	const_reference operator()(size_type i, size_type j) const {
		ASSERT_MSG(i < DIM1 && j < DIM2, "out of range");
		return elems[i][j];
	}
	reference at(size_type i, size_type j) {
		ASSERT_MSG(i < DIM1 && j < DIM2, "out of range");
		return elems[i][j];
	}
	const_reference at(size_type i, size_type j) const {
		ASSERT_MSG(i < DIM1 && j < DIM2, "out of range");
		return elems[i][j];
	}
	TYPE get(size_type i, size_type j) {
		ASSERT_MSG(i < DIM1 && j < DIM2, "out of range");
		return elems[i][j];
	}
	TYPE* getpValue(size_type i, size_type j) {
		ASSERT_MSG(i < DIM1 && j < DIM2, "out of range");
		return elems[i][j];
	}
	void set(size_type i, size_type j, const TYPE& value) {
		ASSERT_MSG(i < DIM1 && j < DIM2, "out of range");
		elems[i][j] = value;
	}
	void set_row(size_type i, const TYPE& value) {
		ASSERT_MSG(i < DIM1, "out of range");
		typedef size_type st;
		for (st j = 0; j < DIM2; j++) {
			this->elems[i][j] = value;
		}
	}
	void set_col(size_type j, const TYPE& value) {
		ASSERT_MSG(j < DIM2, "out of range");
		typedef size_type st;
		for (st i = 0; i < DIM2; i++) {
			this->elems[i][j] = value;
		}
	}
	void assign(const TYPE& value) {
		typedef size_type st;
		for (st i = 0; i < DIM1; i++) {
			for (st j = 0; j < DIM2; j++) {
				this->elems[i][j] = value;
			}
		}
	}
	inline void ones() {
		assign(1);
	}
	inline void zeros() {
		assign(0);
	}
	//
	void show() const {
		std::cout << "> MatrixS " << DIM1 << " x " << DIM2 << "\n";
		std::cout << "> ";
		for (int i = 0; i < DIM1; i++) {
			for (int j = 0; j < DIM2; j++) {
				std::cout << std::scientific << this->elems[i][j] << "  ";
			}
			std::cout << std::endl;
			std::cout << "> ";
		}
		std::cout << "< ----------\n";
	}
};

typedef MatrixS<Float, 2, 2> Matrix2x2;
typedef MatrixS<Float, 3, 3> Matrix3x3;
typedef MatrixS<Float, 4, 4> Matrix4x4;

template<typename TYPE, LarusDef::size_type DIM1, LarusDef::size_type DIM2>
static inline void zeros(MatrixS<Float, DIM1, DIM2>& m ) {
	m.assign(0);
}

template<typename TYPE, LarusDef::size_type DIM1, LarusDef::size_type DIM2>
static inline void ones(MatrixS<Float, DIM1, DIM2>& m ) {
	m.assign(1);
}


/*
 * \brief calculate the determinant of a 2x2 matrix.
 *
 *        Adapted from:
 *        Matrix Inversion
 *        by Richard Carling
 *        from "Graphics Gems", Academic Press, 1990
 */
template<typename TYPE>
static inline TYPE det2x2(const TYPE& a, const TYPE& b, const TYPE& c,
		const TYPE& d) {
	return a * d - b * c;
}

template<typename TYPE>
inline TYPE determinant(const MatrixS<TYPE, 2, 2>& m) {
	TYPE a, b, c, d;
	a = m(0, 0);
	b = m(0, 1);
	c = m(1, 0);
	d = m(1, 1);
	return det2x2(a, b, c, d);
}

/*
 * \brief calculate the determinant of a 3x3 matrix
 *        in the form
 *
 *        | a1,  b1,  c1 |
 *        | a2,  b2,  c2 |
 *        | a3,  b3,  c3 |
 *
 *        Adapted from:
 *        Matrix Inversion
 *        by Richard Carling
 *        from "Graphics Gems", Academic Press, 1990
 */
template<typename TYPE>
static inline TYPE det3x3(const TYPE& a1, const TYPE& a2, const TYPE& a3,
		const TYPE& b1, const TYPE& b2, const TYPE& b3, const TYPE& c1,
		const TYPE& c2, const TYPE& c3) {
	TYPE ans3;
	ans3 = a1 * det2x2(b2, b3, c2, c3) - b1 * det2x2(a2, a3, c2, c3)
			+ c1 * det2x2(a2, a3, b2, b3);
	return ans3;
}

template<typename TYPE>
inline TYPE determinant(const MatrixS<TYPE, 3, 3>& m) {
	const TYPE& a1 = m(0, 0), b1 = m(0, 1), c1 = m(0, 2);
	const TYPE& a2 = m(1, 0), b2 = m(1, 1), c2 = m(1, 2);
	const TYPE& a3 = m(2, 0), b3 = m(2, 1), c3 = m(2, 2);
	return det3x3(a1, a2, a3, b1, b2, b3, c1, c2, c3);
}

/**
 * \brief matrix_determinant:
 * \param m Matrix 4 x 4 .
 *
 * \returns: the value of det(\p m).
 */
template<typename TYPE>
inline TYPE determinant(const MatrixS<TYPE, 4, 4>& m) {
	TYPE ans4;
	const TYPE& a1 = m(0, 0), b1 = m(0, 1), c1 = m(0, 2), d1 = m(0, 3);
	const TYPE& a2 = m(1, 0), b2 = m(1, 1), c2 = m(1, 2), d2 = m(1, 3);
	const TYPE& a3 = m(2, 0), b3 = m(2, 1), c3 = m(2, 2), d3 = m(2, 3);
	const TYPE& a4 = m(3, 0), b4 = m(3, 1), c4 = m(3, 2), d4 = m(3, 3);

	ans4 = a1 * det3x3(b2, b3, b4, c2, c3, c4, d2, d3, d4)
			- b1 * det3x3(a2, a3, a4, c2, c3, c4, d2, d3, d4)
			+ c1 * det3x3(a2, a3, a4, b2, b3, b4, d2, d3, d4)
			- d1 * det3x3(a2, a3, a4, b2, b3, b4, c2, c3, c4);

	return ans4;
}

template<typename TYPE>
inline void adjoint(const MatrixS<TYPE, 4, 4>& m, MatrixS<TYPE, 4, 4>& ma) {
	const TYPE& a1 = m(0, 0), b1 = m(0, 1), c1 = m(0, 2), d1 = m(0, 3);
	const TYPE& a2 = m(1, 0), b2 = m(1, 1), c2 = m(1, 2), d2 = m(1, 3);
	const TYPE& a3 = m(2, 0), b3 = m(2, 1), c3 = m(2, 2), d3 = m(2, 3);
	const TYPE& a4 = m(3, 0), b4 = m(3, 1), c4 = m(3, 2), d4 = m(3, 3);
	/* row column labeling reversed since we transpose rows & columns */
	// the matrix without transpose is called cofactor matrix
	// the matrix with    transpose is called adjoint  matrix
	ma(0, 0) = det3x3(b2, b3, b4, c2, c3, c4, d2, d3, d4);
	ma(1, 0) = -det3x3(a2, a3, a4, c2, c3, c4, d2, d3, d4);
	ma(2, 0) = det3x3(a2, a3, a4, b2, b3, b4, d2, d3, d4);
	ma(3, 0) = -det3x3(a2, a3, a4, b2, b3, b4, c2, c3, c4);

	ma(0, 1) = -det3x3(b1, b3, b4, c1, c3, c4, d1, d3, d4);
	ma(1, 1) = det3x3(a1, a3, a4, c1, c3, c4, d1, d3, d4);
	ma(2, 1) = -det3x3(a1, a3, a4, b1, b3, b4, d1, d3, d4);
	ma(3, 1) = det3x3(a1, a3, a4, b1, b3, b4, c1, c3, c4);

	ma(0, 2) = det3x3(b1, b2, b4, c1, c2, c4, d1, d2, d4);
	ma(1, 2) = -det3x3(a1, a2, a4, c1, c2, c4, d1, d2, d4);
	ma(2, 2) = det3x3(a1, a2, a4, b1, b2, b4, d1, d2, d4);
	ma(3, 2) = -det3x3(a1, a2, a4, b1, b2, b4, c1, c2, c4);

	ma(0, 3) = -det3x3(b1, b2, b3, c1, c2, c3, d1, d2, d3);
	ma(1, 3) = det3x3(a1, a2, a3, c1, c2, c3, d1, d2, d3);
	ma(2, 3) = -det3x3(a1, a2, a3, b1, b2, b3, d1, d2, d3);
	ma(3, 3) = det3x3(a1, a2, a3, b1, b2, b3, c1, c2, c3);
}

template<typename TYPE>
inline void adjoint(const MatrixS<TYPE, 3, 3>& m, MatrixS<TYPE, 3, 3>& ma) {
	/* row column labeling reversed since we transpose rows & columns */
	// the matrix without transpose is called cofactor matrix
	// the matrix with    transpose is called adjoint  matrix
	ma(0, 0) = (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1));
	ma(0, 1) = (m(2, 1) * m(0, 2) - m(0, 1) * m(2, 2));
	ma(0, 2) = (m(0, 1) * m(1, 2) - m(1, 1) * m(0, 2));
	ma(1, 0) = (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2));
	ma(1, 1) = (m(0, 0) * m(2, 2) - m(2, 0) * m(0, 2));
	ma(1, 2) = (m(1, 0) * m(0, 2) - m(0, 0) * m(1, 2));
	ma(2, 0) = (m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1));
	ma(2, 1) = (m(2, 0) * m(0, 1) - m(0, 0) * m(2, 1));
	ma(2, 2) = (m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0));
}

template<typename TYPE>
inline int inverse(const MatrixS<TYPE, 4, 4>& m, MatrixS<TYPE, 4, 4>& m_inv) {
	typedef MatrixS<TYPE, 4, 4> matrix;
	TYPE det = determinant(m);
	if (det == 0.) {
		return -1;
	}
	adjoint(m, m_inv);
	for (typename matrix::size_type i = 0; i < 4; ++i) {
		for (typename matrix::size_type j = 0; j < 4; ++j) {
			m_inv(i, j) /= det;
		}
	}
	return 1;
}

template<typename TYPE>
inline int inverse(const MatrixS<TYPE, 3, 3>& m, MatrixS<TYPE, 3, 3>& m_inv) {
	typedef MatrixS<TYPE, 3, 3> matrix;
	TYPE det = determinant(m);
	if (det == 0.) {
		return -1;
	}
	adjoint(m, m_inv);
	for (typename matrix::size_type i = 0; i < m.iLen(); ++i) {
		for (typename matrix::size_type j = 0; j < m.jLen(); ++j) {
			m_inv(i, j) /= det;
		}
	}
	return 1;
}

template<typename TYPE>
void transpose(const MatrixS<TYPE, 4, 4> m, MatrixS<TYPE, 4, 4> mt) {
	mt(0, 0) = m(0, 0);
	mt(1, 0) = m(0, 1);
	mt(2, 0) = m(0, 2);
	mt(3, 0) = m(0, 3);
	mt(0, 1) = m(1, 0);
	mt(1, 1) = m(1, 1);
	mt(2, 1) = m(1, 2);
	mt(3, 1) = m(1, 3);
	mt(0, 2) = m(2, 0);
	mt(1, 2) = m(2, 1);
	mt(2, 2) = m(2, 2);
	mt(3, 2) = m(2, 3);
	mt(0, 3) = m(3, 0);
	mt(1, 3) = m(3, 1);
	mt(2, 3) = m(3, 2);
	mt(3, 3) = m(3, 3);
}

template<typename TYPE>
void transpose(const MatrixS<TYPE, 3, 3> m, MatrixS<TYPE, 3, 3> mt) {
	mt(0, 0) = m(0, 0);
	mt(1, 0) = m(0, 1);
	mt(2, 0) = m(0, 2);
	mt(0, 1) = m(1, 0);
	mt(1, 1) = m(1, 1);
	mt(2, 1) = m(1, 2);
	mt(0, 2) = m(2, 0);
	mt(1, 2) = m(2, 1);
	mt(2, 2) = m(2, 2);
}

template<typename TYPE>
void transpose(const MatrixS<TYPE, 2, 2> m, MatrixS<TYPE, 2, 2> mt) {
	mt(0, 0) = m(0, 0);
	mt(1, 0) = m(0, 1);
	mt(0, 1) = m(1, 0);
	mt(1, 1) = m(1, 1);
}

}	//namespace=====================================

#endif /* MATRIXT_H_ */
