/************************
 //  \file   Cube.h
 //  \brief
 // 
 //  \author czhou
 //  \date   28 janv. 2015 
 ***********************/
#ifndef CUBE_H_
#define CUBE_H_

#include <iostream>
#include <assert.h>
#include "../TypeDef.h"
#include "../Utility/ArrayList.h"
#include "Matrix.h"

namespace Larus {

template<typename TYPE>
class CubeT {
public:
	// type definitions===================
	typedef TYPE value_type;
	typedef TYPE* pointer;
	typedef const TYPE* const_pointer;
	typedef TYPE& reference;
	typedef const TYPE& const_reference;
	typedef LarusDef::size_type size_type;
	typedef LarusDef::size_type difference_type;
private:
	size_type m_iLen;
	size_type m_jLen;
	size_type m_kLen;
	MatrixT<TYPE> *m_mp;

	void assign_samesize(const CubeT<TYPE>& a);
public:

	//constructor==========================
	CubeT();
	CubeT(const CubeT<TYPE>& a);
	CubeT(size_type iLen, size_type jLen, size_type kLen);

	void reconstruct(size_type iLen, size_type jLen, size_type kLen);
	//=============================================
	CubeT<TYPE>& operator=(const CubeT<TYPE> &a);
	//=============================================
	~CubeT();
	//Capacity=====================================
	size_type size() const;
	size_type iLen() const;
	size_type jLen() const;
	size_type kLen() const;
	bool empty() const;
	//Element access===============================
	MatrixT<TYPE>& operator[](size_type index);
	const MatrixT<TYPE>& operator[](size_type index) const;
	reference at(size_type i, size_type j, size_type k);
	const_reference at(size_type i, size_type j, size_type k) const;
	TYPE get(size_type i, size_type j, size_type k);
	TYPE* getpValue(size_type i, size_type j, size_type k);
	void set(size_type i, size_type j, size_type k, const TYPE& value);
	void set_iplane(size_type i, const TYPE& value);
	void set_jplane(size_type i, const TYPE& value);
	void set_kplane(size_type i, const TYPE& value);
	void assign(const TYPE& value);
	//=============================================
	void swapValue( //
			size_type i1, size_type j1, size_type k1, //
			size_type i2, size_type j2, size_type k2);

	size_type countNumEq(TYPE value); //overload ==

	inline bool testIdxI(size_type);
	inline bool testIdxJ(size_type);
	inline bool testIdxK(size_type);
};

template<typename TYPE>
inline bool isSameSize(const CubeT<TYPE>& a, const CubeT<TYPE>& b) {
	return (b.iLen() == a.iLen() && b.jLen() == a.jLen() && b.kLen() == a.kLen());
}
//===============================================

template<typename TYPE>
inline bool CubeT<TYPE>::testIdxI(LarusDef::size_type i) {
	return (i >= 0 && i < m_iLen) ? true : false;
}

template<typename TYPE>
inline bool CubeT<TYPE>::testIdxJ(LarusDef::size_type j) {
	return (j >= 0 && j < m_jLen) ? true : false;
}

template<typename TYPE>
inline bool CubeT<TYPE>::testIdxK(LarusDef::size_type k) {
	return (k >= 0 && k < m_kLen) ? true : false;
}

template<typename TYPE>
void CubeT<TYPE>::assign_samesize(const CubeT<TYPE>& a) {
	for (size_type i = 0; i < m_iLen; i++) {
		for (size_type j = 0; j < m_jLen; j++) {
			for (size_type k = 0; k < m_jLen; k++) {
				m_mp[i][j][k] = a[i][j][k];
			}
		}
	}
}

template<typename TYPE>
CubeT<TYPE>::CubeT() {
	m_iLen = 0;
	m_jLen = 0;
	m_kLen = 0;
	m_mp = NULL_PTR;
}

template<typename TYPE>
CubeT<TYPE>::CubeT(const CubeT<TYPE>& a) {
	m_iLen = a.iLen();
	m_jLen = a.jLen();
	m_kLen = a.kLen();
	m_mp = new MatrixT<TYPE> [m_iLen];
	for (size_type i = 0; i < m_iLen; i++) {
		m_mp[i].reconstruct(m_jLen, m_kLen);
	}
	assign_samesize(a);
}
template<typename TYPE>
CubeT<TYPE>::CubeT( //
		LarusDef::size_type iLen, //
		LarusDef::size_type jLen, //
		LarusDef::size_type kLen) {
	m_iLen = iLen;
	m_jLen = jLen;
	m_kLen = kLen;
	m_mp = new MatrixT<TYPE> [m_iLen];
	for (size_type i = 0; i < m_iLen; i++) {
		m_mp[i].reconstruct(m_jLen, m_kLen);
	}
}
template<typename TYPE>
void CubeT<TYPE>::reconstruct( //
		LarusDef::size_type iLen, //
		LarusDef::size_type jLen, //
		LarusDef::size_type kLen) {
	m_iLen = iLen;
	m_jLen = jLen;
	m_kLen = kLen;
	//m_total = m_iLen * m_jLen;
	if (m_mp != NULL_PTR) {
		delete[] m_mp;
	}
	m_mp = new MatrixT<TYPE> [m_iLen];
	for (size_type i = 0; i < m_iLen; i++) {
		m_mp[i].reconstruct(m_jLen, m_kLen);
	}
}
template<typename TYPE>
CubeT<TYPE>::~CubeT() {
	delete[] m_mp;
}
template<typename TYPE>
CubeT<TYPE>& CubeT<TYPE>::operator=(const CubeT<TYPE> &a) {
	if (this == &a) {
		return *this;
	}
	if (isSameSize((*this), a)) {
		assign_samesize(a);
	} else {
		reconstruct(a.iLen(), a.jLen(), a.kLen());
		assign_samesize(a);
	}
	return *this;
}

template<typename TYPE>
typename CubeT<TYPE>::size_type CubeT<TYPE>::size() const {
	return m_iLen * m_jLen * m_kLen;
}
template<typename TYPE>
typename CubeT<TYPE>::size_type CubeT<TYPE>::iLen() const {
	return m_iLen;
}
template<typename TYPE>
typename CubeT<TYPE>::size_type CubeT<TYPE>::jLen() const {
	return m_jLen;
}
template<typename TYPE>
typename CubeT<TYPE>::size_type CubeT<TYPE>::kLen() const {
	return m_kLen;
}
template<typename TYPE>
bool CubeT<TYPE>::empty() const {
	return (m_mp == NULL_PTR);
}
template<typename TYPE>
const MatrixT<TYPE>& CubeT<TYPE>::operator[](LarusDef::size_type index) const {
	ASSERT(index >= 0 && index < m_iLen);
	return m_mp[index];
}
template<typename TYPE>
MatrixT<TYPE>& CubeT<TYPE>::operator[](typename CubeT<TYPE>::size_type index) {
	assert(index >= 0 && index < m_iLen);
	return m_mp[index];
}
template<typename TYPE>
typename CubeT<TYPE>::reference CubeT<TYPE>::at(
		typename CubeT<TYPE>::size_type i, typename CubeT<TYPE>::size_type j,
		typename CubeT<TYPE>::size_type k) {
	ASSERT(i >= 0 && i < m_iLen);
	ASSERT(j >= 0 && j < m_jLen);
	ASSERT(k >= 0 && k < m_kLen);
	return m_mp[i][j][k];
}
template<typename TYPE>
typename CubeT<TYPE>::const_reference CubeT<TYPE>::at(
		typename CubeT<TYPE>::size_type i, typename CubeT<TYPE>::size_type j,
		typename CubeT<TYPE>::size_type k) const {
	ASSERT(i >= 0 && i < m_iLen);
	ASSERT(j >= 0 && j < m_jLen);
	ASSERT(k >= 0 && k < m_kLen);
	return m_mp[i][j][k];
}
template<typename TYPE>
TYPE CubeT<TYPE>::get(typename CubeT<TYPE>::size_type i,
		typename CubeT<TYPE>::size_type j, typename CubeT<TYPE>::size_type k) {
	ASSERT(i >= 0 && i < m_iLen);
	ASSERT(j >= 0 && j < m_jLen);
	ASSERT(k >= 0 && k < m_kLen);
	return m_mp[i][j][k];
}
template<typename TYPE>
TYPE* CubeT<TYPE>::getpValue(typename CubeT<TYPE>::size_type i,
		typename CubeT<TYPE>::size_type j, typename CubeT<TYPE>::size_type k) {
	assert(i >= 0 && i < m_iLen);
	assert(j >= 0 && j < m_jLen);
	assert(k >= 0 && k < m_kLen);
	return &m_mp[i][j][k];
}
template<typename TYPE>
void CubeT<TYPE>::set(typename CubeT<TYPE>::size_type i,
		typename CubeT<TYPE>::size_type j, typename CubeT<TYPE>::size_type k,
		const TYPE& value) {
	ASSERT(i >= 0 && i < m_iLen);
	ASSERT(j >= 0 && j < m_jLen);
	ASSERT(k >= 0 && k < m_kLen);
	m_mp[i][j][k] = value;
}
template<typename TYPE>
void CubeT<TYPE>::set_iplane(typename CubeT<TYPE>::size_type i,
		const TYPE& value) {
	ASSERT(i >= 0 && i < m_iLen);
	for (size_type j = 0; j < m_jLen; j++) {
		for (size_type k = 0; k < m_jLen; k++) {
			m_mp[i][j][k] = value;
		}
	}
}
template<typename TYPE>
void CubeT<TYPE>::set_jplane(typename CubeT<TYPE>::size_type j,
		const TYPE& value) {
	ASSERT(j >= 0 && j < m_jLen);
	for (size_type i = 0; i < m_iLen; i++) {
		for (size_type k = 0; k < m_jLen; k++) {
			m_mp[i][j][k] = value;
		}
	}
}

template<typename TYPE>
void CubeT<TYPE>::set_kplane(typename CubeT<TYPE>::size_type k,
		const TYPE& value) {
	ASSERT(k >= 0 && k < m_kLen);
	for (size_type i = 0; i < m_iLen; i++) {
		for (size_type j = 0; j < m_jLen; j++) {
			m_mp[i][j][k] = value;
		}
	}
}
template<typename TYPE>
void CubeT<TYPE>::assign(const TYPE& value) {
	for (size_type i = 0; i < m_iLen; i++) {
		for (size_type j = 0; j < m_jLen; j++) {
			for (size_type k = 0; k < m_jLen; k++) {
				m_mp[i][j][k] = value;
			}
		}
	}
}
template<typename TYPE>
void CubeT<TYPE>::swapValue(typename CubeT<TYPE>::size_type i1,
		typename CubeT<TYPE>::size_type j1, typename CubeT<TYPE>::size_type k1,
		typename CubeT<TYPE>::size_type i2, typename CubeT<TYPE>::size_type j2,
		typename CubeT<TYPE>::size_type k2) {
	assert(i1 >= 0 && i1 < m_iLen && i2 >= 0 && i2 < m_iLen);
	assert(j1 >= 0 && j1 < m_jLen && j2 >= 0 && j2 < m_jLen);
	assert(k1 >= 0 && k1 < m_kLen && k2 >= 0 && k2 < m_kLen);
	if (i1 == i2 && j1 == j2 && k1 == k2) {
		return;
	}
	TYPE tmp;
	tmp = m_mp[i1][j1][k1];
	m_mp[i1][j1][k1] = m_mp[i2][j2][k2];
	m_mp[i2][j2][k2] = tmp;
}
template<typename TYPE>
typename CubeT<TYPE>::size_type CubeT<TYPE>::countNumEq(TYPE value) { //overload ==
	size_type num = 0;
	for (size_type i = 0; i < m_iLen; i++) {
		for (size_type j = 0; j < m_jLen; j++) {
			for (size_type k = 0; k < m_jLen; k++) {
				if (m_mp[i][j] == value) {
					num++;
				}
			}
		}
	}
	return num;
}
} //end of namespace

#endif /* CUBE_H_ */
