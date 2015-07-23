/************************
 //  \file   Space.h
 //  \brief
 // 
 //  \author czhou
 //  \date   20 janv. 2015 
 ***********************/
#ifndef SPACE_H_
#define SPACE_H_

#include <iostream>
#include <assert.h>
#include "../TypeDef.h"
#include "../Utility/ArrayList.h"

namespace Larus {

template<typename TYPE, LarusDef::size_type Dim>
class SpaceT {
public:
	// type definitions===================
	typedef TYPE value_type;
	typedef TYPE* pointer;
	typedef const TYPE* const_pointer;
	typedef TYPE& reference;
	typedef const TYPE& const_reference;
	typedef LarusDef::size_type size_type;
	typedef LarusDef::size_type difference_type;
	//static const int DIM = Dim;
private:
	array<size_type, Dim> m_len;
	arrayListT<TYPE> m_mp;
public:
	//constructor==========================
	SpaceT();
	SpaceT(const SpaceT<TYPE, Dim>& a);
	SpaceT(size_type iLen, size_type = 0, size_type = 0);

	void reconstruct(size_type iLen, size_type = 0, size_type = 0);
	//=============================================
	SpaceT<TYPE, Dim>& operator=(const SpaceT<TYPE, Dim>& a);
	//=============================================
	~SpaceT() {
	}

	//Capacity=====================================
	size_type size() const {
		return m_mp.size();
	}
	size_type iLen() const {
		return m_len[0];
	}
	size_type jLen() const {
		return Dim >= 2 ? m_len[1] : 0;
	}
	size_type kLen() const {
		return Dim >= 3 ? m_len[2] : 0;
	}
	bool empty() const {
		return m_mp.empty();
	}
	//Element access===============================
	//arrayListT<TYPE>& operator[](size_type index);
	//const arrayListT<TYPE>& operator[](size_type index) const;
	size_type to_1d_idx(size_type i, size_type = 0, size_type = 0) const;

	reference operator()(size_type i, size_type = 0, size_type = 0);
	const_reference operator()(size_type i, size_type = 0, size_type = 0) const;
	reference at(size_type i, size_type = 0, size_type = 0);
	const_reference at(size_type i, size_type = 0, size_type = 0) const;

	reference at_1d(size_type i);
	const_reference at_1d(size_type i) const;

	TYPE get(size_type i, size_type = 0, size_type = 0);
	void set(const TYPE& value, size_type i, size_type = 0, size_type = 0);
	//void set_row(size_type i, const TYPE& value);
	//void set_col(size_type i, const TYPE& value);
	void assign(const TYPE& value);
	//element access===============================

	TYPE* getpValue(size_type i, size_type = 0, size_type = 0);
	//void swapValue(size_type i1, size_type j1, size_type i2, size_type j2);
	//size_type countNumEq(const TYPE& value); //overload ==

	inline bool testIdx(size_type dim, size_type idx) const {
		ASSERT(dim < Dim);
		if (idx >= 0 && idx < m_len[dim]) {
			return true;
		} else {
			return false;
		}
	}
	inline bool testIdxIJK(size_type i, size_type j, size_type k) const {
		return testIdx(0, i) && (Dim >= 2) ? testIdx(1, j) :
				true && (Dim >= 3) ? testIdx(2, k) : true;
	}

	inline size_type count_equal(const TYPE& nd) const {
		return m_mp.count_equal(nd);
	}
};

template<typename TYPE, LarusDef::size_type Dim>
SpaceT<TYPE, Dim>::SpaceT() {
	m_len.assign(0);
}
template<typename TYPE, LarusDef::size_type Dim>
SpaceT<TYPE, Dim>::SpaceT(const SpaceT<TYPE, Dim>& a) {
	this->m_len = a.m_len;
	this->m_mp = a.m_mp;
}
template<typename TYPE, LarusDef::size_type Dim>
SpaceT<TYPE, Dim>::SpaceT(size_type iLen, size_type jLen, size_type kLen) {
	this->m_len[0] = iLen;
	if (Dim >= 2) {
		ASSERT(iLen > 0 && jLen > 0);
		this->m_len[1] = jLen;
		this->m_mp.reconstruct(iLen * jLen);
	}
	if (Dim >= 3) {
		ASSERT(iLen > 0 && jLen > 0 && kLen > 0);
		this->m_len[2] = kLen;
		this->m_mp.reconstruct(iLen * jLen * kLen);
	}
}
template<typename TYPE, LarusDef::size_type Dim>
void SpaceT<TYPE, Dim>::reconstruct(size_type iLen, size_type jLen,
		size_type kLen) {
	this->m_len[0] = iLen;
	if (Dim >= 2) {
		ASSERT(iLen > 0 && jLen > 0);
		this->m_len[1] = jLen;
		this->m_mp.reconstruct(iLen * jLen);
	}
	if (Dim >= 3) {
		ASSERT(iLen > 0 && jLen > 0 && kLen > 0);
		this->m_len[2] = kLen;
		this->m_mp.reconstruct(iLen * jLen * kLen);
	}
}
template<typename TYPE, LarusDef::size_type Dim>
SpaceT<TYPE, Dim>& SpaceT<TYPE, Dim>::operator=(const SpaceT<TYPE, Dim>& a) {
	if (this == &a) {
		return *this;
	}
	this->m_len = a.m_len;
	this->m_mp = a.m_mp;
	return *this;
}
template<typename TYPE, LarusDef::size_type Dim>
typename SpaceT<TYPE, Dim>::size_type SpaceT<TYPE, Dim>::to_1d_idx(SpaceT<TYPE, Dim>::size_type i,
		SpaceT<TYPE, Dim>::size_type j, SpaceT<TYPE, Dim>::size_type k) const {
	ASSERT(i < this->m_len[0]);
	if (Dim >= 2)
		ASSERT(j < this->m_len[1]);
	if (Dim >= 3)
		ASSERT(k < this->m_len[2]);
	array<size_type, Dim> inp;
	inp[0] = i;
	inp[1] = j;
	if (Dim >= 3) {
		inp[2] = k;
	}
	size_type idx = 0;
	for (size_type ii = 0; ii < Dim; ii++) {
		size_type b = 1;
		for (size_type jj = ii + 1; jj < Dim; jj++) {
			b *= m_len[jj];
		}
		idx += b * inp[ii];
	}
	return idx;
}

template<typename TYPE, LarusDef::size_type Dim>
TYPE& SpaceT<TYPE, Dim>::at(SpaceT<TYPE, Dim>::size_type i,
		SpaceT<TYPE, Dim>::size_type j, SpaceT<TYPE, Dim>::size_type k) {
	size_type idx = to_1d_idx(i,j,k);
	return m_mp[idx];
}
template<typename TYPE, LarusDef::size_type Dim>
const TYPE& SpaceT<TYPE, Dim>::at(SpaceT<TYPE, Dim>::size_type i,
		SpaceT<TYPE, Dim>::size_type j, SpaceT<TYPE, Dim>::size_type k) const {
	size_type idx = to_1d_idx(i,j,k);
	return m_mp[idx];
}
template<typename TYPE, LarusDef::size_type Dim>
typename SpaceT<TYPE, Dim>::reference SpaceT<TYPE, Dim>::at_1d(
		SpaceT<TYPE, Dim>::size_type i) {
	return m_mp[i];
}
template<typename TYPE, LarusDef::size_type Dim>
typename SpaceT<TYPE, Dim>::const_reference SpaceT<TYPE, Dim>::at_1d(
		SpaceT<TYPE, Dim>::size_type i) const {
	return m_mp[i];
}
template<typename TYPE, LarusDef::size_type Dim>
TYPE& SpaceT<TYPE, Dim>::operator()(SpaceT<TYPE, Dim>::size_type i,
		SpaceT<TYPE, Dim>::size_type j, SpaceT<TYPE, Dim>::size_type k) {
	return at(i, j, k);
}
template<typename TYPE, LarusDef::size_type Dim>
const TYPE& SpaceT<TYPE, Dim>::operator()(SpaceT<TYPE, Dim>::size_type i,
		SpaceT<TYPE, Dim>::size_type j, SpaceT<TYPE, Dim>::size_type k) const {
	return at(i, j, k);
}
template<typename TYPE, LarusDef::size_type Dim>
TYPE SpaceT<TYPE, Dim>::get(SpaceT<TYPE, Dim>::size_type i,
		SpaceT<TYPE, Dim>::size_type j, SpaceT<TYPE, Dim>::size_type k) {
	return at(i, j, k);
}
template<typename TYPE, LarusDef::size_type Dim>
void SpaceT<TYPE, Dim>::set(const TYPE& value, SpaceT<TYPE, Dim>::size_type i,
		SpaceT<TYPE, Dim>::size_type j, SpaceT<TYPE, Dim>::size_type k) {
	this->at(i, j, k) = value;
}
template<typename TYPE, LarusDef::size_type Dim>
void SpaceT<TYPE, Dim>::assign(const TYPE& value) {
	m_mp.assign(value);
}
template<typename TYPE, LarusDef::size_type Dim>
TYPE* SpaceT<TYPE, Dim>::getpValue(SpaceT<TYPE, Dim>::size_type i,
		SpaceT<TYPE, Dim>::size_type j, SpaceT<TYPE, Dim>::size_type k) {
	return &at(i, j, k);
}

}

#endif /* SPACE_H_ */
