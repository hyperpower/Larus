//
//  array.hpp
//  LarusX
//
//  Created by zhou on 13/12/14.
//  Copyright (c) 2014 zhou. All rights reserved.
//

#ifndef _ARRAY_h
#define _ARRAY_h

#include <stddef.h>
#include "../TypeDef.h"

namespace Larus {

template<class T, LarusDef::size_type N = 1>
class array: public ObjectBase {
protected:
	T elems[N];    //!< fixed-size array of elements of type T
public:
	// type definitions
	typedef T value_type;
	typedef T* iterator;
	typedef const T* const_iterator;
	typedef T& reference;
	typedef const T& const_reference;
	typedef LarusDef::size_type size_type;
	typedef LarusDef::size_type difference_type;

	// iterator support
	iterator begin() {
		return elems;
	}
	const_iterator begin() const {
		return elems;
	}

	iterator end() {
		return elems + N;
	}
	const_iterator end() const {
		return elems + N;
	}
	/**
	 * array constructor.
	 */
	array() {
	}
	/**
	 * array copy constructor
	 */
	array(const array<T, N>& p) {
		//loop expand
		if (N == 2) {
			this->elems[0] = p.elems[0];
			this->elems[1] = p.elems[1];
		} else if (N == 3) {
			this->elems[0] = p.elems[0];
			this->elems[1] = p.elems[1];
			this->elems[2] = p.elems[2];
		} else if (N == 4) {
			this->elems[0] = p.elems[0];
			this->elems[1] = p.elems[1];
			this->elems[2] = p.elems[2];
			this->elems[3] = p.elems[3];
		} else {
			for (size_type i = 0; i < N; ++i)
				elems[i] = p.elems[i];
		}
	}

	//operator=
	array<T, N>& operator=(const array<T, N>& rhs) {
		if (&rhs == this) {
			return *this;
		} else {
			if (N == 2) {
				this->elems[0] = rhs.elems[0];
				this->elems[1] = rhs.elems[1];
			} else if (N == 3) {
				this->elems[0] = rhs.elems[0];
				this->elems[1] = rhs.elems[1];
				this->elems[2] = rhs.elems[2];
			} else if (N == 4) {
				this->elems[0] = rhs.elems[0];
				this->elems[1] = rhs.elems[1];
				this->elems[2] = rhs.elems[2];
				this->elems[3] = rhs.elems[3];
			} else {
				for (size_type i = 0; i < N; ++i)
					elems[i] = rhs.elems[i];
			}
		}
		return *this;
	}

	// operator[]
	reference operator[](size_type i) {
		ASSERT_MSG(i < N, "out of range");
		return elems[i];
	}

	const_reference operator[](size_type i) const {
		ASSERT_MSG(i < N, "out of range");
		return elems[i];
	}

	// operator()
	reference operator()(size_type i) {
		ASSERT_MSG(i < N, "out of range");
		return elems[i];
	}

	const_reference operator()(size_type i) const {
		ASSERT_MSG(i < N, "out of range");
		return elems[i];
	}

	// at() with range check
	reference at(size_type i) {
		ASSERT_MSG(i < N, "out of range");
		return elems[i];
	}
	const_reference at(size_type i) const {
		ASSERT_MSG(i < N, "out of range");
		return elems[i];
	}

	// front() and back()
	reference front() {
		return elems[0];
	}

	const_reference front() const {
		return elems[0];
	}

	reference back() {
		return elems[ N > 0 ? N-1 : 0];
	}

	const_reference back() const {
		return elems[ N > 0 ? N-1 : 0];
	}

	// assign one value to all elements
	void assign(const T& value) {
		fill(value);
	}    // A synonym for fill
	/**
	 * \brief fill array with
	 *
	 * \param[in] value The \p value set to array.
	 *
     */
	void fill(const T& value) {
		if (N == 2) {
			this->elems[0] = value;
			this->elems[1] = value;
		} else if (N == 3) {
			this->elems[0] = value;
			this->elems[1] = value;
			this->elems[2] = value;
		} else if (N == 4) {
			this->elems[0] = value;
			this->elems[1] = value;
			this->elems[2] = value;
			this->elems[3] = value;
		} else {
			for (size_type i = 0; i < N; ++i)
				elems[i] = value;
		}
	}
	/**
	 * \brief get size of the array
     */
	static size_type size() {
		return N;
	}
	static bool isEmpty() {
		return false;
	}

	inline value_type* data(){
		return elems;
	}
	inline const value_type* data() const{
		return elems;
	}

	void swap(size_type i1, size_type i2) {
		assert(i1 >= 0 && i1 < N && i2 >= 0 && i2 < N);
		if (i1 == i2) {
			return;
		}
		T tmp;
		tmp = elems[i1];
		elems[i1] = elems[i2];
		elems[i2] = tmp;
	}
};

//partial 2 2 2 2 2 2 2 2 2 2 2 2
//===============================
template<class T>
class array_2: public array<T, 2> {
public:
	//empty constructor===================
	array_2() {
		this->elems[0] = T();
		this->elems[1] = T();
	}

	//constructor=========================
	array_2(const T& a, const T& b) {
		this->elems[0] = a;
		this->elems[1] = b;
	}

	//copy constructor====================
	array_2(const array_2<T>& p) {
		this->elems[0] = p.elems[0];
		this->elems[1] = p.elems[1];
	}

	//operator=
	array_2<T>& operator=(const array_2<T>& rhs) {
		if (&rhs == this) {
			return *this;
		} else {
			this->elems[0] = rhs.elems[0];
			this->elems[1] = rhs.elems[1];
		}
		return *this;
	}

	void assign(const T& value) {
		fill(value);
	}    // A synonym for fill
	void fill(const T& value) {
		this->elems[0] = value;
		this->elems[1] = value;
	}
};

//partial 3 3 3 3 3 3 3 3 3 3 3 3
//===============================
template<class T>
class array_3: public array<T, 3> {
public:
	//empty constructor===================
	array_3() {
		this->elems[0] = T();
		this->elems[1] = T();
		this->elems[2] = T();
	}

	//constructor=========================
	array_3(const T& a, const T& b, const T& c) {
		this->elems[0] = a;
		this->elems[1] = b;
		this->elems[2] = c;
	}

	//copy constructor====================
	array_3(const array_3<T>& p) {
		this->elems[0] = p.elems[0];
		this->elems[1] = p.elems[1];
		this->elems[2] = p.elems[2];
	}

	//operator=
	array_3<T>& operator=(const array_3<T>& rhs) {
		if (&rhs == this) {
			return *this;
		} else {
			this->elems[0] = rhs.elems[0];
			this->elems[1] = rhs.elems[1];
			this->elems[2] = rhs.elems[2];
		}
		return *this;
	}

	void assign(const T& value) {
		fill(value);
	}    // A synonym for fill

	void fill(const T& value) {
		this->elems[0] = value;
		this->elems[1] = value;
		this->elems[2] = value;
	}
};

//partial 4 4 4 4 4 4 4 4 4 4 4 4
//===============================
template<class T>
class array_4: public array<T, 4> {
public:
	//empty constructor===================
	array_4() {
		this->elems[0] = T();
		this->elems[1] = T();
		this->elems[2] = T();
		this->elems[3] = T();
	}

	//constructor=========================
	array_4(const T& a, const T& b, const T& c, const T& d) {
		this->elems[0] = a;
		this->elems[1] = b;
		this->elems[2] = c;
		this->elems[3] = d;
	}

	//copy constructor====================
	array_4(const array_4<T>& p) {
		this->elems[0] = p.elems[0];
		this->elems[1] = p.elems[1];
		this->elems[2] = p.elems[2];
		this->elems[3] = p.elems[3];
	}

	//operator=
	array_4<T>& operator=(const array_4<T>& rhs) {
		if (&rhs == this) {
			return *this;
		} else {
			this->elems[0] = rhs.elems[0];
			this->elems[1] = rhs.elems[1];
			this->elems[2] = rhs.elems[2];
			this->elems[3] = rhs.elems[3];
		}
		return *this;
	}

	void assign(const T& value) {
		fill(value);
	}    // A synonym for fill
	void fill(const T& value) {
		this->elems[0] = value;
		this->elems[1] = value;
		this->elems[2] = value;
		this->elems[3] = value;
	}
};

}

#endif
