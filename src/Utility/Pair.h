//
//  Pair.h
//  LarusX_0_1
//
//  Created by zhou on 27/12/14.
//  Copyright (c) 2014 zhou. All rights reserved.
//

#ifndef Pair_h
#define Pair_h
#include <stddef.h>
#include "../TypeDef.h"

namespace Larus {
template<class T, class V>
class Pair: public ObjectBase {
public:
	T first;
	V second;
	Pair();
	Pair(const T&, const V&);
	Pair(const Pair<T, V>&);
	Pair<T, V>& operator=(const Pair<T, V>&);
	void reconstruct(const T&, const V&);
};

template<class T, class V>
Pair<T, V>::Pair() :
		first(), second() {
}

template<class T, class V>
Pair<T, V>::Pair(const T& _t, const V& _v) :
		first(_t), second(_v) {
}

template<class T, class V>
Pair<T, V>::Pair(const Pair<T, V>& p) {
	this->first = p.first;
	this->second = p.second;
}

template<class T, class V>
Pair<T, V>& Pair<T, V>::operator=(const Pair<T, V>& rhs) {
	if (&rhs == this) {
		return *this;
	} else {
		this->first = rhs.first;
		this->second = rhs.second;
	}
	return *this;
}
template<class T, class V>
void Pair<T, V>::reconstruct(const T& t, const V& v) {
	this->first = t;
	this->second = v;
}

/// Returns a pair object with (a,b)
template<typename T1, typename T2>
inline Pair<T1, T2> makePair(const T1& a, const T2& b) {
	return (Pair<T1, T2>(a, b));
}

}

#endif
