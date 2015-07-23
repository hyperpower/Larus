/*
 * ts_define.h
 *
 *  Created on: May 30, 2015
 *      Author: zhou
 */

#ifndef TS_TS_DEFINE_H_
#define TS_TS_DEFINE_H_

#include "stdio.h"
namespace LarusTS {

typedef double Float;



enum Intersect {
	TS_OUT = -1, TS_ON = 0, TS_IN = 1
};

#define _return_val_if_fail(expr,val)  {       \
		if (!(expr))                           \
			return (val);                      \
	    }

#define _return_if_fail(expr)  {               \
		if (!(expr))                           \
			return ;                           \
	    }

template<class V>
inline int SIGN(const V& x) {
	return ((x) > 0. ? 1 : -1);
}
template<class V>
inline int ORIENT1D(const V& a,const V& b){
	return ((a) > (b) ? 1 : (a) < (b) ? -1 : 0);
}

template<class V>
static int sortp(V p[], std::size_t n) {
	int sign = 1;
	std::size_t i, j;

	for (i = 0; i < n - 1; ++i){
		for (j = 0; j < n - 1 - i; ++j){
			if (long(p[j + 1]) < long(p[j])) {
				V tmp = p[j];

				p[j] = p[j + 1];
				p[j + 1] = tmp;
				sign = -sign;
			}
		}
	}
	return sign;
}

}

#endif /* TS_TS_DEFINE_H_ */
