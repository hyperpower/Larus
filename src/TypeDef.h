/*
 * TypeDef.h
 *
 *  Created on: 18 janv. 2014
 *      Author: zhou
 */

#ifndef TYPEDEF_H_
#define TYPEDEF_H_
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>

//#include <omp.h>

//#define _OPENMP_
#ifdef _OPENMP_

#define OPENMP_Num_Threads 4
#endif
namespace LarusDef {

const double PI = 3.141592653589793238;
typedef int size_type;
}

namespace Larus {

#define _USE_SIMD_

#define NULL_PTR (NULL)
#define ASSERT(expr) assert(expr)
#define ASSERT_MSG(expr, msg) assert((expr)&&(msg))
#define CAST(type, p)      ((type)p)
#define _IF_TRUE_RETRUN(expr)  if(expr){return;};
#define _IF_FALSE_RETRUN(expr)  if(false==(expr)){return;};

typedef LarusDef::size_type idx_t;

class ObjectBase {
public:
	virtual ~ObjectBase() {
	}
	;
};

class DataBase {
public:
	virtual ~DataBase() {
	}
	;
};

typedef DataBase* pDataBase;

class GeoBase {
public:
	virtual ~GeoBase() {
	}
	;
};

typedef GeoBase* pGeoBase;

typedef double Float;

typedef void* utPointer;
typedef const void* const_utPointer;

const Float SMALL = 1e-15;
const Float PINF = 1e255;

template<class TYPE>
inline TYPE MAX(TYPE a, TYPE b) {
	return (a >= b) ? a : b;
}
template<class TYPE>
inline TYPE MAX(TYPE a, TYPE b, TYPE c) {
	return MAX(MAX(a, b), c);
}
template<class TYPE>
inline TYPE MAX(TYPE a, TYPE b, TYPE c, TYPE d) {
	return MAX(MAX(a, b), MAX(c, d));
}
template<class TYPE>
inline TYPE MAX(TYPE a, TYPE b, TYPE c, TYPE d,TYPE a2, TYPE b2, TYPE c2, TYPE d2) {
	return MAX(MAX(a, b,c, d), MAX(a2,b2,c2,d2));
}
template<class TYPE>
inline TYPE MIN(TYPE a, TYPE b) {
	return (a <= b) ? a : b;
}
template<class TYPE>
inline TYPE ABS(TYPE a) {
	return (a <= 0.0) ? -a : a;
}

bool isEqual(Float a, Float b);
bool isZero(Float a);

Float get_wall_time();

Float get_cpu_time();

}

#endif /* TYPEDEF_H_ */
