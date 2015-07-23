/************************
 //  \file   Arithmetic.h
 //  \brief
 // 
 //  \author czhou
 //  \date   12 févr. 2014 
 ***********************/
#ifndef ARITHMETIC_H_
#define ARITHMETIC_H_

#include "../TypeDef.h"

#include "../Utility/ArrayList.h"
#include <math.h>
#include <iostream>

namespace Larus {

Float calDiscreteErrorL2(const arrayList &ao, const arrayList &ap);
Float calDiscreteErrorL1(const arrayList &ao, const arrayList &ap);
Float calDiscreteErrorLinf(const arrayList &ao, const arrayList &ap);

Float cal_weighted_arithmetic_mean(const arrayList &a, const arrayList &w);

inline int StepFun(Float x) {
	return (x <= 0) ? 0 : 1;
}

inline int StepFun(int x) {
	return (x <= 0) ? 0 : 1;
}

inline int sign(Float x) {
	if (x < 0.0) {
		return -1;
	} else if (x > 0.0) {
		return 1;
	} else {
		return 0;
	}
}

enum TYPE_Range {
	Range_oo, Range_oc, Range_co, Range_cc,
};

template<typename TYPE>
inline bool isInRange(TYPE down, TYPE value, TYPE up, TYPE_Range range) {
	switch (range) {
	case Range_oo:
		return (down < value && value < up) ? true : false;
	case Range_oc:
		return (down < value && value <= up) ? true : false;
	case Range_co:
		return (down <= value && value < up) ? true : false;
	case Range_cc:
		return (down <= value && value <= up) ? true : false;
	}
	return false;
}

Float randf(Float r1, Float r2);

bool isEqual(Float a, Float b);

bool isZero(Float a);

// this function return a^2+b^2
Float cal_sqare_sum(Float &a, Float &b);

// this function return (a+b)*(a+b)
Float cal_sum_sqare(Float &a, Float &b);

Float roundto(Float, int);

int count_significance_digit(Float a);

Float volume_of_circular_truncated_cone(Float R, Float r, Float h);

template<class TYPE>
inline TYPE max(TYPE a, TYPE b, bool (*Comp_ge)(const TYPE&, const TYPE&)) {
	return Comp_ge(a, b) ? a : b;
}

template<class TYPE>
inline TYPE max(TYPE a, TYPE b, TYPE c,
		bool (*Comp_ge)(const TYPE&, const TYPE&)) {
	TYPE tmp = Comp_ge(a, b) ? a : b;
	return Comp_ge(tmp, c) ? tmp : c;
}

template<class TYPE>
inline TYPE min(TYPE a, TYPE b, bool (*Comp_le)(const TYPE&, const TYPE&)) {
	return Comp_le(a, b) ? a : b;
}

template<class TYPE>
inline TYPE min(TYPE a, TYPE b, TYPE c,
		bool (*Comp_le)(const TYPE&, const TYPE&)) {
	TYPE tmp = Comp_le(a, b) ? a : b;
	return Comp_le(tmp, c) ? tmp : c;
}

template<class TYPE>
inline TYPE mid(TYPE a, TYPE b, TYPE c,
		bool (*Comp_ge)(const TYPE&, const TYPE&)) {
	int idx = Comp_ge(a, b) ? 1 : 2;
	if (idx == 1)
		idx = Comp_ge(a, c) ? 1 : 3;
	else
		idx = Comp_ge(b, c) ? 2 : 3;
	if (idx == 1)
		return Comp_ge(b, c) ? b : c;
	else if (idx == 2)
		return Comp_ge(a, c) ? a : c;
	else
		return Comp_ge(a, b) ? a : b;
}

template<class TYPE>
inline void sort(const TYPE& a, const TYPE& b, const TYPE& c,  //
		bool (*Comp)(const TYPE&, const TYPE&), //
		TYPE& big, TYPE& mid, TYPE& small) {
	int idx = Comp(a, b) ? 1 : 2;
	if (idx == 1)
		idx = Comp(a, c) ? 1 : 3;
	else
		idx = Comp(b, c) ? 2 : 3;
	if (idx == 1) {
		big = a;
		if (Comp(b, c)) {
			mid = b;
			small = c;
		} else {
			mid = c;
			small = b;
		}
		return;
	}
	if (idx == 2) {
		big = b;
		if (Comp(a, c)) {
			mid = a;
			small = c;
		} else {
			mid = c;
			small = a;
		}
		return;
	}
	big = c;
	if (Comp(a, b)) {
		mid = a;
		small = b;
	} else {
		mid = b;
		small = a;
	}
}

template<class TYPE>
bool Comp_great(const TYPE& a, const TYPE& b) {
	return a > b;
}

template<class TYPE>
bool Comp_less(const TYPE& a, const TYPE& b) {
	return a < b;
}

template<class TYPE>
inline void sort(const TYPE& a, const TYPE& b, const TYPE& c,  //
		bool (*Comp)(const TYPE&, const TYPE&), //
		int& big, int& mid, int& small) {
	int idx = Comp(a, b) ? 0 : 1;
	if (idx == 0)
		idx = Comp(a, c) ? 0 : 2;
	else
		idx = Comp(b, c) ? 1 : 2;
	if (idx == 0) {
		big = 0;
		if (Comp(b, c)) {
			mid = 1;
			small = 2;
		} else {
			mid = 2;
			small = 1;
		}
		return;
	}
	if (idx == 1) {
		big = 1;
		if (Comp(a, c)) {
			mid = 0;
			small = 2;
		} else {
			mid = 2;
			small = 0;
		}
		return;
	}
	big = 2;
	if (Comp(a, b)) {
		mid = 0;
		small = 1;
	} else {
		mid = 1;
		small = 0;
	}
}
template<class TYPE>
inline void swap(TYPE& a, TYPE& b) //
		{
	TYPE tmp = a;
	a = b;
	b = tmp;
}

template<class TYPE>
inline void sort_incr(TYPE& a, TYPE& b, TYPE& c) //
		{
	if (b < a) {
		swap(a, b);
	}
	if (c < a) {
		swap(a, c);
	}
	if (c < b) {
		swap(b, c);
	}
}

template<class TYPE>
int _quadratic_discriminant(const TYPE& a, const TYPE& b, const TYPE& c,
		Float& discri) {
	discri = b * b - 4.0 * a * c;
	if (discri == 0) {
		return 1;
	} else if (discri > 0) {
		return 2;
	} else {
		return 0;
	}
}

template<class TYPE>
int solve_quadratic_equation(const TYPE& a, const TYPE& b, const TYPE& c,
		Float& x1, Float& x2) {
	Float discri = 0;
	int numroot = _quadratic_discriminant(a, b, c, discri);
	if (numroot == 2) {
		x1 = (-b - sqrt(discri)) / 2 / a;
		x2 = (-b + sqrt(discri)) / 2 / a;
		return 2;
	} else if (numroot == 1) {
		x1 = -b / 2 / a;
		x2 = x1;
		return 1;
	} else {
		return 0;
	}
}

template<class TYPE>
int solve_cubic_equation(const TYPE& a, const TYPE& b, const TYPE& c,
		const TYPE& d, Float& x1, Float& x2, Float& x3) {
	ASSERT(a != 0);
	Float A = b * b - 3.0 * a * c;
	Float B = b * c - 9.0 * a * d;
	Float C = c * c - 3.0 * b * d;
	Float discri = B * B - 4.0 * A * C;
	//case 1 has three equal real roots
	if (A == 0 && B == 0) {
		x1 = -b / 3.0 / a;
		x2 = x1;
		x3 = x1;
		return 1;
	}
	if (discri > 0) {
		Float Y1 = A * b + 1.5 * a * (-B + sqrt(discri));
		Float Y2 = A * b + 1.5 * a * (-B - sqrt(discri));
		Float cuberY1 = Y1 < 0 ? -pow(-Y1, 1.0 / 3.0) : pow(Y1, 1.0 / 3.0);
		Float cuberY2 = Y2 < 0 ? -pow(-Y2, 1.0 / 3.0) : pow(Y2, 1.0 / 3.0);
		x1 = (-b - cuberY1 - cuberY2) / 3.0 / a;
		//ignore complex roots
		x2 = x1;
		x3 = x1;
		return 2;
	}
	if (discri == 0) {
		Float K = B / A;
		x1 = -b / a + K;
		x2 = K / 2.0;
		x3 = x2;
		sort_incr(x1, x2, x3);
		return 3;
	}
	if (discri < 0) {
		Float T = (2.0 * A * b - 3.0 * a * B) / (2.0 * pow(A, 1.5));
		Float sita3 = acos(T) / 3.0;
		x1 = (-b - 2.0 * sqrt(A) * cos(sita3)) / (3.0 * a);
		x2 = (-b + sqrt(A) * (cos(sita3) + sqrt(3.0) * sin(sita3))) / (3.0 * a);
		x3 = (-b + sqrt(A) * (cos(sita3) - sqrt(3.0) * sin(sita3))) / (3.0 * a);
		sort_incr(x1, x2, x3);
		return 4;
	}
	return -1;
}

} //namespace Larus

#endif /* ARITHMETIC_H_ */
