/************************
 //  \file   Arithmetic.cpp
 //  \brief
 // 
 //  \author zhou
 //  \date   28 avr. 2014 
 ***********************/

#include <stddef.h>
#include <assert.h>
#include <math.h>
#include "Arithmetic.h"
#include <iostream>

namespace Larus {

Float calDiscreteErrorL2(const arrayList &ao, const arrayList &ap) {
	assert(ao.size() == ap.size());
	int n = ao.size();
	Float errsum = 0.0;
	for (int i = 0; i < n; i++) {
		errsum += (ap[i] - ao[i]) * (ap[i] - ao[i]);
	}
	return sqrt(errsum);
}
Float calDiscreteErrorL1(const arrayList &ao, const arrayList &ap) {
	assert(ao.size() == ap.size());
	int n = ao.size();
	Float errsum = 0.0;
	for (int i = 0; i < n; i++) {
		errsum += abs(ap[i] - ao[i]);
	}
	return errsum;
}

Float calDiscreteErrorLinf(const arrayList &ao, const arrayList &ap) {
	assert(ao.size() == ap.size());
	int n = ao.size();
	arrayList aabs(n);
	for (int i = 0; i < n; i++) {
		aabs[i] = ABS(ap[i] - ao[i]);
	}
	return aabs.findMax();
}


Float randf(Float &r1, Float &r2) {
	assert(r1 != r2);
	Float rnum1 = rand() % 100;
	Float rmax = max(r1, r2);
	Float rmin = min(r1, r2);
	return rmin + (rmax - rmin) / 100 * rnum1;
}

Float cal_weighted_arithmetic_mean(const arrayList &a, const arrayList &w) {
	assert(a.size() == w.size());
	int len = a.size();
	double sumup = 0.0;
	double sumdown = 0.0;
	for (int i = 0; i < len; i++) {
		sumup += w[i] * a[i];
		sumdown += w[i];
	}
	return sumup / sumdown;
}

Float cal_sqare_sum(Float &a, Float &b) {
	return a * a + b * b;
}

Float cal_sum_sqare(Float &a, Float &b) {
	return (a + b) * (a + b);
}

// round to n digit of a
// roundto(3.145,1)=3.1    roundto(3.145,2)=3.15
Float roundto(Float a, int n) {
	return round(a * pow(10, n)) / pow(10, n);
}

int count_significance_digit(Float a) {
	for (int i = 0; i < 30; i++) {
		if (roundto(a, i) == a) {
			return i;
		}
	}
	return -1;
}

Float volume_of_circular_truncated_cone(Float R, Float r, Float h) {
	return LarusDef::PI * h * (R * R + r * r + R * r) / 3;
}

}
