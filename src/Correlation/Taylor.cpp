/************************
 //  \file   Taylor.cpp
 //  \brief
 // 
 //  \author czhou
 //  \date   24 nov. 2014 
 ***********************/

#include "Taylor.h"
#include <math.h>

namespace Larus {

Float cal_Fr_200r(Float Eo) {
	return 0.34 / pow((1.0 + 3805.15 / pow(Eo, 3.06)), 0.58);
}

Float Fun_L(Float R, Float A, Float B, Float C, Float G) {
	return A / pow(1 + pow(R / B, C), G);
}

Float cal_Fr_10r200(Float Eo, Float R) {
	const Float a = 0.34;
	const Float b = 14.793;
	const Float c = -3.06;
	const Float d = 0.58;
	const Float e = 31.08;
	const Float f = 29.868;
	const Float g = -1.96;
	const Float h = -0.49;
	const Float i = -1.45;
	const Float j = 24.867;
	const Float k = -9.93;
	const Float l = -0.094;
	const Float m = -1.0295;

	Float A = Fun_L(Eo, a, b, c, d);
	Float B = Fun_L(Eo, e, f, g, h);
	Float C = Fun_L(Eo, i, j, k, l);
	Float G = m / C;

	return Fun_L(R, A, B, C, G);
}

Float cal_Fr_r10(Float Eo, Float R) {
	return 0.009494 / pow(1 + 6197 / pow(Eo, 2.561), 0.5793) * pow(R, 1.026);
}

Float cal_R(Float D, Float g, Float rho_l, Float rho_g, Float mu_l) {
	return sqrt(pow(D, 3) * g * (rho_l - rho_g) * rho_l) / mu_l;
}

Float cal_Nf_from_Eo_Mo(Float Eo,Float Mo){
	return pow(pow(Eo,3)/Mo,0.25);
}
Float cal_Eo_from_Nf_Mo(Float Nf,Float Mo){
	return pow((pow(Nf,4)*Mo),1.0/3.0);
}
Float cal_Mo_from_Nf_Eo(Float Nf,Float Eo){
	return pow(Eo,3)/pow(Nf,4);
}

Float cal_Fr_Viana2003(Float Eo, Float R) {
	//if(R<=10){
	//	return cal_Fr_r10(Eo,R);
	//}else if(R>10 && R<200){
	return cal_Fr_10r200(Eo, R);
	//}else{
	//	return cal_Fr_10r200(Eo,R);
	//return cal_Fr_200r(Eo);
	//}
}

//Wallis (1969)
Float cal_m(Float Re) {
	if (Re > 250) {
		return 10;
	} else if (18 <= Re && Re <= 250) {
		return 69 * pow(Re, -0.35);
	} else {
		return 25;
	}
}

Float cal_Fr_Wallis1969(Float Re, Float Eo) {
	Float m = cal_m(Re);
	return 0.345 * (1 - exp(-0.01 * Re / 0.345)) * (1 - exp((3.37 - Eo) / m));
}

Float Brown_1965(Float sita, Float v, Float g, Float R, Float Utb) {
	return pow((3 * v * Utb * pow((R - sita), 2.0) / (2 * g * (R - sita))),
			1.0 / 3.0) - sita;
}

Float Brown_1965_2(Float sita, Float v, Float g, Float R, Float Utb) {
	return pow(sita, 3) + (3 * v * Utb / 2 / g) * sita
			- (3 * v * Utb / 2 / g) * R;
}

Float d_Brown_1965_2(Float sita, Float v, Float g, Float R, Float Utb) {
	return 3 * pow(sita, 2) + (3 * v * Utb / 2 / g);
}

Float newtonMethod(int n, Float assum, Float eps, Float (*f)(Float),
		Float (*fd)(Float)) {
	//int step = 0;
	Float x = assum;
	//Float Next = 0;
	Float A = f(x);
	Float B = fd(x);
	for (int i=0; i < n; i++) {
		if (abs(f(x)) < eps) {
			break;
		} else {
			x = x - A / B;
		}
	}
	return x;
}



}
