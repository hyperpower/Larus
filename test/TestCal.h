/************************
 //  \file   TestCal.h
 //  \brief
 // 
 //  \author czhou
 //  \date   14 oct. 2014 
 ***********************/
#ifndef TESTCAL_H_
#define TESTCAL_H_

#include "TestCal.h"
#include "../src/Algebra/Arithmetic.h"
#include "../src/Correlation/Dimensionless.h"
#include "../src/Correlation/Taylor.h"
#include <sstream>
#include <math.h>
#include "../src/IO/IO_tablefile.h"

namespace Larus {

//Wallis (1969)

Float cal_Buoyancy_Reynolds(Float D, Float g, Float rhol, Float rhog,
		Float mul) {
	return pow(pow(D, 3) * g * (rhol - rhog) * mul, 0.5) / mul;
}

Float cal_Fr(Float D, Float g, Float Re, Float Eo, Float m) {
	return 0.345 * sqrt(g * D) * (1 - exp(-0.01 * Re / 0.345))
			* (1 - exp((3.37 - Eo) / m));
}

void test_correlation() {
	cout << "test correlation\n";
	Float rho_l = 1000;
	Float rho_g = 1;
	Float mu_l = 10;
	//Float mu_g = 0.1;
	//Float D=1;
	Float g = 1;
	Float moset = 1;
	int bi = 0;
	int be = 1000;
	Float eoset = 0;
	Float eostep = 1;
	arrayList amo(be);
	arrayList aeo(be);

	arrayList afr(be);
	arrayList afr2(be);
	for (int i = bi; i < be; i++) {
		amo[i] = moset;
		aeo[i] = eoset + i * eostep;
	}
	aeo.show();
	for (int i = 0; i < be; i++) {
		Mo_group mog;
		mog.set_group(amo[i], rho_l, rho_g, mu_l, 0.0, g);
		mog.cal_sigma();
		//mog.show();
		Eo_group eog;
		eog.set_group(aeo[i], rho_l, rho_g, 0.0, mog.sigma, g);
		eog.cal_D();
		//eog.show();
		Float R = cal_R(eog.D, g, rho_l, rho_g, mu_l);
		//cout<<R<<endl;
		afr[i] = cal_Fr_Viana2003(eog.Eo, R);
		afr2[i] = cal_Fr_Wallis1969(R, eog.Eo);
	}
	//cout<<cal_Fr_10r200(60, 200)<<endl;
	//cout<<cal_Fr_200r(60)<<endl;
	//cout<<pow(14.793,3.06)<<endl;
	ostringstream sstream;
	sstream.precision(2);
	sstream.width(6);
	sstream << std::scientific << moset;
	string filename = "mo-" + sstream.str() + ".txt";
	TableFile tfile(filename);
	tfile.val.appendCol(amo);
	tfile.val.appendCol(aeo);
	tfile.val.appendCol(afr);
	tfile.val.appendCol(afr2);
	tfile.output();
	//cout << "Property:  \n";
	//cout << "rho_l " << rho_l << endl;
	//cout << "rho_g " << rho_g << endl;
	//cout << "mu_l  " << mu_l << endl;
	//cout << "mu_g  " << mu_g << endl;
	//cout << "Mo    " << mog.cal_Mo() << endl;
	//cout << "R     " << R << endl;
	//cout << eog.cal_Eo() << "  " << FrV << endl;
}

void test_initial_film() {
	cout << "test filmthickness\n";
	Float rho_l = 1000;
	Float rho_g = 1;
	Float mu_l = 0.001;

	Float g = 1;

	Float moset = 225;
	Float nfset = 183.7117307087;

	int bi = 0;
	int be = 1;
	arrayList amo(be);
	arrayList anf(be);
	arrayList aeo(be);

	arrayList afr(be);
	arrayList afr2(be);

	arrayList au(be);
	arrayList au2(be);

	arrayList asita(be);
	arrayList asita2(be);

	Float nfstep = 1;
	for (int i = bi; i < be; i++) {
		amo[i] = moset;

		anf[i] = nfset + i * nfstep;

		aeo[i] = cal_Eo_from_Nf_Mo(anf[i], amo[i]);

		afr[i] = cal_Fr_Viana2003(aeo[i], anf[i]);
		afr2[i] = cal_Fr_Wallis1969(anf[i], aeo[i]);

		Mo_group mog;
		mog.set_group(amo[i], rho_l, rho_g, mu_l, 0.0, g);
		mog.cal_sigma();
		mog.show();

		Eo_group eog;
		eog.set_group(aeo[i], rho_l, rho_g, 0.0, mog.sigma, g);
		eog.cal_D();
		eog.show();
		Fr_group_0 frg;
		frg.set_group(afr[i], 0.0, eog.D, g);
		frg.cal_u();

		Fr_group_0 frg2;
		frg2.set_group(afr2[i], 0.0, eog.D, g);
		frg2.cal_u();
		frg.show();

		Float a = 1;
		Float b = 0;
		Float c = 1.5 * mu_l / rho_l / g * frg.u;
		Float d = -1.5 * mu_l / rho_l / g * frg.u * 0.5 * frg.D;
		Float x1 = 0;
		Float x2 = 0;
		Float x3 = 0;
		//cout << "a = " << a << endl;
		//cout << "b = " << b << endl;
		//cout << "c = " << c << endl;
		//cout << "d = " << d << endl;

		solve_cubic_equation(a, b, c, d, x1, x2, x3);
		//cout << "cas = " << cas << endl;
		//cout << "x1 = " << x1 << endl;
		//cout << "x2 = " << x2 << endl;
		//cout << "x3 = " << x3 << endl;
		Float a2 = 1;
		Float b2 = 0;
		Float c2 = 1.5 * mu_l / rho_l / g * frg2.u;
		Float d2 = -1.5 * mu_l / rho_l / g * frg2.u * 0.5 * frg2.D;
		Float x12 = 0;
		Float x22 = 0;
		Float x32 = 0;
		//cout << "a = " << a << endl;
		//cout << "b = " << b << endl;
		//cout << "c = " << c << endl;
		//cout << "d = " << d << endl;

		solve_cubic_equation(a2, b2, c2, d2, x12, x22, x32);
		//cout << "cas = " << cas << endl;
		//cout << "x1 = " << x1 << endl;
		//cout << "x2 = " << x2 << endl;
		//cout << "x3 = " << x3 << endl;
		asita[i] = x1 / frg.D;
		asita2[i] = x12 / frg2.D;

		cout << "sita/D v = " << asita[i] << endl;
		cout << "sita/D w = " << asita2[i] << endl;
	}
}

void test_filmthick() {
	cout << "test filmthickness\n";
	Float rho_l = 1000;
	Float rho_g = 1;
	Float mu_l = 0.001;

	Float g = 1;

	Float moset = 0.328;
	Float nfset = 1;

	int bi = 0;
	int be = 1000;
	arrayList amo(be);
	arrayList anf(be);
	arrayList aeo(be);

	arrayList afr(be);
	arrayList afr2(be);

	arrayList au(be);
	arrayList au2(be);

	arrayList asita(be);
	arrayList asita2(be);

	Float nfstep = 1;
	for (int i = bi; i < be; i++) {
		amo[i] = moset;

		anf[i] = nfset + i * nfstep;

		aeo[i] = cal_Eo_from_Nf_Mo(anf[i], amo[i]);

		afr[i] = cal_Fr_Viana2003(aeo[i], anf[i]);
		afr2[i] = cal_Fr_Wallis1969(anf[i], aeo[i]);

		Mo_group mog;
		mog.set_group(amo[i], rho_l, rho_g, mu_l, 0.0, g);
		mog.cal_sigma();
		//mog.show();

		Eo_group eog;
		eog.set_group(aeo[i], rho_l, rho_g, 0.0, mog.sigma, g);
		eog.cal_D();
		//eog.show();
		Fr_group_0 frg;
		frg.set_group(afr[i], 0.0, eog.D, g);
		frg.cal_u();

		Fr_group_0 frg2;
		frg2.set_group(afr2[i], 0.0, eog.D, g);
		frg2.cal_u();
		//frg.show();

		Float a = 1;
		Float b = 0;
		Float c = 1.5 * mu_l / rho_l / g * frg.u;
		Float d = -1.5 * mu_l / rho_l / g * frg.u * 0.5 * frg.D;
		Float x1 = 0;
		Float x2 = 0;
		Float x3 = 0;
		//cout << "a = " << a << endl;
		//cout << "b = " << b << endl;
		//cout << "c = " << c << endl;
		//cout << "d = " << d << endl;

		solve_cubic_equation(a, b, c, d, x1, x2, x3);
		//cout << "cas = " << cas << endl;
		//cout << "x1 = " << x1 << endl;
		//cout << "x2 = " << x2 << endl;
		//cout << "x3 = " << x3 << endl;
		Float a2 = 1;
		Float b2 = 0;
		Float c2 = 1.5 * mu_l / rho_l / g * frg2.u;
		Float d2 = -1.5 * mu_l / rho_l / g * frg2.u * 0.5 * frg2.D;
		Float x12 = 0;
		Float x22 = 0;
		Float x32 = 0;
		//cout << "a = " << a << endl;
		//cout << "b = " << b << endl;
		//cout << "c = " << c << endl;
		//cout << "d = " << d << endl;

		solve_cubic_equation(a2, b2, c2, d2, x12, x22, x32);
		//cout << "cas = " << cas << endl;
		//cout << "x1 = " << x1 << endl;
		//cout << "x2 = " << x2 << endl;
		//cout << "x3 = " << x3 << endl;
		asita[i] = x1 / frg.D;
		asita2[i] = x12 / frg2.D;
	}
	ostringstream sstream;
	sstream.precision(2);
	sstream.width(6);
	sstream << std::scientific << moset;
	string filename = "mo-" + sstream.str() + "-film.txt";
	TableFile tfile(filename);
	tfile.val.appendCol(amo);
	tfile.val.appendCol(aeo);
	tfile.val.appendCol(anf);
	tfile.val.appendCol(afr);
	tfile.val.appendCol(afr2);
	tfile.val.appendCol(asita);
	tfile.val.appendCol(asita2);
	tfile.output();
}

}




#endif /* TESTCAL_H_ */
