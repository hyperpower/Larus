/************************
 //  \file   Dimensionless.h
 //  \brief
 // 
 //  \author czhou
 //  \date   24 nov. 2014 
 ***********************/
#ifndef DIMENSIONLESS_H_
#define DIMENSIONLESS_H_

//Eotvos number
#include "../TypeDef.h"
#include <math.h>
#include <iostream>

namespace Larus {

//Eo group
struct Eo_group {
	Float Eo;
	Float rho_l;
	Float rho_g;
	Float D;
	Float sigma;
	Float g;

	void set_group(Float Eo, Float rho_l, Float rho_g, Float D, Float sigma,
			Float g) {
		this->Eo = Eo;
		this->rho_l = rho_l;
		this->rho_g = rho_g;
		this->D = D;
		this->sigma = sigma;
		this->g = g;
	}

	Float cal_Eo() {
		Eo = g * (rho_l - rho_g) * D * D / sigma;
		return Eo;
	}

	Float cal_g() {
		g = Eo * sigma / (rho_l - rho_g) / D / D;
		return g;
	}

	Float cal_sigma() {
		sigma = g * (rho_l - rho_g) * D * D / Eo;
		return sigma;
	}

	Float cal_D() {
		D = sqrt(Eo * sigma / g / (rho_l - rho_g));
		return D;
	}

	void show() {
		std::cout << "    g * D^2 * (rho_l-rho_g)  \n";
		std::cout << "Eo =-------------------------\n";
		std::cout << "             sigma         \n";
		std::cout << "Eo    = " << Eo << std::endl;
		std::cout << "rho_l = " << rho_l << std::endl;
		std::cout << "rho_g = " << rho_g << std::endl;
		std::cout << "D     = " << D << std::endl;
		std::cout << "sigma = " << sigma << std::endl;
		std::cout << "g     = " << g << std::endl;
	}

};

struct Mo_group {
	Float Mo;
	Float rho_l;
	Float rho_g;
	Float mu_l;
	Float sigma;
	Float g;

	void set_group(Float Mo, Float rho_l, Float rho_g, Float mu_l, Float sigma,
			Float g) {
		this->Mo = Mo;
		this->rho_l = rho_l;
		this->rho_g = rho_g;
		this->mu_l = mu_l;
		this->sigma = sigma;
		this->g = g;
	}

	Float cal_Mo() {
		this->Mo = g * mu_l * mu_l * mu_l * mu_l * (rho_l - rho_g) / rho_l
				/ rho_l / pow(sigma, 3.0);
		return this->Mo;
	}

	Float cal_sigma() {
		this->sigma = pow(
				(g * mu_l * mu_l * mu_l * mu_l * (rho_l - rho_g) / rho_l / rho_l
						/ Mo), 1.0 / 3.0);
		return this->sigma;
	}

	void show() {
		std::cout << "    g * mu_l^4 * (rho_l-rho_g)  \n";
		std::cout << "Mo =----------------------------\n";
		std::cout << "       rho_l^2 * sigma^3        \n";
		std::cout << "Mo    = " << Mo << std::endl;
		std::cout << "rho_l = " << rho_l << std::endl;
		std::cout << "rho_g = " << rho_g << std::endl;
		std::cout << "mu_l  = " << mu_l << std::endl;
		std::cout << "sigma = " << sigma << std::endl;
		std::cout << "g     = " << g << std::endl;
	}
};

struct Fr_group_0 {
	Float Fr;
	Float u;
	Float D;
	Float g;

	Float cal_Fr() {
		this->Fr = u / sqrt(g * D);
		return this->Fr;
	}

	void set_group(Float Fr, Float u, Float D, Float g) {
		this->Fr = Fr;
		this->u = u;
		this->D = D;
		this->g = g;
	}

	Float cal_u(){
		this->u=Fr*sqrt(g*D);
		return this->D;
	}

	void show() {
			std::cout << "              u  \n";
			std::cout << "Fr =-----------------------\n";
			std::cout << "       (g * D)^(1/2)       \n";
			std::cout << "Fr    = " << Fr << std::endl;
			std::cout << "u     = " << u << std::endl;
			std::cout << "D     = " << D << std::endl;
			std::cout << "g     = " << g << std::endl;
		}

};

struct Fr_group_1 {
	Float Fr;
	Float rho_l;
	Float rho_g;
	Float u;
	Float D;
	Float g;

	Float cal_Fr() {
		Fr = u / sqrt((rho_l - rho_g) / rho_l * g * D);
		return Fr;
	}

	void set_group(Float Fr, Float rho_l, Float rho_g, Float u, Float D,
			Float g) {
		this->Fr = Fr;
		this->rho_l = rho_l;
		this->rho_g = rho_g;
		this->u = u;
		this->D = D;
		this->g = g;
	}
};

struct Re_group_1 {
	Float Re;
	Float rho_l;
	Float rho_g;
	Float mu_l;
	Float D;
	Float g;

	Float cal_Re() {
		Re = sqrt(pow(D, 3.0) * g * (rho_l - rho_g) * rho_l) / mu_l;
		return Re;
	}
};

inline Float cal_inverse_viscosity_number(Float rho_l, Float mu_l, Float g,
		Float D) {
	return rho_l * sqrt(g * D * D * D) / mu_l;
}

inline Float cal_inverse_viscosity_number(Float Eo, Float Mo) {
	return pow(Eo * Eo * Eo / Mo, 0.25);
}

}

#endif /* DIMENSIONLESS_H_ */
