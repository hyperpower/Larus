/************************
 //  \file   Taylor.h
 //  \brief
 // 
 //  \author czhou
 //  \date   24 nov. 2014 
 ***********************/
#ifndef TAYLOR_H_
#define TAYLOR_H_

#include "../TypeDef.h"
#include "Dimensionless.h"

namespace Larus{
//Flavia Viana, Raimundo Pardo, Rodolfo Yanez, Jose L. Trallero, and
//Daniel D. Joseph. Universal correlation for the rise velocity of long gas
//bubbles in round pipes. Journal of Fluid Mechanics, 494:379–398, 2003.
Float cal_Fr_Viana2003(Float Eo, Float R);

Float cal_R(Float D, Float g, Float rho_l, Float rho_g, Float mu_l);

Float cal_Nf_from_Eo_Mo(Float,Float);
Float cal_Eo_from_Nf_Mo(Float,Float);
Float cal_Mo_from_Nf_Eo(Float,Float);

Float cal_Fr_10r200(Float, Float);
Float cal_Fr_r10(Float Eo, Float R);
Float cal_Fr_200r(Float Eo);

Float cal_Fr_Wallis1969(Float Re, Float Eo);

Float Brown_1965_2(Float sita, Float v, Float g, Float R, Float Utb);
Float d_Brown_1965_2(Float sita, Float v, Float g, Float R, Float Utb);

Float newtonMethod(int n, Float assum, Float eps, Float (*f)(Float),
		Float (*fd)(Float));
}



#endif /* TAYLOR_H_ */
