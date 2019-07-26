/*
 * build_lorentzian.cpp
 *
 *  Created on: 22 Feb 2016
 *      Author: obenomar
 */
#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>

using Eigen::VectorXd;

VectorXd build_l_mode_a1etaa3_simu(const VectorXd x_l, double H_l, double fc_l, double f_s, double eta, double a3, double gamma_l, const int l, VectorXd V){
/*
 * This model includes:
 *      - splitting a1
 *      - a second order effect eta, independent from the frequency (identified to the centrifugal force if no magnetic field effect)
 *      - latitudinal effect a3
*/

	const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), result(Nxl);
	double Qlm, clm;

	result.setZero();
	for(int m=-l; m<=l; m++){
		if(l != 0){
			Qlm=(l*(l+1) - 3*pow(m,2))/((2*l - 1)*(2*l + 3)); // accounting for a2
			if(l == 1){
				clm=m; // a3 for l=1
			}
			if(l == 2){
				clm=(5*pow(m,3) - 17*m)/3.; // a3 for l=2
			}
			if(l == 3){
				clm=0; // a3 NOT YET IMPLEMENTED FOR l=3
			}
			profile=(x_l - tmp.setConstant(fc_l*(1. + eta*Qlm) + m*f_s + clm*a3)).array().square();
			profile=4*profile/pow(gamma_l,2);
		} else{
			profile=(x_l - tmp.setConstant(fc_l)).array().square();
			profile=4*profile/pow(gamma_l,2);
		}
		result=result+ H_l*V(m+l)* ((tmp.setConstant(1) + profile)).cwiseInverse();
	}

return result;
}


VectorXd build_l_mode_act_simu(const VectorXd x_l, double H_l, double fc_l, double f_s, double eta, double a3, double b, double alpha, double asym, double gamma_l, const int l, VectorXd V){
/*
 * This model includes:
 *      - Asymetry of Lorentzian asym
 *      - splitting a1
 *      - centrifugal force effect eta
 *      - latitudinal effect a3
 *      - effect of magnetic field of the form b.nu^alpha
*/
	const long Nxl=x_l.size();
        VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
	double Qlm, clm;

    //std::cout << "           Begin build mode" << std::endl;
    
	result.setZero();
	for(int m=-l; m<=l; m++){
		if(l != 0){
			Qlm=(l*(l+1) - 3*pow(m,2))/((2*l - 1)*(2*l + 3)); // accounting for a2
			if(l == 1){
				clm=m; // a3 for l=1
			}
			if(l == 2){
				clm=(5*pow(m,3) - 17*m)/3.; // a3 for l=2
			}
			if(l == 3){
				clm=0; // a3 NOT YET IMPLEMENTED FOR l=3
			}
		
			profile=(x_l - tmp.setConstant(fc_l + Qlm*(eta*fc_l*pow(f_s,2) + b*pow(fc_l*1e-3,alpha)) + m*f_s + clm*a3) ).array().square();
			profile=4*profile/pow(gamma_l,2);
		} else{
			profile=(x_l - tmp.setConstant(fc_l)).array().square();
			profile=4*profile/pow(gamma_l,2);
		}
		tmp.setConstant(1);
		asymetry=(tmp + asym*(x_l/fc_l - tmp)).array().square() + (tmp2.setConstant(0.5*gamma_l*asym/fc_l)).array().square();
		//std::cout << "asymetry=" << asymetry.transpose() <<std::endl;
		//tmp2=((tmp.setConstant(1) + profile)).cwiseInverse();
		result=result+ H_l*V(m+l)*asymetry.cwiseProduct(((tmp.setConstant(1) + profile)).cwiseInverse());
		//std::cout << "in build_lorentzian" <<std::endl;
	}

    //std::cout << "           Exit build mode" << std::endl;
    
//std::cout << "result=" << result.transpose() <<std::endl;
return result;
}
