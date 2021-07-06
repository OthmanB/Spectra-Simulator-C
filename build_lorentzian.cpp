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
#include "external/integrate/activity.h"

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


VectorXd build_l_mode_a1a2a3(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s, const double a2, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V){
/*
 * This model includes:
 *      - Asymetry of Lorentzian asym
 *      - splitting a1
 *      - an Asphericity parameter eta
 *      - latitudinal effect a3
*/
    const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
    double Qlm, clm, a2_terms;

    result.setZero();
    for(int m=-l; m<=l; m++){
        if(l != 0){
            //Qlm=(l*(l+1) - 3*pow(m,2))/((2*l - 1)*(2*l + 3)); // accounting for eta
            if(l == 1){
                clm=m; // a3 for l=1
                a2_terms=(3*m*m - 2)*a2;  // From Takashi note and Pnl decomposition: c2(n,l) = [3m*m - l(l+1)] / (2l-1)
            }
            if(l == 2){
                clm=(5*pow(m,3) - 17*m)/3.; // a3 for l=2
                a2_terms=(m*m -2)*a2;
            }
            if(l == 3){
                clm=0; // a3 NOT YET IMPLEMENTED FOR l=3
                a2_terms=(3*m*m - 12)*a2/5;
            }
            profile=(x_l - tmp.setConstant(fc_l + m*f_s + a2_terms + clm*a3)).array().square();
            profile=4*profile/pow(gamma_l,2);
        } else{
            profile=(x_l - tmp.setConstant(fc_l)).array().square();
            profile=4*profile/pow(gamma_l,2);
        }
        if(asym == 0){ //Model with no asymetry
            result=result+ H_l*V(m+l)* ((tmp.setConstant(1) + profile)).cwiseInverse();
        } else{
            tmp.setConstant(1);
            asymetry=(tmp + asym*(x_l/fc_l - tmp)).array().square() + (tmp2.setConstant(0.5*gamma_l*asym/fc_l)).array().square();
            result=result+ H_l*V(m+l)*asymetry.cwiseProduct(((tmp.setConstant(1) + profile)).cwiseInverse());
        }
    }

return result;
}


VectorXd build_l_mode_a1etaAlma3(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s, 
    const double eta0, const double epsilon_nl, const VectorXd& thetas, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V){
/*
 * This model includes:
 *      - Asymetry of Lorentzian asym
 *      - splitting a1
 *      - an Asphericity term eta (Centrifugal effect) and due to Active region following Gizon 2002, AN, 323, 251. 
 *			Note that thetas[0] is the mid latitude of the active region and thetas[1] is the angulare size (called delta in the activity.cpp code that contain Alm())
 *      - latitudinal effect a3
 */
    const long Nxl=x_l.size();
    const long double Dnl=0.75;
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
    double Qlm, clm, G, CF_term, AR_term;

    result.setZero();
    for(int m=-l; m<=l; m++){
        //G=Glm(l, m, thetas[0], thetas[1]); // WARNING WARNING WARNING: HERE WE APPLY ALM ONLY TO l>0 BUT MAY BE NOT CORRECT. NOTE: l=0,m=0 might only be shifted
        if(l != 0){
            Qlm=(l*(l+1) - 3*pow(m,2))/((2*l - 1)*(2*l + 3)); // accounting for eta ... be carefull as dnl=2/3 is not accounted for here (see papini-gizon 2020)
            Qlm=Qlm*2./3.; // Adding Dnl=2/3, see Papini&Gizon. Beware: Old version of IDL postMCMC tool is not compatible anymore (as Dnl was added there)
	    if(l == 1){
                clm=m; // a3 for l=1
            }
            if(l == 2){
                clm=(5*pow(m,3) - 17*m)/3.; // a3 for l=2
            }
            if(l == 3){
                clm=(pow(m,3)-7*m)/2; // a3 implemented on 30/04/2021
            }
            CF_term=eta0*Dnl*pow(f_s*1e-6,2)*Qlm; //(4./3.)*pi*Dnl*pow(a1*1e-6,2.)/(rho*G);
            AR_term=epsilon_nl*Alm(l, m, thetas[0]*M_PI/180., thetas[1]*M_PI/180., "gate");
            //std::cout << thetas << std::endl;
            //exit(EXIT_SUCCESS);
            //std::cout << "(" << l << "," << m << ") : " << "d_CF=" << CF_term*fc_l  << "            d_AR=" << AR_term*fc_l  << "         m.a1=" << m*f_s << std::endl;
            profile=(x_l - tmp.setConstant(fc_l*(1. + CF_term + AR_term) + m*f_s + clm*a3)).array().square();
            profile=4*profile/pow(gamma_l,2);
        } else{
            profile=(x_l - tmp.setConstant(fc_l)).array().square();
            profile=4*profile/pow(gamma_l,2);
        }
        if(asym == 0){ //Model with no asymetry
            result=result+ H_l*V(m+l)* ((tmp.setConstant(1) + profile)).cwiseInverse();
        } else{
            tmp.setConstant(1);
            asymetry=(tmp + asym*(x_l/fc_l - tmp)).array().square() + (tmp2.setConstant(0.5*gamma_l*asym/fc_l)).array().square();
            result=result+ H_l*V(m+l)*asymetry.cwiseProduct(((tmp.setConstant(1) + profile)).cwiseInverse());
        }
    }
    /*std::cout << "------" << std::endl;
    std::cout << "XXXXXX" << std::endl;
    std::cout << "------" << std::endl;
    */
return result;
}
