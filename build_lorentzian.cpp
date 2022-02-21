/*
 * build_lorentzian.cpp
 *
 *  Created on: 22 Feb 2016
 *      Author: obenomar
 */
#include <math.h>
#include <Eigen/Dense>
#include "build_lorentzian.h"
#include <iostream>
#include <iomanip>
#include "external/integrate/activity.h"
#include "acoefs.h"
using Eigen::VectorXd;
using Eigen::VectorXi;



VectorXd build_l_mode_a1l_etaa3(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s1, const double f_s2, const double eta0, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V){
    /*
     * This model IS WITHOUT THE ASSUMPTION S11=S22. It includes:
     *      - Asymetry of Lorentzian asym
     *      - splitting a1(l=1) and a1(l=2). ASSUMES a1(l=3) =  (a1(1) + a1(2))/2.
     *      - an Asphericity parameter eta
     *      - latitudinal effect a3(l=2) only. We consider a3(l=3)=0
     */
    const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
    double f_s;
    
    result.setZero();
    for(int m=-l; m<=l; m++){
        if(l != 0){
            if(l == 1){
                f_s=f_s1;
            }
            if(l == 2){
                f_s=f_s2;
            }
            if(l == 3){
                f_s=(f_s1 + f_s2)/2.; // APPROXIMATION
            }
            profile=(x_l - tmp.setConstant(fc_l*(1. + eta0*pow(f_s*1e-6,2)*Qlm(l,m)) + m*f_s + Pslm(3,l,m)*a3)).array().square();
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

VectorXd build_l_mode_a1l_a2a3(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s1, const double f_s2, const double a2, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V){
    /*
     * This model IS WITHOUT THE ASSUMPTION S11=S22. It includes:
     *      - Asymetry of Lorentzian asym
     *      - splitting a1(l=1) and a1(l=2). ASSUMES a1(l=3) =  (a1(1) + a1(2))/2.
     *      - an Asphericity parameter eta
     *      - latitudinal effect a3(l=2) only. We consider a3(l=3)=0
     */
    const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
    double f_s, a2_terms;
    
    result.setZero();
    for(int m=-l; m<=l; m++){
        if(l != 0){
            a2_terms=Pslm(2,l,m)*a2;
            if(l == 1){
                //clm=0; // a3=0 for l=1 BY DEFINITION
                f_s=f_s1;
                //a2_terms=(3*m*m - 2)*a2;  // From Takashi note and Pnl decomposition: c2(n,l) = [3m*m - l(l+1)] / (2l-1)
            }
            if(l == 2){
                //clm=(5*pow(m,3) - 17*m)/3.; // a3 for l=2
                f_s=f_s2;
                //a2_terms=(m*m -2)*a2;
            }
            if(l == 3){
                //clm=(pow(m,3)-7*m)/2; // a3 implemented on 30/04/2021
                f_s=(f_s1 + f_s2)/2.; // APPROXIMATION
                //a2_terms=(3*m*m - 12)*a2/5;
            }
            profile=(x_l - tmp.setConstant(fc_l + m*f_s + a2_terms + Pslm(3,l,m)*a3)).array().square(); // a1=f_s , a2 and a3 coefficients
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

VectorXd build_l_mode_a1etaa3(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s, const double eta0, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V){
/*
 * This model includes:
 *      - Asymetry of Lorentzian asym
 *      - splitting a1
 *      - an Asphericity parameter eta
 *      - latitudinal effect a3
*/
	const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);

	result.setZero();
	for(int m=-l; m<=l; m++){
		if(l != 0){
            /*if(l == 1){
				clm=m; // a3 for l=1
			}
			if(l == 2){
				clm=(5*pow(m,3) - 17*m)/3.; // a3 for l=2
			}
			if(l == 3){
				clm=(pow(m,3)-7*m)/2; // a3 implemented on 30/04/2021
			}
            */
			profile=(x_l - tmp.setConstant(fc_l*(1. + eta0*pow(f_s*1e-6,2)*Qlm(l,m)) + m*f_s + Pslm(3,l,m)*a3)).array().square();
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
 *             Currently we use a hard-coded filter type "gate" which is rough but match the Gizon paper. "gauss" is also available and might be our final choice.
 *             Once we could compare the method adapted from Gizon works on global fits
 *      - latitudinal effect a3
 */
    const std::string filter_type="gate"; // The alternative is also "gauss" to have more smooth edges for the Activity effect
    const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
    double clm, G, CF_term, AR_term;

    result.setZero();
    for(int m=-l; m<=l; m++){
        //G=Glm(l, m, thetas[0], thetas[1]); // WARNING WARNING WARNING: HERE WE APPLY GLM ONLY TO l>0 BUT MAY BE NOT CORRECT. NOTE: l=0,m=0 might only be shifted
        if(l != 0){
            /*
            if(l == 1){
                clm=m; // a3 for l=1 WRONG ! ERROR NOTED ON 18/11/2021
            }
            if(l == 2){
                clm=(5*pow(m,3) - 17*m)/3.; // a3 for l=2
            }
            if(l == 3){
               clm=(pow(m,3)-7*m)/2; // a3 implemented on 30/04/2021
            }
            */
            CF_term=eta0*pow(f_s*1e-6,2)*Qlm(l,m); //(4./3.)*pi*pow(a1*1e-6,2.)/(rho*G);
            AR_term=epsilon_nl*Alm_deg(l, m, thetas[0], thetas[1], filter_type);
            /*std::cout << "(l,m) =  (" << l << "," << m << ")" << std::endl;
            std::cout << "thetas =" << thetas << std::endl;
            std::cout << "Alm =" << Alm_deg(l, m, thetas[0], thetas[1], filter_type) << std::endl;
            std::cout << "(" << l << "," << m << ") : " << "d_CF=" << CF_term*fc_l  << "            d_AR=" << AR_term*fc_l  << "         m.a1=" << m*f_s << std::endl;
            */
            profile=(x_l - tmp.setConstant(fc_l*(1. + CF_term + AR_term) + m*f_s + Pslm(3,l,m)*a3)).array().square();
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

VectorXd build_l_mode_a1a2a3(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s, const double a2, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V){
/*
 * This model includes:
 *      - Asymetry of Lorentzian asym
 *      - splitting a1
 *      - latitudinal effect a3
*/
    const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
    double clm, a2_terms;

    result.setZero();
    for(int m=-l; m<=l; m++){
        if(l != 0){
           clm=Pslm(3,l,m); // Changes made on 18/11/2021 : Use of acoefs.cpp
           a2_terms=Pslm(2,l,m)*a2;
           /*if(l == 1){
                clm=m; // a3 for l=1
                a2_terms=(3*m*m - 2)*a2;  // From Takashi note and Pnl decomposition: c2(n,l) = [3m*m - l(l+1)] / (2l-1)
            }
            if(l == 2){
                clm=(5*pow(m,3) - 17*m)/3.; // a3 for l=2
                a2_terms=(m*m -2)*a2;
            }
            if(l == 3){
                clm=(pow(m,3)-7*m)/2; // a3 implemented on 30/04/2021
                a2_terms=(3*m*m - 12)*a2/5;
            }
            */
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

VectorXd build_l_mode_act_simu(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s, const double eta, const double a3, 
        const double b, const double alpha, const double asym, const double gamma_l, const int l, const VectorXd& V){
/*
 * This model includes:
 *      - Asymetry of Lorentzian asym
 *      - splitting a1
 *      - centrifugal force effect eta
 *      - latitudinal effect a3
 *      - effect of magnetic field of the form b.nu^alpha
*/
    const double a2=0;
    const double a4=0;
    const double a5=0;
    const double a6=0;
    const double eta0=eta*1e-12; // Equivalence valid only due to change of notation overtime. Do not assume it to be true in other functions ! The 1e-12 comes from a1 (in Hz) in build_l_mode_aj()
    VectorXd result;

    std::cout << "Warning: b and alpha are ignored as this function is obselete" << std::endl;
    result=build_l_mode_aj(x_l, H_l, fc_l, f_s, a2, a3, a4, a5, a6, eta0, asym, gamma_l, l, V);
return result;
}

VectorXd build_l_mode_aj(const VectorXd& x_l, const double H_l, const double fc_l, 
        const double a1, const double a2, const double a3, const double a4, const double a5, const double a6, 
        const double eta0, const double asym, const double gamma_l, const int l, const VectorXd& V){
/*
 * This model includes:
 *      - Asymetry of Lorentzian asym
 *      - splittings in the form of a-coefficients aj with j={1,6}
 *      - eta: If eta0 > 0, account for the centrifugal distorsion in the a2 term (a2_CF) ==> This means the measured a2 = a2_AR : It WILL NOT include centrifugal effects
 *                         but the model DO account for it. In that case, eta0 =3./(4.*pi*rho*G) Should be the value given to that function: NOTE THAT Benomar2018 HAS A MISTAKE IN THAT FORMULATION (Eq. S6)
 *             If eta0 <=0 0, set a2_CF(eta) = 0 such that the measured a2 = a2_CF + a2_AR
*/
    const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
    double nu_nlm;

    result.setZero();
    for(int m=-l; m<=l; m++){
        if(l != 0){
            nu_nlm=fc_l + a1*Pslm(1,l,m) + a2*Pslm(2,l,m) + a3*Pslm(3,l,m) + a4*Pslm(4,l,m)+ a5*Pslm(5,l,m) + a6*Pslm(6,l,m);
            if (eta0 > 0){
                nu_nlm = nu_nlm + fc_l*eta0*Qlm(l,m)*pow(a1*1e-6,2);
            } 
            profile=(x_l - tmp.setConstant(nu_nlm)).array().square();
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

VectorXd build_l_mode_a1l_etaa3_v2(const VectorXd& x_l, const VectorXd& H_lm, const double fc_l, const double f_s1, const double f_s2, const double eta0, const double a3, const double asym, const double gamma_l, const int l){
    /*
     * This model IS WITHOUT THE ASSUMPTION S11=S22. It includes:
     *      - Asymetry of Lorentzian asym
     *      - splitting a1(l=1) and a1(l=2). ASSUMES a1(l=3) =  (a1(1) + a1(2))/2.
     *      - an Asphericity parameter eta
     *      - latitudinal effect a3(l=2) only. We consider a3(l=3)=0
     *		- Inclination IS NOT IMPOSED contrary to build_l_mode_a1l_etaa3. INSTEAD H_lm should be of dimension l(l+1) and provide all the heights 
     */
    const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
    double f_s;
    result.setZero();
    for(int m=-l; m<=l; m++){
        if(l != 0){
            switch (l){
                case 1:
                    f_s=f_s1;
                    break;
                case 2:
                    f_s=f_s2;
                    break;
                case 3:
                    f_s=(f_s1 + f_s2)/2.; // APPROXIMATION
                    break;
            }
            profile=(x_l - tmp.setConstant(fc_l*(1. + eta0*pow(f_s*1e-6,2)*Qlm(l,m)) + Pslm(1,l,m)*f_s + Pslm(3,l,m)*a3)).array().square();
            //profile=(x_l - tmp.setConstant(fc_l*(1. + eta*Qlm) + Pslm(1,l,m)*a1 + clm*a3)).array().square();
            profile=4*profile/pow(gamma_l,2);
        } else{
            profile=(x_l - tmp.setConstant(fc_l)).array().square();
            profile=4*profile/pow(gamma_l,2);
        }
        if(asym == 0){ //Model with no asymetry
            result=result+ H_lm(m+l)* ((tmp.setConstant(1) + profile)).cwiseInverse();
        } else{
            tmp.setConstant(1);
            asymetry=(tmp + asym*(x_l/fc_l - tmp)).array().square() + (tmp2.setConstant(0.5*gamma_l*asym/fc_l)).array().square();
            result=result+ H_lm(m+l)*asymetry.cwiseProduct(((tmp.setConstant(1) + profile)).cwiseInverse());
        }
    }
    
    std::cout << "NEED CHECKS in build_l_mode_a1l_etaa3_v2: The function was never verified" << std::endl;
    exit(EXIT_SUCCESS);

    return result;
}

VectorXd build_l_mode_a1etaa3_v2(const VectorXd& x_l, const VectorXd& H_lm, const double fc_l, const double f_s, const double eta0, const double a3, const double asym, const double gamma_l, const int l){
/*
 * This model includes:
 *      - Asymetry of Lorentzian asym
 *      - splitting a1
 *      - an Asphericity parameter eta
 *      - latitudinal effect a3
*/
	const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
	double clm;

	/*std::cout << " ---------- " << std::endl;
	std::cout << " l = " << l << std::endl;
	std::cout << "H_lm =" << H_lm << std::endl;
	std::cout << "fc_l =" << fc_l << std::endl;
	std::cout << "eta =" << eta << std::endl;
	std::cout << "a3 =" << a3 << std::endl;
	std::cout << "asym =" << asym << std::endl;
	std::cout << "gamma_l =" << gamma_l << std::endl;
	std::cout << " ---------- " << std::endl;
*/
	result.setZero();
	for(int m=-l; m<=l; m++){
		if(l != 0){
            //clm=Pslm(3,l,m); // Changes made on 18/11/2021 : Use of acoefs.cpp
            /*
            if(l == 1){
				clm=m; // a3 for l=1
			}
			if(l == 2){
				clm=(5*pow(m,3) - 17*m)/3.; // a3 for l=2
			}
			if(l == 3){
				clm=(pow(m,3)-7*m)/2; // a3 implemented on 30/04/2021
			}
            */
			profile=(x_l - tmp.setConstant(fc_l*(1. + eta0*pow(f_s*1e-6,2)*Qlm(l,m)) + m*f_s + Pslm(3,l,m)*a3)).array().square();
			profile=4*profile/pow(gamma_l,2);
		} else{
			profile=(x_l - tmp.setConstant(fc_l)).array().square();
			profile=4*profile/pow(gamma_l,2);
		}
		if(asym == 0){ //Model with no asymetry
			result=result+ H_lm(m+l)* ((tmp.setConstant(1) + profile)).cwiseInverse();
		} else{
			tmp.setConstant(1);
			asymetry=(tmp + asym*(x_l/fc_l - tmp)).array().square() + (tmp2.setConstant(0.5*gamma_l*asym/fc_l)).array().square();
			result=result+ H_lm(m+l)*asymetry.cwiseProduct(((tmp.setConstant(1) + profile)).cwiseInverse());
		}
	}

return result;
}

VectorXd build_l_mode_asym_act(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s, const double eta, const double a3, const double b, const double alpha, const double asym, const double gamma_l, const int l, const VectorXd& V){
/*
 * ---- OBSELETE FUNCTION ----
 * This model includes:
 *      - Asymetry of Lorentzian asym
 *      - splitting a1
 *      - centrifugal force effect eta.nu.a1^2 (rotation-induced oblateness)
 *      - latitudinal effect a3
 *      - effect of magnetic field of the form b.nu^alpha
 * ---------------------------
*/
    VectorXd result;

    std::cout << "Obselete Function" << std::endl;
    std::cout << "The program will exit now" << std::endl;
    exit(EXIT_SUCCESS);

//std::cout << "result=" << result.transpose() <<std::endl;
return result;
}

VectorXd optimum_lorentzian_calc_a1l_etaa3(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, const double f_s1, const double f_s2, const double eta0, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c){
    /*
     function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
     that contains the lorentzian model.
     BEWARE: USES build_l_mode_a1l_etaa3() ==> Asphericity is a linear term in nu
     */
    //const double c=20.;
    double pmin, pmax, f_s;
    VectorXi ivals;
    VectorXd m0, x_l, y_out(y.size());
    y_out=y;;
    switch(l){
        case 0:
            f_s=0.;
            break;
        case 1:
            f_s=f_s1;
            break;
        case 2:
            f_s=f_s2;
            break;
        case 3:
            f_s=(f_s1 + f_s2)/2.;
            break;
    }
    
    ivals=set_imin_imax(x, l, fc_l, gamma_l, f_s, c, step);    
    x_l=x.segment(ivals[0], ivals[1]-ivals[0]);
    
    m0=build_l_mode_a1l_etaa3(x_l, H_l, fc_l, f_s1, f_s2, eta0, a3, asym, gamma_l, l, V);
    //mall.setZero();
    //mall.segment(imin, imax-imin)=m0;
    y_out.segment(ivals[0], ivals[1]-ivals[0])= y_out.segment(ivals[0], ivals[1]-ivals[0]) + m0;
    return y_out;
}


VectorXd optimum_lorentzian_calc_a1l_a2a3(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, const double f_s1, const double f_s2, const double a2, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c){
    /*
     function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
     that contains the lorentzian model.
    BEWARE: USES build_l_mode_a1a2a3() ==> LINEAR DEPENDENCE OF Asphericity IN NU IS NOT ACCOUNTED FOR... THIS DEPENDENCE MAY BE IMPLEMENTED AT HIGHER LEVEL when calling this function
     */
    double pmin, pmax, f_s;
    VectorXi ivals;
    VectorXd m0, x_l, y_out(y.size());
    y_out=y;;
    switch(l){
        case 0:
            f_s=0.;
            break;
        case 1:
            f_s=f_s1;
            break;
        case 2:
            f_s=f_s2;
            break;
        case 3:
            f_s=(f_s1 + f_s2)/2.;
            break;
    }

    ivals=set_imin_imax(x, l, fc_l, gamma_l, f_s, c, step);    
    x_l=x.segment(ivals[0], ivals[1]-ivals[0]);

    m0=build_l_mode_a1l_a2a3(x_l, H_l, fc_l, f_s1, f_s2, a2, a3, asym, gamma_l, l, V);
    y_out.segment(ivals[0], ivals[1]-ivals[0])= y_out.segment(ivals[0], ivals[1]-ivals[0]) + m0;
    return y_out;
}


VectorXd optimum_lorentzian_calc_a1etaa3(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, const double f_s, const double eta0, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c){
/*
	function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
	that contains the lorentzian model.
	BEWARE: USES build_l_mode_a1etaa3() ==> Asphericity is a linear term in nu
*/
	//const double c=20.;
	double pmin, pmax;
	VectorXi ivals;
	VectorXd m0, x_l, y_out(y.size());
    y_out=y;

    ivals=set_imin_imax(x, l, fc_l, gamma_l, f_s, c, step);    
    x_l=x.segment(ivals[0], ivals[1]-ivals[0]);
 
	m0=build_l_mode_a1etaa3(x_l, H_l, fc_l, f_s, eta0, a3, asym, gamma_l, l, V);
    y_out.segment(ivals[0], ivals[1]-ivals[0])= y_out.segment(ivals[0], ivals[1]-ivals[0]) + m0;
return y_out;
}

VectorXd optimum_lorentzian_calc_a1etaAlma3(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, const double f_s, const double eta0, const double epsilon_nl, const VectorXd& thetas, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c){
/*
    function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
    that contains the lorentzian model.
    BEWARE: USES build_l_mode_a1etaGlma3() ==> Asphericity is a linear term in nu
*/
    //const double c=20.;
    double pmin, pmax;
    VectorXi ivals;
    VectorXd m0, x_l, y_out(y.size());
    y_out=y;
    //VectorXd mall(x.size());

    ivals=set_imin_imax(x, l, fc_l, gamma_l, f_s, c, step);    
    x_l=x.segment(ivals[0], ivals[1]-ivals[0]);
 
    m0=build_l_mode_a1etaAlma3(x_l, H_l, fc_l, f_s, eta0, epsilon_nl, thetas, a3, asym, gamma_l, l, V);
    //mall.setZero();
    //mall.segment(imin, imax-imin)=m0;
    y_out.segment(ivals[0], ivals[1]-ivals[0])= y_out.segment(ivals[0], ivals[1]-ivals[0]) + m0;
return y_out;
}


VectorXd optimum_lorentzian_calc_a1a2a3(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, const double f_s, const double a2, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c){
/*
    function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
    that contains the lorentzian model.
    BEWARE: USES build_l_mode_a1a2a3() ==> LINEAR DEPENDENCE OF Asphericity IN NU IS NOT ACCOUNTED FOR... THIS DEPENDENCE MAY BE IMPLEMENTED AT HIGHER LEVEL when calling this function
*/
    //const double c=20.;
    double pmin, pmax;
    VectorXi ivals;
    VectorXd m0, x_l, y_out(y.size());
    y_out=y;

    ivals=set_imin_imax(x, l, fc_l, gamma_l, f_s, c, step);    
    x_l=x.segment(ivals[0], ivals[1]-ivals[0]);

    m0=build_l_mode_a1a2a3(x_l, H_l, fc_l, f_s, a2, a3, asym, gamma_l, l, V);
    y_out.segment(ivals[0], ivals[1]-ivals[0])= y_out.segment(ivals[0], ivals[1]-ivals[0]) + m0;
return y_out;
}

VectorXd optimum_lorentzian_calc_aj(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, 
        const double a1, const double a2, const double a3, const double a4, const double a5, const double a6, 
        const double eta0, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c){
/*
    function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
    that contains the lorentzian model.
    BEWARE: USES build_l_mode_aj() ==> LINEAR DEPENDENCE OF Asphericity IN NU IS NOT ACCOUNTED FOR... THIS DEPENDENCE MAY BE IMPLEMENTED AT HIGHER LEVEL when calling this function
*/
    double pmin, pmax;
    VectorXi ivals;
    VectorXd m0, x_l, y_out(y.size());
    y_out=y;

    ivals=set_imin_imax(x, l, fc_l, gamma_l, a1, c, step);    
    x_l=x.segment(ivals[0], ivals[1]-ivals[0]);

    m0=build_l_mode_aj(x_l, H_l, fc_l, a1, a2, a3, a4, a5, a6, eta0, asym, gamma_l, l, V);
    y_out.segment(ivals[0], ivals[1]-ivals[0])= y_out.segment(ivals[0], ivals[1]-ivals[0]) + m0;
return y_out;
}


VectorXd optimum_lorentzian_calc_a1etaa3_v2(const VectorXd& x, const VectorXd& y, const VectorXd& H_lm, const double fc_l, const double f_s, const double eta0, const double a3, const double asym, const double gamma_l, const int l, const double step, const double c){
/*
	function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
	that contains the lorentzian model.
	BEWARE: USES build_l_mode_a1etaa3() ==> Asphericity is a linear term in nu
*/
	double pmin, pmax;
    VectorXi ivals;
	VectorXd m0, x_l, y_out(y.size());
    y_out=y;

    ivals=set_imin_imax(x, l, fc_l, gamma_l, f_s, c, step);    
    x_l=x.segment(ivals[0], ivals[1]-ivals[0]);

	m0=build_l_mode_a1etaa3_v2(x_l, H_lm, fc_l, f_s, eta0, a3, asym, gamma_l, l);
    y_out.segment(ivals[0], ivals[1]-ivals[0])= y_out.segment(ivals[0], ivals[1]-ivals[0]) + m0;
return y_out;
}


VectorXd optimum_lorentzian_calc_a1acta3(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, const double f_s, const double eta, const double a3, 
		const double b, const double alpha, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c){
/*
	function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
	that contains the lorentzian model.
	BEWARE: USES build_l_mode_asym_act() ==> Include mode asymetry and effect of activity

*/
    VectorXd y_out;
    
    std::cout << "Obselete Function" << std::endl;
    std::cout << "The program will exit now" << std::endl;
    exit(EXIT_SUCCESS);

return y_out;
}

VectorXd optimum_lorentzian_calc_a1l_etaa3_v2(const VectorXd& x, const VectorXd& y, const VectorXd& H_lm, const double fc_l, const double f_s1, const double f_s2, const double eta0, const double a3, const double asym, const double gamma_l, const int l, const double step, const double c){
    /*
     function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
     that contains the lorentzian model.
     BEWARE: USES build_l_mode_a1l_etaa3() ==> Asphericity is a linear term in nu

     This function differs from optimum_lorentzian_calc_a1l_etaa3 by the fact that it fits directly the (l,m) heights instead of considering H_l and V==ratios
     Thus an l=1 will have H_1m = [ H(m=-1), H(m=0), H(m=1)] components, etc...
     */
    double f_s;
    VectorXd m0, x_l, y_out(y.size());
    VectorXi ivals;
    y_out=y;;
    switch(l){
        case 0:
            f_s=0.;
            break;
        case 1:
            f_s=f_s1;
            break;
        case 2:
            f_s=f_s2;
            break;
        case 3:
            f_s=(f_s1 + f_s2)/2.;
            break;
    }
    ivals=set_imin_imax(x, l, fc_l, gamma_l, f_s, c, step);    
    x_l=x.segment(ivals[0], ivals[1]-ivals[0]);
    m0=build_l_mode_a1l_etaa3_v2(x_l, H_lm, fc_l, f_s1, f_s2, eta0, a3, asym, gamma_l, l);
    //mall.setZero();
    //mall.segment(imin, imax-imin)=m0;
    y_out.segment(ivals[0], ivals[1]-ivals[0])= y_out.segment(ivals[0], ivals[1]-ivals[0]) + m0;
    return y_out;
}


// ------ Common functions -----
double Qlm(const int l, const int m){
    // Qlm term that is used in the centrifugal distortion calculation
    // The final Qlm include the dnl=2/3 term that is sometime considered separately by some authors
    // See papini-gizon 2020 or Gizon 2002 for the definition
    const long double Dnl=2./3;
    double Qlm;
    Qlm=(l*(l+1) - 3*pow(m,2))/((2*l - 1)*(2*l + 3));
    Qlm=Qlm * Dnl;
    return Qlm;
}

VectorXi set_imin_imax(const VectorXd& x, const int l, const double fc_l, const double gamma_l, const double f_s, const double c, const double step){

    VectorXd pvals(2);
    VectorXi ivals(2);
    if(gamma_l >= 1 && f_s >= 1){
        if(l != 0){
            pvals[0]=fc_l - c*(l*f_s + gamma_l);
            pvals[1]=fc_l + c*(l*f_s + gamma_l);
        } else{
            pvals[0]=fc_l - c*gamma_l*2.2;
            pvals[1]=fc_l + c*gamma_l*2.2;
        }
    }
    if(gamma_l <= 1 && f_s >= 1){
        if(l !=0){
            pvals[0]=fc_l -c*(l*f_s + 1);
            pvals[1]=fc_l + c*(l*f_s + 1);
        } else{
            pvals[0]=fc_l - c*2.2;
            pvals[1]=fc_l + c*2.2;
        }
    }
    if(gamma_l >= 1 && f_s <= 1){
        if(l != 0){
            pvals[0]=fc_l - c*(l + gamma_l);
            pvals[1]=fc_l + c*(l + gamma_l);
        } else{
            pvals[0]=fc_l - c*2.2*gamma_l;
            pvals[1]=fc_l + c*2.2*gamma_l;
        }
    }
    if(gamma_l <= 1 && f_s <= 1){
        if(l !=0){
            pvals[0]=fc_l - c*(l+1);
            pvals[1]=fc_l + c*(l+1);
        } else{
            pvals[0]=fc_l -c*2.2;
            pvals[1]=fc_l +c*2.2;
        }
    }
    // ---------- Handling boundaries ----------
    //     case when the proposed values lead to (pmax - step) < min(x)....
    if( (pvals[1] - step) < x.head(1)(0)){
        pvals[1]=x.head(1)(0) +c;  // bundary = first value of x plus c
    }
    //     case when the proposed values lead to (pmin + step) >= max(x)...
    if( (pvals[0] + step) >= x.tail(1)(0)){
        pvals[0]=x.tail(1)(0) - c; // bundary = last value of x minus c
    }

    ivals[0]=floor((pvals[0]-x.head(1)(0))/step); // Here it is assumed that there is a regular grid WITHOUT WHOLES (CANNOT USE a .FILT file)
    ivals[1]=ceil((pvals[1]-x.head(1)(0))/step);
    
    if(ivals[0] < 0){ivals[0]=0;} // ensure that we always stay within the vector x
    if(ivals[1] > x.size()){ivals[1]=x.size();} // ensure that we always stay within the vector x
    if(ivals[1]-ivals[0] <= 0){
        std::cout << "Warning imax -imin <= 0 : imin=" << ivals[0] << "   imax=" << ivals[1] << std::endl;
        std::cout << " - pmin=" << pvals[0] << "   pmax=" << pvals[1] << std::endl;
        std::cout << " - step=" << step << std::endl;
        std::cout << " --------" << std::endl;
        std::cout << " - l=" << l << std::endl;
        std::cout << " - fc_l=" << fc_l << std::endl;
        //std::cout << " - H_lm=" << H_lm << std::endl;
        std::cout << " - gamma_l=" << gamma_l << std::endl;
        std::cout << " - f_s=" << f_s << std::endl;
        //std::cout << " - eta=" << eta << std::endl;
        //std::cout << " - a3=" << a3 << std::endl;
        //std::cout << " - asym=" << asym << std::endl;
        std::cout << " --------" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    //std::cout << "xmin=" << x.head(1) << "(microHz)" << std::endl;
    //std::cout << "xmax=" << x.tail(1) << "(microHz)" << std::endl;
    
    //std::cout << "step=" << step << std::endl;
    //std::cout << "pmin=" << pmin << "(microHz)  pmax=" << pmax << "(microHz)" << std::endl;
    //std::cout << "imin=" << imin << "  imax=" << imax << std::endl;
    //std::cin.ignore();
    
    return ivals;
}
