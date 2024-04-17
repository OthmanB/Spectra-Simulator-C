/**
 * @file decompose_nu.cpp
 * @brief Functions of decompose_nu_nl.
 *
 * This file contains the function to decompose nu(n,l) into the asymptotic terms. 
 *
 * @date 24 Feb 2023
 * @author obenomar
 */

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include "../../linspace.h"
#include "../../linfit.h"
#include "decompose_nu.h"
#include "data.h"
//#include "string_handler.h" // Original from the standalone decomposition program
#include "../../ioproc.h" // Replacement of string_handler.h

using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;

#include <Eigen/Core>

Eigen::VectorXd gradient(const Eigen::VectorXd& x, double (*f)(const Eigen::VectorXd&)) {
    const double h = 1e-4;
    const int n = x.size();
    Eigen::VectorXd grad(n);
    Eigen::VectorXd x_center = x;
    //#pragma omp parallel for
    for (int i = 0; i < n; i++) {
        Eigen::VectorXd x_plus = x;
        Eigen::VectorXd x_minus = x;
        if (x(i) + h < 1.0) {
            x_plus(i) += h;
        }
        if (x(i) - h > 0.0) {
            x_minus(i) -= h;
        }
        double f_plus = f(x_plus);
        double f_minus = f(x_minus);
        double f_center = f(x);
        if (x(i) + h < 1.0 && x(i) - h > 0.0) {
            grad(i) = (f_plus - f_minus) / (2.0 * h) - (1.0 / 6.0) * h * h * (f(x_plus + x_minus - 2.0 * x_center) - 2.0 * f_center + f(x_plus + x_minus - 2.0 * x_center)) / (h * h);
        } else if (x(i) + h < 1.0) {
            grad(i) = (f_plus - f_center) / h - (1.0 / 2.0) * h * f(x_plus + x - 2.0 * x_center) / (h * h);
        } else if (x(i) - h > 0.0) {
            grad(i) = (f_center - f_minus) / h - (1.0 / 2.0) * h * f(x - 2.0 * x_center + x_minus) / (h * h);
        } else {
            grad(i) = 0.0;
        }
    }

    return grad;
}

#include <Eigen/Core>

Eigen::VectorXd gradient_discrete(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
    const int n = x.size();
    Eigen::VectorXd grad(n);

    #pragma omp parallel for
    
    for (int i = 0; i < n; i++) {
        if (i == 0) {
            grad(i) = (y(i+1) - y(i)) / (x(i+1) - x(i));
        } else if (i == n-1) {
            grad(i) = (y(i) - y(i-1)) / (x(i) - x(i-1));
        } else {
            grad(i) = (y(i+1) - y(i-1)) / (x(i+1) - x(i-1));
        }
    }
    
    return grad;
}





//VectorXd where_in_range(const VectorXd& vec, const double value_min, const double value_max, const bool strict){
/*
 * Gives the indexes of values of an array within a specified range 
 *
*/
/*
   int cpt;
   VectorXd index(vec.size());
	
	cpt=0;
	for(int i=0; i<vec.size(); i++){
		if (strict == 1){
			if(vec[i] > value_min && vec[i] < value_max){
				index[cpt]=i;
				cpt=cpt+1;
			}
		} else{
			if(vec[i] >= value_min && vec[i] <= value_max){
				index[cpt]=i;
				cpt=cpt+1;
			}		
		}		
	}
	if(cpt >=1){
		index.conservativeResize(cpt);
	} else{
		index.resize(1);
		index[0]=-1;
	}
	return index;
}
*/

Data_asympt_p decompose_nu_nl(const int l, const VectorXd fl0, const VectorXd fl, const double Cfactor, const bool verbose){
    /*
        Function that takes the p modes frequencies and decompose them
        ie, identifies, Dnu, epsilon, n, d0l and O(2) (higher order) terms
        This allows to extract the different contribution to the modes 
        and use the following properties:
            - Dnu and epsilon are determined by the l=0 frequencies
            - d0l + O(2) and the radial order n are determined by solving analytically the linear fit
                to the frequencies at the fixed pre-determined value of the slope = Dnu
            - O(2) is determined by assuming it to be a 0-mean perturbation around d0l
    */
    // Declaration
	Data_asympt_p output;
	MatrixXd n_all(4,fl.size());
	VectorXd O2_term, fitcoef, func, sol(4);
    VectorXi n_s_optimum;
    VectorXd tmp, n_ref, n_l, n_best, grad;
	double Dnu, epsilon, n0_0, d0l, e0, n0_l;
    output.error_status=false; // No error by default

	// Compute Dnu and epsilon
	n_ref=linspace(0, fl0.size()-1, fl0.size());
	fitcoef=linfit(n_ref, fl0);
    Dnu=fitcoef[0];
    epsilon=modf(fitcoef[1]/fitcoef[0], &n0_0);
    //grad=gradient_discrete(n_ref, fl0); // Alternative way using the gradient in order to get Dnu. But then you don't have epsilon
    //std::cout << "grad (O1)= " << grad << std::endl;
    //std::cout << "                                           " << grad.mean() << std::endl;
    //std::cout << "Dnu = " << Dnu << std::endl;
    // Least square on l=1 frequencies with a fix slope
    //   1. Isolate d0l + O(2)
    tmp.resize(fl.size());   tmp.setConstant((epsilon + 1.*l/2)*Dnu);
    func=fl - tmp;// This is Dnu_ref.n + d0l + O(2)
    //    2. Determine a first guess for n
    e0=modf(func[0]/Dnu, &n0_l);
    n_l=linspace(n0_l, n0_l + func.size()-1, func.size());
    //    3. Make a list of adjacent potential n and find the one that ensure |Dnu_ref| >> |d0l| + O(2) and d0l < 0
    tmp.resize(n_l.size());
    tmp.setConstant(1);
    n_all.row(0)=n_l-tmp;
    n_all.row(1)=n_l;
    n_all.row(2)=n_l+tmp;
    tmp.setConstant(2);
    n_all.row(3)=n_l+tmp;
    for(int n_s=0;n_s<4;n_s++){
        sol[n_s]=func.mean() - Dnu*n_all.row(n_s).mean();
    }
    n_s_optimum=where_in_range(sol, -Dnu*Cfactor, Dnu*Cfactor, false); // Identify solutions close to 0 as |d0l|<<|Dnu_ref|
    if (n_s_optimum.size() > 1){
        std::cout << "Error: multiple solutions of d0l found" << std::endl;
        std::cout << "       Debug required" << std::endl;
        std::cout << "      d0l =" <<  d0l << std::endl;
        exit(EXIT_FAILURE);
	} else{
        if (n_s_optimum[0] == -1){
            std::cout << "Error: Optimum values not found within the table of solutions for n !" << std::endl;
            std::cout << "       Did you put the correct l associated to fl ?" << std::endl;
            std::cout << "       l       = " << l << std::endl;
            std::cout << "       fl      = " << fl.transpose() << std::endl;
            std::cout << "       sol     = " << sol.transpose() << std::endl;
            std::cout << "       Dnu     = " << Dnu << std::endl;
            std::cout << "       Cfactor = " << Cfactor << std::endl;
            output.error_status=true; // Return the structure but flag that an error happened
            return output;
            //exit(EXIT_FAILURE);
        } else{
    	    n_best=n_all.row(n_s_optimum[0]);  //This is the list of identified n following the conditions in the where
		    d0l=sol[n_s_optimum[0]];  // This is the best solution for d0l + O(2)
        }
	}
    //    4. Identify the residual variation by assuming that O(2) is a random term of 0-mean
    tmp.resize(n_best.size());
    tmp.setConstant(d0l);
    O2_term=func- n_best*Dnu - tmp;
    if (verbose == true){
        std::cout << "--"  << std::endl;
        std::cout << "Best list of n matching conditions |Dnu_ref| >> |d0l| + O(2) and d0l < 0:"   << std::endl;
        std::cout << n_best.transpose()   << std::endl;
        std::cout << "---"   << std::endl;
        std::cout << "Identified best solution for d0l:"   << std::endl;
        std::cout << d0l   << std::endl;
        std::cout << "---"  << std::endl;
        std::cout << "Residual of the constrained fit:"  << std::endl;
        std::cout << O2_term.transpose()  << std::endl;
	}
	output.n=n_best;
	output.Dnu=Dnu;
	output.epsilon=epsilon;
	output.d0l=d0l;
	output.O2_term=O2_term;
    return output;
}
