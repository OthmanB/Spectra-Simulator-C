/**
 * @file noise_models.cpp
 * @brief Functions related to noise models.
 *
 * This file contains the Functions related to noise models.
 *
 * @date 24 Feb 2016
 * @author obenomar
 */
#include <math.h>
#include <Eigen/Dense>
#include "noise_models.h"
#include <iostream>
#include <iomanip>

using Eigen::VectorXd;
using Eigen::MatrixXd;


VectorXd harvey_like(const VectorXd& noise_params, VectorXd& x, VectorXd& y, const int Nharvey){
	/* This function calculate a sum of harvey like profile + a white noise and adds 
	   these profiles to a initial input vector y. The function assumes that the 
	   inputs are in the following order: [H0, tc0, p0, ..., Hn, tcn, pn, N0]
	*/

	const long Nx=x.size();
	int cpt=0;
	VectorXd ones(Nx), white_noise(Nx), tmp(Nx); //, y_out(Nx);
	
	white_noise.setConstant(noise_params.tail(1)(0));
	//y_out.setZero();

	for(long i=0; i<Nharvey;i++){
		if(noise_params(cpt+1) != 0){
		  		tmp=((1e-3)*noise_params(cpt+1)*x).array().pow(noise_params(cpt+2)); 
				tmp=noise_params(cpt)*(tmp + ones.setConstant(1)).cwiseInverse();
				y= y + tmp;
		} 	
		cpt=cpt+3;
	}
	y=y + white_noise;

return y;
}


VectorXd harvey_like(const MatrixXd& noise_params, VectorXd x, VectorXd& y){
	/* This function calculate a sum of harvey like profile + a white noise and adds 
	   these profiles to a initial input vector y. The function assumes that the 
	   inputs are in the following order: [H0, tc0, p0, ..., Hn, tcn, pn, N0]
	*/
	const long Nx=x.size();
	int cpt=0;
	VectorXd ones(Nx), white_noise(Nx), tmp(Nx); //, y_out(Nx);
	
	const int Nrows=noise_params.rows();
	const int Nharvey=Nrows-1; // The last line is assumed to always contain the white noise.. Then Nharvey all_lines-1
	white_noise.setConstant(noise_params(Nrows-1, 0)); // check!

	for(long i=0; i<Nharvey;i++){
		if(noise_params(i,0) > 0 && noise_params(i,1) >0 && noise_params(i,2) >0){ // If all values of the line are valid...
		  		tmp=((1e-3)*noise_params(i,1)*x).array().pow(noise_params(i,2)); 
				tmp=noise_params(i,0)*(tmp + ones.setConstant(1)).cwiseInverse();
				y= y + tmp;
		} 	
		cpt=cpt+3;
	}
	y=y + white_noise;
return y;
}

VectorXd harvey_like(const MatrixXd& noise_params, const VectorXd& x){
	/* ALTERNATE FORM OF ABOVE
	   This function calculate a sum of harvey like profile + a white noise and adds 
	   these profiles to a initial input vector y. The function assumes that the 
	   inputs are in the following order: [H0, tc0, p0, ..., Hn, tcn, pn, N0]
	*/

	const long Nx=x.size();
	double H, tau, p;
	VectorXd s_noise(Nx), spec_noise(Nx), ones(Nx);

	spec_noise.setZero();	
	ones.setOnes();
	for(int i=0; i<noise_params.rows(); i++){
		H=noise_params(i,0);
		tau=noise_params(i,1);
		p=noise_params(i,2);
		
		s_noise.setZero();
		if(H>0 && (tau<0 || p<0) ){
			if(tau==-2 && p==-2){
				s_noise.setConstant(H);
				//std::cout << "Condition H ok but others at -2" << std::endl;
				spec_noise=spec_noise + s_noise;
			} else {
				std::cout << "Problem with the definition of the White noise" << std::endl;
				std::cout << "   H is defined positive but tau and/or p are not set to -2" << std::endl;
				std::cout << "   For a white noise, please set tau and p together to -2" << std::endl;
				std::cout << "The program will stop now" << std::endl;
				exit(EXIT_FAILURE);
			}
		} else {
			if(H<0 && (tau>0 || p>0)){
				std::cout << "Problem with the definition of the Height parameter for the noise" << std::endl;
				std::cout << " if tau or p positive, H cannot be negative" << std::endl;
				std::cout << "The program will stop now" << std::endl;
				exit(EXIT_FAILURE);
			}
			if( ((tau < 0) && (p > 0)) || ((tau > 0) && (p < 0)) ){
				std::cout << "Problem with the definition of the Harvey-like profile" << std::endl;
				std::cout << "   H is defined but tau or p are set to -1" << std::endl;
				std::cout << "   For a Harvey-like profile, please set tau and p together" << std::endl;
				std::cout << "The program will stop now" << std::endl;
				exit(EXIT_FAILURE);
			}
			if( (tau >0 && p >0) ){ // If all the parameters of the Harvey noise are defined, then...
				s_noise=((1e-3)*tau*x).array().pow(p); 
				//std::cout << "x[0]=" << x[0] << "   before inverse s_noise[0]=" << s_noise[0] << std::endl;
				s_noise=H*(s_noise + ones).cwiseInverse();
				spec_noise=spec_noise + s_noise;
				//std::cout << "main" << std::endl;
				//std::cout << "x[0]=" << x[0] << "   s_noise[0]=" << s_noise[0] << "   spec_noise[0]=" << spec_noise[0] << std::endl;
			}			
		}		
	}
return spec_noise;
}

VectorXd harvey1985(const VectorXd& noise_params, const VectorXd& x, const VectorXd& y, const int Nharvey){
	/* This function calculate a sum of harvey profile + a white noise and adds 
	 * these profiles to a initial input vector y. The Harvey profile differ from the
	 * Harvey-like by the fact that H0_1985= H0_like * tc is correlated to the timescale tc. 
	 * There is also a 2pi factor in the denominator. The function assumes that the 
	 * inputs are in the following order: [H0, tc0, p0, ..., Hn, tcn, pn, N0]
	*/

	const long double pi = M_PI;
	const long Nx=x.size();
	int cpt=0;
	VectorXd ones(Nx), white_noise(Nx), tmp(Nx), y_out(Nx);
	
	white_noise.setConstant(noise_params.tail(1)(0));
 	y_out=y;

	for(long i=0; i<Nharvey;i++){
		if(noise_params(cpt+1) != 0){
			tmp=((1e-3)*2*pi*noise_params(cpt+1)*x).array().pow(noise_params(cpt+2)); // Denominator
			tmp=noise_params(cpt)*noise_params(cpt+1)*(tmp + ones.setConstant(1)).cwiseInverse(); // Numerator/Denominator
			y_out= y_out + tmp; // Numerator/Denominator + white noise
		}
		cpt=cpt+3;
	}
	y_out=y_out + white_noise;

return y_out;
}




VectorXd harvey_1985(const MatrixXd& noise_params, const VectorXd& x){
	/* ALTERNATE FORM OF ABOVE
	 * This function calculate a sum of harvey profile + a white noise. 
	 * The Harvey profile differ from the
	 * Harvey-like by the fact that H0_1985= H0_like * tc is correlated to the timescale tc. 
	 * There is also a 2pi factor in the denominator. The function assumes that the 
	 * inputs are in the following order: [H0, tc0, p0, ..., Hn, tcn, pn, N0]
	*/
	const long double pi = M_PI;

	const long Nx=x.size();
	double H, tau, p;
	VectorXd s_noise(Nx), spec_noise(Nx), ones(Nx);

	spec_noise.setZero();	
	ones.setOnes();
	for(int i=0; i<noise_params.rows(); i++){
		H=noise_params(i,0);
		tau=noise_params(i,1);
		p=noise_params(i,2);
		
		s_noise.setZero();
		if(H>0 && (tau<0 || p<0) ){
			if(tau==-2 && p==-2){
				s_noise.setConstant(H);
				spec_noise=spec_noise + s_noise;
			} else {
				std::cout << "Problem with the definition of the White noise" << std::endl;
				std::cout << "   H is defined positive but tau and/or p are not set to -2" << std::endl;
				std::cout << "   For a white noise, please set tau and p together to -2" << std::endl;
				std::cout << "The program will stop now" << std::endl;
				exit(EXIT_FAILURE);
			}
		} else {
			if(H<=0 && (tau>0 || p>0)){
				std::cout << "Problem with the definition of the Height parameter for the noise" << std::endl;
				std::cout << " if tau or p positive, H cannot be negative" << std::endl;
				std::cout << "The program will stop now" << std::endl;
				exit(EXIT_FAILURE);
			}
			if( ((tau < 0) && (p > 0)) || ((tau > 0) && (p < 0)) ){
				std::cout << "Problem with the definition of the Harvey-like profile" << std::endl;
				std::cout << "   H is defined but tau or p are set to -1" << std::endl;
				std::cout << "   For a Harvey-like profile, please set tau and p together" << std::endl;
				std::cout << "The program will stop now" << std::endl;
				exit(EXIT_FAILURE);
			}
			if( (tau >0 && p >0) ){ // If all the parameters of the Harvey noise are defined, then...
				s_noise=((1e-3)*2*pi*tau*x).array().pow(p); // Denominator
				s_noise=H*tau* (s_noise + ones).cwiseInverse();
				spec_noise=spec_noise + s_noise;
			}			
		}		
	}
return spec_noise;
}




// 6 Dec 2023
double get_ksinorm(const double b, const double c, const Eigen::VectorXd& x) {
    double integral = 0.0;
    double h = x(1) - x(0); // assuming x is equally spaced
    for (int i = 0; i < x.size(); i++) {
        double term = 1.0 / (1.0 + std::pow(x(i) / b, c));
        
        if (i == 0 || i == x.size() - 1) {
            integral += 0.5 * term;
        } else {
            integral += term;
        }
    }  
    integral *= h;
    double ksi = b/integral;
    return ksi;
}

// 6 Dec 2023
VectorXd eta_squared_Kallinger2014(const VectorXd& x){
	const double x_nyquist=x.maxCoeff();
	VectorXd eta(x.size());
	eta=sin(0.5*M_PI*x.array()/x_nyquist)/(0.5*M_PI*x.array()/x_nyquist);
	if (x[0] == 0){ // This to avoid the Division by 0
		eta[0]=1;
	}
	return eta.array().square();
}

// 6 Dec 2023
VectorXd Kallinger2014_to_harveylike(const double numax, const double mu_numax, const VectorXd& noise_params, const VectorXd& x){
/*
	This function converts the parameters as they are defined into the TAMCMC and into the Spectrum simulator and 
	for the Kallinger+2014 noise implementation, into Harvey-like parameters, compatible with my functions 
	performing Harvey-like fits. 
	Be carefull with a1 and a2: These must here be provided independently to each others while in Kallinger+2014
	is has a1=a2. However, MCMC fit show discrepancies between these two, so we prefer here to separate them.
	It is up to the user to provide a highler level function that will generate a=k.numax^s and then (a1, a2) from it.
	
	Using notations from Table 2 of Kallinger+2014 (https://arxiv.org/pdf/1408.0817.pdf)
	Not that here we assume the instrumental noise to be Pinstrument(nu) = 0
	The noise_params must have parameters in this order:
	    - Granulation noise: 4 parameters to create a0 and b0
		- Noise a1, a2 : ka, sa, t
		- Noise b1: k1, s1, ( and c1, the slope of the SuperLorentzian)
		- Noise b2: k2, s2, ( and c2, the slope of the SuperLorentzian)
	Such that at the end we have: [ka,sa,t,k1,s1,c1, k2,s2,c2, N0]
*/
	VectorXd harvey_params(10);
	// Compute the Leakage effect as a sinc function (Eq 1 of Kallinger+2014)
	//const VectorXd eta_squared=eta_squared_Kallinger2014(x);
	// Compute b1, b2 and a
	const double a0=std::abs(noise_params[0]*std::pow(std::abs(numax),noise_params[1])); // Very Low frequencies
	const double b0=std::abs(noise_params[2]*std::pow(std::abs(numax + mu_numax),noise_params[3]));
	const double c0=std::abs(noise_params[4]);
	const double a1=noise_params[5];
	const double a2=noise_params[6];
	const double b1=std::abs(noise_params[7]*std::pow(std::abs(numax + mu_numax),noise_params[8])); // Intermediate
	const double b2=std::abs(noise_params[10]*std::pow(std::abs(numax + mu_numax),noise_params[11])); // Below numax
	const double c1=std::abs(noise_params[9]);
	const double c2=std::abs(noise_params[12]);
	const double N0=std::abs(noise_params[13]);
	// Compute the normalisation constants ksi1 and ksi2
	const double ksi0=get_ksinorm(b0, c0, x);
	const double ksi1=get_ksinorm(b1, c1, x);
	const double ksi2=get_ksinorm(b2, c2, x);

	// WARNING: eta_squared is not applied here. It is computed as a filter in artificial_spectrum.cpp
	harvey_params[0]=ksi0 * std::pow(a0,2)/b0; // a0^2/b0 :[ppm]^2 / [microHz]  ==> No problem of norm
	harvey_params[1]=1e3/b0; // b0 is in microHz in Kallinger+2014... my Harvey-like function needs mHz. Hence 1e3
	harvey_params[2]=c0;
	harvey_params[3]=ksi1 * std::pow(a1,2)/b1; 
	harvey_params[4]=1e3/b1; // b1 is in microHz in Kallinger+2014... my Harvey-like function needs  Giga-sec = . Hence 1e3
	harvey_params[5]=c1;
	harvey_params[6]=ksi2 * std::pow(a2,2)/b2;
	harvey_params[7]=1e3/b2; // b2 is in microHz in Kallinger+2014... my Harvey-like function needs mHz. Hence 1e3
	harvey_params[8]=c2;
	harvey_params[9]=N0;
	return harvey_params;
}


VectorXd Kallinger2014(const double numax, const double mu_numax, const VectorXd& noise_params,const VectorXd& x, const VectorXd& y){
/*
	Original definition used in TAMCMC for fitting
	Using notations from Table 2 of Kallinger+2014 (https://arxiv.org/pdf/1408.0817.pdf)
	Not that here we assume the instrumental noise to be Pinstrument(nu) = 0
	The noise_params must have parameters in this order:
	    - Granulation noise: 5 parameters to create a0 and b0
		- Noise a1, a2 : ka, sa, t
		- Noise b1: k1, s1, ( and c1, the slope of the SuperLorentzian)
		- Noise b2: k2, s2, ( and c2, the slope of the SuperLorentzian)
	Such that at the end we have: [ka,sa,t,k1,s1,c1, k2,s2,c2, N0]
*/
	const long Nx=x.size();
	VectorXd ones(Nx), white_noise(Nx), tmp0(Nx), tmp1(Nx), tmp2(Nx), y0(Nx), Power(Nx);
	ones.setOnes();
	// Compute the Leakage effect as a sinc function (Eq 1 of Kallinger+2014)
	const VectorXd eta_squared=eta_squared_Kallinger2014(x);
	// Compute b1, b2 and a
	const double a0=std::abs(noise_params[0]*std::pow(std::abs(numax),noise_params[1])); // Very Low frequencies
	const double b0=std::abs(noise_params[2]*std::pow(std::abs(numax + mu_numax),noise_params[3]));
	const double c0=std::abs(noise_params[4]);
	const double a1=noise_params[5];
	const double a2=noise_params[6];
	const double b1=std::abs(noise_params[7]*std::pow(std::abs(numax + mu_numax),noise_params[8])); // Intermediate
	const double b2=std::abs(noise_params[10]*std::pow(std::abs(numax + mu_numax),noise_params[11])); // Below numax
	const double c1=std::abs(noise_params[9]);
	const double c2=std::abs(noise_params[12]);
	const double N0=std::abs(noise_params[13]);
	// Compute the normalisation constants ksi1 and ksi2
	const double ksi0=get_ksinorm(b0, c0, x);
	const double ksi1=get_ksinorm(b1, c1, x);
	const double ksi2=get_ksinorm(b2, c2, x);
	y0=y;
	// White noise first, added to the input y-vector
	white_noise.setConstant(N0);
 	Power=y + white_noise;
	// Granulation SuperLorentzian
	tmp0=(x/b0).array().pow(c0); // Denominator
	tmp1=(eta_squared * ksi0 * std::pow(a0,2)/b0).cwiseProduct((tmp0 + ones).cwiseInverse()); // Numerator/Denominator
	Power=Power + tmp1;
	// First SuperLorentzian
	tmp0=(x/b1).array().pow(c1); // Denominator
	tmp1=(eta_squared * ksi1 * std::pow(a1,2)/b1).cwiseProduct((tmp0 + ones).cwiseInverse()); // Numerator/Denominator
	Power=Power + tmp1;
	// Second SuperLorentzian
	tmp0=(x/b2).array().pow(c2); // Denominator
	tmp2= (eta_squared *ksi2*std::pow(a2,2)/b2).cwiseProduct((tmp0 + ones).cwiseInverse()); // Numerator/Denominator
	Power=Power + tmp2;
	return Power;

}