/*
 * noise_models.cpp
 *
 *  Created on: 24 Feb 2016
 *      Author: obenomar
 */
#include <math.h>
#include <Eigen/Dense>
//#include <string>
#include "noise_models.h"
#include <iostream>
#include <iomanip>

using Eigen::VectorXd;
using Eigen::MatrixXd;


//VectorXd harvey_like(const VectorXd noise_params, VectorXd x, VectorXd y, const int Nharvey);
//VectorXd harvey1985(const VectorXd noise_params, VectorXd x, VectorXd y, const int Nharvey);

VectorXd harvey_like(const VectorXd noise_params, VectorXd x, VectorXd y, const int Nharvey){
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


VectorXd harvey_like(const MatrixXd noise_params, VectorXd x, VectorXd y){
	/* This function calculate a sum of harvey like profile + a white noise and adds 
	   these profiles to a initial input vector y. The function assumes that the 
	   inputs are in the following order: [H0, tc0, p0, ..., Hn, tcn, pn, N0]
	*/

//	std::cout << "noise_models::harvey_like UNDER CONSTRUCTION... NEED CHECKS" << std::endl;
//	exit(EXIT_SUCCESS);

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

VectorXd harvey_like(const MatrixXd noise_params, const VectorXd x){
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
		
		/*
		std::cout << "--------" << std::endl;
		std::cout << "H=" << H << std::endl;
		std::cout << "tau=" << tau << std::endl;
		std::cout << "p=" << p << std::endl;
		std::cout << "--------" << std::endl;
		*/

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
				s_noise=((1e-3)*tau*x).array().pow(p); 
				//std::cout << "x[0]=" << x[0] << "   before inverse s_noise[0]=" << s_noise[0] << std::endl;
				s_noise=H*(s_noise + ones).cwiseInverse();
				spec_noise=spec_noise + s_noise;
				//std::cout << "main" << std::endl;
				//std::cout << "x[0]=" << x[0] << "   s_noise[0]=" << s_noise[0] << "   spec_noise[0]=" << spec_noise[0] << std::endl;
			}			
		}
		//std::cout << "s_noise=" << s_noise << std::endl;
		//std::cout << "s_noise.size()=" << s_noise.size() << std::endl;		
	}

return spec_noise;

}

VectorXd harvey1985(const VectorXd noise_params, VectorXd x, VectorXd y, const int Nharvey){
	/* This function calculate a sum of harvey profile + a white noise and adds 
	 * these profiles to a initial input vector y. The Harvey profile differ from the
	 * Harvey-like by the fact that H0_1985= H0_like * tc is correlated to the timescale tc. 
	 * There is also a 2pi factor in the denominator. The function assumes that the 
	 * inputs are in the following order: [H0, tc0, p0, ..., Hn, tcn, pn, N0]
	*/

	const long double pi = 3.141592653589793238L;
	const long Nx=x.size();
	int cpt=0;
	VectorXd ones(Nx), white_noise(Nx), tmp(Nx); //, y_out(Nx);
	
	white_noise.setConstant(noise_params.tail(1)(0));
 
	for(long i=0; i<Nharvey;i++){
		if(noise_params(cpt+1) != 0){
			tmp=((1e-3)*2*pi*noise_params(cpt+1)*x).array().pow(noise_params(cpt+2)); // Denominator
			tmp=noise_params(cpt)*noise_params(cpt+1)*(tmp + ones.setConstant(1)).cwiseInverse(); // Numerator/Denominator
			y= y + tmp; // Numerator/Denominator + white noise
		}
		cpt=cpt+3;
	}
	y=y + white_noise;

return y;
}


VectorXd harvey_1985(const MatrixXd noise_params, const VectorXd x){
	/* ALTERNATE FORM OF ABOVE
	 * This function calculate a sum of harvey profile + a white noise. 
	 * The Harvey profile differ from the
	 * Harvey-like by the fact that H0_1985= H0_like * tc is correlated to the timescale tc. 
	 * There is also a 2pi factor in the denominator. The function assumes that the 
	 * inputs are in the following order: [H0, tc0, p0, ..., Hn, tcn, pn, N0]
	*/
	const long double pi = 3.141592653589793238L;

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
				//std::cout << "x[0]=" << x[0] << "   before inverse s_noise[0]=" << s_noise[0] << std::endl;
				//s_noise=H*2*pi*tau*(s_noise + ones).cwiseInverse();
				s_noise=H*tau* (s_noise + ones).cwiseInverse();
				spec_noise=spec_noise + s_noise;
				//std::cout << "main" << std::endl;
				//std::cout << "x[0]=" << x[0] << "   s_noise[0]=" << s_noise[0] << "   spec_noise[0]=" << spec_noise[0] << std::endl;
			}			
		}
		//std::cout << "s_noise=" << s_noise << std::endl;
		//std::cout << "s_noise.size()=" << s_noise.size() << std::endl;		
	}

return spec_noise;

}


