/**
 * @file noise_models.cpp
 * @brief Functions for generating noise models.
 *
 * This file contains functions for generating noise models. It includes functions for calculating the sum of a Harvey-like profile and white noise, as well as the sum of a Harvey profile and white noise.
 * @date 24 Feb 2016
 * @author obenomar
 */
#include <math.h>
#include <Eigen/Dense>

using Eigen::VectorXd;

/**
 * @brief Calculates the sum of a Harvey-like profile and white noise.
 *
 * This function calculates the sum of a Harvey-like profile and white noise, and adds these profiles to an initial input vector y. The function assumes that the inputs are in the following order: [H0, tc0, p0, ..., Hn, tcn, pn, N0].
 *
 * @param noise_params The parameters for the noise model.
 * @param x The x-coordinates of the data points.
 * @param y The y-coordinates of the data points.
 * @param Nharvey The number of Harvey profiles to include in the sum.
 * @return VectorXd The resulting vector after adding the Harvey-like profile and white noise.
 */
VectorXd harvey_like(const VectorXd& noise_params, const VectorXd& x, const VectorXd& y, const int Nharvey){
	const long Nx=x.size();
	int cpt=0;
	VectorXd ones(Nx), white_noise(Nx), tmp(Nx), y_out(Nx); //, y_out(Nx);
	
	white_noise.setConstant(noise_params.tail(1)(0));
	y_out=y;

	for(long i=0; i<Nharvey;i++){
		if(noise_params(cpt+1) != 0){
		  		tmp=((1e-3)*noise_params(cpt+1)*x).array().pow(noise_params(cpt+2)); 
				tmp=noise_params(cpt)*(tmp + ones.setConstant(1)).cwiseInverse();
				y_out= y_out + tmp;
		} 	
		cpt=cpt+3;
	}
	y_out=y_out + white_noise;

return y_out;
}


/**
 * @brief Calculates the sum of a Harvey profile and white noise.
 *
 * This function calculates the sum of a Harvey profile and white noise, and adds these profiles to an initial input vector y. The Harvey profile differs from the Harvey-like profile by the fact that H0_1985 = H0_like * tc is correlated to the timescale tc. There is also a 2pi factor in the denominator. The function assumes that the inputs are in the following order: [H0, tc0, p0, ..., Hn, tcn, pn, N0].
 *
 * @param noise_params The parameters for the noise model.
 * @param x The x-coordinates of the data points.
 * @param y The y-coordinates of the data points.
 * @param Nharvey The number of Harvey profiles to include in the sum.
 * @return VectorXd The resulting vector after adding the Harvey profile and white noise.
 */
VectorXd harvey1985(const VectorXd& noise_params, const VectorXd& x, const VectorXd& y, const int Nharvey){
	const long double pi = 3.141592653589793238L;
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

