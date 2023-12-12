/**
 * @file noise_models.h
 * @brief Declarations of functions related to noise models.
 *
 * This file contains the declarations of functions related to noise models.
 *
 * @date 24 Feb 2016
 * @author obenomar
 */
#pragma once
#include <math.h>
#include <Eigen/Dense>

using Eigen::VectorXd;
using Eigen::MatrixXd;


/**
 * @brief Calculate a sum of Harvey-like profile and white noise and add them to an initial input vector.
 *
 * @param noise_params The vector of noise parameters. The parameters are assumed to be in the following order: [H0, tc0, p0, ..., Hn, tcn, pn, N0].
 * @param x The input vector x.
 * @param y The modified vector y.
 * @param Nharvey The number of Harvey profiles to consider.
 * @return The modified vector y.
 *
 * This function calculates a sum of a Harvey-like profile and white noise and adds them to an initial input vector. The function assumes that the inputs are in the following order: [H0, tc0, p0, ..., Hn, tcn, pn, N0].
 *
 * @note The function modifies the vector y and returns it.
 */
VectorXd harvey_like(const VectorXd& noise_params, VectorXd& x, VectorXd& y, const int Nharvey);



/**
 * @brief Calculate a sum of Harvey-like profile and white noise and add them to an initial input vector.
 *
 * @param noise_params The matrix of noise parameters.
 * @param x The input vector x.
 * @param y The initial input vector y.
 * @return The modified vector y.
 *
 * This function calculates a sum of a Harvey-like profile and white noise and adds them to an initial input vector y. The function assumes that the inputs are in the following order: [H0, tc0, p0, ..., Hn, tcn, pn, N0].
 *
 * @note The last row of noise_params is assumed to always contain the white noise.
 */
VectorXd harvey_like(const MatrixXd& noise_params, VectorXd x, VectorXd& y);

/**
 * @brief Calculate a sum of Harvey-like profile and white noise and add them to an initial input vector.
 *
 * @param noise_params The matrix of noise parameters.
 * @param x The input vector x.
 * @return The modified vector spec_noise.
 *
 * This function calculates a sum of a Harvey-like profile and white noise and adds them to an initial input vector. The function assumes that the inputs are in the following order: [H0, tc0, p0, ..., Hn, tcn, pn, N0].
 *
 * @note The function returns the modified vector spec_noise.
 */
VectorXd harvey_like(const MatrixXd& noise_params, const VectorXd& x);

/**
 * @brief Calculate a sum of Harvey profile and white noise and add them to an initial input vector.
 * 
 * @param noise_params The matrix of noise parameters.
 * @param x The input vector x.
 * @return The modified vector spec_noise.
 * 
 * This function calculates a sum of Harvey profile and white noise and adds them to an initial input vector. The function assumes that the inputs are in the following order: [H0, tc0, p0, ..., Hn, tcn, pn, N0].
 * 
 * @note The function returns the modified vector spec_noise.
 */
VectorXd harvey_1985(const MatrixXd& noise_params, const VectorXd& x);

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
VectorXd harvey1985(const VectorXd& noise_params, const VectorXd& x, const VectorXd& y, const int Nharvey);


/**
 * @brief Compute eta squared using definition from Kallinger+2014 (Eq. 1)
 * 
 * 
 * @param x 
 * @date 6 December 2023
 * @return VectorXd containing eta^2
 */
VectorXd eta_squared_Kallinger2014(const VectorXd& x);


/**
 * @brief Calculates the norm ksi associated to the noise profile as per defined by Kallinger+2014.
 *
 * @param b The b parameter of the noise model.
 * @param c The c parameter of the noise model.
 * @param x The x-coordinates of the data points.
 * @date 6 December 2023
 * @return double The norm ksi.
 */
double get_ksinorm(const double b, const double c, const Eigen::VectorXd& x);

/**
 * @brief Calculates the noise background as in Kallinger+2014
 * 
 * This function calculates the noise background as in Kallinger+2014 (https://arxiv.org/pdf/1408.0817.pdf)
 * It uses notations similar to their Table 2
 * Note that here we assume the instrumental noise to be Pinstrument(nu) = 0
 * Also note that one would need to add eta^2 to the Gaussian envelope if this noise is used for fitting
 * The noise_params must be have the parameters provided in this order:
 * 	 Granulation Amplitude a : ka, sa, t
 *   Characteristic frequency ~ 1/timescale b1: k1, s1, ( and c1, the slope of the SuperLorentzian)
 *   Characteristic frequency ~ 1/timescale b2: k2, s2, ( and c2, the slope of the SuperLorentzian)
 * Such that at the end we have: [ka,sa,t,k1,s1,c1, k2,s2,c2, N0]
 * @param numax The frequency at maximum power of the modes 
 * @param mu_numax A bias parameter to describe inacurracy of the relation Power=f(numax)
 * @param noise_params The noise parameters structured as explained in the function help
 * @param x A vector that contains the frequencies of the data points
 * @param y An initial vector of the same size as x and that contains the power. This vector could be initialised to 0 or be initialised using an instrumental noise
 * @return VectorXd of the same size as x and y that contains the power of the noise background added to the initial y vector.
 * @date 6 December 2023
 * @note The function returns the noise background.
 */
VectorXd Kallinger2014(const double numax, const double mu_numax, const VectorXd& noise_params,const VectorXd& x, const VectorXd& y);

/**
 * @brief Convert Kallinger+2014 noise parameters into parameters compatible with my definition of Harvey-like
 * 
 * 	This function converts the parameters as they are defined into the TAMCMC and into the Spectrum simulator and 
 *	for the Kallinger+2014 noise implementation, into Harvey-like parameters, compatible with my functions 
 *	performing Harvey-like fits. 
 *	Be carefull with a1 and a2: These must here be provided independently to each others while in Kallinger+2014
 *	is has a1=a2. However, MCMC fit show discrepancies between these two, so we prefer here to separate them.
 *	It is up to the user to provide a highler level function that will generate a=k.numax^s and then (a1, a2) from it.
 *
 * This function calculates the noise background as in Kallinger+2014 (https://arxiv.org/pdf/1408.0817.pdf)
 * It uses notations similar to their Table 2
 * Note that here we assume the instrumental noise to be Pinstrument(nu) = 0
 * The noise_params must be have the parameters provided in this order:
 *   Granulation Amplitude a : ka, sa, t
 *   Characteristic frequency ~ 1/timescale b1: k1, s1, ( and c1, the slope of the SuperLorentzian)
 *   Characteristic frequency ~ 1/timescale b2: k2, s2, ( and c2, the slope of the SuperLorentzian)
 * Such that at the end we have three Harvey-like compatible values: [a0, b0,p0,a1,b1,p1, a2,b2,p2, N0]
 * @param numax The frequency at maximum power of the modes 
 * @param mu_numax A bias parameter to describe inacurracy of the relation Power=f(numax)
 * @param noise_params The noise parameters structured as explained in the function help
 * @param x A vector that contains the frequencies of the data points
 * @return VectorXd of size 10 that contains the parameters of the Harvey-like parameters in the form [a0, b0,p0,a1,b1,p1, a2,b2,p2, N0]
 * @date 6 December 2023
 * @note The function returns the noise background.
 */
VectorXd Kallinger2014_to_harveylike(const double numax, const double mu_numax, const VectorXd& noise_params, const VectorXd& x);