/**
 * @file acoefs.h
 * @brief Functions that handle the acoefficients
 * @date 16 Nov 2021
 * @author Othman Benomar
 *
 * This file contains functions that handle the acoefficients. It can create nu(n,l,m) from acoeffs for l<=3, j=1,2,3,4,5,6. It can also decompose into aj coefficients. 
 */


#pragma once 
#include <math.h>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

/**
 * @brief Calculates the Hslm term using the Ritzoller1991 formula.
 *
 * @param s The s index of the coefficient.
 * @param l The degree l of the coefficient.
 * @param m The azimuthal order m of the coefficient.
 * @return The calculated Hslm term.
 */
long double Hslm_Ritzoller1991(const int s,const int l, const int m);


/**
 * @brief Calculates the Pslm polynomials.
 *
 * These coefficients take the Ritzoller1991 coefficients and normalize them by Pslm(l)=l as specified in Schou, JCD, Thompson 1994.
 *
 * @param s The s index of the coefficient.
 * @param l The degree l of the coefficient.
 * @param m The azimuthal order m of the coefficient.
 * @return The calculated Polynomial Pslm.
 */
long double Pslm(const int s,const int l,const int m);

/**
 * @brief Symmetric Splitting Tnlm.
 *
 * This function calculates the Tnlm splitting using the given nu_nlm frequencies and the specified l value.
 *
 * @param nu_nlm The nu_nlm coefficients.
 * @param l The degree l of the coefficient.
 * @return The calculated Tnlm coefficients.
 */
VectorXd Tnlm(VectorXd& nu_nlm, const int l);

/**
 * @brief Anti-Symmetric Splitting Snlm.
 *
 * This function calculates the Snlm splitting using the given nu_nlm frequencies and the specified l value.
 *
 * @param nu_nlm The nu_nlm coefficients.
 * @param l The degree of the coefficient.
 * @return The calculated Snlm coefficients.
 */
VectorXd Snlm(VectorXd& nu_nlm, const int l);

/** 
 * @brief Compute frequencies nu_nlm from a series of a-coefficients and provided the central frequency without splitting nunl0
 * 
 * This function computes frequencies nu_nlm from a series of a-coefficients and the central frequency nunl0. It uses the Legendre polynomials Pslm to calculate the nu_nlm.
 * 
 * @param nunl0 The central frequency without splitting.
 * @param l The degree l of the coefficient.
 * @param a1 The a1 coefficient.
 * @param a2 The a2 coefficient.
 * @param a3 The a3 coefficient.
 * @param a4 The a4 coefficient.
 * @param a5 The a5 coefficient.
 * @param a6 The a6 coefficient.
 * @return The calculated nu_nlm frequencies.
 */
VectorXd nunlm_from_acoefs(const long double nunl0, const int l, 
	const long double a1, const long double a2, const long double a3, const long double a4, const long double a5, const long double a6);

/**
 * @brief Evaluate the analytical a-coefficients for a given mode (n,l)
 * 
 * This function takes the splitted frequencies of a given mode (n,l) and determines the analytical a-coefficients from a1 to a6 for l<=3. It uses the Snlm and Tnlm functions to calculate the Snlm and Tnlm coefficients.
 * 
 * @param l The degree l of the mode.
 * @param nu_nls The splitted frequencies of the mode.
 * @return The calculated a-coefficients.
 */
VectorXd eval_acoefs(const int l, VectorXd& nu_nls);
