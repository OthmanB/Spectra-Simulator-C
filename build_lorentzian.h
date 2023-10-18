/**
 * @file build_lorentzian.h
 * @brief Contains the implementation of the build_l_mode_a1l_etaa3 function.
 *
 * This file provides the implementation of the build_l_mode_a1l_etaa3 function, which builds a Lorentzian mode with asymmetry and splitting based on the given parameters.
 *
 * @date 22 Feb 2016
 * @author obenomar
 */
#pragma once
#include <math.h>
#include <Eigen/Dense>
#include "function_rot.h"
# include <iostream>
# include <iomanip>
#include <string>

using Eigen::VectorXd;
using Eigen::VectorXi;

/**
 * @brief Calculate the Qlm term used in the centrifugal distortion calculation.
 *
 * This function calculates the Qlm term used in the centrifugal distortion calculation based on the given parameters.
 * The final Qlm includes the dnl=2/3 term that is sometimes considered separately by some authors.
 * See Papini-Gizon 2020 or Gizon 2002 for the definition.
 *
 * @param l The mode degree.
 * @param m The mode order.
 * @return The calculated Qlm term.
 */
double Qlm(const int l, const int m);

/**
 * @brief Set the indices imin and imax based on the given parameters.
 *
 * This function sets the indices imin and imax based on the given parameters, such as the frequency range, mode degree, central frequency, mode width, splitting frequency, step size, and constant c.
 * The function calculates the proposed values pmin and pmax based on the parameters and handles the boundaries to ensure that the indices stay within the frequency range.
 * The function also checks if imax - imin <= 0 and prints a warning message if this condition is met.
 *
 * @param x The frequency range.
 * @param l The mode degree.
 * @param fc_l The central frequency of the mode.
 * @param gamma_l The mode width.
 * @param f_s The splitting frequency.
 * @param c The truncation constant. It limits the range of the computation. Used to (1) accelerate the fitting and (2) measure the asymetry.
 * @param step The step size.
 * @return The indices imin and imax.
 */
VectorXi set_imin_imax(const VectorXd& x, const int l, const double fc_l, const double gamma_l, const double f_s, const double c, const double step);

/**
 * @brief Build a Lorentzian mode with asymmetry and splitting by dropping the assumption S11=S22
 *
 * This function builds a Lorentzian mode with asymmetry and splitting based on the given parameters. Dropping the asumption for S11=S22 means a different a1 splitting for l=1 and l=2. But we keep ASSUMES a1(l=3) =  (a1(1) + a1(2))/2.
 *
 * @param x_l The frequency range for the mode.
 * @param H_l The height of the mode.
 * @param fc_l The central frequency of the mode.
 * @param f_s1 The splitting frequency for l=1.
 * @param f_s2 The splitting frequency for l=2.
 * @param eta0 The asphericity parameter.
 * @param a3 The latitudinal effect for l=2.
 * @param asym The asymmetry parameter.
 * @param gamma_l The mode width.
 * @param l The mode degree.
 * @param V The mode visibility.
 * @return The Lorentzian mode.
 */
VectorXd build_l_mode_a1l_etaa3(const VectorXd& x_l, const double H_l,  const double fc_l,  const double f_s1,  const double f_s2,  const double eta,  const double a3,  const double asym,  const double gamma_l, const int l,  const VectorXd& V);

/**
 * @brief Build a Lorentzian mode with asymmetry and splitting (model a1l_a2a3)
 *
 * This function builds a Lorentzian mode with asymmetry and splitting based on the given parameters.
 * This model includes:
 * - Asymmetry of Lorentzian asym
 * - Splitting a1(l=1) and a1(l=2). ASSUMES a1(l=3) = (a1(1) + a1(2))/2.
 * - An Asphericity parameter eta
 * - Latitudinal effect a3(l=2) only. We consider a3(l=3)=0
 *
 * @param x_l The frequency range for the mode.
 * @param H_l The height of the mode.
 * @param fc_l The central frequency of the mode.
 * @param f_s1 The splitting frequency for l=1.
 * @param f_s2 The splitting frequency for l=2.
 * @param a2 The a2 coefficient.
 * @param a3 The a3 coefficient.
 * @param asym The asymmetry parameter.
 * @param gamma_l The mode width.
 * @param l The mode degree.
 * @param V The mode visibility.
 * @return The Lorentzian mode.
 */
VectorXd build_l_mode_a1l_a2a3(const VectorXd& x_l, const double H_l,  const double fc_l,  const double f_s1,  const double f_s2,  const double a2,  const double a3,  const double asym,  const double gamma_l, const int l,  const VectorXd& V);

/**
 * @brief Build a Lorentzian mode with asymmetry and splitting (model a1etaa3)
 *
 * This function builds a Lorentzian mode with asymmetry and splitting based on the given parameters.
 * This model includes:
 * - Asymmetry of Lorentzian asym
 * - Splitting a1
 * - An Asphericity parameter eta
 * - Latitudinal effect a3
 *
 * @param x_l The frequency range for the mode.
 * @param H_l The height of the mode.
 * @param fc_l The central frequency of the mode.
 * @param f_s The splitting frequency.
 * @param eta0 The asphericity parameter.
 * @param a3 The a3 coefficient.
 * @param asym The asymmetry parameter.
 * @param gamma_l The mode width.
 * @param l The mode degree.
 * @param V The mode visibility.
 * @return The Lorentzian mode.
 */
VectorXd build_l_mode_a1etaa3(const VectorXd& x_l,  const double H_l,  const double fc_l,  const double f_s,  const double eta, const double a3,  const double asym,  const double gamma_l, const int l, const VectorXd& V);

// OBSELETE FUNCTION
VectorXd build_l_mode_a1acta3(const VectorXd& x_l,  const double H_l,  const double fc_l,  const double f_s,  const double eta, const double a3,  const double b,  const double alpha,  const double asym, double gamma_l, const int l,  const VectorXd& V);

/**
 * @brief Build a Lorentzian mode with asymmetry and splitting (model a1etaa3_v2)
 *
 * This function builds a Lorentzian mode with asymmetry and splitting based on the given parameters.
 * This model includes:
 * - Asymmetry of Lorentzian asym
 * - Splitting a1
 * - An Asphericity parameter eta
 * - Latitudinal effect a3
 *
 * @param x_l The frequency range for the mode.
 * @param H_lm The heights of the mode.
 * @param fc_l The central frequency of the mode.
 * @param f_s The splitting frequency.
 * @param eta0 The eta0 coefficient.
 * @param a3 The a3 coefficient.
 * @param asym The asymmetry parameter.
 * @param gamma_l The mode width.
 * @param l The mode degree.
 * @return The Lorentzian mode.
 */
VectorXd build_l_mode_a1etaa3_v2(const VectorXd& x_l, const VectorXd& H_lm,  const double fc_l,  const double f_s,  const double eta,  const double a3,  const double asym,  const double gamma_l, const int l);

/**
 * @brief Build a Lorentzian mode with asymmetry and splitting (model a1l_etaa3_v2)
 * 
 * This function builds a Lorentzian mode with asymmetry and splitting based on the given parameters.
 * This model includes:
 * - Asymmetry of Lorentzian asym
 * - Splitting a1(l=1) and a1(l=2). ASSUMES a1(l=3) = (a1(1) + a1(2))/2.
 * - An Asphericity parameter eta
 * - Latitudinal effect a3(l=2) only. We consider a3(l=3) = 0
 * - Inclination IS NOT IMPOSED contrary to build_l_mode_a1l_etaa3. INSTEAD H_lm should be of dimension l(l+1) and provide all the heights
 * 
 * @param x_l The frequency range for the mode.
 * @param H_lm The heights of the mode.
 * @param fc_l The central frequency of the mode.
 * @param f_s1 The splitting frequency for l=1.
 * @param f_s2 The splitting frequency for l=2.
 * @param eta0 The eta0 coefficient.
 * @param a3 The a3 coefficient.
 * @param asym The asymmetry parameter.
 * @param gamma_l The mode width.
 * @param l The mode degree.
 * @return The Lorentzian mode.
 */
VectorXd build_l_mode_a1l_etaa3_v2(const VectorXd& x_l, const VectorXd& H_lm,  const double fc_l,  const double f_s1,  const double f_s2,  const double eta,  const double a3,  const double asym, double gamma_l, const int l);

VectorXd build_l_mode_a1a2a3(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s, const double a2, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V);

/**
 * @brief Calculate the optimized Lorentzian model with asymmetry and splitting (model a1l_etaa3)
 *
 * This function calculates the Lorentzian model with asymmetry and splitting based on the given parameters.
 * This model includes:
 * - Asymmetry of Lorentzian asym
 * - Splitting a1(l=1) and a1(l=2). Assumes a1(l=3) = (a1(1) + a1(2))/2.
 * - An Asphericity parameter eta
 * - Latitudinal effect a3(l=2) only. We consider a3(l=3) = 0
 * - Inclination is not imposed, instead H_l should be of dimension l(l+1) and provide all the heights
 *
 * @param x The frequency range.
 * @param y The original vector.
 * @param H_l The heights of the mode.
 * @param fc_l The central frequency of the mode.
 * @param f_s1 The splitting frequency for l=1.
 * @param f_s2 The splitting frequency for l=2.
 * @param eta0 The eta0 coefficient.
 * @param a3 The a3 coefficient.
 * @param asym The asymmetry parameter.
 * @param gamma_l The mode width.
 * @param l The mode degree.
 * @param V The vector V.
 * @param step The step size.
 * @param c The truncation constant. It limits the range of the computation. Used to (1) accelerate the fitting and (2) measure the asymetry.
 * @return The Lorentzian model.
 */
VectorXd optimum_lorentzian_calc_a1l_etaa3(const VectorXd& x, const VectorXd& y, const  double H_l, const  double fc_l,  const double f_s1,  const double f_s2,  const double eta,  const double a3,  const double asym,  const double gamma_l, const int l,  const VectorXd& V,  const double step, const double c);

/**
 * @brief Calculate the optimized Lorentzian model with asymmetry and splitting (model a1etaa3)
 *
 * This function calculates the Lorentzian model with asymmetry and splitting based on the given parameters.
 * This model includes:
 * - Asphericity of Lorentzian a1etaa3
 * - Splitting frequency f_s
 * - Latitudinal effect a3(l=2) only. We consider a3(l=3) = 0
 * - Inclination is not imposed, instead H_l should be of dimension l(l+1) and provide all the heights
 *
 * @param x The frequency range.
 * @param y The original vector.
 * @param H_l The heights of the mode.
 * @param fc_l The central frequency of the mode.
 * @param f_s The splitting frequency.
 * @param eta0 The eta0 coefficient.
 * @param a3 The a3 coefficient.
 * @param asym The asymmetry parameter.
 * @param gamma_l The mode width.
 * @param l The mode degree.
 * @param V The vector V.
 * @param step The step size.
 * @param c The truncation constant. It limits the range of the computation. Used to (1) accelerate the fitting and (2) measure the asymetry.
 * @return The optimized Lorentzian model.
 */
VectorXd optimum_lorentzian_calc_a1etaa3(const VectorXd& x, const VectorXd& y,  const double H_l,  const double fc_l,  const double f_s,  const double eta,  const double a3,  const double asym,  const double gamma_l, const int l,  const VectorXd& V,  const double step, const double c);

// OBSELETE FUNCTION
VectorXd optimum_lorentzian_calc_a1acta3(const VectorXd& x, const VectorXd& y,  const double H_l,  const double fc_l,  const double f_s,  const double eta,  const double a3, 
		 const double b,  const double alpha,  const double asym,  const double gamma_l, const int l,  const VectorXd& V,  const double step, const double c);

/**
 * @brief Calculate the optimized Lorentzian model with asymmetry and splitting (model a1l_etaa3_v2)
 *
 * This function calculates the Lorentzian model with asymmetry and splitting based on the given parameters.
 * This model includes:
 * - Asphericity of Lorentzian a1l_etaa3_v2
 * - Splitting frequency f_s
 * - Inclination is not imposed
 *
 * @param x The frequency range.
 * @param y The original vector.
 * @param H_lm The heights of the mode for each (l,m) component.
 * @param fc_l The central frequency of the mode.
 * @param f_s1 The splitting frequency for l=1.
 * @param f_s2 The splitting frequency for l=2.
 * @param eta0 The eta0 coefficient.
 * @param a3 The a3 coefficient.
 * @param asym The asymmetry parameter.
 * @param gamma_l The mode width.
 * @param l The mode degree.
 * @param step The step size.
 * @param c The truncation constant. It limits the range of the computation. Used to (1) accelerate the fitting and (2) measure the asymetry.
 * @return The optimized Lorentzian model.
 */
VectorXd optimum_lorentzian_calc_a1l_etaa3_v2(const VectorXd& x,  const VectorXd& y,  const VectorXd& H_lm,  const double fc_l,  const double f_s1,  const double f_s2,  const double eta,  const double a3,  const double asym,  const double gamma_l, const int l, double step, const double c);

/**
 * @brief Calculate the optimized Lorentzian model with asymmetry and splitting (model a1etaa3_v2)
 *
 * This function calculates the Lorentzian model with asymmetry and splitting based on the given parameters.
 * This model includes:
 * - Asphericity of Lorentzian a1etaa3
 * - Splitting frequency f_s
 * - Latitudinal effect a3(l=2) only. We consider a3(l=3) = 0
 * - Inclination is not imposed, instead H_l should be of dimension l(l+1) and provide all the heights
 *
 * @param x The frequency range.
 * @param y The original vector.
 * @param H_l The heights of the mode.
 * @param fc_l The central frequency of the mode.
 * @param f_s The splitting frequency.
 * @param eta0 The eta0 coefficient.
 * @param a3 The a3 coefficient.
 * @param asym The asymmetry parameter.
 * @param gamma_l The mode width.
 * @param l The mode degree.
 * @param V The vector V.
 * @param step The step size.
 * @param c The truncation constant. It limits the range of the computation. Used to (1) accelerate the fitting and (2) measure the asymetry.
 * @return The optimized Lorentzian model.
 */
VectorXd optimum_lorentzian_calc_a1etaa3_v2(const VectorXd& x,  const VectorXd& y,  const VectorXd& H_lm,  const double fc_l,  const double f_s,  const double eta,  const double a3,  const double asym,  const double gamma_l, const int l,  const double step, const double c);

/**
 * @brief Calculate the optimized Lorentzian model with asymmetry and splitting (model a1l_a2a3)
 *
 * This function calculates the Lorentzian model with asymmetry and splitting based on the given parameters.
 * This model includes:
 * - Splitting a1(l=1) and a1(l=2). Assumes a1(l=3) = (a1(1) + a1(2))/2.
 * - An Asphericity parameter a2.
 * - Latitudinal effect a3(l=2) only. We consider a3(l=3) = 0.
 * - Inclination is not imposed.
 *
 * @param x The frequency range.
 * @param y The original vector.
 * @param H_l The heights of the mode.
 * @param fc_l The central frequency of the mode.
 * @param f_s1 The splitting frequency for l=1.
 * @param f_s2 The splitting frequency for l=2.
 * @param a2 The a2 coefficient.
 * @param a3 The a3 coefficient.
 * @param asym The asymmetry parameter.
 * @param gamma_l The mode width.
 * @param l The mode degree.
 * @param V The vector V.
 * @param step The step size.
 * @param c The truncation constant. It limits the range of the computation. Used to (1) accelerate the fitting and (2) measure the asymetry.
 * @return The optimized Lorentzian model.
 */
VectorXd optimum_lorentzian_calc_a1l_a2a3(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, const double f_s1, const double f_s2, const double a2, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c);

/**
 * @brief Calculate the optimized Lorentzian model with asymmetry and splitting (model a1a2a3)
 *
 * This function calculates the Lorentzian model with asymmetry and splitting based on the given parameters.
 * This model includes:
 * - Asphericity parameters a1, a2, and a3.
 * - Splitting frequency f_s.
 * - Inclination is not imposed.
 *
 * @param x The frequency range.
 * @param y The original vector.
 * @param H_l The heights of the mode.
 * @param fc_l The central frequency of the mode.
 * @param f_s The splitting frequency.
 * @param a2 The a2 coefficient.
 * @param a3 The a3 coefficient.
 * @param asym The asymmetry parameter.
 * @param gamma_l The mode width.
 * @param l The mode degree.
 * @param V The vector V.
 * @param step The step size.
 * @param c The truncation constant. It limits the range of the computation. Used to (1) accelerate the fitting and (2) measure the asymetry.
 * @return The optimized Lorentzian model.
 */
VectorXd optimum_lorentzian_calc_a1a2a3(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, const double f_s, const double a2, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c);

/**
 * @brief Build a Lorentzian mode with asymmetry and splitting (model a1etaAlma3)
 *
 * This function builds a Lorentzian mode with asymmetry and splitting based on the given parameters.
 * This model includes:
 * - Asymmetry of Lorentzian asym
 * - Splitting a1
 * - An Asphericity term eta (Centrifugal effect) and due to Active region following Gizon 2002, AN, 323, 251.
 *   Currently we use a hard-coded filter type "gate" which is rough but matches the Gizon paper. "gauss" is also available and might be our final choice.
 *   Once we could compare the method adapted from Gizon works on global fits
 * - Latitudinal effect a3
 *
 * @param x_l The frequency range for the mode.
 * @param H_l The height of the mode.
 * @param fc_l The central frequency of the mode.
 * @param f_s The splitting frequency.
 * @param eta0 The asphericity parameter.
 * @param epsilon_nl The epsilon_nl coefficient.
 * @param thetas The thetas coefficients.
 * @param a3 The a3 coefficient.
 * @param asym The asymmetry parameter.
 * @param gamma_l The mode width.
 * @param l The mode degree.
 * @param V The mode visibility.
 * @return The Lorentzian mode.
 */
VectorXd build_l_mode_a1etaAlma3(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s, 
    const double eta0, const double epsilon_nl, const VectorXd& thetas, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V);

/**
 * @brief Calculate the optimized Lorentzian model with asymmetry and splitting (model a1etaAlma3)
 *
 * This function calculates the Lorentzian model with asymmetry and splitting based on the given parameters.
 * This model includes:
 * - Asphericity of Lorentzian a1etaAlma3
 * - Splitting frequency f_s
 * - Latitudinal effect a3(l=2) only. We consider a3(l=3) = 0
 * - Inclination is not imposed, instead H_l should be of dimension l(l+1) and provide all the heights
 *
 * @param x The frequency range.
 * @param y The original vector.
 * @param H_l The heights of the mode.
 * @param fc_l The central frequency of the mode.
 * @param f_s The splitting frequency.
 * @param eta0 The eta0 coefficient.
 * @param epsilon_nl The epsilon_nl coefficient.
 * @param thetas The thetas vector.
 * @param a3 The a3 coefficient.
 * @param asym The asymmetry parameter.
 * @param gamma_l The mode width.
 * @param l The mode degree.
 * @param V The vector V.
 * @param step The step size.
 * @param c The truncation constant. It limits the range of the computation. Used to (1) accelerate the fitting and (2) measure the asymetry.
 * @return The optimized Lorentzian model.
 */
VectorXd optimum_lorentzian_calc_a1etaAlma3(const VectorXd& x, const VectorXd& y,  const double H_l,  const double fc_l,  const double f_s, 
		const double eta0, const double epsilon_nl, const VectorXd& thetas, const double a3,  const double asym,  const double gamma_l, const int l,  const VectorXd& V,  const double step, const double c);

/**
 * @brief Build a Lorentzian mode with asymmetry and splitting (model act_simu)
 *
 * This function builds a Lorentzian mode with asymmetry and splitting based on the given parameters.
 * This model includes:
 * - Asymmetry of Lorentzian asym
 * - Splitting a1
 * - Centrifugal force effect eta
 * - Latitudinal effect a3
 * - Effect of magnetic field of the form b * nu^alpha
 *
 * @param x_l The frequency range for the mode.
 * @param H_l The height of the mode.
 * @param fc_l The central frequency of the mode.
 * @param f_s The splitting frequency.
 * @param eta The eta coefficient.
 * @param a3 The a3 coefficient.
 * @param b The b coefficient.
 * @param alpha The alpha coefficient.
 * @param asym The asymmetry parameter.
 * @param gamma_l The mode width.
 * @param l The mode degree.
 * @param V The mode visibility.
 * @return The Lorentzian mode.
 */
VectorXd build_l_mode_act_simu(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s, const double eta, const double a3, 
        const double b, const double alpha, const double asym, const double gamma_l, const int l, const VectorXd& V); // Only here for compatibiility reasons

/**
 * @brief Build a Lorentzian mode with asymmetry and splitting (model aj)
 *
 * This function builds a Lorentzian mode with asymmetry and splitting based on the given parameters.
 * This model includes:
 * - Asymmetry of Lorentzian asym
 * - Splittings in the form of a-coefficients aj with j=1,6
 * - eta: If eta0 > 0, account for the centrifugal distortion in the a2 term (a2_CF)
 *   This means the measured a2 = a2_AR. It will not include centrifugal effects, but the model does account for it.
 *   In that case, eta0 = 3./(4.*pi*rho*G) should be the value given to that function.
 *   Note that Benomar2018 has a mistake in that formulation (Eq. S6).
 *   If eta0 <= 0, set a2_CF(eta) = 0 such that the measured a2 = a2_CF + a2_AR.
 *
 * @param x_l The frequency range for the mode.
 * @param H_l The height of the mode.
 * @param fc_l The central frequency of the mode.
 * @param a1 The a1 coefficient.
 * @param a2 The a2 coefficient.
 * @param a3 The a3 coefficient.
 * @param a4 The a4 coefficient.
 * @param a5 The a5 coefficient.
 * @param a6 The a6 coefficient.
 * @param eta0 The eta0 coefficient.
 * @param asym The asymmetry parameter.
 * @param gamma_l The mode width.
 * @param l The mode degree.
 * @param V Mode visibilities.
 * */
VectorXd build_l_mode_aj(const VectorXd& x_l, const double H_l, const double fc_l, 
        const double a1, const double a2, const double a3, const double a4, const double a5, const double a6, 
        const double eta, const double asym, const double gamma_l, const int l, const VectorXd& V);


/**
 * @brief Calculate the optimized Lorentzian model with asymmetry and splitting (model aj)
 *
 * This function calculates the Lorentzian model with asymmetry and splitting based on the given parameters.
 * This model includes:
 * - Asphericity parameters a1, a2, a3, a4, a5, and a6.
 * - Eta0 coefficient.
 * - Inclination is not imposed.
 *
 * @param x The frequency range.
 * @param y The original vector.
 * @param H_l The heights of the mode.
 * @param fc_l The central frequency of the mode.
 * @param a1 The a1 coefficient.
 * @param a2 The a2 coefficient.
 * @param a3 The a3 coefficient.
 * @param a4 The a4 coefficient.
 * @param a5 The a5 coefficient.
 * @param a6 The a6 coefficient.
 * @param eta0 The eta0 coefficient.
 * @param asym The asymmetry parameter.
 * @param gamma_l The mode width.
 * @param l The mode degree.
 * @param V The vector V.
 * @param step The step size.
 * @param c The truncation constant. It limits the range of the computation. Used to (1) accelerate the fitting and (2) measure the asymetry.
 * @return The optimized Lorentzian model.
 */
VectorXd optimum_lorentzian_calc_aj(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, 
        const double a1, const double a2, const double a3, const double a4, const double a5, const double a6, 
        const double eta0, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c);
