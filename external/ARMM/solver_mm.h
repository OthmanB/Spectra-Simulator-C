/**
 * @file solver_mm.h
 * @brief Contains functions for solving the asymptotic relation of the mixed modes.
 *
 * This file contains functions that are adapted from the solver_mm.py function. These functions are used to solve the asymptotic relation of the mixed modes. The development of these functions was based on reading papers by Benoit Mosser and Charlotte Gehand.
 *
 * Papers referenced:
 * - Mosser, B., et al. "Mixed modes in red giants: a window on stellar evolution." Astronomy & Astrophysics 540 (2012): A143.
 * - Mosser, B., et al. "Scaling relations for the width of the red-giant branch in the Kepler field." Astronomy & Astrophysics 517 (2010): A22.
 * - Mosser, B., et al. "The universal red-giant oscillation pattern - An automated determination with CoRoT data." Astronomy & Astrophysics 525 (2011): L9.
 * - Mosser, B., et al. "Mixed modes in red giants: a window on stellar evolution." Astronomy & Astrophysics 572 (2014): L5.
 * - Gehan, C. "Étude des modes mixtes dans les étoiles géantes rouges." PhD thesis, Université de Toulouse (2019).
 *
 * Note that the examples and tests functions in this file were built using asymptotic relations in the original Python code. However, these functions should be applicable to any set of values for nu_p(nu), nu_g(nu), Dnu_p(nu), and DPl(nu) (meaning handling glitches).
 */
#pragma once
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>

#include "version_solver.h"
#include "data_solver.h"
#include "string_handler.h"
#include "interpol.h"
#include "derivatives_handler.h"
#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
#endif

using Eigen::MatrixXd;
using Eigen::VectorXi;
using Eigen::VectorXd;

/**
 * @brief Removes duplicate values from a given EigenVectorXd.
 * @param nu_m_all The input EigenVectorXd.
 * @param tol The tolerance for considering two values as duplicates.
 * @return The EigenVectorXd with unique values.
 */
Eigen::VectorXd removeDuplicates(const Eigen::VectorXd& nu_m_all, double tol);

/**
 * @brief Function that detects sign changes in a given vector.
 * If the sign went from + to -, tag it with a -1.
 * If the sign went from - to +, tag it with a +1.
 * If there is no change of sign, tag it with a 0.
 * 0 is dealt as a zone of change of sign as well. eg. if we pass from 0 to 2, then the result is +1.
 * @param x The input vector for which we want to know sign changes.
 * @param return_indices If true (default), the function returns positions at which the sign changed.
 *                       If false, it returns a vector of tags for the sign changes.
 * @return The vector of sign changes or the positions at which the sign changed.
 */
VectorXi sign_change(const VectorXd& x, bool return_indices=true);

/**
 * @brief Calculates the difference between a single value and a given value.
 * @param nu The input value.
 * @param nu_p The value to subtract from the input value.
 * @return The resulting value.
 */
VectorXd pnu_fct(const VectorXd& nu, const long double nu_p);

/**
 * @brief Calculates the difference between a single value and a given value.
 * @param nu The input value.
 * @param nu_p The value to subtract from the input value.
 * @return The resulting value.
 */
long double pnu_fct(const long double nu, const long double nu_p);

/**
 * @brief Calculates the gnu function for each element of the input vector.
 * @param nu The input vector.
 * @param nu_g The value of nu_g.
 * @param Dnu_p The value of Dnu_p.
 * @param DPl The value of DPl.
 * @param q The value of q.
 * @return The resulting vector.
 */
VectorXd gnu_fct(const VectorXd& nu, const long double nu_g, const long double Dnu_p, const long double DPl, const long double q);

/**
 * @brief Calculates the gnu function for a single value.
 * @param nu The input value.
 * @param nu_g The value of nu_g.
 * @param Dnu_p The value of Dnu_p.
 * @param DPl The value of DPl.
 * @param q The value of q.
 * @return The resulting value.
 */
long double gnu_fct(const long double nu, const long double nu_g, const long double Dnu_p, const long double DPl, const long double q);

/**
 * @brief A small function that generates a series of p modes using the asymptotic relation at the second order.
 * @param Dnu_p The value of Dnu_p.
 * @param np The value of np.
 * @param epsilon The value of epsilon.
 * @param l The value of l.
 * @param delta0l The value of delta0l.
 * @param alpha The value of alpha.
 * @param nmax The value of nmax.
 * @param r The extra term to add to Dnu_p (default is 0).
 * @return The resulting value of nu_p.
 */
long double asympt_nu_p(const long double Dnu_p, const int np, const long double epsilon, const int l, 
	const long double delta0l, const long double alpha, const long double nmax, long double r=0);

/**
 * @brief A small function that generates a series of p modes based on a shifting of a series of l=0 modes and on the asymptotic relation.
 * @param nu_l0 The input vector of l=0 modes.
 * @param Dnu_p The value of Dnu_p.
 * @param np The value of np.
 * @param epsilon The value of epsilon.
 * @param l The value of l.
 * @param delta0l The value of delta0l.
 * @param r The extra term to add to Dnu_p (default is 0).
 * @return The resulting value of nu_p.
 */
long double asympt_nu_p_from_l0(const VectorXd& nu_l0, const long double Dnu_p, const int np, const long double epsilon, const int l, 
	const long double delta0l, long double r=0);

/**
 * @brief A small function that generates a series of p modes based on a shifting of a series of l=0 modes and on the asymptotic relation.
 * @param nu_l0 The input vector of l=0 modes.
 * @param Dnu_p The value of Dnu_p.
 * @param l The value of l.
 * @param delta0l The value of delta0l.
 * @param fmin The minimum frequency value (default is -1).
 * @param fmax The maximum frequency value (default is -1).
 * @return The resulting vector of nu_p.
 */
VectorXd asympt_nu_p_from_l0_Xd(const VectorXd& nu_l0, const long double Dnu_p, const int l, const long double delta0l, long double fmin=-1, long double fmax=-1);

/**
 * @brief Compute the asymptotic relation for the g modes.
 * @param DPl The value of DPl.
 * @param ng The value of ng.
 * @param alpha The value of alpha.
 * @param r An optional parameter that can be added to the Period (e.g. a random quantity).
 * @return The resulting value of nu_g.
 */
long double asympt_nu_g(const long double DPl, const int ng, const long double alpha, long double r=0);

/**
 * @brief Solve the mixed mode asymptotic relation between p modes and g modes.
 * @param nu_p The frequency of a given p-mode (in microHz).
 * @param nu_g The frequency of a given g-mode (in microHz).
 * @param Dnu_p The large separation for p modes (in microHz).
 * @param DPl The period spacing for g modes (in seconds).
 * @param q The coupling term (no unit, should be between 0 and 1).
 * @param numin The minimum frequency considered for the solution (in microHz).
 * @param numax The maximum frequency considered for the solution (in microHz).
 * @param resol The base resolution for the interpolated base function.
 * @param returns_axis If true, returns nu, pnu, and gnu.
 * @param verbose If true, print additional information.
 * @param factor Define how fine the new tiny grid used for performing the interpolation will be.
 * @return A structure with the solutions that match p(nu) = g(nu).
 */
Data_coresolver solver_mm(const long double nu_p, const long double nu_g, const long double Dnu_p, const long double DPl, const long double q, 
	const long double numin, const long double numax, const long double resol, const bool returns_axis=false, const bool verbose=false, const long double factor=0.05);

/**
 * @brief Solve the mixed mode asymptotic relation between p modes and g modes.
 * @param Dnu_p The large separation for p modes (in microHz).
 * @param epsilon The value of epsilon.
 * @param el The value of el.
 * @param delta0l The value of delta0l.
 * @param alpha_p The value of alpha_p.
 * @param nmax The value of nmax.
 * @param DPl The period spacing for g modes (in seconds).
 * @param alpha The value of alpha.
 * @param q The coupling term (no unit, should be between 0 and 1).
 * @param sigma_p The value of sigma_p.
 * @param fmin The minimum frequency considered for the solution (in microHz).
 * @param fmax The maximum frequency considered for the solution (in microHz).
 * @param resol The base resolution for the interpolated base function.
 * @param returns_pg_freqs If true, returns nu_m, nu_p, nu_g, dnup, and dPg.
 * @param verbose If true, print additional information.
 * @return A structure with the solutions that match p(nu) = g(nu).
 */
Data_eigensols solve_mm_asymptotic_O2p(const long double Dnu_p, const long double epsilon, const int el, const long double delta0l, const long double alpha_p, 
	const long double nmax, const long double DPl, const long double alpha, const long double q, const long double sigma_p, 
	const long double fmin, const long double fmax, const long double resol, bool returns_pg_freqs=true, bool verbose=false);

/**
 * @brief Solve the mixed mode asymptotic relation between p modes and g modes using l=0 frequencies as input.
 * @param nu_l0_in Frequencies for the l=0 modes. Used to derive nu_l1_p and therefore Dnu and epsilon.
 * @param el Degree of the mode.
 * @param delta0l First order shift related to core structure (and to D0).
 * @param DPl Average Period spacing of the g modes.
 * @param alpha Phase offset for the g modes.
 * @param q Coupling strength.
 * @param sigma_p Standard deviation controlling the randomization of individual p modes. Set it to 0 for no spread.
 * @param resol Control the grid resolution. Might be set to the resolution of the spectrum.
 * @param returns_pg_freqs If true, returns the values for calculated p and g modes.
 * @param verbose If true, print the solution on screen.
 * @param freq_min The minimum frequency considered for the solution (in microHz).
 * @param freq_max The maximum frequency considered for the solution (in microHz).
 * @return A structure with the solutions that match p(nu) = g(nu).
 */
Data_eigensols solve_mm_asymptotic_O2from_l0(const VectorXd& nu_l0_in, const int el, const long double delta0l, 
    const long double DPl, const long double alpha, const long double q, const long double sigma_p, 
	const long double resol, bool returns_pg_freqs=true, bool verbose=false, const long double freq_min=0, const long double freq_max=1e6);

/**
 * @brief Solve the mixed mode asymptotic relation between p modes and g modes using l=0 frequencies as input.
 * @param nu_p_all Frequencies for the l modes.
 * @param el Degree of the mode.
 * @param DPl Average Period spacing of the g modes.
 * @param alpha Phase offset for the g modes.
 * @param q Coupling strength.
 * @param sigma_p Standard deviation controlling the randomization of individual p modes. Set it to 0 for no spread.
 * @param resol Control the grid resolution. Might be set to the resolution of the spectrum.
 * @param returns_pg_freqs If true, returns the values for calculated p and g modes.
 * @param verbose If true, print the solution on screen.
 * @param freq_min The minimum frequency considered for the solution (in microHz).
 * @param freq_max The maximum frequency considered for the solution (in microHz).
 * @return A structure with the solutions that match p(nu) = g(nu).
 */
Data_eigensols solve_mm_asymptotic_O2from_nupl(const VectorXd& nu_p_all, const int el, 
    const long double DPl, const long double alpha, const long double q, const long double sigma_p, 
	const long double resol, bool returns_pg_freqs=true, bool verbose=false, const long double freq_min=0, const long double freq_max=1e6);
