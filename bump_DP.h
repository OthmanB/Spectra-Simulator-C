// -------------------
// ---- Functions adapted from the bump_DP.py function ----
// This contains all the function that describes the Bumped period spacing
// and how they can be used to derived mode amplitudes from inertia ratio
// as they have been developed and tested. This arises from the following publications:
// https://arxiv.org/pdf/1509.06193.pdf (Inertia and ksi relation)
// https://www.aanda.org/articles/aa/pdf/2015/08/aa26449-15.pdf (Eq. 17 for the rotation - splitting relation)
// https://arxiv.org/pdf/1401.3096.pdf (Fig 13 and 14 for the evolution - rotation relation in SG and RGB) 
// https://arxiv.org/pdf/1505.06087.pdf (Eq. 3 for determining log(g) using Teff and numax, used for getting the evolution stage in conjonction with Fig. 13 from above)
// https://iopscience.iop.org/article/10.1088/2041-8205/781/2/L29/pdf (Inertia and Height relation)
// https://arxiv.org/pdf/1707.05989.pdf (Fig.5, distribution of surface rotation for 361 RGB stars)

// ------------------
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "version_solver.h"
#include "data.h"
//#include "string_handler.h" // REPLACED BY IOPROC.h
#include "ioproc.h"
#include "interpol.h"
#include "noise_models.h" // get the harvey_1985 function
#include "solver_mm.h"

/*
# the ksi function as defined in equation 14 of Mosser+2017 (https://arxiv.org/pdf/1509.06193.pdf)
# Inputs: 
# 	nu: The freqency axis (microHz)
# 	nu_p : Frequenc of a p-mode (microHz)
#	nu_g : Frequency of a g-mode (microHz)
#	Dnu_p: Large separation. DOES NOT NEED TO BE A CONSTANT (fonction of np or nu to account for glitches)
#   DPl: Period Spacing (seconds)
#	q : Coupling term (no unit)
*/
VectorXd ksi_fct1(const VectorXd& nu, const long double nu_p, const long double nu_g, const long double Dnu_p, const long double DPl, const long double q);

/*
# Variant of ksi_fct that deals with arrays for nu_p, nu_g, Dnu_p, DPl
# This requires that nu_p and Dnu_p have the same dimension
# Also nu_g and DPl must have the same dimension
# Additional parameter:
#  - norm-method: When set to "fast", normalise by the max of the ksi_pg calculated at the
#				   Frequencies nu given by the user
#				   When set to 'exact', normalise by the max of a heavily interpolated function
#				   of ksi_pg. Allows a much a higher precision, but will be slower
#				   This could be usefull as in case of low ng, the norm is badly estimated in
#				   "fast" mode. Then we need to use a more continuous function to evaluate the norm
*/
VectorXd ksi_fct2(const VectorXd& nu, const VectorXd& nu_p, const VectorXd& nu_g, const VectorXd& Dnu_p, const VectorXd& DPl, const long double q, const std::string norm_method="fast");

VectorXd gamma_l_fct2(const VectorXd& ksi_pg, const VectorXd& nu_m, const VectorXd& nu_p_l0, const VectorXd& width_l0, const VectorXd& hl_h0_ratio, const int el);

VectorXd h_l_rgb(const VectorXd& ksi_pg);

// Put here the code for reading template files that contain heights and width profiles
template_file read_templatefile(const std::string file, const bool ignore_errors=true);

Data_2vectXd width_height_load_rescale(const VectorXd& nu_star, const long double Dnu_star, const long double numax_star, const std::string file);

// Simple way of computing the core rotation from the surface rotation. Used if we want uniform 
// distribution of rotation in the envelope and a uniform population of core-to-envelope ratios 
// 	 (1) rot_envelope: average rotation in the envelope
//	 (2) core2envelope_star: average rotation in the core 
Data_rot2zone rot_2zones_v2(const long double rot_envelope, const long double core2envelope_star, std::string output_file_rot=" ");

// Simple way of computing the core rotation from the surface rotation. Used if we want uniform 
// distribution of rotation in the envelope and a uniform population of core rotation
// 	 (1) rot_envelope: average rotation in the envelope
//	 (2) rot_core: average rotation in the core 
Data_rot2zone rot_2zones_v3(const long double rot_envelope, const long double rot_core, std::string output_file_rot=" ");

// Function that determine the rotation in the envelope, here approximated to be the surface rotation.
// Inspired by the surface rotation from Ceillier et al. 2017 (https://arxiv.org/pdf/1707.05989.pdf), Fig. 5
// They have a skewed distribution going for ~30 days to ~160 days with a peak around 60 days. 
// For simplicity, I just insert a truncated gaussian distribution with rotation between 30  and 90 and a median of 60.
// The truncation happens at sigma. Values are given in days
// Returns: 
//	rot_s: rotation frequency in microHz
long double rot_envelope(long double med=60., long double sigma=3.);

/*
# Splitting of modes assuming a two-zone (or two-zone averaged) rotation profile
# rot_envelope and rot_core must be in Hz units (or micro, nano, etc...)
# ksi_pg: The ksi function that describes the degree of mixture between p and g modes in function of the more frequency
# rot_envelope: average rotation in the envelope. Must be a scalar
# rot_core: average rotation in the core. Must be a scalar
# Returns:
#	dnu_rot: A vector of same size as ksi_pg
*/
VectorXd dnu_rot_2zones(const VectorXd& ksi_pg, const long double rot_envelope, const long double rot_core);

//Assumptions: nu_max is derived from the requested Dnu_star parameter using the relation from Stello+2009. 
//	Dnu ~ 0.263*numax^0.77 (no uncertainty implemented here)
long double numax_from_stello2009(const long double Dnu_star, const long double spread);


/*
# The main function that generate a set of parameters used to generate Lorentzian profiles FOR SIMULATIONS IN THE C++ SIMULATOR
# Assumptions: 
#     - The frequencies of l=0, l=2 and l=3 modes follow exactly the asymtptotic relation for the p mdoes
#	  - The frequencies of the l=1 modes follow exactly the asymptotitc relation for the mixed modes
#	  - Widths of l=0, 2, 3 modes are rescaled using the synthetic relation from Appourchaux+2014 applied to the solar profile
#	  - Widths of l=1 mixed modes are defined using the ksi function, scaled using l=0 modes. Thus l=0 modes fixes the upper
#		limit of the mode width
#	  - Heights of l=0 modes are rescaled using the measured heights of the solar modes
#	  - Heights of l=1 mixed modes are set based on equipartition of energy accounting for damping and excitation assuming no radiative pressure
#		GrosJean+201X showed that this was valid for not-too-evolved RGB and for SG.
#	  - Bolometric visibilities in height are assumed to be 1, 1.5, 0.5, 0.07 for l=0,1,2,3 respectively
#	  - Splitting is implemented in different ways... see the various models
# 	  Warning for widths and heights: if fmin/fmax is too large, you may have an error as an extrapolation 
#									  would be required, which I forbid by code. An fmin/fmax that englobes
#									  10 radial orders should not pose any issue.
# Input: 
#	Dnu_star: Large separation for the p modes
#   epsilon_star: phase offset term for the asymtptotic relation of the p modes
#   delta0l_star, alpha_p_star, nmax_star: Instead of D0_star, these parameters can be used to create 2nd order asymptotic p modes (see Mosser+2018, Eq.22) 
#   DPl_star: The period spacing for l=1 g modes 
#	alpha_star: The phase offset term for the asymptotic relation of the g modes
#   q_star: Coupling coeficient between p and g modes
#	fmin: Minimum frequency for the modes that should be included in the calculations
#   fmax: Maximum frequency for the modes that should be included in the calculations
# Outputs:
#	nu_lx: Frequencies of the l=x modes. x is between 0 and 3
#	nu_p_l1: Base p modes frequencies used to build the frequencies for the l=1 mixed modes
#	nu_g_l1: Base p modes frequencies used to build the frequencies for the l=1 mixed modes
#	width_lx: Widths of the l=x modes. x is between 0 and 3
#   height_lx: Heights of the l=x modes. x is between 0 and 3 
*/

Params_synthetic_star make_synthetic_asymptotic_star(Cfg_synthetic_star cfg_star);

Cfg_synthetic_star test_make_synthetic_asymptotic_star_sg(void);

Cfg_synthetic_star test_make_synthetic_asymptotic_star_rgb(void);
