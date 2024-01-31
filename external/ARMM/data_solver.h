/**
 * @file data_solver.h
 * @brief Header file that contains all kinds of classes/structures used to process and/or encapsulate data.
 *
 * This file contains the declarations for various classes and structures used to process and encapsulate data. It provides functionality for mixed modes calculation, simulations, core solver data, eigen solutions data, rotation to zone data, envelope asphericity data, synthetic star configuration data, and synthetic star parameter data.
 *
 * @date 22 Feb 2016
 * @author obenomar
 */

#pragma once
#include <Eigen/Dense>
#include <string>
#include <vector>
#include "gnuplot-iostream.h"

using Eigen::MatrixXd;
using Eigen::VectorXi;
using Eigen::VectorXd;

// ----------------------------------------
// ----- For mixed modes calculation ------
// ----------------------------------------

/**
 * @struct Data_2vectXd
 * @brief Structure for holding two VectorXd objects.
 *
 * This structure is used for mixed modes calculation and holds two VectorXd objects.
 */
struct Data_2vectXd{
	VectorXd vecXd1; ///< Arbitrary vector.
	VectorXd vecXd2; ///< Another arbitrary vector.
};

/**
 * @struct template_file
 * @brief Structure for simulations to derive height and widths of modes.
 *
 * This structure is used for simulations only and is used to derive height and widths of modes. It contains reference data for simulations.
 */
struct template_file{
	std::string ID_ref; ///< Star Identification Number.
	double numax_ref; ///< Frequency at maximum power of the reference star.
	double Dnu_ref;  ///< Large frequency separation of the reference star.
	double epsilon_ref; ///< epsilon_p (phase offset) for the reference star.
	MatrixXd data_ref; ///< Matrix containing detailled information on modes (frequencies, widths, heights,...) for the reference star
};

/**
 * @struct Data_coresolver
 * @brief Structure for core solver data.
 *
 * This structure is used to hold core solver data, including vectors for nu_m, ysol, nu, pnu, and gnu.
 */
struct Data_coresolver{
	VectorXd nu_m; ///< Frequencies of the mixed modes.
	VectorXd ysol;  ///< Full List of frequencies solutions found by the solver.
	VectorXd nu; ///< Frequencies used to compute the interpolation grid.
	VectorXd pnu; ///< Terms including the p modes in the implicit equation for the p modes.
	VectorXd gnu; ///< Terms including the g modes in the implicit equation for the mixed modes.
};

/**
 * @struct Data_eigensols
 * @brief Structure for eigen solutions data.
 *
 * This structure is used to hold eigen solutions data, including vectors for nu_p, nu_g, nu_m, dnup, and dPg.
 */
struct Data_eigensols{
	VectorXd nu_p; ///< Frequencies of p modes used to compute mixed modes.
	VectorXd nu_g; ///< Frequencies of g modes used to compute mixed modes.
	VectorXd nu_m; ///< Frequencies of mixed modes.
	VectorXd dnup; ///< Frequency derivatives for p modes ~ &Delta&nu;.
	VectorXd dPg;  ///< Frequency derivatives for g modes ~ &Delta&Pi.
};

/**
 * @struct Data_eigensols_all
 * @brief Structure for all eigen solutions data.
 *
 * This structure is used to hold all eigen solutions data, including a vector for nu_l0.
 */
//struct Data_eigensols_all{
//	VectorXd nu_l0;
//};


/**
 * @struct Data_rot2zone
 * @brief Structure for rotation to zone data.
 *
 * This structure is used to hold rotation to zone data, including variables for rot_core and rot_env.
 */
struct Data_rot2zone{
	long double rot_core; ///< Average core rotation rate.
	long double rot_env; ///< Average envelope rotation rate.
};

/**
 * @struct Envelope_lat_dif_rot
 * @brief Structure for envelope lateral differential rotation data.
 *
 * This structure is used to hold envelope lateral differential rotation data, including variables for a3_l2, a3_l3, and a5_l3.
 */
struct Envelope_lat_dif_rot{ // Added on 13 Sep 2023
	long double a3_l2=0; ///< a3(l=2) coefficient.
	long double a3_l3=0; ///< a3(l=3) coefficient.
	long double a5_l3=0; ///< a5(l=3) coefficient.
};

/**
 * @struct Envelope_asphericity
 * @brief Structure for envelope asphericity data.
 *
 * This structure is used to hold envelope asphericity data, including variables for a2_l1, a2_l2, a2_l3, a4_l2, a4_l3, and a6_l3.
 */
struct Envelope_asphericity{
	long double a2_l1; ///< a2(l=1) coefficient.
	long double a2_l2; ///< a2(l=2) coefficient.
	long double a2_l3; ///< a2(l=3) coefficient.
	long double a4_l2; ///< a4(l=2) coefficient.
	long double a4_l3; ///< a4(l=3) coefficient.
	long double a6_l3; ///< a6(l=3) coefficient.
};

/**
 * @struct Cfg_synthetic_star
 * @brief Structure for synthetic star configuration data.
 *
 * This structure is used to hold synthetic star configuration data, including variables for various parameters such as Teff_star, numax_star, Dnu_star, epsilon_star, delta0l_percent_star, beta_p_star, alpha_p_star, nmax_star, DPl_star, alpha_g_star, q_star, fmin, fmax, maxHNR_l0, noise_params_harvey_like, Gamma_max_l0, rot_env_input, rot_ratio_input, rot_core_input, env_lat_dif_rot, env_aspher, output_file_rot, Vl, H0_spread, filetemplate, sigma_p, sigma_m, Hfactor, Wfactor, inclination, nu_nl, Nf_el, and use_nu_nl.
 */
struct Cfg_synthetic_star{
	long double Teff_star;  ///< Effective temperature of the star.
	long double numax_star;  ///< Frequency of maximum power of the star.
	long double Dnu_star; ///< Large frequency separation of the star.
	long double epsilon_star; ///< Epsilon_p (p modes phase offset) of the star.
	long double delta0l_percent_star; ///< Small frequency separation between l=0 and l.
	long double beta_p_star; ///< Normalised 2nd Order p modes Curvature Beta_p of the star (correction from the asymptotic).
	long double alpha_p_star; ///< 2nd Order p modes Curvature Beta_p of the star (correction from the asymptotic) as defined in Mosser+2017.
	long double nmax_star; ///< Fractional adial order at numax of the star.
	long double DPl_star; ///< Period Spacing for g modes.
	long double alpha_g_star; ///< epsilon_g (g modes phase offset).
	long double q_star; ///< Coupling strength between p and g modes.
	long double fmin;  ///< Low bound for frequencies.
	long double fmax; ///< High bound for frequencies.
	long double maxHNR_l0; ///< Maximum Height to Noise ratio (at numax) for l=0 modes.
	VectorXd noise_params_harvey_like; ///< Noise parameters (Usually, Harvey-like parameters).
	long double Gamma_max_l0; ///< Width at numax for l=0.
	long double rot_env_input; ///< average stellar envelope rotation.
	long double rot_ratio_input; ///< average ratio between core and envelope rotation rot_core/rot_env.
	long double rot_core_input; ///< average stellar core rotation .
	Envelope_lat_dif_rot env_lat_dif_rot; ///< Structure with odds a-coeficients j={3,5} (latitudinal differential information).
	Envelope_asphericity env_aspher; ///< Structure with the even a-coeficients j={2,4,6} (asphericity information).
	std::string output_file_rot; ///< Output file for rotation information.
	VectorXd Vl; ///< Bolometric mode visibility for l={1,2,3}.
	long double H0_spread; ///< Gaussian dispersion in % of the height of l=0 modes. Allows to introduce some randomness to the artificially created star.
	std::string filetemplate; ///< Template file describing the reference star parameters and used for constructing the artificial star.
	long double sigma_p; ///< Randomness coefficient in p modes distribution (OBSELETE).
	long double sigma_m; ///< Randomness coefficient in g modes distribution (OBSELETE).
	long double Hfactor; ///< Empirical correction factor to the asymptotic for l=1 heights (mixed modes).
	long double Wfactor; ///< Empirical correction factor to the asymptotic for l=1 widths  (mixed modes).
	long double inclination; ///< Stellar inclination
	MatrixXd nu_nl; ///< Frequencies of the modes, are here if provided by a template (e.g a theoretical model) and handled by the MCMC model
	VectorXi Nf_el; ///< The numbers of modes 
	bool use_nu_nl=false; ///< If set to true, use the nu_nl frequencies instead of computing them from the asymptotic. These must be set
    bool legacynoise=true; ///< Defines if the noise is strictly speaking a sum of Harvey-like + white noise (legacynoise=false) or if it is a single Harvey + white noise scaling quantity with only numax (legacynoise=true)
};


/**
 * @struct Params_synthetic_star
 * @brief Structure for synthetic star parameter data.
 *
 * This structure is used to hold synthetic star parameter data, including vectors for various parameters such as nu_l0, nu_p_l1, nu_g_l1, nu_m_l1, nu_l2, nu_l3, width_l0, width_l1, width_l2, width_l3, height_l0, height_l1, height_l2, height_l3, a1_l1, a1_l2, a1_l3, a2_l1, a2_l2, a2_l3, a3_l2, a3_l3, a4_l2, a4_l3, a5_l3, and a6_l3.
 */
struct Params_synthetic_star{
	bool failed=false; ///< Variable that is set to true if the mixed modes computation somewhat failed
	VectorXd nu_l0; ///< Frequencies of l=0 p modes.
	VectorXd nu_p_l1; ///< Frequencies of l=1 p modes used in the compuation of the mixed modes.
	VectorXd nu_g_l1;  ///< Frequencies of l=1 g modes used in the compuation of the mixed modes.
	VectorXd nu_m_l1;  ///< Frequencies of the mixed modes, solutions of the implicit asymptotic equation for mixed modes.
	VectorXd nu_l2;  ///< Frequencies of l=2 p modes.
	VectorXd nu_l3;  ///< Frequencies of l=3 p modes.
	VectorXd width_l0;  ///< Widths of l=0 p modes.
	VectorXd width_l1;  ///< Widths of l=1 mixed modes.
	VectorXd width_l2;   ///< Widths of l=2 p modes.
	VectorXd width_l3;  ///< Widths of l=3 p modes.
	VectorXd height_l0;  ///< Heights of l=0 p modes.
	VectorXd height_l1;  ///< Heights of l=1 mixed modes.
	VectorXd height_l2;  ///< Heights of l=2 p modes.
	VectorXd height_l3; ///< Heights of l=3 p modes.
	VectorXd a1_l1;  ///< a1(n,l=1) coefficients.
	VectorXd a1_l2;  ///< a1(n,l=2) coefficients.
	VectorXd a1_l3; ///< a1(n,l=3) coefficients.
	VectorXd a2_l1; ///< a2(n,l=1) coefficients.
	VectorXd a2_l2; ///< a2(n,l=2) coefficients.
	VectorXd a2_l3; ///< a2(n,l=3) coefficients.
	VectorXd a3_l2; ///< a3(n,l=2) coefficients.
	VectorXd a3_l3; ///< a3(n,l=3) coefficients.
	VectorXd a4_l2; ///< a4(n,l=2) coefficients.
	VectorXd a4_l3; ///< a4(n,l=3) coefficients.
	VectorXd a5_l3; ///< a5(n,l=3) coefficients.
	VectorXd a6_l3; ///< a6(n,l=3) coefficients.
};