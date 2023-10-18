/** 
 * @file bump_DP.h
 * @brief Functions computing the full profile of the mixed modes as well as artificial asteroseismic parameters for a RGB stars.
 *
 * This file contains all the functions that describe the Bumped period spacing and how they can be used to derive mode amplitudes from inertia ratio. The functions have been developed and tested based on the following publications:
 * - [Inertia and ksi relation](https://arxiv.org/pdf/1509.06193.pdf)
 * - [Eq. 17 for the rotation - splitting relation](https://www.aanda.org/articles/aa/pdf/2015/08/aa26449-15.pdf)
 * - [Fig 13 and 14 for the evolution - rotation relation in SG and RGB](https://arxiv.org/pdf/1401.3096.pdf)
 * - [Eq. 3 for determining log(g) using Teff and numax, used for getting the evolution stage in conjunction with Fig. 13 from above](https://arxiv.org/pdf/1505.06087.pdf)
 * - [Inertia and Height relation](https://iopscience.iop.org/article/10.1088/2041-8205/781/2/L29/pdf)
 * - [Fig.5, distribution of surface rotation for 361 RGB stars](https://arxiv.org/pdf/1707.05989.pdf)
 */
#pragma once 
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "version_solver.h"
#include "data_solver.h"
#include "string_handler.h"
#include "interpol.h"
#include "noise_models.h" // get the harvey_1985 function
#include "solver_mm.h"
#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
#endif

/**
 * @brief Calculates the ksi function as defined in equation 14 of Mosser+2017 (https://arxiv.org/pdf/1509.06193.pdf).
 *
 * This function calculates the zeta function using the given parameters and equations.
 *
 * @param nu The frequency axis in microHz.
 * @param nu_p The frequency of a p-mode in microHz.
 * @param nu_g The frequency of a g-mode in microHz.
 * @param Dnu_p The large separation, which can be a function of np or nu to account for glitches.
 * @param DPl The period spacing in seconds.
 * @param q The coupling term, without any unit.
 * @return The ksi vector which is the local zeta function.
 */
VectorXd ksi_fct1(const VectorXd& nu, const long double nu_p, const long double nu_g, const long double Dnu_p, const long double DPl, const long double q);

/**
 * @brief Calculates the zeta function for mixed modes in precise mode.
 *
 * This function calculates the ksi_pg vector using the ksi_fct1 function for each combination of nu_p and nu_g.
 *
 * @param nu The frequency vector.
 * @param nu_p The list of p modes in a vector form.
 * @param nu_g The list of g modes in a vector form.
 * @param Dnu_p The local values of the large separation, through a vector.
 * @param DPl The local values of the period spacing in the form of a vector.
 * @param q The coupling strength q between p and g modes.
 * @param norm_method The method for normalizing the ksi_pg vector.
 *                    Possible values: "fast" or "exact". In the exact mode, the ksi_fct2_precise function is used.
 * @return The ksi_pg vector which is the zeta function intervening in mixed modes.
 */
VectorXd ksi_fct2(const VectorXd& nu, const VectorXd& nu_p, const VectorXd& nu_g, const VectorXd& Dnu_p, const VectorXd& DPl, const long double q, const std::string norm_method="fast");

/**
 * @brief Calculates the zeta function for mixed modes in precise mode.
 *
 * This function calculates the ksi_pg vector using the ksi_fct1 function for each combination of nu_p and nu_g.
 * It also performs parallel computation using OpenMP for better performance.
 *
 * @param nu The frequency vector.
 * @param nu_p The list of p modes in a  vector form.
 * @param nu_g The list of g modes in a vector form.
 * @param Dnu_p The local values of the large separation, through a vector.
 * @param DPl The local values of the period spacing in the form of a vector.
 * @param q The coupling strength q between p and g modes.
 * @return The ksi_pg vector which is the zeta function intervening in mixed modes.
 */
Eigen::VectorXd ksi_fct2_precise(const Eigen::VectorXd& nu, const Eigen::VectorXd& nu_p, const Eigen::VectorXd& nu_g, const Eigen::VectorXd& Dnu_p, const Eigen::VectorXd& DPl, const long double q);

/**
 * @brief Calculates the mode width of mixed modes based on the given ksi_pg, nu_m, nu_p_l0, width_l0, hl_h0_ratio, el, and a factor.
 *
 * This function calculates the width_l vector by performing interpolation using the given ksi_pg, nu_m, nu_p_l0, width_l0, hl_h0_ratio, el, and factor.
 * The width0_at_l is calculated using linear interpolation between nu_p_l0 and width_l0 for each element in nu_m.
 * The width_l vector is then calculated by multiplying width0_at_l with (1. - factor*ksi_pg[i]) and dividing by the square root of hl_h0_ratio[i].
 *
 * @param ksi_pg The ksi_pg vector.
 * @param nu_m The nu_m vector.
 * @param nu_p_l0 The nu_p_l0 vector.
 * @param width_l0 The width_l0 vector.
 * @param hl_h0_ratio The hl_h0_ratio vector.
 * @param el The el value.
 * @param factor The factor value.
 * @return The width_l vector.
 *
 * @note The size of nu_p_l0 and width_l0 must be the same.
 * @note The size of ksi_pg and nu_m must be the same.
 * @note If the sizes are inconsistent, the function will print an error message and exit the program.
 */
VectorXd gamma_l_fct2(const VectorXd& ksi_pg, const VectorXd& nu_m, const VectorXd& nu_p_l0, const VectorXd& width_l0, const VectorXd& hl_h0_ratio, const int el, const long double factor=1.0);

/**
 * @brief Calculates the H(l)/H(l=0) ratio based on the given ksi_pg and a factor.
 *
 * This function calculates the h_l_rgb vector by subtracting factor*ksi_pg from a vector of ones, and then taking the square root of each element.
 * If any element in the resulting vector is smaller than a given tolerance, it is set to a value close to 0 at the machine precision level to avoid division by 0 or infinity.
 *
 * @param ksi_pg The ksi_pg vector.
 * @param factor The factor value.
 * @return The h_l_rgb vector.
 */
VectorXd h_l_rgb(const VectorXd& ksi_pg,  const long double factor=1.0);

/**
 * @brief Reads a template file and extracts relevant information.
 *
 * This function reads a template file and extracts the ID reference, numax reference, Dnu reference, epsilon reference, and data reference from the file. The file should be formatted with key-value pairs separated by a delimiter. Comments in the header are ignored.
 *
 * @param file The path to the template file.
 * @param ignore_errors Flag indicating whether to ignore errors if the required keywords are not found in the file.
 * @return A template_file object containing the extracted information.
 */
template_file read_templatefile(const std::string file, const bool ignore_errors=true);

/**
 * @brief Loads and rescales data from a file based on reference values.
 *
 * This function loads data from a file and rescales it based on reference values. The rescaling is done using interpolation and the extracted information from the template file. The rescaled data is then returned as a Data_2vectXd object.
 *
 * @param nu_star The reference frequencies.
 * @param Dnu_star The reference Dnu value.
 * @param numax_star The reference numax value.
 * @param file The path to the template file.
 * @return A Data_2vectXd object containing the rescaled data.
 */
Data_2vectXd width_height_load_rescale(const VectorXd& nu_star, const long double Dnu_star, const long double numax_star, const std::string file);

/**
 * @brief Computes the core rotation from the surface rotation.
 *
 * This function computes the core rotation based on the given surface rotation and the average core-to-envelope ratio. It is used to achieve a uniform distribution of rotation in the envelope and a uniform population of core-to-envelope ratios.
 *
 * @param rot_envelope The average rotation in the envelope (in microHz).
 * @param core2envelope_star The average rotation in the core.
 * @param output_file_rot The path to the output file to save the computed rotation values.
 * @return A Data_rot2zone object containing the computed rotation values.
 */
Data_rot2zone rot_2zones_v2(const long double rot_envelope, const long double core2envelope_star, std::string output_file_rot=" ");

/**
 * @brief Computes the core rotation from the surface rotation.
 *
 * This function computes the core rotation based on the given surface rotation and the average core-to-envelope ratio. It is used to achieve a uniform distribution of rotation in the envelope and a uniform population of core-to-envelope ratios.
 *
 * @param rot_envelope The average rotation in the envelope (in microHz).
 * @param rot_core The average rotation in the core (in microHz).
 * @param output_file_rot The path to the output file to save the computed rotation values.
 * @return A Data_rot2zone object containing the computed rotation values.
 */
Data_rot2zone rot_2zones_v3(const long double rot_envelope, const long double rot_core, std::string output_file_rot=" ");

/**
 * @brief Computes the rotation in the envelope based on a truncated Gaussian distribution.
 *
 * This function determines the rotation in the envelope based on a truncated Gaussian distribution. The distribution is inspired by the surface rotation from Ceillier et al. 2017 (https://arxiv.org/pdf/1707.05989.pdf), Fig. 5, which has a skewed distribution ranging from ~30 days to ~160 days with a peak around 60 days. For simplicity, a truncated Gaussian distribution is used with rotation values between 30 and 90 days and a median of 60 days. The truncation happens at sigma. The values are given in days.
 *
 * @param med The median of the distribution (in days).
 * @param sigma The standard deviation of the distribution (in days).
 * @return The rotation frequency in microHz.
 */
long double rot_envelope(long double med=60., long double sigma=3.);

/**
 * @brief Computes the frequency splitting of modes assuming a two-zone rotation profile.
 *
 * This function calculates the frequency splitting of modes based on a two-zone rotation profile. The rotation rates in the envelope and core are given as input parameters. The function also takes the ksi function, which describes the degree of mixture between p and g modes as a function of frequency.
 *
 * @param ksi_pg The ksi function that describes the degree of mixture between p and g modes as a function of frequency.
 * @param rot_envelope The average rotation in the envelope. Must be a scalar.
 * @param rot_core The average rotation in the core. Must be a scalar.
 * @return A VectorXd object containing the frequency splitting values.
 */
VectorXd dnu_rot_2zones(const VectorXd& ksi_pg, const long double rot_envelope, const long double rot_core);

/**
 * @brief Computes the value of numax based on the Dnu_star parameter using the relation from Stello+2009.
 *
 * This function calculates the value of numax based on the Dnu_star parameter using the relation from Stello+2009. The relation is given by Dnu ~ 0.263 * numax^0.77. The function also allows for adding a uniform spread around numax if the spread parameter is provided.
 *
 * @param Dnu_star The Dnu_star parameter.
 * @param spread The spread around numax, given as a fraction (e.g., 5% is 0.05). Default value is 0.
 * @return The value of numax.
 */
long double numax_from_stello2009(const long double Dnu_star, const long double spread);

/**
 * @brief Generates a set of mode parameters for simulating an evolved star based on the asymptotic relations.
 *
 * This function generates a set of mode parameters for simulating an evolved star based on the asymptotic relations. The function takes a structure `cfg_star` as input, which contains all the important parameters for creating the mode parameters. The mode parameters include frequencies, widths, heights, and a-coefficients for each mode of degree l=0, l=1, l=2, and l=3.
 *
 * @param cfg_star A structure that contains all the important parameters for creating the mode parameters, such as (see cfg_star structure itself for the full list), 
 *                 - Dnu_star: Large separation for the p modes.
 *                 - epsilon_star: Phase offset term for the asymptotic relation of the p modes.
 *                 - delta0l_star, alpha_p_star, nmax_star: Parameters for creating 2nd order asymptotic p modes.
 *                 - DPl_star: Period spacing for l=1 g modes.
 *                 - alpha_star: Phase offset term for the asymptotic relation of the g modes.
 *                 - q_star: Coupling coefficient between p and g modes.
 *                 - fmin: Minimum frequency for the modes that should be included in the calculations.
 *                 - fmax: Maximum frequency for the modes that should be included in the calculations.
 * @return A structure `params_out` that contains the mode parameters for simulating the evolved star.
 *         - nu_lx: Frequencies of the l=x modes, where x is between 0 and 3.
 *         - nu_p_l1: Base p mode frequencies used to build the frequencies for the l=1 mixed modes.
 *         - nu_g_l1: Base g mode frequencies used to build the frequencies for the l=1 mixed modes.
 *         - width_lx: Widths of the l=x modes, where x is between 0 and 3.
 *         - height_lx: Heights of the l=x modes, where x is between 0 and 3.
 *         - aj_lx: The a-coefficients of order j for each mode of degree l=x.
 */
Params_synthetic_star make_synthetic_asymptotic_star(Cfg_synthetic_star cfg_star);

/**
 * @brief Displays the values of the parameters in the Cfg_synthetic_star structure.
 *
 * This function displays the values of the parameters in the Cfg_synthetic_star structure. The parameters include Teff_star, numax_star, Dnu_star, epsilon_star, and other parameters related to the synthetic star simulation.
 *
 * @param cfg The Cfg_synthetic_star structure containing the parameters.
 */
void displayCfgSyntheticStar(const Cfg_synthetic_star& cfg);