
/**
 * @file configure_make_star.h
 * @brief Configures the parameters for creating a synthetic star based on the input parameters.
 *
 * This file contains the implementation of the `configure_make_star` function, which takes a map of input parameters and configures the parameters for creating a synthetic star. The input parameters include Dnu_star, DPl_star, q_star, alpha_g_star, epsilon_star, delta0l_percent_star, rot_env, rot_core, rot_ratio, a3_l2, a3_l3, a5_l3, a2_l1, a2_l2, a2_l3, a4_l2, a4_l3, a6_l3, max_HNR, H0_spread, Gamma_max_l0, Hfactor, Wfactor, numax_star, numax_spread, fmin_in_Dnu, fmax_in_Dnu, output_file_rot, file_template, Vl, nmax_star, nmax_spread, beta_p, alpha_p, params_harvey_like, and inclination.
 */
#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "data_solver.h"

/**
 * @brief Main function of the "make_star" binary file. Configures the parameters for creating a synthetic star based on the input parameters.
 *
 * This function takes a map of input parameters and configures the parameters for creating a synthetic star. The input parameters include Dnu_star, DPl_star, q_star, alpha_g_star, epsilon_star, delta0l_percent_star, rot_env, rot_core, rot_ratio, a3_l2, a3_l3, a5_l3, a2_l1, a2_l2, a2_l3, a4_l2, a4_l3, a6_l3, max_HNR, H0_spread, Gamma_max_l0, Hfactor, Wfactor, numax_star, numax_spread, fmin_in_Dnu, fmax_in_Dnu, output_file_rot, file_template, Vl, nmax_star, nmax_spread, beta_p, alpha_p, params_harvey_like, and inclination.
 *
 * @param input_params A map of input parameters for configuring the synthetic star.
 * @return A Cfg_synthetic_star structure containing the configured parameters for creating the synthetic star.
 */
Cfg_synthetic_star configure_make_star(std::unordered_map<std::string, std::string> input_params);
