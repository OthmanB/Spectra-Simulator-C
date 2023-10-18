/*
 * models_database_grid.h
 *
 * Header file that contains all kind of methods
 * used to generate models for the pulsation/noise
 * 
 *  Created on: 20 Apr 2016
 *      Author: obenomar
 */

/**
 * @file models_database_grid.h
 * @brief Header file that contains models (for grid only)
 *
 * Header file that contains all kind of methods used to generate models for the pulsation/noise when in grid mode.
 * 
 *
 * @date 20 Apr 2016
 * @author obenomar
 */

#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/Dense>
#include "io_star_params.h"
#include "ioproc.h"

/**
 * @brief Generate a configuration file from a reference star with scaled heights and widths.
 *
 * This function generates a configuration file for a simulated star based on a reference star. It rescales the heights and widths of the modes in the reference star based on the input parameters. It also calculates the noise level for the simulated star using the Harvey profile noise model. The resulting configuration file contains the mode parameters and the noise parameters.
 *
 * @param input_params A VectorXd object containing the input parameters for the simulation.
 * @param input_model The Model_data object containing the data for the reference star.
 * @param file_ref_star The file path of the reference star data.
 * @param file_out_modes The file path to write the mode parameters of the simulated star.
 * @param file_out_noise The file path to write the noise parameters of the simulated star.
 */
void generate_cfg_from_refstar_HWscaled(VectorXd input_params, Model_data input_model, std::string file_ref_star, std::string file_out_modes, std::string file_out_noise);

/**
 * @brief Generate a configuration file from a reference star with scaled heights and widths and scaled noise parameters.
 *
 * This function generates a configuration file for a simulated star based on a reference star. It rescales the heights and widths of the modes in the reference star based on the input parameters. It also scales the noise parameters of the simulated star using the input noise scaling coefficients. The resulting configuration file contains the mode parameters and the noise parameters.
 *
 * @param input_params A VectorXd object containing the input parameters for the simulation.
 * @param input_model The Model_data object containing the data for the reference star.
 * @param file_ref_star The file path of the reference star data.
 * @param file_out_modes The file path to write the mode parameters of the simulated star.
 * @param file_out_noise The file path to write the noise parameters of the simulated star.
 */
void generate_cfg_from_refstar_HWscaled_GRANscaled(VectorXd input_params, Model_data input_model, std::string file_ref_star, std::string file_out_modes, std::string file_out_noise);

