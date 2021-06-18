/*
 * models_database_grid.h
 *
 * Header file that contains all kind of methods
 * used to generate models for the pulsation/noise
 * 
 *  Created on: 20 Apr 2016
 *      Author: obenomar
 */

#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/Dense>
#include "io_star_params.h"
#include "ioproc.h"

//void generate_cfg_asymptotic_Hgauss(VectorXd input_params, std::string file_out_modes, std::string file_out_noise);
//void generate_cfg_asymptotic_act_asym_Hgauss(VectorXd input_params, std::string file_out_modes, std::string file_out_noise);
void generate_cfg_from_refstar_HWscaled(VectorXd input_params, Model_data input_model, std::string file_ref_star, std::string file_out_modes, std::string file_out_noise);
void generate_cfg_from_refstar_HWscaled_GRANscaled(VectorXd input_params, Model_data input_model, std::string file_ref_star, std::string file_out_modes, std::string file_out_noise);

