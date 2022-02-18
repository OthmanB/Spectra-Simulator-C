/*
 * models_database.h
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
# include <string>
#include "io_star_params.h"

//void generate_cfg_asymptotic_Hgauss(VectorXd input_params, std::string file_out_modes, std::string file_out_noise);
void generate_cfg_asymptotic_act_asym_Hgauss(VectorXd input_params, std::string file_out_modes, std::string file_out_noise);
void generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string extra);
void generate_cfg_from_synthese_file_Wscaled_a1a2a3asymovGamma(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string extra);
void generate_cfg_from_synthese_file_Wscaled_Alm(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string extra);
void generate_cfg_from_synthese_file_Wscaled_aj(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string extra);
void asymptotic_mm_v1(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path, std::string template_file);
void asymptotic_mm_v2(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path, std::string template_file);
void asymptotic_mm_v3(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path, std::string template_file);
void asymptotic_mm_freeDp_numaxspread_curvepmodes_v1(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path, std::string template_file);
void asymptotic_mm_freeDp_numaxspread_curvepmodes_v2(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path, std::string template_file);
void asymptotic_mm_freeDp_numaxspread_curvepmodes_v3(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path, std::string template_file);
double eta0_fct(const VectorXd& fl0_all); // To compute centrifugal term
