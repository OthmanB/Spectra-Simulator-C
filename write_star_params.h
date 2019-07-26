/*
 * write_star_params.h
 *
 * Header file that contains all kind of methods
 * used to process and/or encapsulate data
 * 
 *  Created on: 22 Feb 2016
 *      Author: obenomar
 */

# pragma once
# include <string>
# include <vector>
# include <Eigen/Dense>
# include "data.h" // contains the structure Data
# include <fstream>

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

bool str_to_bool(const std::string str);
long str_to_lng(const std::string str);
double str_to_dbl(const std::string str);
std::vector<double> str_to_arrdbl(const std::string str, const std::string delimiters);
std::vector<double> arrstr_to_arrdbl(const std::vector<std::string> vals_str);
bool file_exists(const std::string& name);
void write_star_params_act_asym(VectorXd spec_params, MatrixXd mode_params, MatrixXd noise_params, std::string file_out, std::string identifier);
void write_star_mode_params_a1a2a3(MatrixXd mode_params, std::string file_out);
void write_star_mode_params_act_asym(MatrixXd mode_params, std::string file_out);
void write_star_noise_params(MatrixXd noise_params, std::string file_out);
void write_spectrum(VectorXd x, VectorXd y, VectorXd z, std::string file_out, bool write_inmodel);
void write_spectrum_v2(const VectorXd x, const VectorXd y, const VectorXd z, const double scoef1, double scoef2, const std::string file_out);
Config_Data read_main_cfg(std::string cfg_file);			
std::string strtrim(const std::string& str);
std::vector<std::string> strsplit(const std::string str, const std::string delimiters);
Data_Nd read_data_ascii_Ncols(const std::string file_in_name, const std::string delimiter, const bool verbose_data);
VectorXd smooth(VectorXd in, double scoef);
