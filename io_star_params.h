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
# include <string>
# include "ioproc.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;


std::string read_allfile(const std::string file); // Read a full file without caring about the content
std::vector<std::string> read_allfile_vect(const std::string file);  // Read a full file without caring about the content But separate each line in a vector
void write_global_info(VectorXd spec_params, std::string file_out, std::string identifier, bool append=false); // Write Identifier, Cadence and observation duration
void write_star_params_a1a2a3asym(VectorXd spec_params, MatrixXd mode_params, MatrixXd noise_params, std::string file_out, std::string identifier);
void write_star_params_act_asym(VectorXd spec_params, MatrixXd mode_params, MatrixXd noise_params, std::string file_out, std::string identifier);
void write_star_mode_params_a1a2a3(MatrixXd mode_params, std::string file_out, bool append=false);
void write_star_mode_params_act_asym(MatrixXd mode_params, std::string file_out, bool append=false);
void write_star_params_Alm(VectorXd spec_params, MatrixXd mode_params, MatrixXd noise_params, std::string file_out, std::string identifier);
void write_star_params_aj(VectorXd spec_params, MatrixXd mode_params, MatrixXd noise_params, std::string file_out, std::string identifier);
void write_star_mode_params_Alm(MatrixXd mode_params, std::string file_out, bool append=false);
void write_star_mode_params_aj(MatrixXd mode_params, std::string file_out, bool append=false);
void write_star_noise_params(MatrixXd noise_params, std::string file_out, bool append=false);
void write_range_modes(Cfg_synthetic_star cfg_star, Params_synthetic_star params, std::string output_file);
void write_spectrum(VectorXd x, VectorXd y, VectorXd z, std::string file_out, bool write_inmodel); // Alias of the write_spectrum() below, but with fmin=-1, fmax=-1 (de facto optional parameterss)
void write_spectrum(const VectorXd x, const VectorXd y, const VectorXd z, const std::string file_out, const bool write_inmodel, const double fmin, const double fmax); 
void write_spectrum(const VectorXd x, const VectorXd y, const VectorXd z, const double scoef1, double scoef2, const std::string file_out);
VectorXd write_star_model(const MatrixXd mode_params, const MatrixXd noise_params, const std::string file_out, 
                          const std::string identifier, const std::string modelname, const std::string common_template_dir);
MatrixXd bumpoutputs_2_MatrixXd(Params_synthetic_star params, double inc);
Config_Data read_main_cfg(std::string cfg_file);			
Data_Nd read_data_ascii_Ncols(const std::string file_in_name, const std::string delimiter, const bool verbose_data);
Star_params read_star_params(const std::string file_in_name);
VectorXd smooth(VectorXd in, double scoef);
void file_read_error(std::string);

void convert_in_to_model(const std::string file_out, const std::string identifier);