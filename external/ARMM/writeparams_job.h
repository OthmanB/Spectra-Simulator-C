
#pragma once 
// Functions adapted or taken from io_star_params.cpp in the Spectrum Simulator
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include "data_solver.h"
#include "string_handler.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;
 
void copy_cfg(const std::string& filename, const std::string& outputFilename, bool append = false);
void write_range_modes(Cfg_synthetic_star cfg_star, Params_synthetic_star params, std::string output_file, bool append=false);
void write_star_mode_params_asympt_model(MatrixXd mode_params, std::string file_out, bool append=false);
void write_star_noise_params(MatrixXd noise_params, std::string file_out, bool append=false);
void write_star_l1_roots(Params_synthetic_star params_out, std::string file_out, bool append=false);
MatrixXd bumpoutputs_2_MatrixXd_with_aj(Params_synthetic_star params, double inc);
void file_read_error(std::string file_out);