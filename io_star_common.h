/*
 * io_star_common.cpp
 *
 * Contains all kind of methods
 * used to process and/or encapsulate data
 * 
 *  Created on: 21 Jun 2021
 *      Author: obenomar
 */
 # pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include "io_star_common.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

//void common_case0(std::ofstream& outfile, const MatrixXd& mode_params, std::string modelname);
void common_model_MS_Global_a1a2a3_HarveyLike(std::ofstream& outfile, const MatrixXd& mode_params);
void common_model_MS_Global_a1etaAlma3_HarveyLike(std::ofstream& outfile, const MatrixXd& mode_params);
void common_model_MS_Global_aj_HarveyLike(std::ofstream& outfile, const MatrixXd& mode_params);
void write_asym_key(std::ofstream& outfile, const bool fix, const double fix_val);
void write_a1_key(std::ofstream& outfile, const MatrixXd& mode_params, const bool fix, const double fix_val);
void write_a3_key(std::ofstream& outfile, const MatrixXd& mode_params, const bool fix, const double fix_val);
void write_aj_keys_2params_default(std::ofstream& outfile, const MatrixXd& mode_params, const bool fix, const double fix_val, const double j, const VectorXi& aj_ind);
void common_use_file_template(std::ofstream& outfile, const std::string template_file);