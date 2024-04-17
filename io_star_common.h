/**
 * @file io_star_common.h
 * @brief Methods to encapsulate data.
 * 
 * Contains all kind of methods used to process and/or encapsulate data.
 * 
 * @date 21 Jun 2021
 * @author obenomar
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

/**
 * @brief Write common model parameters for a specific model type
 *
 * This function writes the common model parameters for a specific model type to the specified output file stream. The model name and other parameters are determined based on the given mode_params.
 *
 * @param outfile The output file stream to write the model parameters to
 * @param mode_params The mode parameters used to determine the model name and other parameters
 */
void common_model_MS_Global_a1a2a3_HarveyLike(std::ofstream& outfile, const MatrixXd& mode_params);

/**
 * @brief Write common model parameters for a specific model type
 *
 * This function writes the common model parameters for a specific model type to the specified output file stream. The model name and other parameters are determined based on the given mode_params.
 *
 * @param outfile The output file stream to write the model parameters to
 * @param mode_params The mode parameters used to determine the model name and other parameters
 */
void common_model_MS_Global_a1etaAlma3_HarveyLike(std::ofstream& outfile, const MatrixXd& mode_params);

/**
 * @brief Write common model parameters for a specific model type
 *
 * This function writes the common model parameters for a specific model type to the specified output file stream. The model name and other parameters are determined based on the given mode_params.
 *
 * @param outfile The output file stream to write the model parameters to
 * @param mode_params The mode parameters used to determine the model name and other parameters
 */
void common_model_MS_Global_aj_HarveyLike(std::ofstream& outfile, const MatrixXd& mode_params);

/**
 * @brief Write the asymmetry key to the output file stream
 *
 * This function writes the asymmetry key to the specified output file stream. The asymmetry value is determined based on the given fix and fix_val parameters.
 *
 * @param outfile The output file stream to write the asymmetry key to
 * @param fix Boolean value indicating whether the asymmetry is fixed or not
 * @param fix_val The fixed value of the asymmetry (if fix is true)
 */
void write_asym_key(std::ofstream& outfile, const bool fix, const double fix_val);

/**
 * @brief Write the a1 key to the output file stream
 *
 * This function writes the a1 key to the specified output file stream. The value of the a1 key is determined based on the given mode_params, fix, and fix_val parameters.
 *
 * @param outfile The output file stream to write the a1 key to
 * @param mode_params The mode parameters used to determine the value of the a1 key
 * @param fix Boolean value indicating whether the a1 key is fixed or not
 * @param fix_val The fixed value of the a1 key (if fix is true)
 */
void write_a1_key(std::ofstream& outfile, const MatrixXd& mode_params, const bool fix, const double fix_val);

/**
 * @brief Write the a3 key to the output file stream
 *
 * This function writes the a3 key to the specified output file stream. The value of the a3 key is determined based on the given mode_params, fix, and fix_val parameters.
 *
 * @param outfile The output file stream to write the a3 key to
 * @param mode_params The mode parameters used to determine the value of the a3 key
 * @param fix Boolean value indicating whether the a3 key is fixed or not
 * @param fix_val The fixed value of the a3 key (if fix is true)
 */
void write_a3_key(std::ofstream& outfile, const MatrixXd& mode_params, const bool fix, const double fix_val);


/**
 * @brief Write the aj keys with 2 parameters and default values to the output file stream
 *
 * This function writes the aj keys with 2 parameters and default values to the specified output file stream. The values of the aj keys are determined based on the given fix, fix_val, j, and aj_ind parameters.
 *
 * @param outfile The output file stream to write the aj keys to
 * @param mode_params The mode parameters used to determine the default values of the aj keys
 * @param fix Boolean value indicating whether the aj key is fixed or not
 * @param fix_val The fixed value of the aj key (if fix is true)
 * @param j The index of the aj key
 * @param aj_ind The indices of the aj keys in the mode_params matrix
 */
void write_aj_keys_2params_default(std::ofstream& outfile, const MatrixXd& mode_params, const bool fix, const double fix_val, const double j, const VectorXi& aj_ind);

/**
 * @brief Read a small template file and write its content to an output file stream
 *
 * This function reads a small template file that contains all the common parameters that need to be written at the end of the .model file. The content of the template file is written to the specified output file stream.
 *
 * @param outfile The output file stream to write the template content to
 * @param template_file The path to the template file
 */
void common_use_file_template(std::ofstream& outfile, const std::string template_file);