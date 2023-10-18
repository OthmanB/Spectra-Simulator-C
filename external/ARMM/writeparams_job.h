/**
 * @file writeparams_job.h
 * @brief Contains functions for writing and manipulating star mode parameters.
 */
#pragma once 
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
 
/**
 * @brief Copies the contents of one file to another.
 *
 * This function copies the contents of the input file to the output file. It removes trailing white space and tabulation from each line and excludes lines starting with "#". The copied contents are written to the output file.
 *
 * @param filename The path of the input file.
 * @param outputFilename The path of the output file.
 * @param append A flag indicating whether to append to an existing file or create a new file.
 */
void copy_cfg(const std::string& filename, const std::string& outputFilename, bool append = false);

/**
 * @brief Converts synthetic star parameters to a matrix of mode parameters with a specific inclination.
 * 
 * This function takes synthetic star parameters and converts them into a matrix of mode parameters. The mode parameters are organized based on their degree and include additional parameters such as asymmetry and inclination.
 * 
 * @param params The synthetic star parameters.
 * @param inc The inclination angle.
 * @return The matrix of mode parameters.
 */
MatrixXd bumpoutputs_2_MatrixXd_with_aj(Params_synthetic_star params, double inc);

/**
 * @brief Writes range modes to a file.
 *
 * This function takes the configuration and parameters of a synthetic star and writes the range modes to a file. The range modes include the numax, minimum and maximum frequencies, and the position of the curvature for the 2nd order equation of frequencies.
 *
 * @param cfg_star The configuration of the synthetic star.
 * @param params The parameters of the synthetic star.
 * @param output_file The output file path.
 * @param append A flag indicating whether to append to an existing file or create a new file.
 */
void write_range_modes(Cfg_synthetic_star cfg_star, Params_synthetic_star params, std::string output_file, bool append=false);

 /**
 * @brief Writes star mode parameters to a file in a specific format.
 * 
 * This function takes a matrix of mode parameters and writes them to a file. The mode parameters are written in a specific format, including headers and labels for each parameter.
 * 
 * @param mode_params The matrix of mode parameters to be written.
 * @param file_out The output file path.
 * @param append A flag indicating whether to append to an existing file or create a new file.
 */
void write_star_mode_params_asympt_model(MatrixXd mode_params, std::string file_out, bool append=false);

/**
 * @brief Writes star noise parameters to a file.
 *
 * This function takes a matrix of noise parameters and writes them to a file. The noise parameters are organized based on their degree and include additional information such as H0, tau_0, p0, H1, tau_1, p1, and N0.
 *
 * @param noise_params The matrix of noise parameters.
 * @param file_out The output file path.
 * @param append A flag indicating whether to append to an existing file or create a new file.
 */
void write_star_noise_params(MatrixXd noise_params, std::string file_out, bool append=false);

/**
 * @brief Writes l=1 p and g modes to a file.
 *
 * This function takes the parameters of a synthetic star and writes the l=1 p and g modes to a file. The p and g modes are used for the mixed mode computation.
 *
 * @param params_out The parameters of the synthetic star.
 * @param file_out The output file path.
 * @param append A flag indicating whether to append to an existing file or create a new file.
 */
void write_star_l1_roots(Params_synthetic_star params_out, std::string file_out, bool append=false);

/**
 * @brief Displays an error message when a file cannot be read.
 *
 * This function is called when there is an error opening a file for reading. It displays an error message and exits the program.
 *
 * @param file_out The file that could not be read.
 */
void file_read_error(std::string file_out);