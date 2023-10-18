/**
 * @file combi.h
 * @brief Function that handles combination generation and read/write into files
 * 
 * This file contains functions for generating combinations and reading/writing them into files.
 * 
 * @date 10 Oct 2017
 * @author obenomar
 */ 

#pragma once

# include <iostream>
# include <iomanip>
#include <fstream>
# include <Eigen/Dense>
# include <vector>
#include <string>
# include "ioproc.h"

using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;


/**
 * @brief Function that generates all possible combinations recursively
 *
 * This function generates all possible combinations of values using a recursive approach.
 *
 * @param Values The matrix of values for each parameter
 * @param Nvalues The vector containing the number of values for each parameter
 * @param Nparams The number of parameters
 * @param allcombi The matrix to store all the combinations
 * @param z The current row index in the allcombi matrix
 * @param current_combi The vector to store the current combination
 * @param current_param The current parameter index
 */
void generate_combinations(const MatrixXd& Values, const VectorXi& Nvalues, const int Nparams, MatrixXd& allcombi, long& z, VectorXd& current_combi, int current_param);

/**
 * @brief Function that generates all possible combinations
 *
 * This function generates all possible combinations of values using the generate_combinations function.
 *
 * @param Values The matrix of values for each parameter
 * @param Nvalues The vector containing the number of values for each parameter
 * @param Nparams The number of parameters
 *
 * @return The matrix containing all the combinations
 */
MatrixXd define_all_combinations(const MatrixXd& Values, const VectorXi& Nvalues, const int Nparams);

/**
 * @brief Read the ID from the last line of a file containing combinations
 *
 * This function reads the last line of a file containing combinations and extracts the ID from it.
 *
 * @param file_combi The path to the file containing combinations
 * @return The ID extracted from the last line of the file
 */
long read_id_allcombi(std::string file_combi);

/**
 * @brief Write all combinations to a file
 *
 * This function writes all combinations to a file, along with other information such as model name and parameter names.
 *
 * @param allcombi The matrix containing all the combinations
 * @param cte_params The vector of constant parameters
 * @param cfg The configuration data
 * @param fileout The path to the output file
 * @param erase_old_file Flag indicating whether to erase the old file or append to it
 * @param iter The current iteration
 * @param id0 The initial ID
 * @param cte_names The names of the constant parameters
 * @param var_names The names of the variable parameters
 * @param param_names The names of all parameters
 * @return The ID extracted from the last line of the file
 */
std::string write_allcombi(MatrixXd& allcombi, VectorXd& cte_params, Config_Data cfg, std::string fileout, bool erase_old_file, long iter, long id0, 
		    std::vector<std::string> cte_names, std::vector<std::string> var_names, std::vector<std::string> param_names);

