/**
 * @file stellar_models.h
 * @brief High level functions that handle read/write of files from stellar models (ADIPLS, MESA)
 *
 * This file contains the declarations of high level functions that handle the read/write operations
 * of files from stellar models, specifically for the ADIPLS and MESA models.
 *
 * @date 10 Oct 2017
 * @author obenomar
 */

#pragma once

#include <math.h>
#include <Eigen/Dense>
# include <iostream>
# include <iomanip>
# include <vector>
# include "data.h"
# include "ioproc.h" // contains the string handlers
# include "format.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;

/**
 * @brief Function that extracts the relevant information from the table of models.
 *
 * This function extracts the requested row from the input table of models. The data should be arranged by rows, with the first column containing the model index. The function also requires the frequency files to be in a specific format: freq.[Model_index]. For example, freq.000101.24.
 *
 * @param table The input table that contains information from the models.
 * @param table_index The column index that is requested to extract.
 * @param dir_freqs The directory that contains all of the frequency files for each model.
 * @return A structure of type Model_data containing the extracted model parameters and frequencies.
 *
 * @note The function assumes that the table is in the format specified in the input conditions.
 * @note The function uses the Model_data structure defined in data.h.
 */
Model_data get_model_param(const Data_Nd table, const long table_index, const std::string dir_freqs);

/**
 * @brief Function to read the frequency file.
 *
 * This function reads the frequency file and returns the data in a Data_Nd structure.
 *
 * @param file_in_name The name of the frequency file.
 * @param verbose_data Flag to control the verbosity of the data reading process.
 * @return A Data_Nd structure containing the data from the frequency file.
 */
Data_Nd read_freq(const std::string file_in_name, const bool verbose_data);
