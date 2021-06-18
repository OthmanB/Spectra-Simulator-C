/*
 * stellar_models.h
 *
 *  High level functions that handle read/write
 *  of files from stellar models (ADIPLS, MESA)
 * 
 *  Created on: 10 Oct 2017
 *      Author: obenomar
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

Model_data get_model_param(const Data_Nd table, const long table_index, const std::string dir_freqs);
Data_Nd read_freq(const std::string file_in_name, const bool verbose_data);
