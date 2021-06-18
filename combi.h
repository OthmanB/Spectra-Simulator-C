/*
 *  combi.h
 *
 *  Function that handle combination generation
 *  and read/write into files
 * 
 *  Created on: 10 Oct 2017
 *      Author: obenomar
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

MatrixXd define_all_combinations(const MatrixXd Values, const VectorXi Nvalues, const int Nparams);
long read_id_allcombi(std::string file_combi);
std::string write_allcombi(MatrixXd allcombi, VectorXd cte_params, Config_Data cfg, std::string fileout, bool erase_old_file, long iter, long id0, 
		    std::vector<std::string> cte_names, std::vector<std::string> var_names, std::vector<std::string> param_names);

