/*
 * format.h
 *
 *  Functions that format long or str into
 *  a suitable format for read/write
 * 
 *  Created on: 10 Oct 2017
 *      Author: obenomar
 */

#pragma once

//#include <math.h>
//#include <Eigen/Dense>
# include <iostream>
# include <iomanip>
//# include <vector>
# include "ioproc.h" // contains the string handlers

using Eigen::VectorXd;
using Eigen::MatrixXd;

std::string identifier2chain(long identifier);
std::string identifier2chain(std::string identifier);
std::string format_freqname(std::string id);
