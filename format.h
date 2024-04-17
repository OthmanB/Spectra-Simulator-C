/**
 * @file format.h
 * @brief Functions that format long or str into  a suitable format for read/write
 * 
 * 
 * @date 10 Oct 2017
 * @author obenomar
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

/**
 * @brief Convert an identifier to a chain of characters
 *
 * This function converts a given identifier to a chain of characters with leading zeros.
 *
 * @param identifier The identifier to convert
 * @return The chain of characters with leading zeros
 */
std::string identifier2chain(long identifier);

/**
 * @brief Convert a string identifier to a chain of characters
 *
 * This function converts a given string identifier to a chain of characters with leading zeros.
 *
 * @param identifier The string identifier to convert
 * @return The chain of characters with leading zeros
 */
std::string identifier2chain(std::string identifier);

/**
 * @brief Format the name of the frequency file produced by ADIPLS
 * 
 * This function formats the name of the frequency file produced by ADIPLS.
 * It adds leading zeros to the identifier to get a fixed length output string.
 * 
 * @param id The identifier to format
 * @return The formatted string with leading zeros
 */
std::string format_freqname(std::string id);
