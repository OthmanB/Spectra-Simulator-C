/**
 * @file linspace.h
 * @brief Header file for the linspace function
 *
 * Header file for the function that generate equally spaced eigen vector of values
 * 
 *
 * @date 20 Apr 2016
 * @author obenomar
 */

#pragma once
#include <Eigen/Dense>
#include <iostream>

using Eigen::VectorXd;

/**
 * @brief Generate a linearly spaced vector.
 *
 * This function generates a linearly spaced vector between the specified start and end values, with the specified number of elements.
 *
 * @param start_in The starting value of the vector.
 * @param end_in The ending value of the vector.
 * @param num_in The number of elements in the vector.
 * @return A VectorXd object representing the linearly spaced vector.
 *         If num_in is 0, the function returns a vector with a single element of -1.
 */
VectorXd linspace(const long double start_in, const long double end_in, const long num_in);