/** 
 * @file linfit.h
 * @brief A simple linear fit algorithm.
 *
 * This file contains the implementation of a simple linear fit algorithm. It provides a function to perform a linear fit on a set of data points.
 * @date 11 Oct 2017
 * @author obenomar
 */

#pragma once
#include <Eigen/Dense>

using Eigen::VectorXd;

/**
 * @brief Performs a linear fit on a set of data points.
 *
 * This function takes two vectors, x and y, representing the x-coordinates and y-coordinates of the data points, and performs a linear fit to find the best-fit line.
 *
 * @param x The x-coordinates of the data points.
 * @param y The y-coordinates of the data points.
 * @return VectorXd A vector containing the slope and intercept of the best-fit line.
 * 
 * @note The size of x and y must be the same.
 */
VectorXd linfit(const VectorXd& x, const VectorXd& y); // A more optimal way to perform a linear fit