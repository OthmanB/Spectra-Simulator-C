/** 
 * @file linfit.cpp
 * @brief A simple linear fit algorithm.
 *
 * This file contains the implementation of a simple linear fit algorithm. It provides a function to perform a linear fit on a set of data points.
 * @date 11 Oct 2017
 * @author obenomar
 */

# include <Eigen/Dense>
# include <iostream>
# include <iomanip>
# include "linfit.h"

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
VectorXd linfit(const VectorXd& x, const VectorXd& y) {
    if (x.size() != y.size()) {
        std::cout << "x and y do not have the same size!" << std::endl;
        std::cout << "Cannot proceed. The program will exit now." << std::endl;
        exit(EXIT_SUCCESS);
    }
    double sx = x.sum();
    double sy = y.sum();
    double n = x.size();
    double mean_x = sx / n;
    
    VectorXd t = x.array() - mean_x;
    VectorXd tmp = t.array() * y.array();
    double out0 = tmp.sum() / (t.array() * t.array()).sum();
    double out1 = (sy - sx * out0) / n;
    VectorXd out(2);
    out << out0, out1;
    return out;
}



