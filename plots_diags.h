/**
 * @file plots_diags.h
 * @brief Plot the model using Gnuplot.
 *
 * Handling plots using gnuplot
 *
 * @date 06 Oct 2016
 * @author obenomar
 */

#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include "gnuplot-iostream.h"
#include "data.h"
#include "ioproc.h"


using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * @brief Plot the model using Gnuplot.
 *
 * @param x The x values.
 * @param y The y values.
 * @param model The model values.
 * @param scoef1 The first smoothing coefficient.
 * @param scoef2 The second smoothing coefficient.
 * @param file_model The output file name for the plot.
 *
 * This function handles the plotting of the model using Gnuplot. It takes the x values, y values, model values, and smoothing coefficients as input. The output plot is saved in the specified file_model.
 *
 * @note This function uses the Gnuplot library to generate the plot.
 */
void gnuplt_model(VectorXd x, VectorXd y, VectorXd model, double scoef1, double scoef2, std::string file_model);


