/**
 * @file decompose_nu.h
 * @brief Declarations of functions of decompose_nu_nl.
 *
 * This file contains the declarations of the function to decompose nu(n,l) into the asymptotic terms. 
 *
 * @date 24 Feb 2023
 * @author obenomar
 */

#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include "../../linspace.h"
#include "../../linfit.h"
#include "data.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;

/**
 * @brief Decomposes the frequency by identifying the asymptotic elements.
 *
 * This function decomposes the frequency by identifying the asymptotic elements, including nl and O2_l.
 *
 * @param l The angular degree of the frequency.
 * @param fl0 The frequencies for l=0.
 * @param fl The frequencies for the given l.
 * @param Cfactor The correction factor.
 * @param verbose Flag indicating if verbose output is enabled.
 * @return The decomposed asymptotic elements.
 */
Data_asympt_p decompose_nu_nl(const int l, const VectorXd fl0, const VectorXd fl, const double Cfactor, const bool verbose);
