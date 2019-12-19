/*
 * noise_models.cpp
 *
 *  Created on: 24 Feb 2016
 *      Author: obenomar
 */
#pragma once
#include <math.h>
#include <Eigen/Dense>
#include <string>

using Eigen::VectorXd;
using Eigen::MatrixXd;

VectorXd harvey_like(const VectorXd noise_params, VectorXd x, VectorXd y, const int Nharvey);
VectorXd harvey1985(const VectorXd noise_params, VectorXd x, VectorXd y, const int Nharvey);
VectorXd harvey_like(const MatrixXd noise_params, VectorXd x, VectorXd y);