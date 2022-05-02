/*
 * noise_models.cpp
 *
 *  Created on: 24 Feb 2016
 *      Author: obenomar
 */
#pragma once
#include <math.h>
#include <Eigen/Dense>

using Eigen::VectorXd;

VectorXd harvey_like(const VectorXd& noise_params, const VectorXd& x, const VectorXd& y, const int Nharvey);
VectorXd harvey1985(const VectorXd& noise_params, const VectorXd& x, const VectorXd& y, const int Nharvey);
