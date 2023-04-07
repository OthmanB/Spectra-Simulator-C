
#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include "../../linspace.h"
#include "../../linfit.h"
#include "data.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;

Data_asympt_p decompose_nu_nl(const int l, const VectorXd fl0, const VectorXd fl, const double Cfactor, const bool verbose);
//VectorXd where_in_range(const VectorXd& vec, const double value_min, const double value_max, const bool strict);