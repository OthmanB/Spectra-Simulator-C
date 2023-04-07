#pragma once
#include <Eigen/Dense>
#include <cmath>
#include "../../linspace.h"
#include "../../linfit.h"
#include "decompose_nu.h"
#include "data.h"

using Eigen::VectorXd;


Freq_modes rescale_freqs(const double Dnu_star, const double epsilon_star, const Freq_modes freqs_ref, const VectorXd& d0l_star);