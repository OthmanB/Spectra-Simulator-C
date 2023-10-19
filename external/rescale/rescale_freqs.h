/**
 * @file rescale_freqs.h
 * @brief Declarations of functions related to rescaling frequencies.
 *
 * This file contains the declarations of functions related to rescaling frequencies.
 *
 * @date 24 Feb 2023
 * @author obenomar
 */

#pragma once
#include <Eigen/Dense>
#include <cmath>
#include "../../linspace.h"
#include "../../linfit.h"
#include "decompose_nu.h"
#include "data.h"

using Eigen::VectorXd;


/**
 * @brief Rescales the frequencies using the given parameters.
 *
 * This function rescales the frequencies by identifying all of the terms of the asymptotic
 * in order to isolate the residual error term. Then it uses this residual error (with a proper rescaling) to generate a new set 
 * of frequencies following the asymptotic with the desired parameters.
 *
 * @param Dnu_star The desired Dnu value for rescaling.
 * @param epsilon_star The desired epsilon value for rescaling.
 * @param freqs_ref The reference frequencies to be rescaled.
 * @param d0l_star The desired dl values for rescaling.
 * @return The rescaled frequencies.
 */
Freq_modes rescale_freqs(const double Dnu_star, const double epsilon_star, const Freq_modes freqs_ref, const VectorXd& d0l_star);