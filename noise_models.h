/**
 * @file noise_models.h
 * @brief Declarations of functions related to noise models.
 *
 * This file contains the declarations of functions related to noise models.
 *
 * @date 24 Feb 2016
 * @author obenomar
 */
#pragma once
#include <math.h>
#include <Eigen/Dense>

using Eigen::VectorXd;
using Eigen::MatrixXd;


/**
 * @brief Calculate a sum of Harvey-like profile and white noise and add them to an initial input vector.
 *
 * @param noise_params The vector of noise parameters. The parameters are assumed to be in the following order: [H0, tc0, p0, ..., Hn, tcn, pn, N0].
 * @param x The input vector x.
 * @param y The modified vector y.
 * @param Nharvey The number of Harvey profiles to consider.
 * @return The modified vector y.
 *
 * This function calculates a sum of a Harvey-like profile and white noise and adds them to an initial input vector. The function assumes that the inputs are in the following order: [H0, tc0, p0, ..., Hn, tcn, pn, N0].
 *
 * @note The function modifies the vector y and returns it.
 */
VectorXd harvey_like(const VectorXd& noise_params, VectorXd& x, VectorXd& y, const int Nharvey);



/**
 * @brief Calculate a sum of Harvey-like profile and white noise and add them to an initial input vector.
 *
 * @param noise_params The matrix of noise parameters.
 * @param x The input vector x.
 * @param y The initial input vector y.
 * @return The modified vector y.
 *
 * This function calculates a sum of a Harvey-like profile and white noise and adds them to an initial input vector y. The function assumes that the inputs are in the following order: [H0, tc0, p0, ..., Hn, tcn, pn, N0].
 *
 * @note The last row of noise_params is assumed to always contain the white noise.
 */
VectorXd harvey_like(const MatrixXd& noise_params, VectorXd x, VectorXd& y);

/**
 * @brief Calculate a sum of Harvey-like profile and white noise and add them to an initial input vector.
 *
 * @param noise_params The matrix of noise parameters.
 * @param x The input vector x.
 * @return The modified vector spec_noise.
 *
 * This function calculates a sum of a Harvey-like profile and white noise and adds them to an initial input vector. The function assumes that the inputs are in the following order: [H0, tc0, p0, ..., Hn, tcn, pn, N0].
 *
 * @note The function returns the modified vector spec_noise.
 */
VectorXd harvey_like(const MatrixXd& noise_params, const VectorXd& x);

/**
 * @brief Calculate a sum of Harvey profile and white noise and add them to an initial input vector.
 * 
 * @param noise_params The matrix of noise parameters.
 * @param x The input vector x.
 * @return The modified vector spec_noise.
 * 
 * This function calculates a sum of Harvey profile and white noise and adds them to an initial input vector. The function assumes that the inputs are in the following order: [H0, tc0, p0, ..., Hn, tcn, pn, N0].
 * 
 * @note The function returns the modified vector spec_noise.
 */
VectorXd harvey_1985(const MatrixXd& noise_params, const VectorXd& x);

/**
 * @brief Calculates the sum of a Harvey profile and white noise.
 *
 * This function calculates the sum of a Harvey profile and white noise, and adds these profiles to an initial input vector y. The Harvey profile differs from the Harvey-like profile by the fact that H0_1985 = H0_like * tc is correlated to the timescale tc. There is also a 2pi factor in the denominator. The function assumes that the inputs are in the following order: [H0, tc0, p0, ..., Hn, tcn, pn, N0].
 *
 * @param noise_params The parameters for the noise model.
 * @param x The x-coordinates of the data points.
 * @param y The y-coordinates of the data points.
 * @param Nharvey The number of Harvey profiles to include in the sum.
 * @return VectorXd The resulting vector after adding the Harvey profile and white noise.
 */
VectorXd harvey1985(const VectorXd& noise_params, const VectorXd& x, const VectorXd& y, const int Nharvey);