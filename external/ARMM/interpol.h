/**
 * @file interpol.h
 * @brief A collection of functions for interpolation and resampling of arrays.
 */

#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <iomanip>

using Eigen::VectorXd;

/**
 * @brief Linearly interpolates a curve (x, y) at the position x_int.
 * 
 * @param x The x-coordinates of the curve.
 * @param y The y-coordinates of the curve.
 * @param x_int The x-coordinate at which to interpolate.
 * @return double The interpolated y-coordinate.
 */
double lin_interpol(const VectorXd& x, const VectorXd& y, const double x_int);

/**
 * @brief Quadratically interpolates the array a and resamples it to size m.
 * 
 * @param a The array to be interpolated.
 * @param m The desired size of the resampled array.
 * @return VectorXd The resampled array.
 */
VectorXd quad_interpol(const VectorXd& a, const int m);

/**
 * @brief Linearly interpolates a value x in an array a.
 * 
 * @param x The value to be interpolated.
 * @param a The array to interpolate.
 * @return long double The interpolated value.
 */
long double interp1(const long double x, const VectorXd& a);

/**
 * @brief Linearly interpolates an array a and resamples it to size m.
 * 
 * @param a The array to be interpolated.
 * @param m The desired size of the resampled array.
 * @return VectorXd The resampled array.
 */
VectorXd inter1parray(const VectorXd& a, const int m);

/**
 * @brief Calculates the value of a parabola passing through three points.
 * 
 * @param x The x-coordinate at which to evaluate the parabola.
 * @param f_1 The y-coordinate of the point with x-coordinate -1.
 * @param f0 The y-coordinate of the point with x-coordinate 0.
 * @param f1 The y-coordinate of the point with x-coordinate 1.
 * @return long double The value of the parabola at x.
 */
long double parabola(const long double x, const long double f_1, const long double f0, const long double f1);

/**
 * @brief Quadratically interpolates a value x in an array a.
 * 
 * @param x The value to be interpolated.
 * @param a The array to interpolate.
 * @return long double The interpolated value.
 */
long double interp2(const long double x, const VectorXd& a);
