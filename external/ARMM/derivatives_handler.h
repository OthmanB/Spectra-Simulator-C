/**
 * @file derivatives_handler.h
 * @brief A set of function to compute derivation
 * 
 * This is a set of function that can compute the first and second derivatives for any discrete function
 */
#pragma once
#include <iostream>
#include <Eigen/Dense>
#include "../../data.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;


// ---------------------------------------------------------
// -------------- **** 1st Derivative ***** ----------------
// ---------------------------------------------------------

/**
 * @brief Computes the first order backward derivative with 1st order precision.
 * 
 * @param y A vector containing 2 points: y at x-1 (y[0]) and y at x (y[1]).
 * @return Deriv_out A structure containing the computed derivative.
 * 
 * @note The vector y must have a size of at least 2.
 *       If the size is less than 2, the program will exit with an error message.
 */
Deriv_out Fstder_backward_1err_reggrid(const VectorXd& y);

/**
 * @brief Computes the first order backward derivative with 1st order precision.
 *
 * @param y A vector containing 2 points: y at x-1 (y[0]) and y at x (y[1]).
 * @param x A vector containing the corresponding x values for the points in y.
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vectors y and x must have a size of at least 2.
 * If the size is less than 2, the program will exit with an error message.
 */
Deriv_out Fstder_backward_1err_reggrid(const VectorXd& y, const VectorXd& x);

/**
 * @brief Computes the first order forward derivative with 1st order precision.
 *
 * @param y A vector containing 2 points: y at x (y[0]) and y at x+1 (y[1]).
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vector y must have a size of at least 2.
 * If the size is less than 2, the program will exit with an error message.
 */
Deriv_out Fstder_forward_1err_reggrid(const VectorXd& y);

/**
 * @brief Computes the first order forward derivative with 1st order precision.
 *
 * @param y A vector containing 2 points: y at x (y[0]) and y at x+1 (y[1]).
 * @param x A vector containing the corresponding x values for the points in y.
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vectors y and x must have a size of at least 2.
 * If the size is less than 2, the program will exit with an error message.
 */
Deriv_out Fstder_forward_1err_reggrid(const VectorXd& y, const VectorXd& x);

/**
 * @brief Computes the first order centered derivative with 2nd order precision.
 *
 * @param y A vector containing 3 points: y at x-1 (y[0]), y at x (y[1]), and y at x+1 (y[2]).
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vector y must have a size of at least 3.
 * If the size is less than 3, the program will exit with an error message.
 */
Deriv_out Fstder_centered_2err_reggrid(const VectorXd& y);

/**
 * @brief Computes the first order centered derivative with 2nd order precision.
 *
 * @param y A vector containing 3 points: the position y at x-1 (y[0]), y at x (y[1]), and y at x+1 (y[2]).
 * @param x A vector containing the corresponding x values for the points in y.
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vectors y and x must have a size of at least 3.
 * If the size is less than 3, the program will exit with an error message.
 */
Deriv_out Fstder_centered_2err_reggrid(const VectorXd& y, const VectorXd& x);

/**
 * @brief Computes the first order centered derivative with 4th order precision.
 *
 * @param y A vector containing 5 points: the position y at x-2 (y[0]) to y at x+2 (y[4]).
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vector y must have a size of at least 5.
 * If the size is less than 5, the program will exit with an error message.
 */
Deriv_out Fstder_centered_4err_reggrid(const VectorXd y);

/**
 * @brief Computes the first order centered derivative with 4th order precision.
 *
 * @param y A vector containing 5 points: the position y at x-2 (y[0]) to y at x+2 (y[4]).
 * @param x A vector containing the corresponding x values for the points in y.
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vectors y and x must have a size of at least 5.
 * If the size is less than 5, the program will exit with an error message.
 */
Deriv_out Fstder_centered_4err_reggrid(const VectorXd& y, const VectorXd& x);

/**
 * @brief Computes the second order backward derivative with 1st order precision.
 *
 * @param y A vector containing 3 points: the position y at x-2 (y[0]) to y at x (y[2]).
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vector y must have a size of at least 3.
 * If the size is less than 3, the program will exit with an error message.
 */
Deriv_out Scndder_backward_1err_reggrid(const VectorXd& y);

/**
 * @brief Computes the second order backward derivative with 1st order precision.
 *
 * @param y A vector containing 3 points: the position y at x-2 (y[0]) to y at x (y[2]).
 * @param x A vector containing the corresponding x values for the points in y.
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vectors y and x must have a size of at least 3.
 * If the size is less than 3, the program will exit with an error message.
 */
Deriv_out Scndder_backward_1err_reggrid(const VectorXd& y, const VectorXd& x);

/**
 * @brief Computes  the second order forward derivative with 1st order precision.
 *
 * @param y A vector containing 3 points: the position y at x (y[0]) to y at x+2 (y[2]).
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vector y must have a size of at least 3.
 * If the size is less than 3, the program will exit with an error message.
 */
Deriv_out Scndder_forward_1err_reggrid(const VectorXd& y);

/**
 * @brief Computes the second order forward derivative with 1st order precision.
 *
 * @param y A vector containing 3 points: the position y at x (y[0]) to y at x+2 (y[2]).
 * @param x A vector containing the corresponding x values for the points in y.
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vectors y and x must have a size of at least 3.
 * If the size is less than 3, the program will exit with an error message.
 */
Deriv_out Scndder_forward_1err_reggrid(const VectorXd& y, const VectorXd& x);

/**
 * @brief Computes the second order centered derivative with second order precision.
 *
 * @param y A vector containing 3 points: the position y at x-1 (y[0]) to y at x+1 (y[2]).
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vector y must have a size of at least 3.
 * If the size is less than 3, the program will exit with an error message.
 */
Deriv_out Scndder_centered_2err_reggrid(const VectorXd& y);

/**
 * @brief Computes the second order centered derivative with second order precision.
 *
 * @param y A vector containing 3 points: the position y at x-1 (y[0]) to y at x+1 (y[2]).
 * @param x A vector containing the corresponding x values for the points in y.
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vectors y and x must have a size of at least 3.
 * If the size is less than 3, the program will exit with an error message.
 */
Deriv_out Scndder_centered_2err_reggrid(const VectorXd& y, const VectorXd& x);

// Main Functions

/**
 * @brief Computes the first derivative of a list of values with adaptive precision.
 * 
 * The computation is adaptive in the sense that the precision order of the derivative is automatically adjusted depending on the region (lower edge, center, upper edge) in which the derivative is computed.
 * 
 * @param y A vector containing the list of values.
 * @return Deriv_out A structure containing the computed derivative.
 * 
 * @note The vector y must have a size greater than or equal to 2.
 * If the size is less than 2, the behavior is undefined.
 */
Deriv_out Frstder_adaptive_reggrid(const VectorXd& y);

/**
 * @brief Computes the first derivative of a list of values with adaptive precision.
 *
 * The computation is adaptive in the sense that the precision order of the derivative is automatically adjusted depending on the region (lower edge, center, upper edge) in which the derivative is computed.
 *
 * @param y A vector containing the list of values.
 * @param x A vector containing the corresponding x values for the points in y.
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vectors y and x must have the same size and have a size greater than or equal to 2.
 * If the size is less than 2 or the sizes of y and x are not equal, the behavior is undefined.
 */
Deriv_out Frstder_adaptive_reggrid(const VectorXd& y, const VectorXd& x);

/**
 * @brief Computes the second derivative of a list of values with adaptive precision.
 *
 * The computation is adaptive in the sense that the precision order of the derivative is automatically adjusted depending on the region (lower edge, center, upper edge) in which the derivative is computed.
 *
 * @param y A vector containing the list of values.
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vector y must have a size greater than or equal to 3.
 * If the size is less than 3, the behavior is undefined.
 */
Deriv_out Scndder_adaptive_reggrid(const VectorXd& y);

/**
 * @brief Computes the second derivative of a list of values with adaptive precision.
 *
 * The computation is adaptive in the sense that the precision order of the derivative is automatically adjusted depending on the region (lower edge, center, upper edge) in which the derivative is computed.
 *
 * @param y A vector containing the list of values.
 * @param x A vector containing the corresponding x values for the points in y.
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vectors y and x must have the same size and have a size greater than or equal to 3.
 * If the size is less than 3 or the sizes of y and x are not equal, the behavior is undefined.
 */
Deriv_out Scndder_adaptive_reggrid(const VectorXd& y, const VectorXd& x);
