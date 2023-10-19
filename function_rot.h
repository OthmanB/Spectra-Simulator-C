/**
 * @file function_rot.h
 * @brief Functions that compute the mode heights from the stellar inclination.
 * 
 * 
 * @date 11 Feb 2016
 * @author obenomar
 */ 

#include <math.h>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;


/**
 * @brief Calculate the combination coefficient
 *
 * This function calculates the combination coefficient (n choose r) based on the given n and r values.
 *
 * @param n The n value
 * @param r The r value
 * @return The calculated combination coefficient
 *
 * Dependencies: This function requires the following dependencies:
 * - factorial function
 */
double combi(int n, int r);

/**
 * @brief Calculate the dmm coefficient for mode visibilities
 *
 * This function calculates the dmm coefficient for mode visibilities based on the given degree of the mode, m1 and m2 values, and stellar inclination.
 *
 * @param l The degree of the mode
 * @param m1 The m1 value
 * @param m2 The m2 value
 * @param beta The stellar inclination in radians
 * @return The calculated dmm coefficient
 *
 * Dependencies: This function requires the following dependencies:
 * - math.h
 * - combi function
 * - factorial function
 */
double dmm(int l, int m1, int m2, double beta);

/**
 * @brief Calculate the factorial of a number
 *
 * This function calculates the factorial of a given number.
 *
 * @param n The number to calculate the factorial of
 * @return The calculated factorial
 */
int factorial(int n);

/**
 * @brief Calculate the rotation matrix for mode visibilities
 *
 * This function calculates the rotation matrix for mode visibilities based on the given degree of the mode and stellar inclination.
 *
 * @param l The degree of the mode
 * @param beta The stellar inclination in degrees
 * @return A MatrixXd of size 2*l+1 x 2*l+1 containing the rotation matrix
 */
MatrixXd function_rot( int l, double beta);

/**
 * @brief Calculate the amplitude ratio of mode visibilities
 *
 * This function calculates the amplitude ratio of mode visibilities based on the given degree of the mode and stellar inclination.
 *
 * @param l The degree of the mode
 * @param beta The stellar inclination in degrees
 * @return A VectorXd of size 2*l+1 containing the visibilities of the m components
 *
 * Dependencies: This function requires the following dependencies:
 *   - math.h
 *   - Eigen library (version 3.1)
 */
VectorXd amplitude_ratio(int l, double beta);

