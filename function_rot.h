/*
 * function_rot.h
 *
 *  Created on: 11 Feb 2016
 *      Author: obenomar
 */
#include <math.h>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

double combi(int n, int r);
double dmm(int l, int m1, int m2, double beta);
int factorial(int n);
MatrixXd function_rot( int l, double beta);
VectorXd amplitude_ratio(int l, double beta);

