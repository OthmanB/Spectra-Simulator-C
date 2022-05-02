#pragma once
#include <Eigen/Dense>
#include <cmath>

using Eigen::VectorXd;

double lin_interpol(const VectorXd& x, const VectorXd& y, const double x_int);
VectorXd quad_interpol( const VectorXd& a, const int m );
const long double interp1( const long double x,  const VectorXd& a);
long double parabola( const long double x, const long double f_1, const long double f0, const long double f1 );
long double interp2( const long double x, const VectorXd& a);
