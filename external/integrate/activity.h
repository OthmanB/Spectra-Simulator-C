#pragma once
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <Eigen/Dense>
//#include "libIntegrate/src/libIntegrate/_2D/GaussianQuadratures/GaussLegendre.hpp"
#include "GaussLegendre2D.hpp"
using namespace std;


using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

long double sph_norm(const double theta, const double phi, const int l, const int m);
long double sph_norm2(const double theta, const double phi, const int l, const int m, const int dummy1, const int dummy2);
long double Glm(const int l, const int m, const long double theta_min, const long double theta_max);
long double Alm(const int l, const int m, const long double theta0, const long double delta, std::string ftype);
long double Alm_norm_gate(const double theta, const double phi, const int l, const int m, const long double theta0, const long double delta);
long double Alm_norm_gauss(const double theta, const double phi, const int l, const int m, const long double theta0, const long double delta);
VectorXd gauss_filter(const VectorXd theta, const long double theta0, const long double delta);
VectorXd gate_filter(const VectorXd theta, const long double theta0, const long double delta);
long double gauss_filter_cte(const long double theta0, const long double delta);
