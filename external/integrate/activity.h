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

long double sph_norm(const long double theta, const long double phi, const int l, const int m);
long double sph_norm2(const long double theta, const long double phi, const int l, const int m, const int dummy1, const int dummy2);
//long double Glm(const int l, const int m, const long double theta_min, const long double theta_max);
long double integrate_Alm_gate(const int l, const int m, const long double theta0, const long double delta);
long double Alm(const int l, const int m, const long double theta0, const long double delta, std::string ftype); // theta0 and  delta provided in radian here
long double Alm_deg(const int l, const int m, const long double theta0, const long double delta, std::string ftype); // theta0 and delta provided in degress here
long double Alm_norm_gate(const long double theta, const long double phi, const int l, const int m, const long double theta0, const long double delta);
long double Alm_norm_gauss(const long double theta, const long double phi, const int l, const int m, const long double theta0, const long double delta);
VectorXd gauss_filter(const VectorXd theta, const long double theta0, const long double delta);
VectorXd gate_filter(const VectorXd theta, const long double theta0, const long double delta);
long double gauss_filter_cte(const long double theta0, const long double delta);

long double Alm_norm_gate_2pi(const long double theta, const long double phi, const int l, const int m, const long double theta0, const long double delta);
long double Alm_norm_gauss_2pi(const long double theta, const long double phi, const int l, const int m, const long double theta0, const long double delta);
VectorXd gauss_filter_2pi(const VectorXd theta, const long double theta0, const long double delta);
VectorXd gate_filter_2pi(const VectorXd theta, const long double theta0, const long double delta);
long double gauss_filter_cte_2pi(const long double theta0, const long double delta);
