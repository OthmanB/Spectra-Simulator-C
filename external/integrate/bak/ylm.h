#pragma once
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include "libIntegrate/src/libIntegrate/_2D/GaussianQuadratures/GaussLegendre.hpp"
//#include <GaussLegendre.hpp>
using namespace std;

long double sph_norm(const int l, const int m, const double theta, const double phi);
long double Glm(const int l, const int m, const long double theta_min, const long double theta_max);

long double sph_norm_00(const double theta, const double phi);
// Alias for l=1, m=+/-1
long double sph_norm_11(const double theta, const double phi);
// Alias for l=1, m=0
long double sph_norm_10(const double theta, const double phi);
// Alias for l=2, m=+/-2
long double sph_norm_22(const double theta, const double phi);
// Alias for l=2, m=+/-1
long double sph_norm_21(const double theta, const double phi);
// Alias for l=2, m=0
long double sph_norm_20(const double theta, const double phi);
// Alias for l=3, m=+/-3
long double sph_norm_33(const double theta, const double phi);
// Alias for l=3, m=+/-2
long double sph_norm_32(const double theta, const double phi);
// Alias for l=3, m=+/-1
long double sph_norm_31(const double theta, const double phi);
// Alias for l=3, m=0
long double sph_norm_30(const double theta, const double phi);