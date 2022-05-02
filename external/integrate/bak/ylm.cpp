#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include "libIntegrate/src/libIntegrate/_2D/GaussianQuadratures/GaussLegendre.hpp"
//#include <GaussLegendre.hpp>
using namespace std;

//#include <gsl/gsl_integration.h>

long double sph_norm(const int l, const int m, const double theta, const double phi){
//	complex<long double> r;
	long double Re, Im;
	//long double Glm;
	Re=boost::math::spherical_harmonic_r(l, m, theta, phi);
    Im=boost::math::spherical_harmonic_i(l, m, theta, phi);
	return (Re*Re + Im*Im)*std::sin(theta);
}

// Alias for l=0, m=0
long double sph_norm_00(const double theta, const double phi){
	return sph_norm(0, 0, theta, phi);
}
// Alias for l=1, m=+/-1
long double sph_norm_11(const double theta, const double phi){
	return sph_norm(1, 1, theta, phi);
}
// Alias for l=1, m=0
long double sph_norm_10(const double theta, const double phi){
	return sph_norm(1, 0, theta, phi);
}
// Alias for l=2, m=+/-2
long double sph_norm_22(const double theta, const double phi){
	return sph_norm(2, 2, theta, phi);
}
// Alias for l=2, m=+/-1
long double sph_norm_21(const double theta, const double phi){
	return sph_norm(2, 1, theta, phi);
}
// Alias for l=2, m=0
long double sph_norm_20(const double theta, const double phi){
	return sph_norm(2, 0, theta, phi);
}
// Alias for l=3, m=+/-3
long double sph_norm_33(const double theta, const double phi){
	return sph_norm(3, 3, theta, phi);
}
// Alias for l=3, m=+/-2
long double sph_norm_32(const double theta, const double phi){
	return sph_norm(3, 2, theta, phi);
}
// Alias for l=3, m=+/-1
long double sph_norm_31(const double theta, const double phi){
	return sph_norm(3, 1, theta, phi);
}
// Alias for l=3, m=0
long double sph_norm_30(const double theta, const double phi){
	return sph_norm(3, 0, theta, phi);
}

long double Glm(const int l, const int m, const long double theta_min, const long double theta_max){

	//_2D::GQ::GaussLegendreQuadrature<double,64> integrate;
	_2D::GQ::GaussLegendreQuadrature<double,16> integrate;

	const long double phi_min=0;
	const long double phi_max=2.*M_PI;
	long double r;

	switch (l){
		case 0:
			switch (abs(m)){
				case 0:
					r=integrate(sph_norm_00, theta_min, theta_max, phi_min, phi_max);
					break;
				default:
					std::cout << "Error: -l<m<l not respected. Will return -10" << std::endl;
					r=-10.;
					break;
			}
		break;
		case 1:
			switch (abs(m)){
				case 0:
					r=integrate(sph_norm_10, theta_min, theta_max, phi_min, phi_max);
					break;
				case 1:
					r=integrate(sph_norm_11, theta_min, theta_max, phi_min, phi_max);
					break;
				default:
					r=-10;
					std::cout << "Error: -l<m<l not respected. Will return -10" << std::endl;
					break;
			}
		break;
		case 2:
			switch (abs(m)){
				case 0:
					r=integrate(sph_norm_20, theta_min, theta_max, phi_min, phi_max);
					break;
				case 1:
					r=integrate(sph_norm_21, theta_min, theta_max, phi_min, phi_max);
					break;
				case 2:
					r=integrate(sph_norm_22, theta_min, theta_max, phi_min, phi_max);
					break;
				default:
					r=-10;
					std::cout << "Error: -l<m<l not respected. Will return -10" << std::endl;
					break;
			}
		break;
		case 3:
			switch (abs(m)){
				case 0:
					r=integrate(sph_norm_30, theta_min, theta_max, phi_min, phi_max);
					break;
				case 1:
					r=integrate(sph_norm_31, theta_min, theta_max, phi_min, phi_max);
					break;
				case 2:
					r=integrate(sph_norm_32, theta_min, theta_max, phi_min, phi_max);
					break;
				case 3:
					r=integrate(sph_norm_33, theta_min, theta_max, phi_min, phi_max);
					break;
				default:
					r=-10.;
					std::cout << "Error: -l<m<l not respected. Will return -10" << std::endl;
					break;
			}
		break;
		default:
			std::cout << "Error: l>3 not supported" << std::endl;
	}
	return r;
}

/* // For Debug 
int main(void){
	long double theta_min=0.;
	long double theta_max=M_PI/4.;
	const int lmax=3;
	int l, m;
	long double r;
	std::cout << "Integral:" << std::endl;
	std::cout << "    - Polar band of activity:" << std::endl;
	for (l=0; l<=lmax; l++){
		for (m=-l; m<=l;m++){
			r=0;
			r=Glm(l, m, theta_min, theta_max);
			std::cout << "(" << l << "," << m << ") :" << r << std::endl;
		}
		std::cout << " --- " << std::endl;
	}

	theta_min=0. + M_PI/2 - M_PI/4/2;
	theta_max=0. + M_PI/2 + M_PI/4/2;
	std::cout << "    - Equatorial band of activity:" << std::endl;
	for (l=0; l<=lmax; l++){
		for (m=-l; m<=l;m++){
			r=0;
			r=Glm(l, m, theta_min, theta_max);
			std::cout << "(" << l << "," << m << ") :" << r << std::endl;
		}
		std::cout << " --- " << std::endl;
	}

}
*/



