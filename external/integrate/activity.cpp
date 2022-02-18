#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
//#include "libIntegrate/src/libIntegrate/_2D/GaussianQuadratures/GaussLegendre.hpp"
#include "GaussLegendre2D.hpp"
#include <Eigen/Dense>
#include "activity.h"
#include "linspace.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

//#include <GaussLegendre.hpp>
using namespace std;

//#include <gsl/gsl_integration.h>

long double sph_norm(const long double theta, const long double phi, const int l, const int m){
	long double Re, Im;
	Re=boost::math::spherical_harmonic_r(l, m, theta, phi);
    Im=boost::math::spherical_harmonic_i(l, m, theta, phi);
	return (Re*Re + Im*Im)*std::sin(theta);
}


long double sph_norm2(const long double theta, const long double phi, const int l, const int m, const int dummy1, const int dummy2){
	long double Re, Im;
	Re=boost::math::spherical_harmonic_r(l, m, theta, phi);
    Im=boost::math::spherical_harmonic_i(l, m, theta, phi);
	return (Re*Re + Im*Im)*std::sin(theta);
}

//long double Alm_norm_gate(const double theta, const double phi, const int l, const int m, const long double theta0, const long double delta){
long double Alm_norm_gate(const long double theta, const long double phi, const int l, const int m, const long double theta0, const long double delta){
	long double sph, Fmax;
  	VectorXd tmp(1),F;
	sph=sph_norm(theta, phi, l,m);
  	tmp[0]=theta;
    F=gate_filter(tmp, theta0, delta);
    return sph*F[0];
}

//long double Alm_norm_gate(const double theta, const double phi, const int l, const int m, const long double theta0, const long double delta){
long double Alm_norm_gate_2pi(const long double theta, const long double phi, const int l, const int m, const long double theta0, const long double delta){
	long double sph, Fmax;
  	VectorXd tmp(1),F;
	sph=sph_norm(theta, phi, l, m);
  	tmp[0]=theta;
    F=gate_filter_2pi(tmp, theta0, delta);
    return sph*F[0];
}

long double Alm_norm_gauss(const long double theta, const long double phi, const int l, const int m, const long double theta0, const long double delta){
	long double sph, Fmax;
	VectorXd tmp(1),F;
	sph=sph_norm(theta, phi, l, m);
	tmp[0]=theta;
	F=gauss_filter(tmp, theta0, delta);
	return sph*F[0]; // Note on the normalisation: The normalisation using gauss_filter_cte() happens in main.cpp for optimisation purpose 
	//Fmax=gauss_filter_cte(theta0, delta);
	//F=F/Fmax;
}

long double Alm_norm_gauss_2pi(const long double theta, const long double phi, const int l, const int m, const long double theta0, const long double delta){
	long double sph, Fmax;
	VectorXd tmp(1),F;
	sph=sph_norm(theta, phi, l, m);
	tmp[0]=theta;
	F=gauss_filter_2pi(tmp, theta0, delta);
	return sph*F[0]; // Note on the normalisation: The normalisation using gauss_filter_cte() happens in main.cpp for optimisation purpose 
	//Fmax=gauss_filter_cte(theta0, delta);
	//F=F/Fmax;	
}

// A funcion used to perform a Gaussian filtering on the integral term for the 
// 
VectorXd gauss_filter(const VectorXd theta, const long double theta0, const long double delta){
	long double a0;
	long double a1;
	VectorXd F(theta.size());
	for (int i=0; i<theta.size(); i++){
		a0=std::pow(theta[i] - theta0,2)/(2*std::pow(delta,2));
		a1=std::pow(theta[i] - M_PI + theta0, 2)/(2*std::pow(delta,2));	
		F[i]=std::exp(-a0) + std::exp(-a1);
	} 
	return F;
}

// A funcion used to perform a Gaussian filtering on the integral term for the 
// 
VectorXd gauss_filter_2pi(const VectorXd theta, const long double theta0, const long double delta){
	long double a0;
	long double a1, a2;
	VectorXd F(theta.size());
	for (int i=0; i<theta.size(); i++){
		a0=std::pow(theta[i] - theta0,2)/(2*std::pow(delta,2));
		a1=std::pow(theta[i] - M_PI + theta0, 2)/(2*std::pow(delta,2));	
		a2=std::pow(theta[i] - 2*M_PI + theta0, 2)/(2*std::pow(delta,2));	
		F[i]=std::exp(-a0) + std::exp(-a1) + std::exp(-a2);
	} 
	return F;
}

VectorXd gate_filter(const VectorXd theta, const long double theta0, const long double delta){
	VectorXd F(theta.size());
	F.setZero();
	const long double distance_critic=0.4;
	long double distance=M_PI - 2*theta0 - delta;
	for (int i=0; i<theta.size(); i++){
		if ((theta[i] >= (theta0 - delta/2) && theta[i] <= (theta0 + delta/2)) || 
		   (theta[i] >= (M_PI - theta0 - delta/2) && theta[i] <= (M_PI - theta0 + delta/2)) ) {
		   	F[i]=1;
		   }
		 //if (distance<=distance_critic && theta[i] >=(theta0+delta/2) && theta[i] <=(M_PI-theta0-delta/2)){
		 //	std::cout << "distance: " << distance << std::endl;
		 //	std::cout << "(theta0+delta/2) =" << (theta0+delta/2) << std::endl;
		 //	std::cout << "(M_PI-theta0-delta/2) =" << M_PI-theta0-delta/2<< std::endl;
		 //	std::cout << "theta[i] =" << theta[i] << std::endl;
		 //	std::cout << "---------------" << std::endl; 
		 //}
	}
	return F;
}

VectorXd gate_filter_2pi(const VectorXd theta, const long double theta0, const long double delta){
	VectorXd F(theta.size());
	F.setZero();
	for (int i=0; i<theta.size(); i++){
		if ((theta[i] >= (theta0 - delta/2) && theta[i] <= (theta0 + delta/2)) || 
		   (theta[i] >= (M_PI - theta0 - delta/2) && theta[i] <= (M_PI - theta0 + delta/2)) ||
		   (theta[i] >= (2*M_PI - theta0 - delta/2) && theta[i] <= (2*M_PI - theta0 + delta/2)) ){
		   	F[i]=1;
		   }
	}
	return F;
}

long double gauss_filter_cte(const long double theta0, const long double delta){
	/*
		A function that calculate what is the max value of a double-gaussian filter
		that is symetrical towards pi/2 and bounded between [0, pi] and of width delta
                     -                      -
		          -    -                  -   -
		         -       -              -      -
		        -         -            -        -
		       -            -         -          - 
		      -              -       -            -
		    -                  -   -                -
		  -                      -                    -
		--+----------+-----------+---------------------+----> theta
		  0         theta0      pi/2    pi-theta0      pi
	*/

	VectorXd theta= linspace(0, M_PI, 200);
	VectorXd F;
	F=gauss_filter(theta, theta0, delta);
	return F.maxCoeff();
}

long double gauss_filter_cte_2pi(const long double theta0, const long double delta){
	/*
		A function that calculate what is the max value of a triple-gaussian filter
		that is symetrical towards pi/2 and 3pi/2 and bounded between [0, 2pi] and of width delta
          |    -                                         -    
		  |  -    -                                    -   -                                        |
		  | -       -                                 -      -                                      |
		  |-         -                               -        -                                     |
		  -            -                            -          -                                   -|
		  |             -                          -            -                                 - |
		  |              -                       -                -                              -  |
		  |                -                    -                   -                           -   |
		--+-----+----------------+----------------------+-----+-------------------------------------+------> theta
		  0   theta0           pi/2               pi-theta0    pi                                  2pi
	*/
	VectorXd theta= linspace(0, M_PI, 200);
	VectorXd F;
	F=gauss_filter_2pi(theta, theta0, delta);
	return F.maxCoeff();
}

/*
long double Glm(const int l, const int m, const long double theta_min, const long double theta_max){
	_2D::GQ::GaussLegendreQuadrature<double,64> integrate;
	//_2D::GQ::GaussLegendreQuadrature<double,16> integrate;

	const long double phi_min=0;
	const long double phi_max=2.*M_PI;
	const int dummy1=0;
	const int dummy2=0;
	long double r;
	if (std::abs(m)<=l){
		r=integrate(sph_norm2, theta_min, theta_max, phi_min, phi_max, l, m, dummy1, dummy2);
	} else{
		r=-10;
		std::cout << "Glm Error: -l<m<l not respected. Will return -10" << std::endl;
	}
	return r;
}
*/

long double integrate_Alm_gate(const int l, const int m, const long double theta0, const long double delta){
	_2D::GQ::GaussLegendreQuadrature<double,64> integrate;
	const long double theta_min=theta0-delta/2; // Default for ftype='gate'
	const long double theta_max=theta0+delta/2;
	const long double phi_min=0;
	const long double phi_max=2.*M_PI;
	const int dummy =0;
	long double r;
	if (delta != 0){
			r = integrate(sph_norm2, theta_min, theta_max, phi_min, phi_max, l, m, dummy, dummy);
	} else{
		r=0; // When delta is 0, obviously the result is 0
	}
	return r;
}

long double Alm(const int l, const int m, const long double theta0, const long double delta, std::string ftype){

	_2D::GQ::GaussLegendreQuadrature<double,64> integrate;
	//_2D::GQ::GaussLegendreQuadrature<double,16> integrate;
	const long double theta_min=0; // Default for ftype='gate'
	const long double theta_max=M_PI;
	const long double phi_min=0;
	const long double phi_max=2.*M_PI;

	long double r;
	if (std::abs(m)<=l){
		if (ftype == "gate"){
			r=integrate_Alm_gate(l, m, theta0, delta); //A function that integrate directly using Ylm (no filtering, which is faster)
			r=2*r; // Accouting for the North + South Emisphere
			//std::cout << " North + South Emisphere: " << r << std::endl; // FOR DEBUG ONLY
			//std::cout << " South Emisphere: "<<integrate_Alm_gate(l, m, M_PI-theta0, delta) << std::endl; // FOR DEBUG ONLY
			//r=integrate(Alm_norm_gate, theta_min, theta_max, phi_min, phi_max, l, m, theta0, delta);
			//std::cout << " Use of Alm_norm_gate: " << r << std::endl;
			//exit(EXIT_FAILURE);
		}
		if (ftype == "gauss"){
			//theta_min=0;
			//theta_max=M_PI;
			r=integrate(Alm_norm_gauss, theta_min, theta_max, phi_min, phi_max, l, m, theta0, delta);
		}
	} else{
		r=-10;
		std::cout << "Alm Error: -l<m<l not respected. Will return -10" << std::endl;
	}
	return r;
}

/*
// Mini test to compare filtering Alm and Glm 
int main(void){
	double r, Alm_norm;
	const double theta0=M_PI/2;
	const double delta=M_PI/6;
	const double theta_min=theta0-delta/2;
	const double theta_max=theta0+delta/2;

	Alm_norm=gauss_filter_cte(theta0, delta);
	std::cout << "Glm and Alm: " << std::endl;
	for (int l=1; l<=2; l++){
		for (int m=-l;m<=l; m++){
			r=Glm(l, m, theta_min, theta_max);
			std::cout << "Glm(l=" << l << ", m=" << m << ") = " << r << std::endl;
			r=Alm(l, m, theta0, delta, "gate");
			std::cout << "Alm_gate(l=" << l << ", m=" << m << ") = " << r << std::endl;
			r=Alm(l, m, theta0, delta, "gauss");
			std::cout << "Alm_gauss(l=" << l << ", m=" << m << ") = " << r/Alm_norm << std::endl;
			std::cout << " --- " << std::endl;
		}	
	}
}
*/

/*
int main(void){
	const double delta = M_PI/6;
	const double theta0= M_PI/2;
	VectorXd theta(1);
	VectorXd F;
	double Fmax;
	Fmax=gauss_filter_cte(theta0, delta);
	for (int i=0; i<=20; i++){
		theta[0]=i*M_PI/20;
		//F=gate_filter(theta, theta0, delta);
		F=gauss_filter(theta, theta0, delta);
		std::cout << "theta : " << theta << "     F=" << F/Fmax << std::endl;
	}
}
*/

 // For Debug 
/*
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
			//r=Glm(l, m, theta_min, theta_max);
			std::cout << "(" << l << "," << m << ") :" << r << std::endl;
		}
		std::cout << " --- " << std::endl;
	}

}
*/

