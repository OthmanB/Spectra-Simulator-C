#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include "GaussLegendre2D.hpp"
#include <Eigen/Dense>
#include "activity.h"
#include "linspace.h"

void usage(int argc, char* argv[]);
int  options(int argc, char* argv[]);

int main(int argc, char* argv[]){
	const long double delta_limit=0.001;

	//const bool use2pi=0;
	int l;
	double r, Alm_norm;
	double theta0, delta, theta0_limit;
	//const double theta0=M_PI/2;
	//const double delta=M_PI/6;
	std::string ftype;

	int msg_code;
	
	msg_code=options(argc, argv);
	if(msg_code == -1){
		std::cout << "Error detected in options. Cannot proceed. Debug required." << std::endl;
		std::cout << "The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	} else{
//		for (int i=0;i<5; i++){
//			std::cout << argv[i] << std::endl;
//		}
//		std::cout << " ---- " << std::endl;
		std::istringstream(argv[1]) >> l;
		std::istringstream(argv[2]) >> theta0;	
		std::istringstream(argv[3]) >> delta;
		std::istringstream(argv[4]) >> ftype;
		//std::istringstream(argv[5]) >> use2pi;
	}
	theta0_limit=M_PI;
	/*
	if (use2pi == 0){
		theta0_limit=M_PI;
	}
	else{
		theta0_limit=2*M_PI;
	}
	*/
	if (ftype != "gauss" && ftype != "gate"){
		std::cout << "Error for ftype: Only 'gate' and 'gauss' is allowed" << std::endl;
		exit(EXIT_FAILURE);
	}
	std::cout << "#---------------" << std::endl;
	std::cout << "#Configuration: " << std::endl;
	std::cout << "#l=" << l << std::endl;
	std::cout << "#theta0 =" << theta0 << std::endl;
	std::cout << "#delta = " << delta << std::endl;
	std::cout << "#ftype = " << ftype << std::endl;
	//std::cout << "#use2pi = " << use2pi << std::endl;
	std::cout << "#--------------"  << std::endl;
	std::cout << "#l     m      Alm" << std::endl;
	if (delta >=delta_limit && theta0>=0 && theta0< theta0_limit){
		if (ftype == "gauss"){
			Alm_norm=gauss_filter_cte(theta0, delta);
		} else{
			Alm_norm=1;
		}
		//std::cout << "Alm_norm =" << Alm_norm << std::endl;
		for (int m=-l;m<=l; m++){
			r=Alm(l, m, theta0, delta, ftype);
			//std::cout << "r =" << r << std::endl;
			std::cout << l << "   " << m << "  " << r/Alm_norm << std::endl;
		}	
	}
	if (delta <delta_limit && theta0>=0  && theta0< theta0_limit){
		for (int m=-l;m<=l; m++){
			std::cout << l << "   " << m << "  " << 0 << std::endl;
		}
	}
	if (theta0<0 || theta0 > theta0_limit){
		for (int m=-l;m<=l; m++){
			std::cout << l << "   " << m << "  " << -9999 << std::endl;
		}
	}	
}


int options(int argc, char* argv[]){
	int val=-1;	
	
	if(argc == 5){
		val=1; 
	} 
	if (val == -1){ // Error code
		usage(argc, argv);
	} 
//	if (val == 0){ // Version code
//		showversion(); 
//		exit(EXIT_SUCCESS);
//	}
	if (val > 0 ){
		return val; // Execution code val
	} else{
		return -1; // Default value is to return an error code
	}
}


void usage(int argc, char* argv[]){

			std::cout << "Unrecognized argument" << std::endl;
			//std::cout << "     - To execute: " << argv[0] << " <mode_degree> <active_region_colatitude_theta0_in_rad>  <active_region_width_in_rad> <filter_type> <use2pi>" << std::endl;
			std::cout << "     - To execute: " << argv[0] << " <mode_degree> <active_region_colatitude_theta0_in_rad>  <active_region_width_in_rad> <filter_type>" << std::endl;
			//std::cout << "     - To show version: " << argv[0] << " version" << std::endl;
			exit(EXIT_FAILURE);
}