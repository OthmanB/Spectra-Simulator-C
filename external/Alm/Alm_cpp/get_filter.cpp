#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include "GaussLegendre2D.hpp"
#include <Eigen/Dense>
#include "activity.h"
//#include "linspace.h"
#include "string_handler.h"

void usage(int argc, char* argv[]);
int  options(int argc, char* argv[]);

int main(int argc, char* argv[]){
    const int Npts=1000;
	const long double delta_limit=0.001;
	const std::string fileout="../data/outputs/filter.ascii";
	int msg_code;
	int l;
	double Alm_norm;
	double theta0, delta, theta0_limit;
	std::string ftype, tmp;
    VectorXd r;
    VectorXd theta;
	VectorXd range(2);
	std::ofstream fileout_stream;

	msg_code=options(argc, argv);
	if(msg_code == -1){
		std::cout << "Error detected in options. Cannot proceed. Debug required." << std::endl;
		std::cout << "The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	} else{
//		for (int i=0;i<4; i++){
//			std::cout << argv[i] << std::endl;
//		}
//		std::cout << " ---- " << std::endl;
		std::istringstream(argv[1]) >> theta0;	
		std::istringstream(argv[2]) >> delta;
		std::istringstream(argv[3]) >> ftype;
		std::istringstream(argv[4]) >> tmp;
		range=str_to_Xdarr(tmp, ",");
		if(range[0] == -1){
			range[0]=0;
		}
		if(range[1] == -1){
			range[1]=M_PI/2;
		}
	}
	//theta=linspace(range[0], range[1], Npts);
	theta=Eigen::VectorXd::LinSpaced(Npts, range[0], range[1]);
	theta0_limit=M_PI;
	if (ftype != "gauss" && ftype != "gate" && ftype != "triangle"){
		std::cout << "Error for ftype: Only 'gate', 'gauss' or 'triangle' is allowed" << std::endl;
		exit(EXIT_FAILURE);
	}
	//std::cout << "#use2pi = " << use2pi << std::endl;
	if (delta >=delta_limit && theta0>=0 && theta0< theta0_limit){
		if (ftype == "gauss"){
			Alm_norm=gauss_filter_cte(theta0, delta);
		} else{
			Alm_norm=1;
		}
        if (ftype == "triangle"){
            r=triangle_filter(theta, theta0, delta);
        }
        if (ftype == "gauss"){
            r=gauss_filter(theta, theta0, delta);
        }
        if (ftype == "gate"){
            r=gate_filter(theta, theta0, delta);
        }
	}
/*
	std::cout << "# ftype = " << ftype << std::endl;
    std::cout << "# theta0=" << theta0 << std::endl;
    std::cout << "# delta=" << delta << std::endl;
    std::cout << "# theta   /  F(theta, theta0, delta)" << std::endl;
    for (int k=0; k<Npts;k++){
        std::cout << theta[k]  << setw(20) << std::setprecision(12) << r[k] << std::endl; 
    }
*/
	fileout_stream.open(fileout.c_str());
    	if(fileout_stream.is_open()){
			fileout_stream << "# ftype = " << ftype << std::endl;
    		fileout_stream << "# theta0=" << theta0 << std::endl;
    		fileout_stream << "# delta=" << delta << std::endl;
			fileout_stream << "# theta   /  F(theta, theta0, delta)" << std::endl;
		    for (int k=0; k<Npts;k++){
				fileout_stream << theta[k]  << setw(20) << std::setprecision(10) << r[k] << std::endl;
			}
    	} else{
			std::cout << " Unable to open the binary data file " << fileout.c_str() << std::endl;	
			std::cout << " Check that the full path exists" << std::endl;
			std::cout << " The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
    	}
    fileout_stream.close();
	std::cout << " File written succesfully in: " << fileout.c_str() << std::endl;
}


int options(int argc, char* argv[]){
	int val=-1;	
	
	if(argc == 5){
		val=1; 
	} 
	if (val == -1){ // Error code
		usage(argc, argv);
	} 
	if (val > 0 ){
		return val; // Execution code val
	} else{
		return -1; // Default value is to return an error code
	}
}


void usage(int argc, char* argv[]){
            std::cout << "Compute the Filter term of Alm and write a table for all theta between <range>=min,max (no space allowed between the , separator) in a file in data/filter.ascii" << std::endl;
			std::cout << "   <range> can be set to -1,-1 if one wants a range [0, Pi/2] (default calculation range)" << std::endl;
			std::cout << "Unrecognized arguments" << std::endl;
			std::cout << "     - To execute: " << argv[0] << " <active_region_colatitude_theta0_in_rad>  <active_region_width_in_rad> <filter_type>  <range>" << std::endl;
			//std::cout << "     - To show version: " << argv[0] << " version" << std::endl;
			exit(EXIT_FAILURE);
}
