/*
 * io_star_common.cpp
 *
 * Contains all kind of methods
 * used to process and/or encapsulate data
 * 
 *  Created on: 21 Jun 2021
 *      Author: obenomar
 */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include "io_star_common.h"

using Eigen::MatrixXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

// This is somewhat the old original code that was there before version 1.0.0
// Just kept in this function in case it becomes useful at some point 
/*
void common_case0(std::ofstream& outfile, const MatrixXd& mode_params){//, std::string modelname){
	const int Nchars=5;
	const int precision=2;
	const bool fix_a1=0;
	const bool fix_a3=0;

	// Writing the splitting switches and value of splitting and stellar inclination
	outfile << "# Configuration for all common parameters. Check the program for a list of known keywords" << std::endl;
	//outfile << "           model_fullname             " <<  modelname << std::endl;
	write_a1_key(outfile, mode_params, fix_a1, 0); // Last value is irrelevant if fix_a1=0
	outfile << "Asphericity_eta    Fix_Auto       1" << std::endl;
	write_a3_key(outfile, mode_params, fix_a3, 0);
	write_asym_key(outfile, 0, 0); 

	outfile << "Inclination       Uniform    " << mode_params.col(10).mean() << "   0   90" << std::endl;
	outfile << "Visibility_l1     Gaussian       1.5       1.5      0.05"  << std::endl;
	outfile << "Visibility_l2     Gaussian       0.53      0.53     0.05"  << std::endl;
	outfile << "Visibility_l3     Gaussian       0.08      0.08     0.02"  << std::endl;
}
*/

void common_model_MS_Global_a1a2a3_HarveyLike(std::ofstream& outfile, const MatrixXd& mode_params){
	const std::string modelname="model_MS_Global_a1a2a3_HarveyLike";
	double maxHeight=1000.;

	if (mode_params.col(2).maxCoeff()>1000.){
		maxHeight=mode_params.col(2).maxCoeff()*10;
	}

	outfile << "# Configuration for all common parameters. Check the program for a list of known keywords" << std::endl;
	outfile << "           model_fullname             " <<  modelname << std::endl;
	outfile << " freq_smoothness          bool          1.000000          2.000000" << std::endl;
    write_a1_key(outfile, mode_params,  0, 0); // Last value is irrelevant if fix_a1=0
    outfile << "                a2_0                 Fix          0.000000" << std::endl;
    outfile << "                a2_1                 Fix          0.000000" << std::endl;
    outfile << "                a2_2                 Fix          0.000000" << std::endl;
    write_a3_key(outfile, mode_params, 1, 0);
    write_asym_key(outfile, 0, 0);
	outfile << "Inclination       Uniform    " << mode_params.col(10).mean() << "   0   90" << std::endl;
	outfile << "Visibility_l1     Gaussian       1.5       1.5      0.05"  << std::endl;
	outfile << "Visibility_l2     Gaussian       0.53      0.53     0.05"  << std::endl;
	outfile << "Visibility_l3     Gaussian       0.08      0.08     0.02"  << std::endl;
	outfile << "Height     Jeffreys       1.000000          "<< maxHeight  << std::endl;
	outfile << "Width      Fix_Auto       1"  << std::endl;
	outfile << "trunc_c      50" << std::endl;
}

void common_model_MS_Global_a1etaAlma3_HarveyLike(std::ofstream& outfile, const MatrixXd& mode_params){
	const std::string modelname="model_MS_Global_a1etaAlma3_HarveyLike";
	double maxHeight=1000., min_eps_1=-0.01;

	if (mode_params.col(2).maxCoeff()>1000.){
		maxHeight=mode_params.col(2).maxCoeff()*10;
	}
	if (mode_params.col(6).mean()>0){
		std::cout << "Error: epsilon_1 was found to exceed 0... The simulations are currently only designed for the measure of a mean(epsilon) = epsilon_1" << std::endl;
		std::cout << "       In that situation epsilon_1 must be negative to account for the magnetic darkening of the disk" << std::endl;
		std::cout << "       The program will exit now" << std::endl;
		exit(EXIT_SUCCESS);
	}
	if (mode_params.col(6).mean()<-0.01){
		std::cout << "Warning: The value of epsilon_1 is less than -0.01. This is quite strong. Ensure that your setup is correct" << std::endl;
		std::cout << "         Adjusting default priors to accomodate of the input..." << std::endl;
		min_eps_1=4*mode_params.col(6).mean();
	}
	outfile << "# Configuration for all common parameters. Check the program for a list of known keywords" << std::endl;
	outfile << "           model_fullname             " <<  modelname << std::endl;
	outfile << " freq_smoothness                bool          1.000000          2.000000" << std::endl;
    write_a1_key(outfile, mode_params, 0, 0); // Last value is irrelevant if fix_a1=0
    outfile << " epsilon_0             Uniform          " << mode_params.col(6).mean() << "        " << std::setprecision(7)<<  min_eps_1 << "          0.000000" << std::endl; // We only test here the detection of a mean epsilon, which must be negative
    outfile << " epsilon_1             Fix              0.000000         " << std::endl; // This would be to test a slope in epsilon (set to 0 in here)
    outfile << " epsilon_2             Fix              0.000000         " << std::endl; // This would be to test a 2nd Order polynomial in epsilon (set to 0 in here)
    outfile << " theta0                Uniform          " << mode_params.col(7).mean()*M_PI/180. << "         0.000000          3.141592653589793" << std::endl; // Location of the active regions
    outfile << " Dtheta                Uniform          " << mode_params.col(8).mean()*M_PI/180. << "         0.000000         1.047197551196598" << std::endl; // Extent of the active region. Max is pi/3
    write_a3_key(outfile, mode_params, 1, 0);
    //write_asym_key(outfile, 0, 0); // We introduce the asymetry as it can interfer with a2 distorsion measurements
    write_asym_key(outfile, 1, 0); // We introduce the asymetry as it can interfer with a2 distorsion measurements
    outfile << " Inclination           Uniform          " << mode_params.col(11).mean() <<  "         0.000000          90.000000" << std::endl;
	outfile << " Visibility_l1         Gaussian         1.5       1.5      0.05"  << std::endl;
	outfile << " Visibility_l2         Gaussian         0.53      0.53     0.05"  << std::endl;
	outfile << " Visibility_l3         Gaussian         0.08      0.08     0.02"  << std::endl;
	outfile << " Height                Jeffreys         1.000000          "<< std::setprecision(4) << maxHeight  << std::endl;
	outfile << " Width                 Fix_Auto         1"  << std::endl;
	outfile << " trunc_c               30" << std::endl;

}

void write_asym_key(std::ofstream& outfile, const bool fix, const double fix_val){
	const int Nchars=5;
	const int precision=2;
	if (fix == 0){
	    outfile << " Asymetry              Jeffreys_abs      5.0          0.1         200" << std::endl;
	} else{
		outfile << " Asymetry              Fix               0"  << std::endl;
	}
}

void write_a1_key(std::ofstream& outfile, const MatrixXd& mode_params, const bool fix, const double fix_val){
	const int Nchars=5;
	const int precision=2;
	if (fix == 0){
		if (mode_params.col(4).mean() < 0.6){
			outfile << " Splitting_a1          Uniform         " << std::setw(Nchars) << std::setprecision(precision)  << mode_params.col(4).mean() << "          0.0          1.2" << std::endl;
		}
		if (mode_params.col(4).mean() < 2.0 && mode_params.col(4).mean() >= 0.6){
			outfile << " Splitting_a1          Uniform         " << std::setw(Nchars) << std::setprecision(precision)  << mode_params.col(4).mean() << "          0.0          5.0" << std::endl;
		}
		if (mode_params.col(4).mean() < 5.0 && mode_params.col(4).mean() >= 2.0){
			outfile << " Splitting_a1          Uniform         " << std::setw(Nchars) << std::setprecision(precision)  << mode_params.col(4).mean() << "          0.0          7.0" << std::endl;
		}
	} else{
		 outfile << " Splitting_a1          Fix             "  << std::setw(Nchars) << std::setprecision(precision)  << fix_val << std::endl;
	}
}


void write_a3_key(std::ofstream& outfile, const MatrixXd& mode_params, const bool fix, const double fix_val){
	const int Nchars=5;
	const int precision=2;
	const int minval=-0.2*mode_params.col(4).mean(); // 20% of a1 max
	const int maxval=0.2*mode_params.col(4).mean(); // 20% of a1 max
	if (fix == 0){
		outfile << " Splitting_a3          Uniform         " << std::setw(Nchars) << std::setprecision(precision)  <<  mode_params.col(4).mean()*0.01 << "      " << minval << "       " << maxval << std::endl;
	} else{
		outfile << " Splitting_a3          Fix             "  << std::setw(Nchars) << std::setprecision(precision)  << fix_val << std::endl;
	}
}
