/*
 * models_database.cpp
 *
 * Header file that contains all kind of methods
 * used to generate models for the pulsation/noise
 * 
 *  Created on: 20 Apr 2016
 *      Author: obenomar
 */
# include <iostream>
# include <iomanip>
#include <fstream>
# include <Eigen/Dense>
#include "models_database.h"
#include "noise_models.h"
#include "string_handler.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

double *r8vec_normal_01 ( int n, int *seed );

void asymptotic_mm_v1(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path){

	int seed=(unsigned)time(NULL);
	srand(seed);

	// ----- Constants ------
	const double PI = 3.141592653589793238462643;
	const double G=6.667e-8;
	const double Teff_sun=5777;
	const double Dnu_sun=135.1;
	const double numax_sun=3150;
	const double R_sun=6.96342e5;
	const double M_sun=1.98855e30;
	const double rho_sun=M_sun*1e3/(4*PI*pow(R_sun*1e5,3)/3);
	// ----------------------

	//std::cout << "input_params=" << input_params.transpose() << std::endl;
	//std::cout << "input_params.size()=" << input_params.size() << std::endl;
	// ------- Deploy the parameters ------
	double Teff=input_params[0];
	double Dnu=input_params[1];
	double epsilon=input_params[2];
	double alpha=input_params[3];
	double q=input_params[4];
	double hnr_l0=input_params[5];
	double l0_width_at_numax=input_params[6];
	
	double D0=Dnu/100;
	double lmax=3;
	double Nmax_pm=6; // Number of radial order to inject on each side of numax
	double N0=1.; 
	double Hmax_l0=hnr_l0*N0;
	
	double a1;
	
	MatrixXd noise_params(3,3);
	VectorXd input_noise(9);
	input_noise[0]=-1;
	input_noise[1]=-1;
	input_noise[2]=-1;
	input_noise[3]=-1;
	input_noise[4]=-1;
	input_noise[5]=-1; // CHECK THOSE -1 ...
	input_noise[6]=N0;
	input_noise[7]=-2;
	input_noise[8]=-2;
	// ------------------------------------

	// ---- Evaluation of DP ----
	// Super rought estimate derived by visual inspection of the Mosser+2015, Fig.1
	const double c=36.8222;
	const double b=2.63897;
	const double a=0.0168202;
	double DP;
	
	DP=a*pow(Dnu, 2) + b*Dnu + c;
	double *r = r8vec_normal_01 ( 1, &seed );
	DP=DP +  *r*DP*2.5/100.; // Inject a gaussian random error of 2.5%
	// -----------

	// ----------- PYTHON EXTERNAL FUNCTION -------------
	// a. Generate the configuration file for the python function
	int Nchars_spec = 20;
	int precision_spec = 5;
	int sizes;
	
	std::string line0;
	std::ofstream rwfile;
	
	std::cout << "                     - Attempting to write cfg on file " << file_cfg_mm << "..." << std::endl;
	rwfile.open(file_cfg_mm.c_str());
	if(rwfile.is_open()){
		// ---------------------
		rwfile << "# First line: Teff / Dnu / epsilon / D0. Second line: DP1 / alpha / q. Third line coupling / how many l=0 freq on left&right of numax / hmax / width at numax for l=0 / max uniform spread on numax (% or <=0 if off)" << std::endl;
		rwfile << Teff;
		rwfile << std::setw(Nchars_spec) << std::setprecision(precision_spec) << Dnu;
		rwfile << std::setw(Nchars_spec) << std::setprecision(precision_spec) << epsilon << std::setw(Nchars_spec) << D0 << std::endl;
		rwfile << DP << std::setw(Nchars_spec) << alpha << std::endl;
		rwfile << q << std::setw(Nchars_spec) << Nmax_pm <<  std::setw(Nchars_spec)  << Hmax_l0 <<  std::setw(Nchars_spec)  << l0_width_at_numax << std::setw(Nchars_spec) << "-1 " << std::endl;
		rwfile.close();
		std::cout << "Success... starting python3 external program..." << std::endl;
	} else{
		std::cout << "Error! Could not write the configuration file the Python external routine!" << std::endl;
		std::cout << "The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	// b. Call the external function
	const std::string str="python3 -c \"import bump_DP; bump_DP.main_star_generator(config_file='external/ARMM-solver/star_params.global', output_file='" + file_out_modes + "', output_file_range='external/ARMM-solver/star_params.range'" + ", output_file_rot='external/ARMM-solver/star_params.rot')\" ";
	const char *command = str.c_str(); 
	std::cout << "Executing command line: " << std::endl;
	std::cout << "    "  << str << std::endl;
	
	if (Dnu <= 15){
		std::cout << "    Model with small Dnu ==> many mixed modes. This might be long to find the solutions..." << std::endl; 
	}
	system(command);
	// Defining the noise profile parameters
	// Note about the noise: -1 means that it is ignored. -2 mean that the value is irrelevant
	for(int e=0; e<3; e++){
		for(int k=0; k<3; k++){
			noise_params(e, k)=input_noise(3*e + k);
		}
	}
	// A FUNCTION THAT WRITES THE Noise
	write_star_noise_params(noise_params, file_out_noise);

}


void asymptotic_mm_v2(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path){

	int seed=(unsigned)time(NULL);
	srand(seed);
	
	// ----- Constants ------
	const double PI = 3.141592653589793238462643;
	const double G=6.667e-8;
	const double Teff_sun=5777;
	const double Dnu_sun=135.1;
	const double numax_sun=3150;
	const double R_sun=6.96342e5;
	const double M_sun=1.98855e30;
	const double rho_sun=M_sun*1e3/(4*PI*pow(R_sun*1e5,3)/3);
	// ----------------------

	//std::cout << "input_params=" << input_params.transpose() << std::endl;
	//std::cout << "input_params.size()=" << input_params.size() << std::endl;
	// ------- Deploy the parameters ------
	double rot_env=input_params[0];
	double rot_ratio=input_params[1];	
	double Dnu=input_params[2];
	double epsilon=input_params[3];
	double alpha=input_params[4];
	double q=input_params[5];
	double hnr_l0=input_params[6];
	double l0_width_at_numax=input_params[7];
	
	double D0=Dnu/100;
	double lmax=3;
	double Nmax_pm=6; // Number of radial order to inject on each side of numax
	double N0=1.; 
	double Hmax_l0=hnr_l0*N0;
	
	double a1;
	
	MatrixXd noise_params(3,3);
	VectorXd input_noise(9);
	input_noise[0]=-1;
	input_noise[1]=-1;
	input_noise[2]=-1;
	input_noise[3]=-1;
	input_noise[4]=-1;
	input_noise[5]=-1; // CHECK THOSE -1 ...
	input_noise[6]=N0;
	input_noise[7]=-2;
	input_noise[8]=-2;
	// ------------------------------------

	// ---- Evaluation of DP ----
	// Super rought estimate derived by visual inspection of the Mosser+2015, Fig.1
	const double c=36.8222;
	const double b=2.63897;
	const double a=0.0168202;
	double DP;
	
	DP=a*pow(Dnu, 2) + b*Dnu + c;
	double *r = r8vec_normal_01 ( 1, &seed );
	DP=DP +  *r*DP*2.5/100.; // Inject a gaussian random error of 2.5%
	// -----------

	// ----------- PYTHON EXTERNAL FUNCTION -------------
	// a. Generate the configuration file for the python function
	int Nchars_spec = 20;
	int precision_spec = 5;
	int sizes;
	
	std::string line0;
	std::ofstream rwfile;
	
	std::cout << "                     - Attempting to write cfg on file " << file_cfg_mm << "..." << std::endl;
	rwfile.open(file_cfg_mm.c_str());
	if(rwfile.is_open()){
		// ---------------------
		rwfile << "# First line: rot_env / rot_ratio / Dnu / epsilon / D0. Second line: DP1 / alpha / q. Third line coupling / how many l=0 freq on left&right of numax / hmax / width at numax for l=0 / max uniform spread on numax (% or <=0 if off)" << std::endl;
		rwfile << rot_env << std::setw(Nchars_spec) << std::setprecision(precision_spec) << rot_ratio;
		rwfile << std::setw(Nchars_spec) << std::setprecision(precision_spec) << Dnu;
		rwfile << std::setw(Nchars_spec) << std::setprecision(precision_spec) << epsilon << std::setw(Nchars_spec) << D0 << std::endl;
		rwfile << DP << std::setw(Nchars_spec) << alpha << std::endl;
		rwfile << q << std::setw(Nchars_spec) << Nmax_pm <<  std::setw(Nchars_spec)  << Hmax_l0 <<  std::setw(Nchars_spec)  << l0_width_at_numax << std::setw(Nchars_spec) << "-1 " << std::endl;
		rwfile.close();
		std::cout << "Success... starting python3 external program..." << std::endl;
	} else{
		std::cout << "Error! Could not write the configuration file the Python external routine!" << std::endl;
		std::cout << "The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	// b. Call the external function
	const std::string str="python3 -c \"import bump_DP; bump_DP.main_star_generator(config_file='external/ARMM-solver/star_params.global', output_file='" + file_out_modes + "', output_file_range='external/ARMM-solver/star_params.range'" + ", output_file_rot='external/ARMM-solver/star_params.rot', version=2)\" ";
	const char *command = str.c_str(); 
	std::cout << "Executing command line: " << std::endl;
	std::cout << "    "  << str << std::endl;
	
	if (Dnu <= 15){
		std::cout << "    Model with small Dnu ==> many mixed modes. This might be long to find the solutions..." << std::endl; 
	}
	system(command);
	// Defining the noise profile parameters
	// Note about the noise: -1 means that it is ignored. -2 mean that the value is irrelevant
	for(int e=0; e<3; e++){
		for(int k=0; k<3; k++){
			noise_params(e, k)=input_noise(3*e + k);
		}
	}
	// A FUNCTION THAT WRITES THE Noise
	write_star_noise_params(noise_params, file_out_noise);

}

void asymptotic_mm_v3(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path){

	int seed=(unsigned)time(NULL);
	srand(seed);
	
	// ----- Constants ------
	const double PI = 3.141592653589793238462643;
	const double G=6.667e-8;
	const double Teff_sun=5777;
	const double Dnu_sun=135.1;
	const double numax_sun=3150;
	const double R_sun=6.96342e5;
	const double M_sun=1.98855e30;
	const double rho_sun=M_sun*1e3/(4*PI*pow(R_sun*1e5,3)/3);
	// ----------------------

	//std::cout << "input_params=" << input_params.transpose() << std::endl;
	//std::cout << "input_params.size()=" << input_params.size() << std::endl;
	// ------- Deploy the parameters ------
	double rot_env=input_params[0];
	double rot_core=input_params[1];	
	double Dnu=input_params[2];
	double epsilon=input_params[3];
	double alpha=input_params[4];
	double q=input_params[5];
	double hnr_l0=input_params[6];
	double l0_width_at_numax=input_params[7];
	
	double D0=Dnu/100;
	double lmax=3;
	double Nmax_pm=6; // Number of radial order to inject on each side of numax
	double N0=1.; 
	double Hmax_l0=hnr_l0*N0;
	
	double a1;
	
	MatrixXd noise_params(3,3);
	VectorXd input_noise(9);
	input_noise[0]=-1;
	input_noise[1]=-1;
	input_noise[2]=-1;
	input_noise[3]=-1;
	input_noise[4]=-1;
	input_noise[5]=-1; // CHECK THOSE -1 ...
	input_noise[6]=N0;
	input_noise[7]=-2;
	input_noise[8]=-2;
	// ------------------------------------

	// ---- Evaluation of DP ----
	// Super rought estimate derived by visual inspection of the Mosser+2015, Fig.1
	const double c=36.8222;
	const double b=2.63897;
	const double a=0.0168202;
	double DP;
	
	DP=a*pow(Dnu, 2) + b*Dnu + c;
	double *r = r8vec_normal_01 ( 1, &seed );
	DP=DP +  *r*DP*2.5/100.; // Inject a gaussian random error of 2.5%
	// -----------

	// ----------- PYTHON EXTERNAL FUNCTION -------------
	// a. Generate the configuration file for the python function
	int Nchars_spec = 20;
	int precision_spec = 5;
	int sizes;
	
	std::string line0;
	std::ofstream rwfile;
	
	std::cout << "                     - Attempting to write cfg on file " << file_cfg_mm << "..." << std::endl;
	rwfile.open(file_cfg_mm.c_str());
	if(rwfile.is_open()){
		// ---------------------
		rwfile << "# First line: rot_env / rot_core / Dnu / epsilon / D0. Second line: DP1 / alpha. Third line coupling q / how many l=0 freq on left&right of numax / hmax / width at numax for l=0 / max uniform spread on numax (% or <=0 if off)" << std::endl;
		rwfile << rot_env << std::setw(Nchars_spec) << std::setprecision(precision_spec) << rot_core;
		rwfile << std::setw(Nchars_spec) << std::setprecision(precision_spec) << Dnu;
		rwfile << std::setw(Nchars_spec) << std::setprecision(precision_spec) << epsilon << std::setw(Nchars_spec) << D0 << std::endl;
		rwfile << DP << std::setw(Nchars_spec) << alpha << std::endl;
		rwfile << q << std::setw(Nchars_spec) << Nmax_pm <<  std::setw(Nchars_spec)  << Hmax_l0 <<  std::setw(Nchars_spec)  << l0_width_at_numax << std::setw(Nchars_spec) << "-1 " << std::endl;
		rwfile.close();
		std::cout << "Success... starting python3 external program..." << std::endl;
	} else{
		std::cout << "Error! Could not write the configuration file the Python external routine!" << std::endl;
		std::cout << "The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	// b. Call the external function
	const std::string str="python3 -c \"import bump_DP; bump_DP.main_star_generator(config_file='external/ARMM-solver/star_params.global', output_file='" + file_out_modes + "', output_file_range='external/ARMM-solver/star_params.range'" + ", output_file_rot='external/ARMM-solver/star_params.rot', version=3)\" ";
	const char *command = str.c_str(); 
	std::cout << "Executing command line: " << std::endl;
	std::cout << "    "  << str << std::endl;
	
	if (Dnu <= 15){
		std::cout << "    Model with small Dnu ==> many mixed modes. This might be long to find the solutions..." << std::endl; 
	}
	system(command);
	// Defining the noise profile parameters
	// Note about the noise: -1 means that it is ignored. -2 mean that the value is irrelevant
	for(int e=0; e<3; e++){
		for(int k=0; k<3; k++){
			noise_params(e, k)=input_noise(3*e + k);
		}
	}
	// A FUNCTION THAT WRITES THE Noise
	write_star_noise_params(noise_params, file_out_noise);

}

void asymptotic_mm_freeDp_numaxspread_curvepmodes_v1(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path){

	int seed=(unsigned)time(NULL);
	srand(seed);

	// ----- Constants ------
	const double PI = 3.141592653589793238462643;
	const double G=6.667e-8;
	const double Teff_sun=5777;
	const double Dnu_sun=135.1;
	const double numax_sun=3150;
	const double R_sun=6.96342e5;
	const double M_sun=1.98855e30;
	const double rho_sun=M_sun*1e3/(4*PI*pow(R_sun*1e5,3)/3);
	// ----------------------

	//std::cout << "input_params=" << input_params.transpose() << std::endl;
	//std::cout << "input_params.size()=" << input_params.size() << std::endl;
	// ------- Deploy the parameters ------
	double Teff=input_params[0];
	double Dnu=input_params[1];
	double epsilon=input_params[2];
	double DP=input_params[3];
	double alpha=input_params[4];
	double q=input_params[5];
	double hnr_l0=input_params[6];
	double l0_width_at_numax=input_params[7];
	double numax_spread=input_params[8];	

	
	double D0=Dnu/100;
	double lmax=3;
	double Nmax_pm=6; // Number of radial order to inject on each side of numax
	double N0=1.; 
	double Hmax_l0=hnr_l0*N0;
	
	double a1;
	
	MatrixXd noise_params(3,3);
	VectorXd input_noise(9);
	input_noise[0]=-1;
	input_noise[1]=-1;
	input_noise[2]=-1;
	input_noise[3]=-1;
	input_noise[4]=-1;
	input_noise[5]=-1; // CHECK THOSE -1 ...
	input_noise[6]=N0;
	input_noise[7]=-2;
	input_noise[8]=-2;
	// ------------------------------------

	// ----------- PYTHON EXTERNAL FUNCTION -------------
	// a. Generate the configuration file for the python function
	int Nchars_spec = 20;
	int precision_spec = 5;
	int sizes;
	
	std::string line0;
	std::ofstream rwfile;
	
	std::cout << "                     - Attempting to write cfg on file " << file_cfg_mm << "..." << std::endl;
	rwfile.open(file_cfg_mm.c_str());
	if(rwfile.is_open()){
		// ---------------------
		rwfile << "# First line: Teff / Dnu / epsilon / D0. Second line: DP1 / alpha / q. Third line coupling / how many l=0 freq on left&right of numax / hmax / width at numax for l=0 / max uniform spread on numax (% or <=0 if off)" << std::endl;
		rwfile << Teff;
		rwfile << std::setw(Nchars_spec) << std::setprecision(precision_spec) << Dnu;
		rwfile << std::setw(Nchars_spec) << std::setprecision(precision_spec) << epsilon << std::setw(Nchars_spec) << D0 << std::endl;
		rwfile << DP << std::setw(Nchars_spec) << alpha << std::endl;
		rwfile << q << std::setw(Nchars_spec) << Nmax_pm <<  std::setw(Nchars_spec)  << Hmax_l0 <<  std::setw(Nchars_spec)  << l0_width_at_numax << std::setw(Nchars_spec) << numax_spread << std::endl;
		rwfile.close();
		std::cout << "Success... starting python3 external program..." << std::endl;
	} else{
		std::cout << "Error! Could not write the configuration file the Python external routine!" << std::endl;
		std::cout << "The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	// b. Call the external function
	const std::string str="python3 -c \"import bump_DP; bump_DP.main_star_generator(config_file='external/ARMM-solver/star_params.global', output_file='" + file_out_modes + "', output_file_range='external/ARMM-solver/star_params.range'" + ", output_file_rot='external/ARMM-solver/star_params.rot')\" ";
	const char *command = str.c_str(); 
	std::cout << "Executing command line: " << std::endl;
	std::cout << "    "  << str << std::endl;
	
	if (Dnu <= 15){
		std::cout << "    Model with small Dnu ==> many mixed modes. This might be long to find the solutions..." << std::endl; 
	}
	system(command);
	// Defining the noise profile parameters
	// Note about the noise: -1 means that it is ignored. -2 mean that the value is irrelevant
	for(int e=0; e<3; e++){
		for(int k=0; k<3; k++){
			noise_params(e, k)=input_noise(3*e + k);
		}
	}
	// A FUNCTION THAT WRITES THE Noise
	write_star_noise_params(noise_params, file_out_noise);

}

void asymptotic_mm_freeDp_numaxspread_curvepmodes_v2(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path){

	int seed=(unsigned)time(NULL);
	srand(seed);
	
	// ----- Constants ------
	const double PI = 3.141592653589793238462643;
	const double G=6.667e-8;
	const double Teff_sun=5777;
	const double Dnu_sun=135.1;
	const double numax_sun=3150;
	const double R_sun=6.96342e5;
	const double M_sun=1.98855e30;
	const double rho_sun=M_sun*1e3/(4*PI*pow(R_sun*1e5,3)/3);
	// ----------------------

	//std::cout << "input_params=" << input_params.transpose() << std::endl;
	//std::cout << "input_params.size()=" << input_params.size() << std::endl;
	// ------- Deploy the parameters ------
	double rot_env=input_params[0];
	double rot_ratio=input_params[1];	
	double Dnu=input_params[2];
	double epsilon=input_params[3];
	double DP=input_params[4];
	double alpha=input_params[5];
	double q=input_params[6];
	double hnr_l0=input_params[7];
	double l0_width_at_numax=input_params[8];
	double numax_spread=input_params[9];	
	
	double D0=Dnu/100;
	double lmax=3;
	double Nmax_pm=6; // Number of radial order to inject on each side of numax
	double N0=1.; 
	double Hmax_l0=hnr_l0*N0;
	
	double a1;
	
	MatrixXd noise_params(3,3);
	VectorXd input_noise(9);
	input_noise[0]=-1;
	input_noise[1]=-1;
	input_noise[2]=-1;
	input_noise[3]=-1;
	input_noise[4]=-1;
	input_noise[5]=-1; // CHECK THOSE -1 ...
	input_noise[6]=N0;
	input_noise[7]=-2;
	input_noise[8]=-2;
	// ------------------------------------

	// ----------- PYTHON EXTERNAL FUNCTION -------------
	// a. Generate the configuration file for the python function
	int Nchars_spec = 20;
	int precision_spec = 5;
	int sizes;
	
	std::string line0;
	std::ofstream rwfile;
	
	std::cout << "                     - Attempting to write cfg on file " << file_cfg_mm << "..." << std::endl;
	rwfile.open(file_cfg_mm.c_str());
	if(rwfile.is_open()){
		// ---------------------
		rwfile << "# First line: rot_env / rot_ratio / Dnu / epsilon / D0. Second line: DP1 / alpha / q. Third line coupling / how many l=0 freq on left&right of numax / hmax / width at numax for l=0 / max uniform spread on numax (% or <=0 if off)" << std::endl;
		rwfile << rot_env << std::setw(Nchars_spec) << std::setprecision(precision_spec) << rot_ratio;
		rwfile << std::setw(Nchars_spec) << std::setprecision(precision_spec) << Dnu;
		rwfile << std::setw(Nchars_spec) << std::setprecision(precision_spec) << epsilon << std::setw(Nchars_spec) << D0 << std::endl;
		rwfile << DP << std::setw(Nchars_spec) << alpha << std::endl;
		rwfile << q << std::setw(Nchars_spec) << Nmax_pm <<  std::setw(Nchars_spec)  << Hmax_l0;
		rwfile <<  std::setw(Nchars_spec)  << l0_width_at_numax <<  std::setw(Nchars_spec)  << numax_spread << std::endl;
		rwfile.close();
		std::cout << "Success... starting python3 external program..." << std::endl;
	} else{
		std::cout << "Error! Could not write the configuration file the Python external routine!" << std::endl;
		std::cout << "The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	// b. Call the external function
	const std::string str="python3 -c \"import bump_DP; bump_DP.main_star_generator(config_file='external/ARMM-solver/star_params.global', output_file='" + file_out_modes + "', output_file_range='external/ARMM-solver/star_params.range'" + ", output_file_rot='external/ARMM-solver/star_params.rot', version=2)\" ";
	const char *command = str.c_str(); 
	std::cout << "Executing command line: " << std::endl;
	std::cout << "    "  << str << std::endl;
	
	if (Dnu <= 15){
		std::cout << "    Model with small Dnu ==> many mixed modes. This might be long to find the solutions..." << std::endl; 
	}
	system(command);
	// Defining the noise profile parameters
	// Note about the noise: -1 means that it is ignored. -2 mean that the value is irrelevant
	for(int e=0; e<3; e++){
		for(int k=0; k<3; k++){
			noise_params(e, k)=input_noise(3*e + k);
		}
	}
	// A FUNCTION THAT WRITES THE Noise
	write_star_noise_params(noise_params, file_out_noise);

}

void asymptotic_mm_freeDp_numaxspread_curvepmodes_v3(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path){

	int seed=(unsigned)time(NULL);
	srand(seed);
	
	// ----- Constants ------
	const double PI = 3.141592653589793238462643;
	const double G=6.667e-8;
	const double Teff_sun=5777;
	const double Dnu_sun=135.1;
	const double numax_sun=3150;
	const double R_sun=6.96342e5;
	const double M_sun=1.98855e30;
	const double rho_sun=M_sun*1e3/(4*PI*pow(R_sun*1e5,3)/3);
	// ----------------------

	//std::cout << "input_params=" << input_params.transpose() << std::endl;
	//std::cout << "input_params.size()=" << input_params.size() << std::endl;
	// ------- Deploy the parameters ------
	double rot_env=input_params[0];
	double rot_core=input_params[1];	
	double Dnu=input_params[2];
	double epsilon=input_params[3];
	double DP=input_params[4];
	double alpha=input_params[5];
	double q=input_params[6];
	double hnr_l0=input_params[7];
	double l0_width_at_numax=input_params[8];
	double numax_spread=input_params[9];
		
	double D0=Dnu/100;
	double lmax=3;
	double Nmax_pm=6; // Number of radial order to inject on each side of numax
	double N0=1.; 
	double Hmax_l0=hnr_l0*N0;
	
	double a1;
	
	MatrixXd noise_params(3,3);
	VectorXd input_noise(9);
	input_noise[0]=-1;
	input_noise[1]=-1;
	input_noise[2]=-1;
	input_noise[3]=-1;
	input_noise[4]=-1;
	input_noise[5]=-1; // CHECK THOSE -1 ...
	input_noise[6]=N0;
	input_noise[7]=-2;
	input_noise[8]=-2;
	// ------------------------------------

	// ----------- PYTHON EXTERNAL FUNCTION -------------
	// a. Generate the configuration file for the python function
	int Nchars_spec = 20;
	int precision_spec = 5;
	int sizes;
	
	std::string line0;
	std::ofstream rwfile;
	
	std::cout << "                     - Attempting to write cfg on file " << file_cfg_mm << "..." << std::endl;
	rwfile.open(file_cfg_mm.c_str());
	if(rwfile.is_open()){
		// ---------------------
		rwfile << "# First line: rot_env / rot_core / Dnu / epsilon / D0. Second line: DP1 / alpha. Third line coupling q / how many l=0 freq on left&right of numax / hmax / width at numax for l=0 / max uniform spread on numax (% or <=0 if off)" << std::endl;
		rwfile << rot_env << std::setw(Nchars_spec) << std::setprecision(precision_spec) << rot_core;
		rwfile << std::setw(Nchars_spec) << std::setprecision(precision_spec) << Dnu;
		rwfile << std::setw(Nchars_spec) << std::setprecision(precision_spec) << epsilon << std::setw(Nchars_spec) << D0 << std::endl;
		rwfile << DP << std::setw(Nchars_spec) << alpha << std::endl;
		rwfile << q << std::setw(Nchars_spec) << Nmax_pm <<  std::setw(Nchars_spec)  << Hmax_l0 <<  std::setw(Nchars_spec)  << l0_width_at_numax << std::setw(Nchars_spec);
		rwfile <<  std::setw(Nchars_spec)  << numax_spread << std::endl;
		rwfile.close();
		std::cout << "Success... starting python3 external program..." << std::endl;
	} else{
		std::cout << "Error! Could not write the configuration file the Python external routine!" << std::endl;
		std::cout << "The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	// b. Call the external function
	const std::string str="python3 -c \"import bump_DP; bump_DP.main_star_generator(config_file='external/ARMM-solver/star_params.global', output_file='" + file_out_modes + "', output_file_range='external/ARMM-solver/star_params.range'" + ", output_file_rot='external/ARMM-solver/star_params.rot', version=3)\" ";
	const char *command = str.c_str(); 
	std::cout << "Executing command line: " << std::endl;
	std::cout << "    "  << str << std::endl;
	
	if (Dnu <= 15){
		std::cout << "    Model with small Dnu ==> many mixed modes. This might be long to find the solutions..." << std::endl; 
	}
	system(command);
	// Defining the noise profile parameters
	// Note about the noise: -1 means that it is ignored. -2 mean that the value is irrelevant
	for(int e=0; e<3; e++){
		for(int k=0; k<3; k++){
			noise_params(e, k)=input_noise(3*e + k);
		}
	}
	// A FUNCTION THAT WRITES THE Noise
	write_star_noise_params(noise_params, file_out_noise);

}





// ------------------------------------
// ------- MAIN SEQUENCE MODELS ------
// ------------------------------------

void generate_cfg_asymptotic_act_asym_Hgauss(VectorXd input_params, std::string file_out_modes, std::string file_out_noise){

	// ----- Constants ------
	const double PI = 3.141592653589793238462643;
	const double G=6.667e-8;
	const double Teff_sun=5777;
	const double Dnu_sun=135.1;
	const double numax_sun=3150;
	const double R_sun=6.96342e5;
	const double M_sun=1.98855e30;
	const double rho_sun=M_sun*1e3/(4*PI*pow(R_sun*1e5,3)/3);

    std::cout << "      generate_cfg_asymptotic_act_asym_Hgauss" << std::endl;

	VectorXd Visibilities(4);
	Visibilities << 1, 1.5, 0.5, 0.08;
	// ----------------------

	double ks=2; // controls the width of the gaussian for heights

	//std::cout << "input_params=" << input_params.transpose() << std::endl;
	//std::cout << "input_params.size()=" << input_params.size() << std::endl;
	// ------- Deploy the parameters ------
	double numax=input_params[0];
	double Dnu=input_params[1];
	double epsilon=input_params[2];
	double D0=input_params[3];
	double Max_Height=input_params[4];
	double Width=input_params[5];
	double lmax=input_params[6];
	int    Nmax=input_params[7];
	double a1=input_params[8];
	double a3=input_params[9];
	double b=input_params[10];
	double alfa=input_params[11];
	double beta_asym=input_params[12];
	double inc=input_params[13];	
	VectorXd input_noise(9);
	if((input_params.size() - 14) == 9){
		input_noise=input_params.segment(14, 9);
	} else{ // If we did not provide 3 lines of 3 elements, then this means that we have 2 lines of 3 elements + 1 line with one elements (7 parameters)
		input_noise.segment(0, 7)=input_params.segment(14, 7);
		input_noise[7]=-2;
		input_noise[8]=-2;
	}
	// ------------------------------------

    //std::cout << "           - Variables..." << std::endl;

	// --------- Variables ---------
	int k;
	double el, en, n_at_numax, height, eta;
	VectorXd en_list(Nmax);
	MatrixXd nu(int(lmax+1), Nmax), h(int(lmax+1), Nmax), w(int(lmax+1), Nmax), s_a1(int(lmax+1), Nmax), s_eta(int(lmax+1), Nmax), 
		 s_a3(int(lmax+1), Nmax), s_asym(int(lmax+1), Nmax), s_b(int(lmax+1), Nmax), s_alfa(int(lmax+1), Nmax), i(int(lmax+1), Nmax), 
		 mode_params(int(lmax+1)*Nmax, 11), noise_params(3,3);
	// -----------------------------

	if((Nmax % 2) !=0){
		std::cout << "Nmax must be a odd number. Please change that value accordingly" << std::endl;
		std::cout << "The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

	k=0;
	n_at_numax=numax/Dnu - epsilon + el/2 + el*(el+1)*D0/Dnu;
	for(en=-Nmax/2 + 1; en<=Nmax/2; en++){
		en_list[k]=floor(n_at_numax) + en;
		k=k+1;
	}	

	// Define the centrifugal force effect	
	eta=(4./3)*PI * pow( a1*1e-6 ,2) / (G * rho_sun) * pow(Dnu_sun/Dnu,2);
	//std::cout << "Fixing centrifugal term eta = " << eta << std::endl;
	
	// Create a list of frequencies, Height, Width, Splitting, Centrifugal terms, latitudinal terms and stellar inclination

    //std::cout << "           - Parameters..." << std::endl;
    
	for(el=0; el<=lmax; el++){
		for(en=0; en<Nmax; en++){
			if(el == 0 || el == 1){
				nu(el,en)=(en_list[en] + epsilon + el/2)*Dnu - el*(el + 1)*D0;
			}
			if(el == 2 || el == 3){
				nu(el,en)=(en_list[en]-1 + epsilon + el/2)*Dnu - el*(el + 1)*D0;
			}
			if(el > 3){
				std::cout << "Generating frequencies with degree higher than l=3 is not implemented" << std::endl;
				std::cout << "Please set lmax<4" << std::endl;
				std::cout << "The program will exit now" << std::endl;
				exit(EXIT_FAILURE);
			}
			
			height=Max_Height*exp(-0.5 * pow( (nu(el, en) - numax)/(ks*Dnu), 2));
			h(el, en)=height*Visibilities[el];

			w(el,en)=Width;

			s_a1(el,en)=a1;
			s_eta(el,en)=eta;
			s_a3(el,en)=a3;
			s_b(el,en)=b;
			s_alfa(el,en)=alfa;
			s_asym(el,en)=beta_asym;
			i(el,en)=inc;
			/*			
			std::cout << "en_list[en]=" << en_list[en] << std::endl;
			std::cout << "nu(el,en)=" << nu(el,en) << std::endl;
			std::cout << "s_a1(el,en)=" << s_a1(el,en) << std::endl;
			std::cout << "s_eta(el,en)=" << s_eta(el,en) << std::endl;
			std::cout << "s_a3(el,en)=" << s_a3(el,en) << std::endl;
			std::cout << "s_b(el,en)=" << s_b(el,en) << std::endl;
			std::cout << "s_alfa(el,en)=" << s_alfa(el,en) << std::endl;
			std::cout << "s_asym(el,en)=" << s_asym(el,en) << std::endl;
			std::cout << "i(el,en)=" << i(el,en) << std::endl;
			std::cout << "--------------" << std::endl;
			*/
		}
	}

    //std::cout << "           - mode_params preparation..." << std::endl;

	// Summarizing the information into a parameter table
	for(el=0; el<=lmax; el++){
		for(k=el*Nmax; k<(el+1)*Nmax; k++){
			mode_params(k,0)=el;
			mode_params(k,1)=nu(el , k-el*Nmax);
			mode_params(k,2)=h(el , k-el*Nmax);
			mode_params(k,3)=w(el , k-el*Nmax);
			mode_params(k,4)=s_a1(el , k-el*Nmax);
			mode_params(k,5)=s_eta(el , k-el*Nmax);
			mode_params(k,6)=s_a3(el , k-el*Nmax);
			mode_params(k,7)=s_b(el , k-el*Nmax);
			mode_params(k,8)=s_alfa(el , k-el*Nmax);
			mode_params(k,9)=s_asym(el , k-el*Nmax);
			mode_params(k,10)=i(el , k-el*Nmax);
		}
	}
	

	// A FUNCTION THAT WRITES THE PARAMETERS
	write_star_mode_params_act_asym(mode_params, file_out_modes);

    //std::cout << "           - Noise..." << std::endl;

	// Defining the noise profile parameters
	// Note about the noise: -1 means that it is ignored. -2 mean that the value is irrelevant
	for(int e=0; e<3; e++){
		for(int k=0; k<3; k++){
			noise_params(e, k)=input_noise(3*e + k);
		}
	}
	// A FUNCTION THAT WRITES THE Noise
	write_star_noise_params(noise_params, file_out_noise);

    //std::cout << "           - Exit" << std::endl;


}

/* same as generate_cfg_asymptotic but Height follow a gaussian.
 The Height input is therefore the maximum height at numax
 The Width is rescaled from the relation given by Appourchaux et al. 2014 and given for the Sun 
*/
void  generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string extra){

	std::cout << "generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma UNDER CONSTRUCTION" << std::endl;
	exit(EXIT_SUCCESS);

	int lmax_ref Nmax_synthese, lmax, pos_max;
	double HNR, a1_ov_Gamma, a3, beta_asym, inc, HNRref, Height_factor, Gamma_at_numax, a1; 
	VectorXd tmp;
	VectorXd nu_star, h_star, gamma_star, s_a1_star, s_a3_star, s_asym_star, inc_star;
	star_params ref_star;
		
	ref_star=read_star_params(extra); // THIS FUNCTION HAS YET TO BE WRITTEN
	lmax_ref=ref_star.mode_params.col(0).maxCoeff(); // The first column contains the degrees

	// Defining the noise profile parameters
	// Note about the noise: -1 means that it is ignored. -2 mean that the value is irrelevant
	tmp.resize(ref_star.mode_params.col(0); 	
	tmp.setConstant(0);
	local_noise=harvey_like(noise_params, ref_star.mode_params.col(1), tmp); // Generate a list of local noise values for each frequencies
	
// ---- Deploy the parameters -----
	HNR=input_params[0];
	a1_ov_Gamma=input_params[1];
	lmax=input_params[2];
	a3=input_params[3];
	beta_asym=input_params[4];
	inc=input_params[5];
// ---------------------------------

	HNRref=(ref_star.mode_params.col(2)/local_noise).maxCoef(); // This is the maximum HNR of the reference data
	Height_factor=HNR/HNRref;  // compute the coeficient required to get the proper max(HNR)

	pos_max=where_dbl(ref_star.mode_params.col(2), ref_star.mode_params.col(2).maxCoef(), 0.001);
	Gamma_at_numax=ref_star.mode_params(3, pos_max);
		
	a1=a1_ov_Gamma*Gamma_at_numax; // We do no vary the Width. We change the splitting in order to get the wished a1/Gamma0

	if (local_noise.sum()!= 0){
		N0=1; // Imposing the the Noise background is 1
		std::cout <<  "Using N0=" << N0 << " (white noise)" << std::endl;
		std::cout << "HNR of all modes (degree / freq / HNR):" << std::endl;
		for(i =0; i<ref_star.mode_params.rows(); i++){
			std::cout << ref_star.mode_params.col(0) << "  " << ref_star.mode_params.col(1) << "  " << Height_factor * ref_star.mode_params.col(2)/local_noise << std::endl;
		}
	} else{
		std::cout << 'Warning: bruit_local from the stat_synthese file is 0 ==> Cannot compute N0=mean(local_noisr)' << std::endl;
		std::cout << '         The program will stop now' << std::endl;
		exit(EXIT_FAILURE);
	}
	
	// Defining the final size for all of the outptus
	h_star.resize(ref_star.mode_params.rows());
	gamma_star.resize(ref_star.mode_params.rows());
	s_a1_star.resize(ref_star.mode_params.rows());
	s_a3_star.resize(ref_star.mode_params.rows());	
	s_asym_star.resize(ref_star.mode_params.rows());
	inc_star.resize(ref_star.mode_params.rows());

	// Refactoring the heights
	h_star=Height_factor * ref_star.mode_params.col(2)* N0/local_noise; 

	gamma_star=ref_star.mode_params.col(3); // In IDL, AN INTERPOLATION WAS DONE FOR l>0. HERE WE ASSUME THE .in file is whatever the true model should be (no interpolation)
	
	if (a1 > 0){
		s_a1_star.setConstant(a1);
	} else{
		s_a1_star=ref_star.mode_params.col(4); 
		}
	}
	
	if (a3 > 0){
		s_a3_star.setConstant(a3);	
	} else{
		s_a3_star=ref_star.mode_params.col(6);	
	}
	if (inc > 0){ 
		inc_star.setConstant(inc);
	} else{
		inc_star=ref_star.mode_params.col(10);
	}
	
	
	// CONTINUE FROM HERE....
	std::cout << "Implemented until here..."<< std::endl;
	exit(EXIT_SUCCESS);
	mode_params=dblarr((lmax+1)*Nmax, 11)

	for el=0, lmax do begin
			mode_params[el*Nmax:(el+1)*Nmax-1, 0]=el
			mode_params[el*Nmax:(el+1)*Nmax-1, 1]=nu[el,*]
			mode_params[el*Nmax:(el+1)*Nmax-1, 2]=
			mode_params[el*Nmax:(el+1)*Nmax-1, 3]=w[el,*]  
			mode_params[el*Nmax:(el+1)*Nmax-1, 4]=s_a1[el,*] 
			mode_params[el*Nmax:(el+1)*Nmax-1, 5]=0  
			mode_params[el*Nmax:(el+1)*Nmax-1, 6]=s_a3[el,*]
			mode_params[el*Nmax:(el+1)*Nmax-1, 7]=0
			mode_params[el*Nmax:(el+1)*Nmax-1, 8]=0
			mode_params[el*Nmax:(el+1)*Nmax-1, 9]=beta_asym
			mode_params[el*Nmax:(el+1)*Nmax-1, 10]=i[el,*]  
	endfor


	write_star_mode_params_act_asym, mode_params, file_out_modes
	

	noise_params=dblarr(3, 3) ; [ [H0, tau0, p0], [H1, tau1, p1], N0]
	noise_params[0, *]=-1
	noise_params[1, *]=-1
	noise_params[2, 0]=N0 ; We put only a white noise
	noise_params[2, 1:2]=-2
	write_star_noise_params, noise_params, file_out_noise
		
}


/*
// small test program
int main(){

	VectorXd input_params;
	std::string file_out_modes;
	std::string file_out_noise;

	file_out_modes="file_params.txt";
	file_out_noise="file_noise.txt";

	input_params.setZero(14 +9);
	input_params << 2000., 90., 0.2, 1.35, 10., 5., 3, 16, 1.0, 0.025, -0.05, 3., 50., 90, -1, -1, -1, -1, -1, -1, 0.1, -2, -2;  
	//generate_cfg_asymptotic_Hgauss(input_params, file_out_modes, file_out_noise);
	generate_cfg_asymptotic_act_asym_Hgauss(input_params, file_out_modes, file_out_noise);
	std::cout << "End" << std::endl;
}
*/
