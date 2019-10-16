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
		rwfile << "# First line: Teff / Dnu / epsilon / D0. Second line: DP1 / alpha / q. Third line coupling / how many l=0 freq on left&right of numax / hmax / width at numax for l=0" << std::endl;
		rwfile << Teff;
		rwfile << std::setw(Nchars_spec) << std::setprecision(precision_spec) << Dnu;
		rwfile << std::setw(Nchars_spec) << std::setprecision(precision_spec) << epsilon << std::setw(Nchars_spec) << D0 << std::endl;
		rwfile << DP << std::setw(Nchars_spec) << alpha << std::endl;
		rwfile << q << std::setw(Nchars_spec) << Nmax_pm <<  std::setw(Nchars_spec)  << Hmax_l0 <<  std::setw(Nchars_spec)  << l0_width_at_numax << std::endl;
		rwfile.close();
		std::cout << "Success... starting python3 external program..." << std::endl;
	} else{
		std::cout << "Error! Could not write the configuration file the Python external routine!" << std::endl;
		std::cout << "The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	// b. Call the external function
	//const std::string str="python3 -c \"import bump_DP; bump_DP.main_star_generator(config_file='external/ARMM-solver/star_params.global', output_file='external/ARMM-solver/star_params.modes')\" ";
	//const std::string str="python3 -c \"import bump_DP; bump_DP.main_star_generator(config_file='external/ARMM-solver/star_params.global', output_file='" + file_out_modes + "')\" ";
	const std::string str="python3 -c \"import bump_DP; bump_DP.main_star_generator(config_file='external/ARMM-solver/star_params.global', output_file='" + file_out_modes + "', output_file_range='external/ARMM-solver/star_params.range')\" ";
	const char *command = str.c_str(); 
	std::cout << "Executing command line: " << std::endl;
	std::cout << "    "  << str << std::endl;
	
	if (Dnu <= 15){
		std::cout << "    Model with small Dnu ==> many mixed modes. This might be long to find the solutions..." << std::endl; 
	}
	system(command);
	//std::cout << external_path << std::endl;
	//exit(EXIT_SUCCESS);
		

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



/*
void generate_cfg_asymptotic_Hgauss(VectorXd input_params, std::string file_out_modes, std::string file_out_noise){

	// ----- Constants ------
	const double PI = 3.141592653589793238462643;
	const double G=6.667e-8;
	const double Teff_sun=5777;
	const double Dnu_sun=135.1;
	const double numax_sun=3150;
	const double R_sun=6.96342e5;
	const double M_sun=1.98855e30;
	const double rho_sun=M_sun*1e3/(4*PI*pow(R_sun*1e5,3)/3);

	VectorXd Visibilities(4);
	Visibilities << 1, 1.5, 0.5, 0.07;
	// ----------------------

	double ks=2; // controls the width of the gaussian for heights

	// ------- Deploy the parameters ------
	double numax=input_params[0];
	double Dnu=input_params[1];
	double epsilon=input_params[2];
	double D0=input_params[3];
	double Max_Height=input_params[4];
	double Width=input_params[5];
	double    lmax=input_params[6];
	int    Nmax=input_params[7];
	double a1=input_params[8];
	double coef_for_a2=input_params[9];
	double a3=input_params[10];
	double inc=input_params[11];
	const VectorXd input_noise=input_params.segment(12, 3*3);
	//double N0=input_params[12];
	// ------------------------------------


	// --------- Variables ---------
	int k;
	double el, en, n_at_numax, height, a2;
	VectorXd en_list(Nmax);
	MatrixXd nu(int(lmax+1), Nmax), h(int(lmax+1), Nmax), w(int(lmax+1), Nmax), s_a1(int(lmax+1), Nmax), s_a2(int(lmax+1), Nmax), 
		 s_a3(int(lmax+1), Nmax), i(int(lmax+1), Nmax), mode_params(int(lmax+1)*Nmax, 8), noise_params(3,3);
	// -----------------------------

	if((Nmax % 2) !=0){
		std::cout << "Nmax must be a odd number. Please change that value accordingly" << std::endl;
		std::cout << "The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

	n_at_numax=numax/Dnu - epsilon + el/2 + el*(el+1)*D0/Dnu;
	for(en=-Nmax/2 + 1; en<=Nmax/2; en++){
		en_list[k]=floor(n_at_numax) + en;
		k=k+1;
	}	


	// Define the centrifugal force effect	
	a2=coef_for_a2 * (4./3)*PI * pow( a1*1e-6 ,2) / (G * rho_sun) * pow(Dnu_sun/Dnu,2);
	std::cout << "Using coef_for_a2 = " << coef_for_a2 << std::endl;
	std::cout << "   ===> Resulting coeficient a2 = " << a2 << std::endl;
	
	// Create a list of frequencies, Height, Width, Splitting, Centrifugal terms, latitudinal terms and stellar inclination
	for(el=0; el<=lmax; el++){
		for(en=0; en<Nmax; en++){
			if(el == 0 || el == 1){
				nu(el,en)=(en_list[k] + epsilon + el/2)*Dnu - el*(el + 1)*D0;
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
			s_a2(el,en)=a2;
			s_a3(el,en)=a3;
			i(el,en)=inc;
		}
	}

	// Summarizing the information into a parameter table
	for(el=0; el<=lmax; el++){
		for(k=el*Nmax; k<(el+1)*Nmax; k++){
			mode_params(k,0)=el;
			mode_params(k,1)=nu(el , k-el*Nmax);
			mode_params(k,2)=h(el , k-el*Nmax);
			mode_params(k,3)=w(el , k-el*Nmax);
			mode_params(k,4)=s_a1(el , k-el*Nmax);
			mode_params(k,5)=s_a2(el , k-el*Nmax);
			mode_params(k,6)=s_a3(el , k-el*Nmax);
			mode_params(k,7)=i(el , k-el*Nmax);
		}
	}
	// A FUNCTION THAT WRITES THE PARAMETERS
	// ......

	// Defining the noise profile parameters
	// Note about the noise: -1 means that it is ignored. -2 mean that the value is irrelevant
	for(int e=0; e<3; e++){
		for(int k=0; k<3; k++){
			noise_params(e, k)=input_noise(3*e + k);
		}
	}
	// A FUNCTION THAT WRITES THE Noise
	// ......


	std::cout << "WARNING: THIS FUNCTION IS NOT COMPLETELY WRITTEN" << std::endl;
	std::cout << "The program will stop now" << std::endl;
	exit(EXIT_FAILURE);


}
*/

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
