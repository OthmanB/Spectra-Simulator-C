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
# include <string>
# include <Eigen/Dense>
#include <random>
#include "models_database.h"
#include "noise_models.h"
//#include "string_handler.h" // Replaced by ioproc.h on 17/06/2021
#include "ioproc.h"
#include "external/ARMM/bump_DP.h"
#include "linspace.h"
#include "linfit.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

double *r8vec_normal_01 ( int n, int *seed );

void asymptotic_mm_v1(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path, std::string template_file){

	int seed=(unsigned)time(NULL);
	srand(seed);

	std::random_device rd;
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<double> distrib(0 , 1);

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
	const std::string cpath=getcwd(NULL, 0);
	
	const std::string file_range=cpath + "/external/ARMM-solver/star_params.range";
	const double lmax=3;
	const double Nmax_pm=6; // Number of radial order to inject on each side of numax

	// ------- Deploy the parameters ------
	//double Teff=input_params[0];
	//double Dnu=input_params[1];
	//double epsilon=input_params[2];
	//double delta0l_percent=input_params[3];
	//double beta_p=input_params[4];
	double numax_spread=0;
	double nmax_spread=input_params[5];
	double DP_var_percent=input_params[6];
	//double alpha=input_params[7];
	//double q=input_params[8];
	//double hnr_l0=input_params[9];
	//double l0_width_at_numax=input_params[10];
	//double Vl1=input_params[11];
	//double Vl2=input_params[12];
	//double Vl3=input_params[13];
	//double H0_spread=input_params[14];
	double inc_rad, inc_star, inc_y, inc_c;
	double xmin, xmax, a1;
	double H, tau, p, N0;
	MatrixXd mode_params, noise_params(3,3);
	
	Cfg_synthetic_star cfg_star;
	Params_synthetic_star params;
	// ------------------------------------

	// ---- Evaluation of DP ----
	// Super rought estimate derived by visual inspection of the Mosser+2015, Fig.1
	const double c=36.8222;
	const double b=2.63897;
	const double a=0.0168202;
	double DP;
	
	cfg_star.Teff_star=input_params[0];
	cfg_star.Dnu_star=input_params[1];
	cfg_star.epsilon_star=input_params[2];
	cfg_star.delta0l_percent_star=input_params[3];
	cfg_star.beta_p_star=input_params[4];
	
	DP=a*pow(cfg_star.Dnu_star, 2) + b*cfg_star.Dnu_star + c;
	double *r = r8vec_normal_01 ( 1, &seed );
	DP=DP +  *r*DP*DP_var_percent/100.; // Inject a gaussian random error of DP_var_percent
	// -----------
	cfg_star.DPl_star=DP;                
	cfg_star.alpha_g_star=input_params[7];
	cfg_star.q_star=input_params[8];
	cfg_star.maxHNR_l0=input_params[9];
	cfg_star.H0_spread=input_params[14];
	cfg_star.Gamma_max_l0=input_params[10];
	cfg_star.Hfactor=input_params[23];
	cfg_star.Wfactor=input_params[24];
	cfg_star.rot_env_input=-1;
	cfg_star.rot_ratio_input=-1;
	cfg_star.rot_core_input=-1;
	cfg_star.noise_params_harvey_like.resize(8);
	cfg_star.noise_params_harvey_like <<  input_params[15], input_params[16] , input_params[17] , input_params[18] , input_params[19] , input_params[20] , input_params[21]  , input_params[22];    //[A_Pgran ,  B_Pgran , C_Pgran   ,  A_taugran ,  B_taugran  , C_taugran    , p      N0]
	cfg_star.numax_star=numax_from_stello2009(cfg_star.Dnu_star, numax_spread); // Second argument is the random spread on numax
	cfg_star.fmin=cfg_star.numax_star -Nmax_pm*cfg_star.Dnu_star;
	cfg_star.fmax=cfg_star.numax_star +(Nmax_pm+2)*cfg_star.Dnu_star;
	cfg_star.output_file_rot=cpath + "/external/ARMM-solver/star_params.rot";
	cfg_star.Vl.resize(3);
	cfg_star.Vl << input_params[11], input_params[12], input_params[13];
	cfg_star.filetemplate = template_file;
	cfg_star.nmax_star=cfg_star.numax_star/cfg_star.Dnu_star - cfg_star.epsilon_star;
	cfg_star.alpha_p_star=cfg_star.beta_p_star/cfg_star.nmax_star;
	cfg_star.sigma_m=0;
	cfg_star.sigma_p=0;

	if (std::abs(nmax_spread) > 0)
	{
		try
		{
			xmin=cfg_star.nmax_star*(1. - std::abs(nmax_spread)/100.);
			xmax=cfg_star.nmax_star*(1. + std::abs(nmax_spread)/100.);
			cfg_star.nmax_star=xmin + (xmax-xmin)*distrib(gen);
			
		}
		catch(...)
		{
			std::cout << "Error debug info:" << std::endl;
			std::cout << "cfg_star.nmax: " << cfg_star.nmax_star << std::endl;
			std::cout << "nmax_spread: " << nmax_spread << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	// Determination of an isotropic inclination
	inc_rad=distrib(gen)*PI/2;
	inc_y=distrib(gen);
	inc_c=std::cos(inc_rad);
	while (inc_y <=inc_c){
		inc_rad=distrib(gen)*PI/2;
		inc_y=distrib(gen);
		inc_c=std::cos(inc_rad);
	}
	inc_star=inc_rad*180./PI;
	// b. Generate the mode profiles and frequencies
	params=make_synthetic_asymptotic_star(cfg_star);
	mode_params=bumpoutputs_2_MatrixXd(params, inc_star); // get the output in a format that can be written with the writting function

	//el, nu, h, w, a1, eta, a3, b, alfa, asym, inc
	write_star_mode_params_act_asym(mode_params, file_out_modes);
	write_range_modes(cfg_star, params, file_range);

	if (cfg_star.Dnu_star <= 15){
		std::cout << "    Model with small Dnu ==> many mixed modes. This might be long to find the solutions..." << std::endl; 
	}
	
	tau=input_params[18] * pow(cfg_star.numax_star*1e-6,input_params[19]) + input_params[20]; // Granulation timescale (in seconds)
	H=input_params[15] * pow(cfg_star.numax_star*1e-6,input_params[16]) + input_params[17]; // Granulation Amplitude
	H=H/tau ; //This is due to the used definition for the Harvey profile (conversion from Hz to microHz)
	tau=tau/1000. ; //conversion in ksec
	p=input_params[21];// power law:  MUST BE CLOSE TO 2
	N0=input_params[22];
	noise_params(0,0)=-1;
	noise_params(0,1)=-1;
	noise_params(0,2)=-1; 
	noise_params(1,0)=H;
	noise_params(1,1)=tau;
	noise_params(1,2)=p; 
	noise_params(2, 0)=N0; // White noise
	noise_params(2, 1)=-2;
	noise_params(2, 2)=-2;

	// A FUNCTION THAT WRITES THE Noise
	write_star_noise_params(noise_params, file_out_noise);
}


void asymptotic_mm_v2(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path, std::string template_file){

	int seed=(unsigned)time(NULL);
	srand(seed);

	std::random_device rd;
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<double> distrib(0 , 1);
	
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
	const std::string cpath=getcwd(NULL, 0);
	
	const std::string file_range=cpath + "/external/ARMM-solver/star_params.range";
	const double lmax=3;
	const double Nmax_pm=6; // Number of radial order to inject on each side of numax

	//std::cout << "input_params=" << input_params.transpose() << std::endl;
	//std::cout << "input_params.size()=" << input_params.size() << std::endl;
	// ------- Deploy the parameters ------
	//double rot_env=input_params[0];
	//double rot_ratio=input_params[1];	
	//double Dnu=input_params[2];
	//double epsilon=input_params[3];
	//double delta0l_percent=input_params[4];
	//double beta_p=input_params[5];
	double numax_spread=0;
	double nmax_spread=input_params[6];
	double DP_var_percent=input_params[7];
	//double alpha=input_params[8];
	//double q=input_params[9];
	//double hnr_l0=input_params[10];
	//double l0_width_at_numax=input_params[11];
	//double Vl1=input_params[12];
	//double Vl2=input_params[13];
	//double Vl3=input_params[14];
	//double H0_spread=input_params[15];
	
	double inc_rad, inc_star, inc_y, inc_c;
	double xmin, xmax, a1;
	double H, tau, p, N0;
	MatrixXd mode_params, noise_params(3,3);
	
	Cfg_synthetic_star cfg_star;
	Params_synthetic_star params;
	// ------------------------------------

	// ---- Evaluation of DP ----
	// Super rought estimate derived by visual inspection of the Mosser+2015, Fig.1
	const double c=36.8222;
	const double b=2.63897;
	const double a=0.0168202;
	double DP;
	
		// -----------
	cfg_star.Teff_star=-1;
	cfg_star.Dnu_star=input_params[2];
	cfg_star.epsilon_star=input_params[3];
	cfg_star.delta0l_percent_star=input_params[4];
	cfg_star.beta_p_star=input_params[5];
	
	DP=a*pow(cfg_star.Dnu_star, 2) + b*cfg_star.Dnu_star + c;
	double *r = r8vec_normal_01 ( 1, &seed );
	DP=DP +  *r*DP*DP_var_percent/100.; // Inject a gaussian random error of DP_var_percent
	// -----------
	cfg_star.DPl_star=DP;                
	cfg_star.alpha_g_star=input_params[8];
	cfg_star.q_star=input_params[9];
	cfg_star.maxHNR_l0=input_params[10];
	cfg_star.H0_spread=input_params[15];
	cfg_star.Gamma_max_l0=input_params[11];
	cfg_star.Hfactor=input_params[24];
	cfg_star.Wfactor=input_params[25];
	cfg_star.rot_env_input=input_params[0];
	cfg_star.rot_ratio_input=input_params[1];
	cfg_star.rot_core_input=-1;
	cfg_star.noise_params_harvey_like.resize(8);
	cfg_star.noise_params_harvey_like <<  input_params[16], input_params[17] , input_params[18] , input_params[19] , input_params[20] , input_params[21] , input_params[22]  , input_params[23];    //[A_Pgran ,  B_Pgran , C_Pgran   ,  A_taugran ,  B_taugran  , C_taugran    , p      N0]
	cfg_star.numax_star=numax_from_stello2009(cfg_star.Dnu_star, numax_spread); // Second argument is the random spread on numax
	cfg_star.fmin=cfg_star.numax_star -Nmax_pm*cfg_star.Dnu_star;
	cfg_star.fmax=cfg_star.numax_star +(Nmax_pm+2)*cfg_star.Dnu_star;
	cfg_star.output_file_rot=cpath + "/external/ARMM-solver/star_params.rot";
	cfg_star.Vl.resize(3);
	cfg_star.Vl << input_params[12], input_params[13], input_params[14];
	cfg_star.filetemplate = template_file;
	cfg_star.nmax_star=cfg_star.numax_star/cfg_star.Dnu_star - cfg_star.epsilon_star;
	cfg_star.alpha_p_star=cfg_star.beta_p_star/cfg_star.nmax_star;
	cfg_star.sigma_m=0;
	cfg_star.sigma_p=0;

	if (std::abs(nmax_spread) > 0)
	{
		try
		{
			xmin=cfg_star.nmax_star*(1. - std::abs(nmax_spread)/100.);
			xmax=cfg_star.nmax_star*(1. + std::abs(nmax_spread)/100.);
			cfg_star.nmax_star=xmin + (xmax-xmin)*distrib(gen);
			
		}
		catch(...)
		{
			std::cout << "Error debug info:" << std::endl;
			std::cout << "cfg_star.nmax: " << cfg_star.nmax_star << std::endl;
			std::cout << "nmax_spread: " << nmax_spread << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	// Determination of an isotropic inclination
	inc_rad=distrib(gen)*PI/2;
	inc_y=distrib(gen);
	inc_c=std::cos(inc_rad);
	while (inc_y <=inc_c){
		inc_rad=distrib(gen)*PI/2;
		inc_y=distrib(gen);
		inc_c=std::cos(inc_rad);
	}
	inc_star=inc_rad*180./PI;
	
	// b. Generate the mode profiles and frequencies
	params=make_synthetic_asymptotic_star(cfg_star);
	mode_params=bumpoutputs_2_MatrixXd(params, inc_star); // get the output in a format that can be written with the writting function

	//el, nu, h, w, a1, eta, a3, b, alfa, asym, inc
	write_star_mode_params_act_asym(mode_params, file_out_modes);
	write_range_modes(cfg_star, params, file_range);

	if (cfg_star.Dnu_star <= 15){
		std::cout << "    Model with small Dnu ==> many mixed modes. This might be long to find the solutions..." << std::endl; 
	}

	tau=input_params[19] * pow(cfg_star.numax_star*1e-6,input_params[20]) + input_params[21]; // Granulation timescale (in seconds)
	H=input_params[16] * pow(cfg_star.numax_star*1e-6,input_params[17]) + input_params[18]; // Granulation Amplitude
	H=H/tau ; //This is due to the used definition for the Harvey profile (conversion from Hz to microHz)
	tau=tau/1000. ; //conversion in ksec
	p=input_params[22];// power law:  MUST BE CLOSE TO 2
	N0=input_params[23];
	noise_params(0,0)=-1;
	noise_params(0,1)=-1;
	noise_params(0,2)=-1; 
	noise_params(1,0)=H;
	noise_params(1,1)=tau;
	noise_params(1,2)=p; 
	noise_params(2, 0)=N0; // White noise
	noise_params(2, 1)=-2;
	noise_params(2, 2)=-2;

	// A FUNCTION THAT WRITES THE Noise
	write_star_noise_params(noise_params, file_out_noise);

	//std::cout << "Model finished" << std::endl;
}

void asymptotic_mm_v3(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path, std::string template_file){

	int seed=(unsigned)time(NULL);
	srand(seed);

	std::random_device rd;
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<double> distrib(0 , 1);
	
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
	const std::string cpath=getcwd(NULL, 0);
	
	const std::string file_range=cpath + "/external/ARMM-solver/star_params.range";
	const double lmax=3;
	const double Nmax_pm=6; // Number of radial order to inject on each side of numax

	//std::cout << "input_params=" << input_params.transpose() << std::endl;
	//std::cout << "input_params.size()=" << input_params.size() << std::endl;
	// ------- Deploy the parameters ------
	//double rot_env=input_params[0];
	//double rot_core=input_params[1];	
	//double Dnu=input_params[2];
	//double epsilon=input_params[3];
	//double delta0l_percent=input_params[4];
	//double beta_p=input_params[5];
	double numax_spread=0;
	double nmax_spread=input_params[6];
	double DP_var_percent=input_params[7];
	//double alpha=input_params[8];
	//double q=input_params[9];
	//double hnr_l0=input_params[10];
	//double l0_width_at_numax=input_params[11];
	//double Vl1=input_params[12];
	//double Vl2=input_params[13];
	//double Vl3=input_params[14];
	//double H0_spread=input_params[15];

	double inc_rad, inc_star, inc_y, inc_c;
	double xmin, xmax, a1;
	double H, tau, p, N0;
	MatrixXd mode_params, noise_params(3,3);
	
	Cfg_synthetic_star cfg_star;
	Params_synthetic_star params;
	// ------------------------------------

	// ---- Evaluation of DP ----
	// Super rought estimate derived by visual inspection of the Mosser+2015, Fig.1
	const double c=36.8222;
	const double b=2.63897;
	const double a=0.0168202;
	double DP;
	
		// -----------
	cfg_star.Teff_star=-1;
	cfg_star.Dnu_star=input_params[2];
	cfg_star.epsilon_star=input_params[3];
	cfg_star.delta0l_percent_star=input_params[4];
	cfg_star.beta_p_star=input_params[5];
	
	DP=a*pow(cfg_star.Dnu_star, 2) + b*cfg_star.Dnu_star + c;
	double *r = r8vec_normal_01 ( 1, &seed );
	DP=DP +  *r*DP*DP_var_percent/100.; // Inject a gaussian random error of DP_var_percent
	// -----------
	cfg_star.DPl_star=DP;                
	cfg_star.alpha_g_star=input_params[8];
	cfg_star.q_star=input_params[9];
	cfg_star.maxHNR_l0=input_params[10];
	cfg_star.H0_spread=input_params[15];
	cfg_star.Gamma_max_l0=input_params[11];
	cfg_star.Hfactor=input_params[24];
	cfg_star.Wfactor=input_params[25];
	cfg_star.rot_env_input=input_params[0];
	cfg_star.rot_ratio_input=-1;
	cfg_star.rot_core_input=input_params[1];
	cfg_star.noise_params_harvey_like.resize(8);
	cfg_star.noise_params_harvey_like <<  input_params[16], input_params[17] , input_params[18] , input_params[19] , input_params[20] , input_params[21] , input_params[22]  , input_params[23];    //[A_Pgran ,  B_Pgran , C_Pgran   ,  A_taugran ,  B_taugran  , C_taugran    , p      N0]
	cfg_star.numax_star=numax_from_stello2009(cfg_star.Dnu_star, numax_spread); // Second argument is the random spread on numax
	cfg_star.fmin=cfg_star.numax_star -Nmax_pm*cfg_star.Dnu_star;
	cfg_star.fmax=cfg_star.numax_star +(Nmax_pm+2)*cfg_star.Dnu_star;
	cfg_star.output_file_rot=cpath + "/external/ARMM-solver/star_params.rot";
	cfg_star.Vl.resize(3);
	cfg_star.Vl << input_params[12], input_params[13], input_params[14];
	cfg_star.filetemplate = template_file;
	cfg_star.nmax_star=cfg_star.numax_star/cfg_star.Dnu_star - cfg_star.epsilon_star;
	cfg_star.alpha_p_star=cfg_star.beta_p_star/cfg_star.nmax_star;
	cfg_star.sigma_m=0;
	cfg_star.sigma_p=0;

	if (std::abs(nmax_spread) > 0)
	{
		try
		{
			xmin=cfg_star.nmax_star*(1. - std::abs(nmax_spread)/100.);
			xmax=cfg_star.nmax_star*(1. + std::abs(nmax_spread)/100.);
			cfg_star.nmax_star=xmin + (xmax-xmin)*distrib(gen);
			
		}
		catch(...)
		{
			std::cout << "Error debug info:" << std::endl;
			std::cout << "cfg_star.nmax: " << cfg_star.nmax_star << std::endl;
			std::cout << "nmax_spread: " << nmax_spread << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	// Determination of an isotropic inclination
	inc_rad=distrib(gen)*PI/2;
	inc_y=distrib(gen);
	inc_c=std::cos(inc_rad);
	while (inc_y <=inc_c){
		inc_rad=distrib(gen)*PI/2;
		inc_y=distrib(gen);
		inc_c=std::cos(inc_rad);
	}
	inc_star=inc_rad*180./PI;
	
	// b. Generate the mode profiles and frequencies
	params=make_synthetic_asymptotic_star(cfg_star);
	mode_params=bumpoutputs_2_MatrixXd(params, inc_star); // get the output in a format that can be written with the writting function

	//el, nu, h, w, a1, eta, a3, b, alfa, asym, inc
	write_star_mode_params_act_asym(mode_params, file_out_modes);
	write_range_modes(cfg_star, params, file_range);

	if (cfg_star.Dnu_star <= 15){
		std::cout << "    Model with small Dnu ==> many mixed modes. This might be long to find the solutions..." << std::endl; 
	}
	tau=input_params[19] * pow(cfg_star.numax_star*1e-6,input_params[20]) + input_params[21]; // Granulation timescale (in seconds)
	H=input_params[16] * pow(cfg_star.numax_star*1e-6,input_params[17]) + input_params[18]; // Granulation Amplitude
	H=H/tau ; //This is due to the used definition for the Harvey profile (conversion from Hz to microHz)
	tau=tau/1000. ; //conversion in ksec
	p=input_params[22];// power law:  MUST BE CLOSE TO 2
	N0=input_params[23];
	noise_params(0,0)=-1;
	noise_params(0,1)=-1;
	noise_params(0,2)=-1; 
	noise_params(1,0)=H;
	noise_params(1,1)=tau;
	noise_params(1,2)=p; 
	noise_params(2, 0)=N0; // White noise
	noise_params(2, 1)=-2;
	noise_params(2, 2)=-2;
	// A FUNCTION THAT WRITES THE Noise
	write_star_noise_params(noise_params, file_out_noise);

	//exit(EXIT_SUCCESS);
}

void asymptotic_mm_freeDp_numaxspread_curvepmodes_v1(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path, std::string template_file){

	int seed=(unsigned)time(NULL);
	srand(seed);

	std::random_device rd;
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<double> distrib(0 , 1);
	
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
	const std::string cpath=getcwd(NULL, 0);
	
	const std::string file_range=cpath + "/external/ARMM-solver/star_params.range";
	const double lmax=3;
	const double Nmax_pm=6; // Number of radial order to inject on each side of numax

	//std::cout << "input_params=" << input_params.transpose() << std::endl;
	//std::cout << "input_params.size()=" << input_params.size() << std::endl;
	// ------- Deploy the parameters ------
	//double Teff=input_params[0];
	//double Dnu=input_params[1];
	//double epsilon=input_params[2];
	//double delta0l_percent=input_params[3];
	//double beta_p=input_params[4];
	double nmax_spread=input_params[5];
	//double DP=input_params[6];
	//double alpha=input_params[7];
	//double q=input_params[8];
	//double hnr_l0=input_params[9];
	//double l0_width_at_numax=input_params[10];
	double numax_spread=input_params[11]/100; // Spread in %
	//double Vl1=input_params[12];
	//double Vl2=input_params[13];
	//double Vl3=input_params[14];
	//double H0_spread=input_params[15];

	double inc_rad, inc_star, inc_y, inc_c;
	double xmin, xmax, a1;
	double H, tau, p, N0;
	MatrixXd mode_params, noise_params(3,3);
	
	Cfg_synthetic_star cfg_star;
	Params_synthetic_star params;
	// ------------------------------------

	// ---- Evaluation of DP ----
	// Super rought estimate derived by visual inspection of the Mosser+2015, Fig.1
	const double c=36.8222;
	const double b=2.63897;
	const double a=0.0168202;
	
	// -----------
	cfg_star.Teff_star=input_params[0];
	cfg_star.Dnu_star=input_params[1];
	cfg_star.epsilon_star=input_params[2];
	cfg_star.delta0l_percent_star=input_params[3];
	cfg_star.beta_p_star=input_params[4];
	
	cfg_star.DPl_star=input_params[6];                
	cfg_star.alpha_g_star=input_params[7];
	cfg_star.q_star=input_params[8];
	cfg_star.maxHNR_l0=input_params[9];
	cfg_star.H0_spread=input_params[15];
	cfg_star.Gamma_max_l0=input_params[10];
	cfg_star.Hfactor=input_params[24];
	cfg_star.Wfactor=input_params[25];
	cfg_star.rot_env_input=-1;
	cfg_star.rot_ratio_input=-1;
	cfg_star.rot_core_input=-1;
	cfg_star.noise_params_harvey_like.resize(8);
	cfg_star.noise_params_harvey_like <<  input_params[16], input_params[17] , input_params[18] , input_params[19] , input_params[20] , input_params[21] , input_params[22]  , input_params[23];    //[A_Pgran ,  B_Pgran , C_Pgran   ,  A_taugran ,  B_taugran  , C_taugran    , p      N0]
	cfg_star.numax_star=numax_from_stello2009(cfg_star.Dnu_star, numax_spread); // Second argument is the random spread on numax
	cfg_star.fmin=cfg_star.numax_star -Nmax_pm*cfg_star.Dnu_star;
	cfg_star.fmax=cfg_star.numax_star +(Nmax_pm+2)*cfg_star.Dnu_star;
	cfg_star.output_file_rot=cpath + "/external/ARMM-solver/star_params.rot";
	cfg_star.Vl.resize(3);
	cfg_star.Vl << input_params[12], input_params[13], input_params[14];
	cfg_star.filetemplate = template_file;
	cfg_star.nmax_star=cfg_star.numax_star/cfg_star.Dnu_star - cfg_star.epsilon_star;
	cfg_star.alpha_p_star=cfg_star.beta_p_star/cfg_star.nmax_star;
	cfg_star.sigma_m=0;
	cfg_star.sigma_p=0;

	if (std::abs(nmax_spread) > 0)
	{
		try
		{
			xmin=cfg_star.nmax_star*(1. - std::abs(nmax_spread)/100.);
			xmax=cfg_star.nmax_star*(1. + std::abs(nmax_spread)/100.);
			cfg_star.nmax_star=xmin + (xmax-xmin)*distrib(gen);
			
		}
		catch(...)
		{
			std::cout << "Error debug info:" << std::endl;
			std::cout << "cfg_star.nmax: " << cfg_star.nmax_star << std::endl;
			std::cout << "nmax_spread: " << nmax_spread << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	// Determination of an isotropic inclination
	inc_rad=distrib(gen)*PI/2;
	inc_y=distrib(gen);
	inc_c=std::cos(inc_rad);
	while (inc_y <=inc_c){
		inc_rad=distrib(gen)*PI/2;
		inc_y=distrib(gen);
		inc_c=std::cos(inc_rad);
	}
	inc_star=inc_rad*180./PI;
	
	// b. Generate the mode profiles and frequencies
	params=make_synthetic_asymptotic_star(cfg_star);
	mode_params=bumpoutputs_2_MatrixXd(params, inc_star); // get the output in a format that can be written with the writting function

	//el, nu, h, w, a1, eta, a3, b, alfa, asym, inc
	write_star_mode_params_act_asym(mode_params, file_out_modes);
	write_range_modes(cfg_star, params, file_range);

	if (cfg_star.Dnu_star <= 15){
		std::cout << "    Model with small Dnu ==> many mixed modes. This might be long to find the solutions..." << std::endl; 
	}

	tau=input_params[19] * pow(cfg_star.numax_star*1e-6,input_params[20]) + input_params[21]; // Granulation timescale (in seconds)
	H=input_params[16] * pow(cfg_star.numax_star*1e-6,input_params[17]) + input_params[18]; // Granulation Amplitude
	H=H/tau ; //This is due to the used definition for the Harvey profile (conversion from Hz to microHz)
	tau=tau/1000. ; //conversion in ksec
	p=input_params[22];// power law:  MUST BE CLOSE TO 2
	N0=input_params[23];
	noise_params(0,0)=-1;
	noise_params(0,1)=-1;
	noise_params(0,2)=-1; 
	noise_params(1,0)=H;
	noise_params(1,1)=tau;
	noise_params(1,2)=p; 
	noise_params(2, 0)=N0; // White noise
	noise_params(2, 1)=-2;
	noise_params(2, 2)=-2;
	// A FUNCTION THAT WRITES THE Noise
	write_star_noise_params(noise_params, file_out_noise);

}

void asymptotic_mm_freeDp_numaxspread_curvepmodes_v2(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path, std::string template_file){
	
	int seed=(unsigned)time(NULL);
	srand(seed);

	std::random_device rd;
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<double> distrib(0 , 1);
	
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
	const std::string cpath=getcwd(NULL, 0);
	
	const std::string file_range=cpath + "/external/ARMM-solver/star_params.range";
	const double lmax=3;
	const double Nmax_pm=6; // Number of radial order to inject on each side of numax

	//std::cout << "input_params=" << input_params.transpose() << std::endl;
	//std::cout << "input_params.size()=" << input_params.size() << std::endl;
	// ------- Deploy the parameters ------
	//double rot_env=input_params[0];
	//double rot_ratio=input_params[1];	
	//double Dnu=input_params[2];
	//double epsilon=input_params[3];
	//double delta0l_percent=input_params[4];
	//double beta_p=input_params[5];
	double nmax_spread=input_params[6];
	//double DP=input_params[7];
	//double alpha=input_params[8];
	//double q=input_params[9];
	//double hnr_l0=input_params[10];
	//double l0_width_at_numax=input_params[11];
	double numax_spread=input_params[12]/100.;	
	//double Vl1=input_params[13];
	//double Vl2=input_params[14];
	//double Vl3=input_params[15];
	//double H0_spread=input_params[16];
	double inc_rad, inc_star, inc_y, inc_c;
	double xmin, xmax, a1;
	double H, tau, p, N0;
	MatrixXd mode_params, noise_params(3,3);
	
	Cfg_synthetic_star cfg_star;
	Params_synthetic_star params;
	// ------------------------------------

	// ---- Evaluation of DP ----
	// Super rought estimate derived by visual inspection of the Mosser+2015, Fig.1
	const double c=36.8222;
	const double b=2.63897;
	const double a=0.0168202;
	
	// -----------
	cfg_star.Teff_star=-1;
	cfg_star.Dnu_star=input_params[2];
	cfg_star.epsilon_star=input_params[3];
	cfg_star.delta0l_percent_star=input_params[4];
	cfg_star.beta_p_star=input_params[5];
	
	cfg_star.DPl_star=input_params[7];                
	cfg_star.alpha_g_star=input_params[8];
	cfg_star.q_star=input_params[9];
	cfg_star.maxHNR_l0=input_params[10];
	cfg_star.H0_spread=input_params[16];
	cfg_star.Gamma_max_l0=input_params[11];
	cfg_star.Hfactor=input_params[25];
	cfg_star.Wfactor=input_params[26];
	cfg_star.rot_env_input=input_params[0];
	cfg_star.rot_ratio_input=input_params[1];
	cfg_star.rot_core_input=-1;
	cfg_star.noise_params_harvey_like.resize(8);
	cfg_star.noise_params_harvey_like <<  input_params[17], input_params[18] , input_params[19] , input_params[20] , input_params[21] , input_params[22] , input_params[23]  , input_params[24];    //[A_Pgran ,  B_Pgran , C_Pgran   ,  A_taugran ,  B_taugran  , C_taugran    , p      N0]
	cfg_star.numax_star=numax_from_stello2009(cfg_star.Dnu_star, numax_spread); // Second argument is the random spread on numax
	cfg_star.fmin=cfg_star.numax_star -Nmax_pm*cfg_star.Dnu_star;
	cfg_star.fmax=cfg_star.numax_star +(Nmax_pm+2)*cfg_star.Dnu_star;
	cfg_star.output_file_rot=cpath + "/external/ARMM-solver/star_params.rot";
	cfg_star.Vl.resize(3);
	cfg_star.Vl << input_params[13], input_params[14], input_params[15];
	cfg_star.filetemplate = template_file;
	cfg_star.nmax_star=cfg_star.numax_star/cfg_star.Dnu_star - cfg_star.epsilon_star;
	cfg_star.alpha_p_star=cfg_star.beta_p_star/cfg_star.nmax_star;
	cfg_star.sigma_m=0;
	cfg_star.sigma_p=0;

	if (std::abs(nmax_spread) > 0)
	{
		try
		{
			xmin=cfg_star.nmax_star*(1. - std::abs(nmax_spread)/100.);
			xmax=cfg_star.nmax_star*(1. + std::abs(nmax_spread)/100.);
			cfg_star.nmax_star=xmin + (xmax-xmin)*distrib(gen);
			
		}
		catch(...)
		{
			std::cout << "Error debug info:" << std::endl;
			std::cout << "cfg_star.nmax: " << cfg_star.nmax_star << std::endl;
			std::cout << "nmax_spread: " << nmax_spread << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	// Determination of an isotropic inclination
	inc_rad=distrib(gen)*PI/2;
	inc_y=distrib(gen);
	inc_c=std::cos(inc_rad);
	while (inc_y <=inc_c){
		inc_rad=distrib(gen)*PI/2;
		inc_y=distrib(gen);
		inc_c=std::cos(inc_rad);
	}
	inc_star=inc_rad*180./PI;
	
	// b. Generate the mode profiles and frequencies
	params=make_synthetic_asymptotic_star(cfg_star);
	mode_params=bumpoutputs_2_MatrixXd(params, inc_star); // get the output in a format that can be written with the writting function

	//el, nu, h, w, a1, eta, a3, b, alfa, asym, inc
	write_star_mode_params_act_asym(mode_params, file_out_modes);
	write_range_modes(cfg_star, params, file_range);

	if (cfg_star.Dnu_star <= 15){
		std::cout << "    Model with small Dnu ==> many mixed modes. This might be long to find the solutions..." << std::endl; 
	}

	tau=input_params[20] * pow(cfg_star.numax_star*1e-6,input_params[21]) + input_params[22]; // Granulation timescale (in seconds)
	H=input_params[17] * pow(cfg_star.numax_star*1e-6,input_params[18]) + input_params[19]; // Granulation Amplitude
	H=H/tau ; //This is due to the used definition for the Harvey profile (conversion from Hz to microHz)
	tau=tau/1000. ; //conversion in ksec
	p=input_params[23];// power law:  MUST BE CLOSE TO 2
	N0=input_params[24];
	noise_params(0,0)=-1;
	noise_params(0,1)=-1;
	noise_params(0,2)=-1; 
	noise_params(1,0)=H;
	noise_params(1,1)=tau;
	noise_params(1,2)=p; 
	noise_params(2, 0)=N0; // White noise
	noise_params(2, 1)=-2;
	noise_params(2, 2)=-2;
	// A FUNCTION THAT WRITES THE Noise
	write_star_noise_params(noise_params, file_out_noise);

}

void asymptotic_mm_freeDp_numaxspread_curvepmodes_v3(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path, std::string template_file){

	int seed=(unsigned)time(NULL);
	srand(seed);

	std::random_device rd;
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<double> distrib(0 , 1);
	
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
	const std::string cpath=getcwd(NULL, 0);
	
	const std::string file_range=cpath + "/external/ARMM-solver/star_params.range";
	const double lmax=3;
	const double Nmax_pm=6; // Number of radial order to inject on each side of numax

	//std::cout << "input_params=" << input_params.transpose() << std::endl;
	//std::cout << "input_params.size()=" << input_params.size() << std::endl;
	// ------- Deploy the parameters ------
	//double rot_env=input_params[0];
	//double rot_core=input_params[1];	
	//double Dnu=input_params[2];
	//double epsilon=input_params[3];	
	//double delta0l_percent=input_params[4];
	//double beta_p=input_params[5];
	double nmax_spread=input_params[6];
	//double DP=input_params[7];
	//double alpha=input_params[8];
	//double q=input_params[9];
	//double hnr_l0=input_params[10];
	//double l0_width_at_numax=input_params[11];
	double numax_spread=input_params[12]/100.;
	//double Vl1=input_params[13];
	//double Vl2=input_params[14];
	//double Vl3=input_params[15];
	//double H0_spread=input_params[16];

	double inc_rad, inc_star, inc_y, inc_c;
	double xmin, xmax, a1;
	double H, tau, p, N0;
	MatrixXd mode_params, noise_params(3,3);
	
	Cfg_synthetic_star cfg_star;
	Params_synthetic_star params;
	// ------------------------------------

	// ---- Evaluation of DP ----
	// Super rought estimate derived by visual inspection of the Mosser+2015, Fig.1
	const double c=36.8222;
	const double b=2.63897;
	const double a=0.0168202;
	try{ // Attempting to read file_cfg_mm, if exist. 
		std::cout << "Trying to find and read the mixed modes configuration file" << file_cfg_mm << std::endl;
		cfg_star=read_theoretical_freqs(file_cfg_mm);
		if (cfg_star.use_nu_nl == true){
			std::cout << "    " << file_cfg_mm << " read! and use_nu_nl = true" << std::endl;
			std::cout << "    The frequencies of l=0,1,2,3 listed in the file will be used! " << std::endl;
		}
	} catch(...){
		std::cout << "     " << file_cfg_mm << " not found. Pursuing assuming cfg_star.use_nu_nl=false" << std::endl;
		cfg_star.use_nu_nl=false;
	}
	// ----------
	if (cfg_star.use_nu_nl == false){ // we use the parameters defined in the main.cfg
		cfg_star.Dnu_star=input_params[2];
		cfg_star.DPl_star=input_params[7];                
		cfg_star.q_star=input_params[9];
		cfg_star.alpha_g_star=input_params[8];
		cfg_star.epsilon_star=input_params[3];
	} else{
		std::cout << "Using:" << std::endl;
		std::cout << "    cfg_star.Dnu_star = " << cfg_star.Dnu_star << std::endl;
		std::cout << "    cfg_star.DPl_star = " << cfg_star.DPl_star << std::endl;
		std::cout << "    cfg_star.q_star   = " << cfg_star.q_star << std::endl;
		std::cout << "    cfg_star.alpha_g_star   = " << cfg_star.alpha_g_star << std::endl;
		std::cout << "    epsilon_p_star will be calculated using the provided list of frequencies and Dnu_star" << std::endl;
		std::cout << "    Note that delta0l_percent_star is taken from the main.cfg configuration " << std::endl;
	}
	cfg_star.Teff_star=-1;
	cfg_star.delta0l_percent_star=input_params[4];
	cfg_star.beta_p_star=input_params[5];
	
	cfg_star.maxHNR_l0=input_params[10];
	cfg_star.H0_spread=input_params[16];
	cfg_star.Gamma_max_l0=input_params[11];
	cfg_star.Hfactor=input_params[25];
	cfg_star.Wfactor=input_params[26];
	cfg_star.rot_env_input=input_params[0];
	cfg_star.rot_ratio_input=-1;
	cfg_star.rot_core_input=input_params[1];
	cfg_star.noise_params_harvey_like.resize(8);
	cfg_star.noise_params_harvey_like <<  input_params[17], input_params[18] , input_params[19] , input_params[20] , input_params[21] , input_params[22] , input_params[23]  , input_params[24];    //[A_Pgran ,  B_Pgran , C_Pgran   ,  A_taugran ,  B_taugran  , C_taugran    , p      N0]
	cfg_star.numax_star=numax_from_stello2009(cfg_star.Dnu_star, numax_spread); // Second argument is the random spread on numax
	cfg_star.fmin=cfg_star.numax_star -Nmax_pm*cfg_star.Dnu_star;
	cfg_star.fmax=cfg_star.numax_star +(Nmax_pm+2)*cfg_star.Dnu_star;
	cfg_star.output_file_rot=cpath + "/external/ARMM-solver/star_params.rot";
	cfg_star.Vl.resize(3);
	cfg_star.Vl << input_params[13], input_params[14], input_params[15];
	cfg_star.filetemplate = template_file;
	cfg_star.nmax_star=cfg_star.numax_star/cfg_star.Dnu_star - cfg_star.epsilon_star;
	cfg_star.alpha_p_star=cfg_star.beta_p_star/cfg_star.nmax_star;
	cfg_star.sigma_m=0;
	cfg_star.sigma_p=0;

	if (std::abs(nmax_spread) > 0)
	{
		try
		{
			xmin=cfg_star.nmax_star*(1. - std::abs(nmax_spread)/100.);
			xmax=cfg_star.nmax_star*(1. + std::abs(nmax_spread)/100.);
			cfg_star.nmax_star=xmin + (xmax-xmin)*distrib(gen);
			
		}
		catch(...)
		{
			std::cout << "Error debug info:" << std::endl;
			std::cout << "cfg_star.nmax: " << cfg_star.nmax_star << std::endl;
			std::cout << "nmax_spread: " << nmax_spread << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	// Determination of an isotropic inclination
	inc_rad=distrib(gen)*PI/2;
	inc_y=distrib(gen);
	inc_c=std::cos(inc_rad);
	while (inc_y <=inc_c){
		inc_rad=distrib(gen)*PI/2;
		inc_y=distrib(gen);
		inc_c=std::cos(inc_rad);
	}
	inc_star=inc_rad*180./PI;
	
	// b. Generate the mode profiles and frequencies
	params=make_synthetic_asymptotic_star(cfg_star);
	mode_params=bumpoutputs_2_MatrixXd(params, inc_star); // get the output in a format that can be written with the writting function

	//el, nu, h, w, a1, eta, a3, b, alfa, asym, inc
	write_star_mode_params_act_asym(mode_params, file_out_modes);
	write_range_modes(cfg_star, params, file_range);

	if (cfg_star.Dnu_star <= 15){
		std::cout << "    Model with small Dnu ==> many mixed modes. This might be long to find the solutions..." << std::endl; 
	}

	tau=input_params[20] * pow(cfg_star.numax_star*1e-6,input_params[21]) + input_params[22]; // Granulation timescale (in seconds)
	H=input_params[17] * pow(cfg_star.numax_star*1e-6,input_params[18]) + input_params[19]; // Granulation Amplitude
	//std::cout << "pass" << std::endl;
	H=H/tau ; //This is due to the used definition for the Harvey profile (conversion from Hz to microHz)
	tau=tau/1000. ; //conversion in ksec
	p=input_params[23];// power law:  MUST BE CLOSE TO 2
	N0=input_params[24];
	noise_params(0,0)=-1;
	noise_params(0,1)=-1;
	noise_params(0,2)=-1; 
	noise_params(1,0)=H;
	noise_params(1,1)=tau;
	noise_params(1,2)=p; 
	noise_params(2, 0)=N0; // White noise
	noise_params(2, 1)=-2;
	noise_params(2, 2)=-2;
	//std::cout << "pass" << std::endl;
	//exit(EXIT_SUCCESS);

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
	int el, en, k;
	double n_at_numax, height, eta;
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
	//eta=(4./3)*PI * pow( a1*1e-6 ,2) / (G * rho_sun) * pow(Dnu_sun/Dnu,2);
	eta=pow(a1*1e-6,2)*3.*M_PI/(pow(Dnu/Dnu_sun,2.) * rho_sun*G);  // Corrected on 3 March 2022
 
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
	This function use a reference star as a template to generate frequencies and Width, Height profiles
	can be rescaled so that you can modify the HNR but keep the same height profile
	Note that the user here provides a target a1/Width so that a1 is automatically adjusted to match the 
	requested a1/Width. The code will not change the Width so that code is not adapted to test blending between adjacent l modes,
	such as the l=0 and l=2 mode blending. 
*/
void generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string extra){

	int i;
	double N0, HNR, a1_ov_Gamma, a3, beta_asym, inc, HNRmaxref, Height_factor, Gamma_at_numax, a1, Gamma_coef; 
	VectorXi pos_max;
	VectorXd tmp;
	VectorXd HNRref, local_noise,h_star, gamma_star, s_a1_star, s_a3_star, s_asym_star, inc_star;
	Star_params ref_star;
	MatrixXd mode_params, noise_params;

	ref_star=read_star_params(extra); 

	// Defining the noise profile parameters
	// Note about the noise: -1 means that it is ignored. -2 mean that the value is irrelevant
	tmp.resize(ref_star.mode_params.rows()); 	
	tmp.setConstant(0);
	local_noise=harvey_like(ref_star.noise_params, ref_star.mode_params.col(1), tmp); // Generate a list of local noise values for each frequencies

// ---- Deploy the parameters -----
	HNR=input_params[0];
	a1_ov_Gamma=input_params[1];
	Gamma_at_numax=input_params[2];
	a3=input_params[3];
	beta_asym=input_params[4];
	inc=input_params[5];
// ---------------------------------
	HNRref=ref_star.mode_params.col(2);
	HNRref=HNRref.cwiseProduct(local_noise.cwiseInverse());
	pos_max=where_dbl(ref_star.mode_params.col(0), 0, 1e-3);
	//std::cout << "ref_star.mode_params.col(0) =" << ref_star.mode_params.col(0) << std::endl;
	//std::cout <<" pos_max=" << pos_max << std::endl;
	HNRmaxref=0;
	for (int n=0; n<pos_max.size();n++){
		if (HNRmaxref < HNRref[pos_max[n]]){
			HNRmaxref=HNRref[pos_max[n]];
		}
	}
	//HNRmaxref=HNRref.maxCoeff(); // This is the maximum HNR of the reference data
	Height_factor=HNR/HNRmaxref;  // compute the coeficient required to get the proper max(HNR)
	pos_max=where_dbl(HNRref, HNRmaxref, 0.001);
	if (pos_max[0] >= 0){
		Gamma_coef=Gamma_at_numax/ref_star.mode_params(pos_max[0], 3); // Correction coeficient to apply on Gamma(nu) in order to ensure that we have Gamma(nu=numax) = Gamma_at_numax
	} else{
		std::cout << "Error! could not find the max position for the mode Widths profile" << std::endl;
		std::cout << "Code debug required" << std::endl;
		exit(EXIT_FAILURE);
	}

	a1=a1_ov_Gamma*Gamma_at_numax; // We can vary the Width and splitting. But we need to change the splitting in order to get the wished a1/Gamma0

	// Defining the final size for all of the outptus
	gamma_star.resize(ref_star.mode_params.rows());
	gamma_star=Gamma_coef*ref_star.mode_params.col(3); // In IDL, AN INTERPOLATION WAS DONE FOR l>0. HERE WE ASSUME THE .in file is whatever the true model should be (no interpolation)
	// Refactoring the heights
	h_star.resize(ref_star.mode_params.rows());
	//h_star=Height_factor * HNRref * N0; 
	if (local_noise.sum()!= 0){
		//N0=1; // Imposing the the Noise background is 1
		//std::cout <<  "Using N0=" << N0 << " (white noise)" << std::endl;
		std::cout << "Using the Noise profile of the Reference star" << std::endl;
		std::cout << "HNR of all modes (degree / freq / HNR  /  Height /  local_noise):" << std::endl;
		for(i =0; i<ref_star.mode_params.rows(); i++){
			h_star[i]=Height_factor * HNRref[i] * local_noise[i];
			std::cout << "     " << ref_star.mode_params(i,0) << "  " << ref_star.mode_params(i,1) << "  " << Height_factor * HNRref[i]  << "  " << h_star[i] << "  "  << local_noise[i] << std::endl;
		}
	} else{
		std::cout << "Warning: bruit_local from the stat_synthese file is 0 ==> Cannot compute N0=mean(local_noise)" << std::endl;
		std::cout << "         The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
	}

	// Defining the final size for all of the outptus
	//h_star.resize(ref_star.mode_params.rows());
	//gamma_star.resize(ref_star.mode_params.rows());
	s_a1_star.resize(ref_star.mode_params.rows());
	s_a3_star.resize(ref_star.mode_params.rows());	
	s_asym_star.resize(ref_star.mode_params.rows());
	inc_star.resize(ref_star.mode_params.rows());

	//gamma_star=Gamma_coef*ref_star.mode_params.col(3); // In IDL, AN INTERPOLATION WAS DONE FOR l>0. HERE WE ASSUME THE .in file is whatever the true model should be (no interpolation)
	
	if (a1 >= 0){
		s_a1_star.setConstant(a1);
	} else{
		s_a1_star=ref_star.mode_params.col(4); 
	}
	if (a3 >= 0){
		s_a3_star.setConstant(a3);	
	} else{
		s_a3_star=ref_star.mode_params.col(6);	
	}
	if (beta_asym >= 0){
		s_asym_star.setConstant(beta_asym);	
	} else{
		s_asym_star=ref_star.mode_params.col(9);	
	}	
	if (inc >= 0){ 
		inc_star.setConstant(inc);
	} else{
		inc_star=ref_star.mode_params.col(10);
	}
	
	mode_params.setZero(ref_star.mode_params.rows(), 11);

	mode_params.col(0)=ref_star.mode_params.col(0); // List of els
	mode_params.col(1)=ref_star.mode_params.col(1); // List of frequencies
	mode_params.col(2)=h_star;
	mode_params.col(3)=gamma_star; 
	mode_params.col(4)=s_a1_star; 
	mode_params.col(5).setConstant(0); //ref_star.mode_params.col(5); // asphericity coefficient depends on a1, cannot use the one from the ref_star
	mode_params.col(6)=s_a3_star;
	//mode_params.col(7])=0;   // Already to 0 due to initialisation
	//mode_params.col(8)=0;   // Already to 0 due to initialisation
	mode_params.col(9)=s_asym_star;
	mode_params.col(10)=inc_star;  

	// A FUNCTION THAT WRITES THE PARAMETERS
	write_star_mode_params_act_asym(mode_params, file_out_modes);

	// A FUNCTION THAT WRITES THE Noise
	write_star_noise_params(ref_star.noise_params, file_out_noise);
}

/* 
	This function use a reference star as a template to generate frequencies and Width, Height profiles
	can be rescaled so that you can modify the HNR but keep the same height profile
	Note that the user here provides a target a1/Width so that a1 is automatically adjusted to match the 
	requested a1/Width. The code will not change the Width so that code is not adapted to test blending between adjacent l modes,
	such as the l=0 and l=2 mode blending
	It handles a2 and a3 as well as free parameters. Thosse can be described by 2nd Order polynomial function of frequency
	and you can add a random quantity (set in nHz) that will be added as a random 'error' ==> scatter ==> test robustness 
*/
void generate_cfg_from_synthese_file_Wscaled_a1a2a3asymovGamma(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string extra){

	int i;
	double HNR, a1_ov_Gamma, a3, beta_asym, inc, HNRmaxref, Height_factor, Gamma_at_numax, a1, Gamma_coef, fl, a2; 
	VectorXi pos_max;
	VectorXd tmp;
	VectorXd HNRref, local_noise,h_star, gamma_star, s_a1_star, s_a2_star, s_a3_star, s_asym_star, inc_star, a2_terms, a3_terms;
	Star_params ref_star;
	MatrixXd mode_params, noise_params;

	ref_star=read_star_params(extra); 

	// Defining the noise profile parameters
	// Note about the noise: -1 means that it is ignored. -2 mean that the value is irrelevant
	tmp.resize(ref_star.mode_params.rows()); 	
	tmp.setConstant(0);
	local_noise=harvey_like(ref_star.noise_params, ref_star.mode_params.col(1), tmp); // Generate a list of local noise values for each frequencies

// ---- Deploy the parameters -----
	HNR=input_params[0];
	a1_ov_Gamma=input_params[1];
	Gamma_at_numax=input_params[2];
	a2_terms=input_params.segment(3,3);
	a3_terms=input_params.segment(6,3);
	beta_asym=input_params[9];
	inc=input_params[10];
// ---------------------------------
	HNRref=ref_star.mode_params.col(2);
	HNRref=HNRref.cwiseProduct(local_noise.cwiseInverse());
	pos_max=where_dbl(ref_star.mode_params.col(0), 0, 1e-3);
	//std::cout << "ref_star.mode_params.col(0) =" << ref_star.mode_params.col(0) << std::endl;
	//std::cout <<" pos_max=" << pos_max << std::endl;
	HNRmaxref=0;
	for (int n=0; n<pos_max.size();n++){
		if (HNRmaxref < HNRref[pos_max[n]]){
			HNRmaxref=HNRref[pos_max[n]];
		}
	}
	//HNRmaxref=HNRref.maxCoeff(); // This is the maximum HNR of the reference data
	Height_factor=HNR/HNRmaxref;  // compute the coeficient required to get the proper max(HNR)
	pos_max=where_dbl(HNRref, HNRmaxref, 0.001);
	if (pos_max[0] >= 0){
		Gamma_coef=Gamma_at_numax/ref_star.mode_params(pos_max[0], 3); // Correction coeficient to apply on Gamma(nu) in order to ensure that we have Gamma(nu=numax) = Gamma_at_numax
	} else{
		std::cout << "Error! could not find the max position for the mode Widths profile" << std::endl;
		std::cout << "Code debug required" << std::endl;
		exit(EXIT_FAILURE);
	}

	a1=a1_ov_Gamma*Gamma_at_numax; // We can vary the Width and splitting. But we need to change the splitting in order to get the wished a1/Gamma0

	// Defining the final size for all of the outptus
	gamma_star.resize(ref_star.mode_params.rows());
	gamma_star=Gamma_coef*ref_star.mode_params.col(3); // In IDL, AN INTERPOLATION WAS DONE FOR l>0. HERE WE ASSUME THE .in file is whatever the true model should be (no interpolation)
	// Refactoring the heights
	h_star.resize(ref_star.mode_params.rows());
	//h_star=Height_factor * HNRref * N0; 
	if (local_noise.sum()!= 0){
		//N0=1; // Imposing the the Noise background is 1
		//std::cout <<  "Using N0=" << N0 << " (white noise)" << std::endl;
		std::cout << "Using the Noise profile of the Reference star" << std::endl;
		std::cout << "HNR of all modes (degree / freq / HNR  /  Height /  local_noise):" << std::endl;
		for(i =0; i<ref_star.mode_params.rows(); i++){
			h_star[i]=Height_factor * HNRref[i] * local_noise[i];
			std::cout << "     " << ref_star.mode_params(i,0) << "  " << ref_star.mode_params(i,1) << "  " << Height_factor * HNRref[i]  << "  " << h_star[i] << "  "  << local_noise[i] << std::endl;
		}
	} else{
		std::cout << "Warning: bruit_local from the stat_synthese file is 0 ==> Cannot compute N0=mean(local_noise)" << std::endl;
		std::cout << "         The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
	}
	// Defining the final size for all of the outptus
	//h_star.resize(ref_star.mode_params.rows());
	//gamma_star.resize(ref_star.mode_params.rows());
	s_a1_star.resize(ref_star.mode_params.rows());
	s_a2_star.resize(ref_star.mode_params.rows());
	s_a3_star.resize(ref_star.mode_params.rows());	
	s_asym_star.resize(ref_star.mode_params.rows());
	inc_star.resize(ref_star.mode_params.rows());

	// Refactoring the heights
	//h_star=Height_factor * HNRref * N0; 
	//gamma_star=Gamma_coef*ref_star.mode_params.col(3); // In IDL, AN INTERPOLATION WAS DONE FOR l>0. HERE WE ASSUME THE .in file is whatever the true model should be (no interpolation)

	if (a1 >= 0){
		s_a1_star.setConstant(a1);
	} else{
		s_a1_star=ref_star.mode_params.col(4); 
	}
	if (a2_terms[0] != 0 || a2_terms[1] != 0 || a2_terms[2] != 0){
		for (int i=0; i<ref_star.mode_params.rows();i++){ 
			fl=ref_star.mode_params(i,1); // Get the frequencies and iterate over them
			a2=a2_terms[0] + a2_terms[1]*(fl*1e-3) + a2_terms[2]*(fl*fl*1e-6); //three terms: one constant term + one linear in nu + one quadratic in nu [BE CAREFUL WITH UNITS]
			s_a2_star[i]=a2;
		}
		//s_a2_star.setConstant(a2);	
	} else{
		s_a2_star.setConstant(0); // We cannot use the asphericity parameter as in the synthese.in reference files as those consider the asumption of asphericity through the beta parameter, not a2	
	}
	if (a3_terms[0] != 0 || a3_terms[1] != 0 || a3_terms[2] != 0){
		for (int i=0; i<ref_star.mode_params.rows();i++){ 
			fl=ref_star.mode_params(i,1); // Get the frequencies and iterate over them
			a3=a3_terms[0] + a3_terms[1]*(fl*1e-3) + a3_terms[2]*(fl*fl*1e-6); //three terms: one constant term + one linear in nu + one quadratic in nu [BE CAREFUL WITH UNITS]
			s_a3_star[i]=a3;
		}

	} else{
		s_a3_star=ref_star.mode_params.col(6);	
	}
	if (beta_asym >= 0){
		s_asym_star.setConstant(beta_asym);	
	} else{
		s_asym_star=ref_star.mode_params.col(9);	
	}	
	if (inc >= 0){ 
		inc_star.setConstant(inc);
	} else{
		inc_star=ref_star.mode_params.col(10);
	}
	
	mode_params.setZero(ref_star.mode_params.rows(), 11);

	mode_params.col(0)=ref_star.mode_params.col(0); // List of els
	mode_params.col(1)=ref_star.mode_params.col(1); // List of frequencies
	mode_params.col(2)=h_star;
	mode_params.col(3)=gamma_star; 
	mode_params.col(4)=s_a1_star; 
	mode_params.col(5)=s_a2_star; //ref_star.mode_params.col(5); // asphericity coefficient depends on a1, cannot use the one from the ref_star
	mode_params.col(6)=s_a3_star;
	//mode_params.col(7])=0;   // Already to 0 due to initialisation
	//mode_params.col(8)=0;   // Already to 0 due to initialisation
	mode_params.col(9)=s_asym_star;
	mode_params.col(10)=inc_star;  

	// A FUNCTION THAT WRITES THE PARAMETERS
	write_star_mode_params_act_asym(mode_params, file_out_modes);

	// A FUNCTION THAT WRITES THE Noise
	write_star_noise_params(ref_star.noise_params, file_out_noise);
}



/* 
	This function use a reference star as a template to generate frequencies and Width, Height profiles
	can be rescaled so that you can modify the HNR but keep the same height profile
	Note that the user here provides a target a1/Width so that a1 is automatically adjusted to match the 
	requested a1/Width. The code will not change the Width so that code is not adapted to test blending between adjacent l modes,
	such as the l=0 and l=2 mode blending
	It handles aj coeficient up to j=6 as well as free parameters and consider an activity term Alm ~ {a2, a4, a6, ...} following Gizon2002 idea.
	a3 can be given as a polynomial of the frequency and you can add a random quantity (set in nHz) that will be added as a random 'error' ==> scatter ==> test robustness 
*/
void generate_cfg_from_synthese_file_Wscaled_Alm(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string extra){

	int i;
	double HNR, a1_ov_Gamma, a2, a3, a4, a5, a6, beta_asym, inc, HNRmaxref, Height_factor, Gamma_at_numax, a1, Gamma_coef, fl,rho, eta0; 
	VectorXi pos_max;
	VectorXd tmp, xfit, rfit;
	VectorXd HNRref, local_noise,h_star, gamma_star, s_a1_star, s_a2_star, s_a3_star, s_a4_star,s_a5_star,s_a6_star,
		s_eta0_star, s_epsilon_star, s_theta0_star, s_delta_star, s_asym_star, inc_star, activity_terms;
	Star_params ref_star;
	MatrixXd mode_params, noise_params;

	ref_star=read_star_params(extra); 

	eta0=eta0_fct(ref_star.mode_params.col(1));

	// Defining the noise profile parameters
	// Note about the noise: -1 means that it is ignored. -2 mean that the value is irrelevant
	tmp.resize(ref_star.mode_params.rows()); 	
	tmp.setConstant(0);
	local_noise=harvey_like(ref_star.noise_params, ref_star.mode_params.col(1), tmp); // Generate a list of local noise values for each frequencies

// ---- Deploy the parameters -----
	HNR=input_params[0];
	a1_ov_Gamma=input_params[1];
	Gamma_at_numax=input_params[2];
	activity_terms=input_params.segment(3,3); // epsilon_nl , theta0, delta.... used for computing dnu_nlm_activity = nu_nl.epsilon_nl.Alm(theta0, delta)
	a3=input_params[6];
	beta_asym=input_params[7];
	inc=input_params[8];
// ---------------------------------
	HNRref=ref_star.mode_params.col(2);
	HNRref=HNRref.cwiseProduct(local_noise.cwiseInverse());
	pos_max=where_dbl(ref_star.mode_params.col(0), 0, 1e-3);
	//std::cout << "ref_star.mode_params.col(0) =" << ref_star.mode_params.col(0) << std::endl;
	//std::cout <<" pos_max=" << pos_max << std::endl;
	HNRmaxref=0;
	for (int n=0; n<pos_max.size();n++){
		if (HNRmaxref < HNRref[pos_max[n]]){
			HNRmaxref=HNRref[pos_max[n]];
		}
	}
	//HNRmaxref=HNRref.maxCoeff(); // This is the maximum HNR of the reference data
	Height_factor=HNR/HNRmaxref;  // compute the coeficient required to get the proper max(HNR)
	pos_max=where_dbl(HNRref, HNRmaxref, 0.001);
	if (pos_max[0] >= 0){
		Gamma_coef=Gamma_at_numax/ref_star.mode_params(pos_max[0], 3); // Correction coeficient to apply on Gamma(nu) in order to ensure that we have Gamma(nu=numax) = Gamma_at_numax
	} else{
		std::cout << "Error! could not find the max position for the mode Widths profile" << std::endl;
		std::cout << "Code debug required" << std::endl;
		exit(EXIT_FAILURE);
	}

	a1=a1_ov_Gamma*Gamma_at_numax; // We can vary the Width and splitting. But we need to change the splitting in order to get the wished a1/Gamma0

	// Defining the final size for all of the outptus
	gamma_star.resize(ref_star.mode_params.rows());
	gamma_star=Gamma_coef*ref_star.mode_params.col(3); // In IDL, AN INTERPOLATION WAS DONE FOR l>0. HERE WE ASSUME THE .in file is whatever the true model should be (no interpolation)
	// Refactoring the heights
	h_star.resize(ref_star.mode_params.rows());
	//h_star=Height_factor * HNRref * N0; 
	if (local_noise.sum()!= 0){
		//N0=1; // Imposing the the Noise background is 1
		//std::cout <<  "Using N0=" << N0 << " (white noise)" << std::endl;
		std::cout << "Using the Noise profile of the Reference star" << std::endl;
		std::cout << "HNR of all modes (degree / freq / HNR  /  Height /  local_noise):" << std::endl;
		for(i =0; i<ref_star.mode_params.rows(); i++){
			h_star[i]=Height_factor * HNRref[i] * local_noise[i];
			std::cout << "     " << ref_star.mode_params(i,0) << "  " << ref_star.mode_params(i,1) << "  " << Height_factor * HNRref[i]  << "  " << h_star[i] << "  "  << local_noise[i] << std::endl;
		}
	} else{
		std::cout << "Warning: bruit_local from the stat_synthese file is 0 ==> Cannot compute N0=mean(local_noise)" << std::endl;
		std::cout << "         The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
	}

	if (a1 < 0){
		std::cout << "Error: The a1 provided by the user is negative" << std::endl;
		std::cout << "       Only positive values are valid" << std::endl;
		exit(EXIT_FAILURE);
	}
	if (inc < 0){ 
		std::cout << "Error: The inclination provided by the user is negative" << std::endl;
		std::cout << "       Only positive values are valid" << std::endl;
		exit(EXIT_FAILURE);
	}
	mode_params.setZero(ref_star.mode_params.rows(), 12);

	mode_params.col(0)=ref_star.mode_params.col(0); // List of els
	mode_params.col(1)=ref_star.mode_params.col(1); // List of frequencies
	mode_params.col(2)=h_star;
	mode_params.col(3)=gamma_star; 
	mode_params.col(4).setConstant(a1);//=s_a1_star; 
	mode_params.col(5).setConstant(eta0);//=s_eta0_star; 
	mode_params.col(6).setConstant(activity_terms[0]);//=s_epsilon_star;
	mode_params.col(7).setConstant(activity_terms[1]);//=s_theta0_star; 
	mode_params.col(8).setConstant(activity_terms[2]);//=s_delta_star;
	mode_params.col(9).setConstant(a3);//=s_a3_star;
	mode_params.col(10).setConstant(beta_asym);//=s_asym_star;
	mode_params.col(11).setConstant(inc);//=inc_star;  

	// A FUNCTION THAT WRITES THE PARAMETERS
	// THIS MIGHT NEED TO BE WRITTEN IN ANOTHER FORMAT
	write_star_mode_params_Alm(mode_params, file_out_modes);
	// A FUNCTION THAT WRITES THE Noise
	write_star_noise_params(ref_star.noise_params, file_out_noise);
}


/* 
	This function use a reference star as a template to generate frequencies and Width, Height profiles
	can be rescaled so that you can modify the HNR but keep the same height profile
	Note that the user here provides a target a1/Width so that a1 is automatically adjusted to match the 
	requested a1/Width. The code will not change the Width so that code is not adapted to test blending between adjacent l modes,
	such as the l=0 and l=2 mode blending
	It handles aj coeficient up to j=6 as well as free parameters and consider an activity term Alm ~ {a2, a4, a6, ...} following Gizon2002 idea.
	a3 can be given as a polynomial of the frequency and you can add a random quantity (set in nHz) that will be added as a random 'error' ==> scatter ==> test robustness 
*/
void generate_cfg_from_synthese_file_Wscaled_aj(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string extra){

	int i;
	double HNR, a1_ov_Gamma, a2, a3, a4, a5, a6, beta_asym, inc, HNRmaxref, Height_factor, Gamma_at_numax, a1, Gamma_coef, fl,rho, eta0; 
	VectorXi pos_max;
	VectorXd tmp, xfit, rfit;
	VectorXd HNRref, local_noise,h_star, gamma_star, s_a1_star, s_a2_star, s_a3_star, s_a4_star,s_a5_star,s_a6_star,
		s_eta0_star, s_epsilon_star, s_theta0_star, s_delta_star, s_asym_star, inc_star, activity_terms;
	Star_params ref_star;
	MatrixXd mode_params, noise_params;

	ref_star=read_star_params(extra); 

	// Defining the noise profile parameters
	// Note about the noise: -1 means that it is ignored. -2 mean that the value is irrelevant
	tmp.resize(ref_star.mode_params.rows()); 	
	tmp.setConstant(0);
	local_noise=harvey_like(ref_star.noise_params, ref_star.mode_params.col(1), tmp); // Generate a list of local noise values for each frequencies

// ---- Deploy the parameters -----
	HNR=input_params[0];
	a1_ov_Gamma=input_params[1];
	Gamma_at_numax=input_params[2];
	a2=input_params[3];
	a3=input_params[4];
	a4=input_params[5];
	a5=input_params[6];
	a6=input_params[7];
	beta_asym=input_params[8];
	inc=input_params[9];
// ---------------------------------
	HNRref=ref_star.mode_params.col(2);
	HNRref=HNRref.cwiseProduct(local_noise.cwiseInverse());
	pos_max=where_dbl(ref_star.mode_params.col(0), 0, 1e-3);
	//std::cout << "ref_star.mode_params.col(0) =" << ref_star.mode_params.col(0) << std::endl;
	//std::cout <<" pos_max=" << pos_max << std::endl;
	HNRmaxref=0;
	for (int n=0; n<pos_max.size();n++){
		if (HNRmaxref < HNRref[pos_max[n]]){
			HNRmaxref=HNRref[pos_max[n]];
		}
	}
	//HNRmaxref=HNRref.maxCoeff(); // This is the maximum HNR of the reference data
	Height_factor=HNR/HNRmaxref;  // compute the coeficient required to get the proper max(HNR)
	pos_max=where_dbl(HNRref, HNRmaxref, 0.001);
	if (pos_max[0] >= 0){
		Gamma_coef=Gamma_at_numax/ref_star.mode_params(pos_max[0], 3); // Correction coeficient to apply on Gamma(nu) in order to ensure that we have Gamma(nu=numax) = Gamma_at_numax
	} else{
		std::cout << "Error! could not find the max position for the mode Widths profile" << std::endl;
		std::cout << "Code debug required" << std::endl;
		exit(EXIT_FAILURE);
	}

	a1=a1_ov_Gamma*Gamma_at_numax; // We can vary the Width and splitting. But we need to change the splitting in order to get the wished a1/Gamma0

	// Defining the final size for all of the outptus
	gamma_star.resize(ref_star.mode_params.rows());
	gamma_star=Gamma_coef*ref_star.mode_params.col(3); // In IDL, AN INTERPOLATION WAS DONE FOR l>0. HERE WE ASSUME THE .in file is whatever the true model should be (no interpolation)
	// Refactoring the heights
	h_star.resize(ref_star.mode_params.rows());
	//h_star=Height_factor * HNRref * N0; 
	if (local_noise.sum()!= 0){
		//N0=1; // Imposing the the Noise background is 1
		//std::cout <<  "Using N0=" << N0 << " (white noise)" << std::endl;
		std::cout << "Using the Noise profile of the Reference star" << std::endl;
		std::cout << "HNR of all modes (degree / freq / HNR  /  Height /  local_noise):" << std::endl;
		for(i =0; i<ref_star.mode_params.rows(); i++){
			h_star[i]=Height_factor * HNRref[i] * local_noise[i];
			std::cout << "     " << ref_star.mode_params(i,0) << "  " << ref_star.mode_params(i,1) << "  " << Height_factor * HNRref[i]  << "  " << h_star[i] << "  "  << local_noise[i] << std::endl;
		}
	} else{
		std::cout << "Warning: bruit_local from the stat_synthese file is 0 ==> Cannot compute N0=mean(local_noise)" << std::endl;
		std::cout << "         The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
	}
	if (a1 < 0){
		std::cout << "Error: The a1 provided by the user is negative" << std::endl;
		std::cout << "       Only positive values are valid" << std::endl;
		exit(EXIT_FAILURE);
	}
	if (inc < 0){ 
		std::cout << "Error: The inclination provided by the user is negative" << std::endl;
		std::cout << "       Only positive values are valid" << std::endl;
		exit(EXIT_FAILURE);
	}
	mode_params.setZero(ref_star.mode_params.rows(), 16);

	mode_params.col(0)=ref_star.mode_params.col(0); // List of els
	mode_params.col(1)=ref_star.mode_params.col(1); // List of frequencies
	mode_params.col(2)=h_star;
	mode_params.col(3)=gamma_star; 
	mode_params.col(4).setConstant(a1);//=s_a1_star; 
	mode_params.col(5).setConstant(a2);//=s_a2_star;
	mode_params.col(6).setConstant(a3);//=s_a3_star;
	mode_params.col(7).setConstant(a4);//=s_a4_star;
	mode_params.col(8).setConstant(a5);//=s_a5_star;
	mode_params.col(9).setConstant(a6);//=s_a6_star;
	mode_params.col(10).setConstant(beta_asym);//=s_asym_star;
	mode_params.col(11).setConstant(inc);//=inc_star;  
	// A FUNCTION THAT WRITES THE PARAMETERS
	// THIS MIGHT NEED TO BE WRITTEN IN ANOTHER FORMAT
	write_star_mode_params_aj(mode_params, file_out_modes);

	// A FUNCTION THAT WRITES THE Noise
	write_star_noise_params(ref_star.noise_params, file_out_noise);
}


// ------ Common ------
double eta0_fct(const VectorXd& fl0_all){
	const double G=6.667e-8;
    const double Dnu_sun=135.1;
    const double R_sun=6.96342e5; //in km
    const double M_sun=1.98855e30; //in kg
    const double rho_sun=M_sun*1e3/(4*M_PI*std::pow(R_sun*1e5,3)/3); //in g.cm-3
    double rho, eta0;
    VectorXd xfit, rfit;
    xfit=linspace(0, fl0_all.size()-1, fl0_all.size());
    rfit=linfit(xfit, fl0_all); // rfit[0] = Dnu 
    rho=pow(rfit[0]/Dnu_sun,2.) * rho_sun;
    //eta0=3./(4.*M_PI*rho*G); // WRONG BECAUSE I WRONGLY CONSIDERED OMEGA ~ a1. It should be OMEGA ~ 2.pi.a1 of course
    eta0=3.*M_PI/(rho*G); 
    return eta0;
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
