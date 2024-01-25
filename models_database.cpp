/**
 * @file models_database.cpp
 * @brief Header file that contains models
 *
 * File that contains all kind of methods used to generate models for the pulsation/noise
 * 
 *
 * @date 20 Apr 2016
 * @author obenomar
 */

# include <iostream>
# include <iomanip>
#include <fstream>
# include <string>
# include <Eigen/Dense>
#include <random>
#include "models_database.h"
#include "noise_models.h"
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

	// ------- Deploy the parameters ------
	//double rot_env=input_params[0];
	//double rot_ratio=input_params[1];	
	//double a2_l1_core=input_params[2];	
	//double a2_l1_env=input_params[3];	
	//double a2_l2_env=input_params[4];	
	//double a2_l3_env=input_params[5];	
	//double a3_l2_env=input_params[6];	
	//double a3_l3_env=input_params[7];	
	//double a4_l2_env=input_params[8];	
	//double a4_l3_env=input_params[9];		
	//double a5_l3_env=input_params[10];	
	//double a6_l3_env=input_params[11];	
	//double Dnu=input_params[12];
	//double epsilon=input_params[13];
	//double delta0l_percent=input_params[14];
	//double beta_p=input_params[15];
	double nmax_spread=input_params[16];
	//double DP=input_params[17];
	//double alpha=input_params[18];
	//double q=input_params[19];
	//double hnr_l0=input_params[20];
	//double l0_width_at_numax=input_params[21];
	double numax_spread=input_params[22]/100.;	
	//double Vl1=input_params[23];
	//double Vl2=input_params[24];
	//double Vl3=input_params[25];
	//double H0_spread=input_params[26];
	//double A_Pgran=input_params[27];
	//double B_Pgran=input_params[28];
	//double C_Pgran=input_params[29];
	//double A_taugran=input_params[30];
	//double B_taugran=input_params[31];
	//double C_taugran=input_params[32];
	//double P=input_params[33];
	//double N0=input_params[34];
	//double Hfactor=input_params[35];
	//double Wfactor=input_params[36];

	double inc_rad, inc_star, inc_y, inc_c;
	double xmin, xmax, a1;
	double H, tau, p, N0;
	MatrixXd mode_params, noise_params(3,3);
	
	Cfg_synthetic_star cfg_star;
	Params_synthetic_star params;
	// ------------------------------------
	
	// -----------
	cfg_star.Teff_star=-1;
	cfg_star.Dnu_star=input_params[12];
	cfg_star.epsilon_star=input_params[13];
	cfg_star.delta0l_percent_star=input_params[14];
	cfg_star.beta_p_star=input_params[15];
	
	cfg_star.DPl_star=input_params[17];                
	cfg_star.alpha_g_star=input_params[18];
	cfg_star.q_star=input_params[19];
	cfg_star.maxHNR_l0=input_params[20];
	cfg_star.H0_spread=input_params[26];
	cfg_star.Gamma_max_l0=input_params[21];
	cfg_star.Hfactor=input_params[35];
	cfg_star.Wfactor=input_params[36];
	cfg_star.rot_env_input=input_params[0];
	cfg_star.rot_ratio_input=input_params[1];
	cfg_star.rot_core_input=-1;
	cfg_star.env_aspher.a2_l1=0; // SET TO 0 for l=1 modes. Ideally, would need a2_l1_core mixed with a2_l1_env
	cfg_star.env_aspher.a2_l2=input_params[4];
	if(input_params[5] <= -9999){ // If the user want l2 and l3 having the same aj coefficient, they need to put l3 to -9999 or smaller
		cfg_star.env_aspher.a2_l3=cfg_star.env_aspher.a2_l2;
	} else{
		cfg_star.env_aspher.a2_l3=input_params[5];
	}	
	cfg_star.env_aspher.a4_l2=input_params[8];
	if(input_params[5] <= -9999){ // If the user want l2 and l3 having the same aj coefficient, they need to put l3 to -9999 or smaller
		cfg_star.env_aspher.a4_l3=cfg_star.env_aspher.a4_l2;
	} else{
		cfg_star.env_aspher.a4_l3=input_params[9];
	}
	cfg_star.env_aspher.a6_l3=input_params[11];
	cfg_star.env_lat_dif_rot.a3_l2=input_params[6];
	if(input_params[5] <= -9999){ // If the user want l2 and l3 having the same aj coefficient, they need to put l3 to -9999 or smaller
		cfg_star.env_lat_dif_rot.a3_l3=cfg_star.env_lat_dif_rot.a3_l2;
	} else{
		cfg_star.env_lat_dif_rot.a3_l3=input_params[7];
	}
	cfg_star.env_lat_dif_rot.a5_l3=input_params[10];
	
	
	cfg_star.noise_params_harvey_like.resize(8);
	cfg_star.noise_params_harvey_like <<  input_params[17], input_params[18] , input_params[19] , input_params[20] , input_params[21] , input_params[22] , input_params[23]  , input_params[24];    //[A_Pgran ,  B_Pgran , C_Pgran   ,  A_taugran ,  B_taugran  , C_taugran    , p      N0]
	cfg_star.numax_star=numax_from_stello2009(cfg_star.Dnu_star, numax_spread); // Second argument is the random spread on numax
	cfg_star.fmin=cfg_star.numax_star -Nmax_pm*cfg_star.Dnu_star;
	cfg_star.fmax=cfg_star.numax_star +(Nmax_pm+2)*cfg_star.Dnu_star;
	cfg_star.output_file_rot=cpath + "/external/ARMM-solver/star_params.rot";
	cfg_star.Vl.resize(3);
	cfg_star.Vl << input_params[23], input_params[24], input_params[25];
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

	write_star_mode_params_aj(mode_params, file_out_modes);
	write_range_modes(cfg_star, params, file_range);

	if (cfg_star.Dnu_star <= 15){
		std::cout << "    Model with small Dnu ==> many mixed modes. This might be long to find the solutions..." << std::endl; 
	}

	tau=input_params[30] * pow(cfg_star.numax_star*1e-6,input_params[31]) + input_params[32]; // Granulation timescale (in seconds)
	H=input_params[27] * pow(cfg_star.numax_star*1e-6,input_params[28]) + input_params[29]; // Granulation Amplitude
	H=H/tau ; //This is due to the used definition for the Harvey profile (conversion from Hz to microHz)
	tau=tau/1000. ; //conversion in ksec
	p=input_params[33];// power law:  MUST BE CLOSE TO 2
	N0=input_params[34];
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
	//double a2_l1_core=input_params[2];	
	//double a2_l1_env=input_params[3];	
	//double a2_l2_env=input_params[4];	
	//double a2_l3_env=input_params[5];	
	//double a3_l2_env=input_params[6];	
	//double a3_l3_env=input_params[7];	
	//double a4_l2_env=input_params[8];	
	//double a4_l3_env=input_params[9];		
	//double a5_l3_env=input_params[10];	
	//double a6_l3_env=input_params[11];	
	//double Dnu=input_params[12];
	//double epsilon=input_params[13];
	//double delta0l_percent=input_params[14];
	//double beta_p=input_params[15];
	double nmax_spread=input_params[16];
	//double DP=input_params[17];
	//double alpha=input_params[18];
	//double q=input_params[19];
	//double hnr_l0=input_params[20];
	//double l0_width_at_numax=input_params[21];
	double numax_spread=input_params[22]/100.;	
	//double Vl1=input_params[23];
	//double Vl2=input_params[24];
	//double Vl3=input_params[25];
	//double H0_spread=input_params[26];
	//double A_Pgran=input_params[27];
	//double B_Pgran=input_params[28];
	//double C_Pgran=input_params[29];
	//double A_taugran=input_params[30];
	//double B_taugran=input_params[31];
	//double C_taugran=input_params[32];
	//double P=input_params[33];
	//double N0=input_params[34];
	//double Hfactor=input_params[35];
	//double Wfactor=input_params[36];

	double inc_rad, inc_star, inc_y, inc_c;
	double xmin, xmax, a1;
	double H, tau, p, N0;
	MatrixXd mode_params, noise_params(3,3);
	
	Cfg_synthetic_star cfg_star;
	Params_synthetic_star params;
	// ------------------------------------

	// Attempting to read file_cfg_mm, if exist. 
	std::cout << "Trying to find and read the mixed modes configuration file" << file_cfg_mm << std::endl;
	cfg_star=read_theoretical_freqs(file_cfg_mm, false); // try to read with critical = false
	if (cfg_star.use_nu_nl == true){
		std::cout << "    " << file_cfg_mm << " read! " << std::endl;
		std::cout << "    Found use_nu_nl = true ==> " << std::endl;
		std::cout << "    The frequencies of l=0,1,2,3 listed in the file will be used! " << std::endl;
	} else{
		std::cout << "     " << file_cfg_mm << " not found or cfg_star.use_nu_nl is set to false in the file. Pursuing..." << std::endl;
	}
	// ----------
	if (cfg_star.use_nu_nl == false){ // we use the parameters defined in the main.cfg
		cfg_star.Dnu_star=input_params[12];
		cfg_star.DPl_star=input_params[17];                
		cfg_star.q_star=input_params[19];
		cfg_star.alpha_g_star=input_params[18];
		cfg_star.epsilon_star=input_params[13];
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
	cfg_star.delta0l_percent_star=input_params[14];
	cfg_star.beta_p_star=input_params[15];

	cfg_star.maxHNR_l0=input_params[20];
	cfg_star.H0_spread=input_params[26];
	cfg_star.Gamma_max_l0=input_params[21];
	cfg_star.Hfactor=input_params[35];
	cfg_star.Wfactor=input_params[36];

	cfg_star.rot_env_input=input_params[0];
	cfg_star.rot_ratio_input=-1;
	cfg_star.rot_core_input=input_params[1];
	cfg_star.env_aspher.a2_l1=0; // SET TO 0 for l=1 modes. Ideally, would need a2_l1_core mixed with a2_l1_env
	cfg_star.env_aspher.a2_l2=input_params[4];
/*
	0 : nurot_env  
	1 : nurot_core  
	2 : a2_l1_core    
	3 : a2_l1_env    
	4 : a2_l2_env   
	5 : a2_l3_env   
	6 : a3_l2_env    
	7 : a3_l3_env     
	8 : a4_l2_env    
	9 : a4_l3_env    
	10: a5_l3_env    
	11: a6_l3_env     
	12: Dnu    
	13: epsilon    
	14: delta0l_percent  
	15: beta_p_star  
	16: nmax_spread  
	17: DP1    
	18: alpha     
	19: q       
	20: SNR    
	21: maxGamma   
	22: numax_spread        
	23: Vl1     
	24: Vl2    
	25: Vl3   
	26: H0_spread      
	27: A_Pgran   
	28: B_Pgran  
	29: C_Pgran    
	30: A_taugran   
	31: B_taugran    
	32: C_taugran    
	33: P      
	34: N0    
	35: Hfactor    
	36: Wfactor
*/
	if(input_params[5] <= -9999){ // If the user want l2 and l3 having the same aj coefficient, they need to put l3 to -9999 or smaller
		cfg_star.env_aspher.a2_l3=cfg_star.env_aspher.a2_l2;
	} else{
		cfg_star.env_aspher.a2_l3=input_params[5];
	}	
	cfg_star.env_aspher.a4_l2=input_params[8];
	if(input_params[9] <= -9999){ // If the user want l2 and l3 having the same aj coefficient, they need to put l3 to -9999 or smaller
		cfg_star.env_aspher.a4_l3=cfg_star.env_aspher.a4_l2;
	} else{
		cfg_star.env_aspher.a4_l3=input_params[9];
	}
	cfg_star.env_aspher.a6_l3=input_params[11];
	cfg_star.env_lat_dif_rot.a3_l2=input_params[6];
	if(input_params[7] <= -9999){ // If the user want l2 and l3 having the same aj coefficient, they need to put l3 to -9999 or smaller
		cfg_star.env_lat_dif_rot.a3_l3=cfg_star.env_lat_dif_rot.a3_l2;
	} else{
		cfg_star.env_lat_dif_rot.a3_l3=input_params[7];
	}
	cfg_star.env_lat_dif_rot.a5_l3=input_params[10];

	cfg_star.noise_params_harvey_like.resize(8);
	cfg_star.noise_params_harvey_like <<  input_params[27], input_params[28] , input_params[29] , input_params[30] , input_params[31] , input_params[32] , input_params[33]  , input_params[34];    //[A_Pgran ,  B_Pgran , C_Pgran   ,  A_taugran ,  B_taugran  , C_taugran    , p      N0]
	cfg_star.numax_star=numax_from_stello2009(cfg_star.Dnu_star, numax_spread); // Second argument is the random spread on numax
	cfg_star.fmin=cfg_star.numax_star -Nmax_pm*cfg_star.Dnu_star;
	cfg_star.fmax=cfg_star.numax_star +(Nmax_pm+2)*cfg_star.Dnu_star;
	cfg_star.output_file_rot=cpath + "/external/ARMM-solver/star_params.rot";
	cfg_star.Vl.resize(3);
	cfg_star.Vl << input_params[23], input_params[24], input_params[25];
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

	write_star_mode_params_aj(mode_params, file_out_modes);
	write_range_modes(cfg_star, params, file_range);

	if (cfg_star.Dnu_star <= 15){
		std::cout << "    Model with small Dnu ==> many mixed modes. This might be long to find the solutions..." << std::endl; 
	}

	tau=input_params[30] * pow(cfg_star.numax_star*1e-6,input_params[31]) + input_params[32]; // Granulation timescale (in seconds)
	H=input_params[27] * pow(cfg_star.numax_star*1e-6,input_params[28]) + input_params[29]; // Granulation Amplitude
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
	//exit(EXIT_SUCCESS);
}


void asymptotic_mm_freeDp_numaxspread_curvepmodes_v3_GRANscaled_Kallinger2014(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path, std::string template_file){

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
	double nmax_spread=input_params[16];
	double numax_spread=input_params[22]/100.;	

	double inc_rad, inc_star, inc_y, inc_c;
	double xmin, xmax, a1;
	double H, tau, p, N0;
	MatrixXd mode_params, noise_params(4,3);
	
	Cfg_synthetic_star cfg_star;
	Params_synthetic_star params;
	// ------------------------------------

	// Attempting to read file_cfg_mm, if exist. 
	std::cout << "Trying to find and read the mixed modes configuration file" << file_cfg_mm << std::endl;
	cfg_star=read_theoretical_freqs(file_cfg_mm, false); // try to read with critical = false
	if (cfg_star.use_nu_nl == true){
		std::cout << "    " << file_cfg_mm << " read! " << std::endl;
		std::cout << "    Found use_nu_nl = true ==> " << std::endl;
		std::cout << "    The frequencies of l=0,1,2,3 listed in the file will be used! " << std::endl;
	} else{
		std::cout << "     " << file_cfg_mm << " not found or cfg_star.use_nu_nl is set to false in the file. Pursuing..." << std::endl;
	}
	// ----------
	if (cfg_star.use_nu_nl == false){ // we use the parameters defined in the main.cfg
		cfg_star.Dnu_star=input_params[12];
		cfg_star.DPl_star=input_params[17];                
		cfg_star.q_star=input_params[19];
		cfg_star.alpha_g_star=input_params[18];
		cfg_star.epsilon_star=input_params[13];
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
	cfg_star.delta0l_percent_star=input_params[14];
	cfg_star.beta_p_star=input_params[15];

	cfg_star.maxHNR_l0=input_params[20];
	cfg_star.H0_spread=input_params[26];
	cfg_star.Gamma_max_l0=input_params[21];
	cfg_star.Hfactor=input_params[41];
	cfg_star.Wfactor=input_params[42];

	cfg_star.rot_env_input=input_params[0];
	cfg_star.rot_ratio_input=-1;
	cfg_star.rot_core_input=input_params[1];
	cfg_star.env_aspher.a2_l1=0; // SET TO 0 for l=1 modes. Ideally, would need a2_l1_core mixed with a2_l1_env
	cfg_star.env_aspher.a2_l2=input_params[4];
/*
	0 : nurot_env  
	1 : nurot_core  
	2 : a2_l1_core    
	3 : a2_l1_env    
	4 : a2_l2_env   
	5 : a2_l3_env   
	6 : a3_l2_env    
	7 : a3_l3_env     
	8 : a4_l2_env    
	9 : a4_l3_env    
	10: a5_l3_env    
	11: a6_l3_env     
	12: Dnu    
	13: epsilon    
	14: delta0l_percent  
	15: beta_p_star  
	16: nmax_spread  
	17: DP1    
	18: alpha     
	19: q       
	20: SNR    
	21: maxGamma   
	22: numax_spread        
	23: Vl1     
	24: Vl2    
	25: Vl3   
	26: H0_spread      
	---> 14 Noise parameters from Kallinger2014 follow after
	41: Hfactor
	42: Wfactor
	43: Tobs
	44: Cadence
*/

	if(input_params[5] <= -9999){ // If the user want l2 and l3 having the same aj coefficient, they need to put l3 to -9999 or smaller
		cfg_star.env_aspher.a2_l3=cfg_star.env_aspher.a2_l2;
	} else{
		cfg_star.env_aspher.a2_l3=input_params[5];
	}	
	cfg_star.env_aspher.a4_l2=input_params[8];
	if(input_params[9] <= -9999){ // If the user want l2 and l3 having the same aj coefficient, they need to put l3 to -9999 or smaller
		cfg_star.env_aspher.a4_l3=cfg_star.env_aspher.a4_l2;
	} else{
		cfg_star.env_aspher.a4_l3=input_params[9];
	}
	cfg_star.env_aspher.a6_l3=input_params[11];
	cfg_star.env_lat_dif_rot.a3_l2=input_params[6];
	if(input_params[7] <= -9999){ // If the user want l2 and l3 having the same aj coefficient, they need to put l3 to -9999 or smaller
		cfg_star.env_lat_dif_rot.a3_l3=cfg_star.env_lat_dif_rot.a3_l2;
	} else{
		cfg_star.env_lat_dif_rot.a3_l3=input_params[7];
	}
	cfg_star.env_lat_dif_rot.a5_l3=input_params[10];
	// This must precede the noise definition
	cfg_star.numax_star=numax_from_stello2009(cfg_star.Dnu_star, numax_spread); // Second argument is the random spread on numax
	const double numax_star_nospread=numax_from_stello2009(cfg_star.Dnu_star, 0); // Used by Kallinger2014 relations
	//
	// --- Deploy noise parameters ----
	// Raw noise params start at index=16
	//Raw noise params: k_Agran         s_Agran         k_taugran       s_taugran       c0              ka              ks              k1              s1              c1              k2              s2              c2              N0
	//Need conversion to: k_Agran         s_Agran         k_taugran       s_taugran  c0        Na1      Na2     k1     s1    c1      k2    s2      c1
	//a1,a2 but be derived first here....
	const double ka=input_params[32];
	const double ks=input_params[33];
	const double Na1=ka*std::pow(cfg_star.numax_star, ks);
	const double Na2=Na1;
	VectorXd noise_params_Kallinger(14);
	noise_params_Kallinger.segment(0,6)=input_params.segment(27,6);
	noise_params_Kallinger[5]=Na1;
	noise_params_Kallinger[6]=Na2;
	noise_params_Kallinger.segment(7,7)=input_params.segment(27+6+2-1,7);
	//std::cout << "ka = " << ka << std::endl;
	//std::cout << "ks = " << ks << std::endl;
	//for(int k=0;k<noise_params_Kallinger.size();k++){
	//	std::cout << "noise_params_Kallinger[" << k << "] = " << noise_params_Kallinger[k] << std::endl;
	//}
	//std::cout << "Verify that all these affectation are correct... " << std::endl;
	//exit(EXIT_SUCCESS);
	// Defining the new noise harvey parameters. These will be used to compute the 
	// local noise of the new star once we have the rescaled frequencies
	// The noise params here follows Kallinger et al. 2014
	VectorXd x;
	const double Tobs=input_params[43];
	const double Cadence=input_params[44];
	const double df=1e6/(Tobs * 86400.);
	const double Delta=1e6/Cadence; // /2
	const int Ndata=Delta/df;
	//std::cout << "Tobs=input_params[41] = " << Tobs << std::endl;
	//std::cout << "Cadence=input_params[42] = " << Cadence << std::endl;
	//std::cout << "Delta = " << Delta << std::endl;
	x.setLinSpaced(Ndata, 0, Delta);
	const double mu_numax=0; // This parameter is redundant with numax_spread as this already introduce some jitter in the mode position. 
	const VectorXd noise_params_harvey=Kallinger2014_to_harveylike(numax_star_nospread, mu_numax, noise_params_Kallinger, x);
	cfg_star.noise_params_harvey_like=noise_params_harvey; 
	// ------------
	// ------------ 
	cfg_star.fmin=cfg_star.numax_star -Nmax_pm*cfg_star.Dnu_star;
	cfg_star.fmax=cfg_star.numax_star +(Nmax_pm+2)*cfg_star.Dnu_star;
	cfg_star.output_file_rot=cpath + "/external/ARMM-solver/star_params.rot";
	cfg_star.Vl.resize(3);
	cfg_star.Vl << input_params[23], input_params[24], input_params[25];
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

	write_star_mode_params_aj(mode_params, file_out_modes);
	write_range_modes(cfg_star, params, file_range);

	noise_params(0,0)=noise_params_harvey[0];
	noise_params(0,1)=noise_params_harvey[1];
	noise_params(0,2)=noise_params_harvey[2];
	noise_params(1,0)=noise_params_harvey[3];
	noise_params(1,1)=noise_params_harvey[4];
	noise_params(1,2)=noise_params_harvey[5]; 
	noise_params(2, 0)=noise_params_harvey[6];
	noise_params(2, 1)=noise_params_harvey[7];
	noise_params(2, 2)=noise_params_harvey[8];
	noise_params(3, 0)=noise_params_harvey[9]; // White noise
	noise_params(3, 1)=-2;
	noise_params(3, 2)=-2;
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

	std::random_device rd;
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<double> distrib(0 , 1);

	int i;
	double HNR, a1_ov_Gamma, a2, a3, a4, a5, a6, beta_asym, inc, 
			HNRmaxref, Height_factor, Gamma_at_numax, a1, Gamma_coef, fl,rho, eta0,
			Dnu, epsilon, delta0l_percent, numax_spread; 
	double xmin, xmax;
	VectorXi pos;
	VectorXd tmp, xfit, rfit, d0l(3), f_rescaled_lin;
	VectorXd HNRref, local_noise,h_star, gamma_star, s_a1_star, s_a2_star, s_a3_star, s_a4_star,s_a5_star,s_a6_star,
		s_eta0_star, s_epsilon_star, s_theta0_star, s_delta_star, s_asym_star, inc_star, activity_terms;
	Star_params ref_star;
	Freq_modes f_rescaled, f_ref;
	MatrixXd mode_params, noise_params;

	ref_star=read_star_params(extra); 
	mode_params.setZero(ref_star.mode_params.rows(), 16); // Final table of values that define a star

	// Defining the noise profile parameters
	// Note about the noise: -1 means that it is ignored. -2 mean that the value is irrelevant
	tmp.resize(ref_star.mode_params.rows()); 	
	tmp.setConstant(0);
	local_noise=harvey_like(ref_star.noise_params, ref_star.mode_params.col(1), tmp); // Generate a list of local noise values for each frequencies
// ---- Deploy the parameters -----
	Dnu=input_params[0];
	epsilon=input_params[1];
	delta0l_percent=input_params[2];
	HNR=input_params[3];
	a1_ov_Gamma=input_params[4];
	Gamma_at_numax=input_params[5];
	a2=input_params[6];
	a3=input_params[7];
	a4=input_params[8];
	a5=input_params[9];
	a6=input_params[10];
	beta_asym=input_params[11];
	inc=input_params[12];
	const double H_spread=input_params[13];
	const double nu_spread=input_params[14];
	const double gamma_spread=input_params[15];
	const double do_flat_noise=input_params[16]; // if <= 0 : Do not make flat noise. Any positive value will make noise. Flat noise may not be be fix

// ---------------------------------
	// --- Perform rescaling of frequencies  ---
	d0l << delta0l_percent*Dnu/100., delta0l_percent*Dnu/100., delta0l_percent*Dnu/100.; // small separation l=1,2,3
	pos=where_dbl(ref_star.mode_params.col(0), 0, 1e-3); // pick l=0
	if(pos[0] != -1){
		f_ref.fl0.resize(pos.size());
		for(int n=0; n<pos.size();n++){
			f_ref.fl0[n]=ref_star.mode_params(pos[n],1);
		}
	} else{
		std::cerr << "Error in generate_cfg_from_synthese_file_Wscaled_aj : You must have at least l=0 to perform a rescale" << std::endl;
		exit(EXIT_FAILURE);
	}
	if(pos[0] != -1){
		pos=where_dbl(ref_star.mode_params.col(0), 1, 1e-3); // pick l=1
		f_ref.fl1.resize(pos.size());
		for(int n=0; n<pos.size();n++){
			f_ref.fl1[n]=ref_star.mode_params(pos[n],1);
		}
	}
	pos=where_dbl(ref_star.mode_params.col(0), 2, 1e-3); // pick l=2
	if(pos[0] != -1){
		f_ref.fl2.resize(pos.size());
		for(int n=0; n<pos.size();n++){
			f_ref.fl2[n]=ref_star.mode_params(pos[n],1);
		}
	}
	pos=where_dbl(ref_star.mode_params.col(0), 3, 1e-3); // pick l=3
	if(pos[0] != -1){
		f_ref.fl3.resize(pos.size());
		for(int n=0; n<pos.size();n++){
			f_ref.fl3[n]=ref_star.mode_params(pos[n],1);
		}
	}

	f_rescaled=rescale_freqs(Dnu, epsilon, f_ref, d0l);
	if (f_rescaled.error_status == true){
		std::cerr << "Error while rescaling: There is likely an issue in frequency tagging." << std::endl;
		std::cerr << "                       Debug in generate_cfg_from_synthese_file_Wscaled_aj required " << std::endl;
		exit(EXIT_FAILURE);
	}
	f_rescaled_lin.resize(f_ref.fl0.size()+ f_ref.fl1.size()+f_ref.fl2.size()+ f_ref.fl3.size());
	for(int n=0; n<f_ref.fl0.size(); n++){
		mode_params(n,0)=0;
		f_rescaled_lin[n]=f_rescaled.fl0[n];
		// Adding to f_rescaled_lin[n] a uniform random quantity that is bounded by xmin=f_rescaled_lin[n]*(1 - nu_spread) and xmax=f_rescaled_lin[n]*(1+nu_spread)
		if (nu_spread > 0)
		{
				xmin=f_rescaled_lin[n]*(1. - nu_spread/100.);
				xmax=f_rescaled_lin[n]*(1. + nu_spread/100.);
				f_rescaled_lin[n]=xmin + (xmax-xmin)*distrib(gen);
		}
	}
	for(int n=0; n<f_ref.fl1.size(); n++){
		mode_params(n+f_ref.fl0.size(),0)=1;
		f_rescaled_lin[n+f_ref.fl0.size()]=f_rescaled.fl1[n];
		if (nu_spread > 0)
		{
				xmin=f_rescaled_lin[n+f_ref.fl0.size()]*(1. - nu_spread/100.);
				xmax=f_rescaled_lin[n+f_ref.fl0.size()]*(1. + nu_spread/100.);
				f_rescaled_lin[n+f_ref.fl0.size()]=xmin + (xmax-xmin)*distrib(gen);
		}
	}
	for(int n=0; n<f_ref.fl2.size(); n++){
		mode_params(n+f_ref.fl0.size()+f_ref.fl1.size(),0)=2;
		f_rescaled_lin[n+f_ref.fl0.size()+f_ref.fl1.size()]=f_rescaled.fl2[n];
		if (nu_spread > 0)
		{
				xmin=f_rescaled_lin[n+f_ref.fl0.size()+f_ref.fl1.size()]*(1. - nu_spread/100.);
				xmax=f_rescaled_lin[n+f_ref.fl0.size()+f_ref.fl1.size()]*(1. + nu_spread/100.);
				f_rescaled_lin[n+f_ref.fl0.size()+f_ref.fl1.size()]=xmin + (xmax-xmin)*distrib(gen);
		}
	}
	for(int n=0; n<f_ref.fl3.size(); n++){
		mode_params(n+f_ref.fl0.size()+f_ref.fl1.size()+f_ref.fl2.size(),0)=3;
		f_rescaled_lin[n+f_ref.fl0.size()+f_ref.fl1.size()+f_ref.fl2.size()]=f_rescaled.fl3[n];
		if (nu_spread > 0)
		{
				xmin=f_rescaled_lin[n+f_ref.fl0.size()+f_ref.fl1.size()+f_ref.fl2.size()]*(1. - nu_spread/100.);
				xmax=f_rescaled_lin[n+f_ref.fl0.size()+f_ref.fl1.size()+f_ref.fl2.size()]*(1. + nu_spread/100.);
				f_rescaled_lin[n+f_ref.fl0.size()+f_ref.fl1.size()+f_ref.fl2.size()]=xmin + (xmax-xmin)*distrib(gen);
		}
	}

	// ---- Rescaling Height and Width profiles ----
	HNRref=ref_star.mode_params.col(2);
	HNRref=HNRref.cwiseProduct(local_noise.cwiseInverse());
	pos=where_dbl(ref_star.mode_params.col(0), 0, 1e-3);
	//std::cout << "ref_star.mode_params.col(0) =" << ref_star.mode_params.col(0) << std::endl;
	//std::cout <<" pos=" << pos << std::endl;
	HNRmaxref=0;
	for (int n=0; n<pos.size();n++){
		if (HNRmaxref < HNRref[pos[n]]){
			HNRmaxref=HNRref[pos[n]];
		}
	}
	//HNRmaxref=HNRref.maxCoeff(); // This is the maximum HNR of the reference data
	Height_factor=HNR/HNRmaxref;  // compute the coeficient required to get the proper max(HNR)
	pos=where_dbl(HNRref, HNRmaxref, 0.001);
	if (pos[0] >= 0){
		Gamma_coef=Gamma_at_numax/ref_star.mode_params(pos[0], 3); // Correction coeficient to apply on Gamma(nu) in order to ensure that we have Gamma(nu=numax) = Gamma_at_numax
	} else{
		std::cout << "Error! could not find the max position for the mode Widths profile" << std::endl;
		std::cout << "Code debug required" << std::endl;
		exit(EXIT_FAILURE);
	}

	a1=a1_ov_Gamma*Gamma_at_numax; // We can vary the Width and splitting. But we need to change the splitting in order to get the wished a1/Gamma0

	// Defining the final size for all of the outptus
	gamma_star.resize(ref_star.mode_params.rows());
	gamma_star=Gamma_coef*ref_star.mode_params.col(3); 
	if (gamma_spread > 0)
	{
		xmin=gamma_star[i]*(1. - gamma_spread/100.);
		xmax=gamma_star[i]*(1. + gamma_spread/100.);
		gamma_star[i]=xmin + (xmax-xmin)*distrib(gen);
	}
	
	// Refactoring the heights
	h_star.resize(ref_star.mode_params.rows());
	//h_star=Height_factor * HNRref * N0; 
	if (local_noise.sum()!= 0){
		//N0=1; // Imposing the the Noise background is 1
		//std::cout <<  "Using N0=" << N0 << " (white noise)" << std::endl;
		std::cout << "Using the Noise profile of the Reference star" << std::endl;
		std::cout << "HNR of all modes (degree / freq_template  / freq_rescaled / HNR  /  Height /  local_noise):" << std::endl;
		for(i =0; i<ref_star.mode_params.rows(); i++){
			if (do_flat_noise <=0){
				h_star[i]=Height_factor * HNRref[i] * local_noise[i];
			} else{
				h_star[i]=Height_factor * HNRref[i] * do_flat_noise; // Here do_flat_noise will contain the wished N0 value
			}
			// Adding a unifomly random value, with maximum range 2*H_spread
			if (H_spread > 0)
			{
					xmin=h_star[i]*(1. - H_spread/100.);
					xmax=h_star[i]*(1. + H_spread/100.);
					h_star[i]=xmin + (xmax-xmin)*distrib(gen);
			}
			std::cout << "     " << ref_star.mode_params(i,0) << "  " << ref_star.mode_params(i,1)  << "  "  << f_rescaled_lin[i] << "  " << Height_factor * HNRref[i]  << "  " << h_star[i] << "  "  << local_noise[i] << std::endl;
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


	// Filling the output Matrix. 
	//   Note on assumptions here: We assume that the template has monotonically increasing frequency, for each l. 
	//                             It is also assumed that block of l are in the strict order l=0,1,2,3

	//mode_params.col(0)=ref_star.mode_params.col(0); // List of els
	//mode_params.col(1)=ref_star.mode_params.col(1); // List of frequencies
	mode_params.col(1)=f_rescaled_lin; // Frequencies in a flat array
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
	if (do_flat_noise <=0){
		write_star_noise_params(ref_star.noise_params, file_out_noise);
	} else{
		MatrixXd noise_params(2,3);
		noise_params(0,0)=0;
		noise_params(0,1)=1;
		noise_params(0,2)=1;
		noise_params(1,0)=do_flat_noise;
		noise_params(1,1)=-2;
		noise_params(1,2)=-2;
		write_star_noise_params(noise_params, file_out_noise);
	}

}


/* 
	This function use a reference star as a template to generate frequencies and Width, Height profiles
	can be rescaled so that you can modify the HNR but keep the same height profile
	Note that the user here provides a target a1/Width so that a1 is automatically adjusted to match the 
	requested a1/Width. The code will not change the Width so that code is not adapted to test blending between adjacent l modes,
	such as the l=0 and l=2 mode blending
	It handles aj coeficient up to j=6 as well as free parameters and consider an activity term Alm ~ {a2, a4, a6, ...} following Gizon2002 idea.
	a3 can be given as a polynomial of the frequency and you can add a random quantity (set in nHz) that will be added as a random 'error' ==> scatter ==> test robustness 
	The noise background is scaled using numax, itself derived from the relation between Dnu and numax modulo a random term to add some dispersion.
    It Implement A Single Harvey profile noise which scales with numax (+ a White Noise). The definition of the Harvey parameters are as defined by Karoff et al. 2010
                Recommended coefficients for the scaling are Pgran = A numax^B + C with A=10^-4 and B=-2, C=0. t_gran = A numax^B + C with A=1 and B=-1 and C=0
				The formulation is also compatible with Kallinger et al. 2014.
*/
void generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string extra){

	std::random_device rd;
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<double> distrib(0 , 1);

	int i;
	double HNR, a1_ov_Gamma, a2, a3, a4, a5, a6, beta_asym, inc, 
			HNRmaxref, Height_factor, Gamma_at_numax, a1, Gamma_coef, fl,rho, eta0,
			Dnu, epsilon, delta0l_percent, numax_star, numax_spread; 
	double xmin, xmax;
	VectorXi pos;
	VectorXd tmp, xfit, rfit, d0l(3), f_rescaled_lin;
	VectorXd HNRref, local_noise,local_noise_new, h_star, gamma_star, s_a1_star, s_a2_star, s_a3_star, s_a4_star,s_a5_star,s_a6_star,
		s_eta0_star, s_epsilon_star, s_theta0_star, s_delta_star, s_asym_star, inc_star, noise_params_harvey1985(4);
	Star_params ref_star;
	Freq_modes f_rescaled, f_ref;
	MatrixXd mode_params, noise_params(3,3);

// ---- Deploy the parameters -----
	Dnu=input_params[0];
	epsilon=input_params[1];
	delta0l_percent=input_params[2];
	HNR=input_params[3];
	a1_ov_Gamma=input_params[4];
	Gamma_at_numax=input_params[5];
	a2=input_params[6];
	a3=input_params[7];
	a4=input_params[8];
	a5=input_params[9];
	a6=input_params[10];
	beta_asym=input_params[11];
	inc=input_params[12];
	numax_spread=input_params[21];
	const double H_spread=input_params[22];
	const double nu_spread=input_params[23];
// ---------------------------------

	ref_star=read_star_params(extra); 
	mode_params.setZero(ref_star.mode_params.rows(), 16); // Final table of values that define a star

	// Defining the noise profile parameters
	// Note about the noise: -1 means that it is ignored. -2 mean that the value is irrelevant
	tmp.resize(ref_star.mode_params.rows()); 	
	tmp.setConstant(0);
	local_noise=harvey_like(ref_star.noise_params, ref_star.mode_params.col(1), tmp); // Generate a list of local noise values for each frequencies

	// Defining the new noise harvey parameters. These will be used to compute the 
	// local noise of the new star once we have the rescaled frequencies
	numax_star=numax_from_stello2009(Dnu, numax_spread); // Second argument is the random spread on numax
	std::cout << "   - numax = " << numax_star << std::endl;
	// The noise params here follow Karoff 2010. The formulation is also compatible with Kallinger et al. 2014
	//noise_params_harvey_like=[A_Pgran ,  B_Pgran , C_Pgran   ,  A_taugran ,  B_taugran  , C_taugran    , p      N0] // 
	noise_params_harvey1985[0] = input_params[13] * std::pow(numax_star*1e-6,input_params[14]) + input_params[15]; // Granulation Amplitude
	noise_params_harvey1985[1] = input_params[16] * std::pow(numax_star*1e-6, input_params[17]) + input_params[18]; // Granulation timescale (in seconds)
	noise_params_harvey1985[0] = noise_params_harvey1985[0]/noise_params_harvey1985[1];
	noise_params_harvey1985[1]= noise_params_harvey1985[1]/1000.; // timescale
	noise_params_harvey1985[2]=input_params[19]; // slope
	noise_params_harvey1985[3]=input_params[20]; // white noise

	// --- Perform rescaling of frequencies  ---
	d0l << delta0l_percent*Dnu/100., delta0l_percent*Dnu/100., delta0l_percent*Dnu/100.; // small separation l=1,2,3
	pos=where_dbl(ref_star.mode_params.col(0), 0, 1e-3); // pick l=0
	if(pos[0] != -1){
		f_ref.fl0.resize(pos.size());
		for(int n=0; n<pos.size();n++){
			f_ref.fl0[n]=ref_star.mode_params(pos[n],1);
		}
	} else{
		std::cerr << "Error in generate_cfg_from_synthese_file_Wscaled_aj : You must have at least l=0 to perform a rescale" << std::endl;
		exit(EXIT_FAILURE);
	}
	if(pos[0] != -1){
		pos=where_dbl(ref_star.mode_params.col(0), 1, 1e-3); // pick l=1
		f_ref.fl1.resize(pos.size());
		for(int n=0; n<pos.size();n++){
			f_ref.fl1[n]=ref_star.mode_params(pos[n],1);
		}
	}
	pos=where_dbl(ref_star.mode_params.col(0), 2, 1e-3); // pick l=2
	if(pos[0] != -1){
		f_ref.fl2.resize(pos.size());
		for(int n=0; n<pos.size();n++){
			f_ref.fl2[n]=ref_star.mode_params(pos[n],1);
		}
	}
	pos=where_dbl(ref_star.mode_params.col(0), 3, 1e-3); // pick l=3
	if(pos[0] != -1){
		f_ref.fl3.resize(pos.size());
		for(int n=0; n<pos.size();n++){
			f_ref.fl3[n]=ref_star.mode_params(pos[n],1);
		}
	}

	f_rescaled=rescale_freqs(Dnu, epsilon, f_ref, d0l);
	if (f_rescaled.error_status == true){
		std::cerr << "Error while rescaling: There is likely an issue in frequency tagging." << std::endl;
		std::cerr << "                       Debug in generate_cfg_from_synthese_file_Wscaled_aj required " << std::endl;
		exit(EXIT_FAILURE);
	}
	f_rescaled_lin.resize(f_ref.fl0.size()+ f_ref.fl1.size()+f_ref.fl2.size()+ f_ref.fl3.size());
	xmin=-nu_spread;
	xmax=nu_spread;
	for(int n=0; n<f_ref.fl0.size(); n++){
		mode_params(n,0)=0;
		f_rescaled_lin[n]=f_rescaled.fl0[n];
		// Adding to f_rescaled_lin[n] a uniform random quantity that is bounded by xmi and xmax
		if (nu_spread > 0)
		{
				f_rescaled_lin[n]=f_rescaled_lin[n] + (xmax-xmin)*distrib(gen);
		}
	}
	for(int n=0; n<f_ref.fl1.size(); n++){
		mode_params(n+f_ref.fl0.size(),0)=1;
		f_rescaled_lin[n+f_ref.fl0.size()]=f_rescaled.fl1[n];
		if (nu_spread > 0)
		{
				f_rescaled_lin[n+f_ref.fl0.size()]=f_rescaled_lin[n+f_ref.fl0.size()] + (xmax-xmin)*distrib(gen);
		}
	}
	for(int n=0; n<f_ref.fl2.size(); n++){
		mode_params(n+f_ref.fl0.size()+f_ref.fl1.size(),0)=2;
		f_rescaled_lin[n+f_ref.fl0.size()+f_ref.fl1.size()]=f_rescaled.fl2[n];
		if (nu_spread > 0)
		{
				f_rescaled_lin[n+f_ref.fl0.size()+f_ref.fl1.size()]=f_rescaled_lin[n+f_ref.fl0.size()+f_ref.fl1.size()] + (xmax-xmin)*distrib(gen);
		}
	}
	for(int n=0; n<f_ref.fl3.size(); n++){
		mode_params(n+f_ref.fl0.size()+f_ref.fl1.size()+f_ref.fl2.size(),0)=3;
		f_rescaled_lin[n+f_ref.fl0.size()+f_ref.fl1.size()+f_ref.fl2.size()]=f_rescaled.fl3[n];
		if (nu_spread > 0)
		{
				f_rescaled_lin[n+f_ref.fl0.size()+f_ref.fl1.size()+f_ref.fl2.size()]=f_rescaled_lin[n+f_ref.fl0.size()+f_ref.fl1.size()+f_ref.fl2.size()] + (xmax-xmin)*distrib(gen);
		}
	}

	// Compute the local noise for the new star
	local_noise_new.resize(f_rescaled_lin.size());
	local_noise_new.setZero();
	local_noise_new=harvey1985(noise_params_harvey1985, f_rescaled_lin, local_noise_new, 1); // Iterate on Noise_l0 to update it by putting the noise profile with one harvey profile

	// ---- Rescaling Height and Width profiles ----
	HNRref=ref_star.mode_params.col(2);
	HNRref=HNRref.cwiseProduct(local_noise.cwiseInverse());
	pos=where_dbl(ref_star.mode_params.col(0), 0, 1e-3);

	HNRmaxref=0;
	for (int n=0; n<pos.size();n++){
		if (HNRmaxref < HNRref[pos[n]]){
			HNRmaxref=HNRref[pos[n]];
		}
	}

	//HNRmaxref=HNRref.maxCoeff(); // This is the maximum HNR of the reference data
	Height_factor=HNR/HNRmaxref;  // compute the coeficient required to get the proper max(HNR)
	pos=where_dbl(HNRref, HNRmaxref, 0.001);
	if (pos[0] >= 0){
		Gamma_coef=Gamma_at_numax/ref_star.mode_params(pos[0], 3); // Correction coeficient to apply on Gamma(nu) in order to ensure that we have Gamma(nu=numax) = Gamma_at_numax
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
	if (local_noise_new.sum()!= 0){
		//N0=1; // Imposing the the Noise background is 1
		//std::cout <<  "Using N0=" << N0 << " (white noise)" << std::endl;
		std::cout << "Using the Noise profile of the New star, using provided noise parameters and Karoff et al. 2010 prescription" << std::endl;
		std::cout << "HNR of all modes (degree / freq_template  / freq_rescaled / HNR  /  Height /  local_noise(noise params)):" << std::endl;
		for(i =0; i<ref_star.mode_params.rows(); i++){
			h_star[i]=Height_factor * HNRref[i] * local_noise_new[i];
			// Adding a unifomly random value, with maximum range 2*H_spread
			if (H_spread > 0)
			{
					xmin=h_star[i]*(1. - H_spread/100.);
					xmax=h_star[i]*(1. + H_spread/100.);
					h_star[i]=xmin + (xmax-xmin)*distrib(gen);
			}
			std::cout << "     " << ref_star.mode_params(i,0) << "  " << ref_star.mode_params(i,1)  << "  "  << f_rescaled_lin[i] << "  " << Height_factor * HNRref[i]  << "  " << h_star[i] << "  "  << local_noise_new[i] << std::endl;
		}
	} else{
		std::cerr << "Error: local_noise_new has no valid value. Debug required." << std::endl;
		std::cerr << "         The program will stop now" << std::endl;
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
	
	// Filling the output Matrix. 
	//   Note on assumptions here: We assume that the template has monotonically increasing frequency, for each l. 
	//                             It is also assumed that block of l are in the strict order l=0,1,2,3

	//mode_params.col(0)=ref_star.mode_params.col(0); // List of els
	//mode_params.col(1)=ref_star.mode_params.col(1); // List of frequencies
	mode_params.col(1)=f_rescaled_lin; // Frequencies in a flat array
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

	noise_params(0,0)=-1;
	noise_params(0,1)=-1;
	noise_params(0,2)=-1; 
	noise_params(1,0)=noise_params_harvey1985[0];
	noise_params(1,1)=noise_params_harvey1985[1];
	noise_params(1,2)=noise_params_harvey1985[2]; 
	noise_params(2, 0)=noise_params_harvey1985[3]; // White noise
	noise_params(2, 1)=-2;
	noise_params(2, 2)=-2;

	// A FUNCTION THAT WRITES THE Noise
	write_star_noise_params(noise_params, file_out_noise);
	//exit(EXIT_SUCCESS);
}


/* 
	This function use a reference star as a template to generate frequencies and Width, Height profiles
	can be rescaled so that you can modify the HNR but keep the same height profile
	Note that the user here provides a target a1/Width so that a1 is automatically adjusted to match the 
	requested a1/Width. The code will not change the Width so that code is not adapted to test blending between adjacent l modes,
	such as the l=0 and l=2 mode blending
	It handles aj coeficient up to j=6 as well as free parameters and consider an activity term Alm ~ {a2, a4, a6, ...} following Gizon2002 idea.
	a3 can be given as a polynomial of the frequency and you can add a random quantity (set in nHz) that will be added as a random 'error' ==> scatter ==> test robustness 
	The noise background is scaled using numax, itself derived from the relation between Dnu and numax modulo a random term to add some dispersion.
    It Implement closely the Kallinger+2014 noise model, composed of 3 Harvey-like profiles. See Table 2 of Kallinger et al. 2014 and its implementation in noise_models.cpp.
*/
void generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled_Kallinger2014(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string extra){

	std::random_device rd;
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<double> distrib(0 , 1);

	int i;
	double HNR, a1_ov_Gamma, a2, a3, a4, a5, a6, beta_asym, inc, 
			HNRmaxref, Height_factor, Gamma_at_numax, a1, Gamma_coef, fl,rho, eta0,
			Dnu, epsilon, delta0l_percent, numax_star, numax_spread; 
	double xmin, xmax;
	VectorXi pos;
	VectorXd tmp, xfit, rfit, d0l(3), f_rescaled_lin;
	VectorXd HNRref, local_noise,local_noise_new, h_star, gamma_star, s_a1_star, s_a2_star, s_a3_star, s_a4_star,s_a5_star,s_a6_star,
		s_eta0_star, s_epsilon_star, s_theta0_star, s_delta_star, s_asym_star, inc_star, noise_params_harvey1985(4);
	Star_params ref_star;
	Freq_modes f_rescaled, f_ref;
	MatrixXd mode_params, noise_params(4,3);

// ---- Deploy the mode parameters -----
	Dnu=input_params[0];
	epsilon=input_params[1];
	delta0l_percent=input_params[2];
	HNR=input_params[3];
	a1_ov_Gamma=input_params[4];
	Gamma_at_numax=input_params[5];
	a2=input_params[6];
	a3=input_params[7];
	a4=input_params[8];
	a5=input_params[9];
	a6=input_params[10];
	beta_asym=input_params[11];
	inc=input_params[12];
	numax_spread=input_params[13];
	const double H_spread=input_params[14];
	const double nu_spread=input_params[15];
	
	numax_star=numax_from_stello2009(Dnu, numax_spread); // Second argument is the random spread on numax
	const double numax_star_nospread=numax_from_stello2009(Dnu, 0); // Used by Kallinger2014 relations
	std::cout << "   - numax = " << numax_star << std::endl;
// ---------------------------------
// --- Deploy noise parameters ----
	// Raw noise params start at index=16
	//Raw noise params: k_Agran         s_Agran         k_taugran       s_taugran       c0              ka              ks              k1              s1              c1              k2              s2              c2              N0
	//Need conversion to: k_Agran         s_Agran         k_taugran       s_taugran  c0        Na1      Na2     k1     s1    c1      k2    s2      c1
	//a1,a2 but be derived first here....
	const double ka=input_params[21];
	const double ks=input_params[22];
	const double Na1=ka*std::pow(numax_star, ks);
	const double Na2=Na1;
	VectorXd noise_params_Kallinger(14);
	noise_params_Kallinger.segment(0,6)=input_params.segment(16,6);
	noise_params_Kallinger[5]=Na1;
	noise_params_Kallinger[6]=Na2;
	noise_params_Kallinger.segment(7,7)=input_params.segment(16+6+2-1,7);
	//std::cout << "ka = " << ka << std::endl;
	//std::cout << "ks = " << ks << std::endl;
	//for(int k=0;k<noise_params_Kallinger.size();k++){
	//	std::cout << "noise_params_Kallinger[" << k << "] = " << noise_params_Kallinger[k] << std::endl;
	//}
	//std::cout << "Verify that all these affectation are correct... " << std::endl;
	//exit(EXIT_SUCCESS);

	ref_star=read_star_params(extra, false);  // no verbose here
	mode_params.setZero(ref_star.mode_params.rows(), 16); // Final table of values that define a star

	// Defining the noise profile parameters
	// Note about the noise: -1 means that it is ignored. -2 mean that the value is irrelevant
	tmp.resize(ref_star.mode_params.rows()); 	
	tmp.setConstant(0);
	if (ref_star.noise_model == "Harvey-like"){
		local_noise=harvey_like(ref_star.noise_params, ref_star.mode_params.col(1), tmp); // Generate a list of local noise values for each frequencies
	} else{
		std::cerr << "Noise model not supported for the reference star. " << std::endl;
		std::cerr << "Debug required in models_database.cpp...." << std::endl;
		exit(EXIT_FAILURE);
	}
	// Defining the new noise harvey parameters. These will be used to compute the 
	// local noise of the new star once we have the rescaled frequencies
	// The noise params here follows Kallinger et al. 2014
	VectorXd x;
	const double Tobs=input_params[30];
	const double Cadence=input_params[31];
	const double df=1e6/(Tobs * 86400.);
	const double Delta=1e6/Cadence; // /2
	const int Ndata=Delta/df;
	//std::cout << "Tobs=input_params[30] = " << Tobs << std::endl;
	//std::cout << "Cadence=input_params[31] = " << Cadence << std::endl;
	//std::cout << "Delta = " << Delta << std::endl;
	x.setLinSpaced(Ndata, 0, Delta);
	const double mu_numax=0; // This parameter is redundant with numax_spread as this already introduce some jitter in the mode position. 
	const VectorXd noise_params_harvey=Kallinger2014_to_harveylike(numax_star_nospread, mu_numax, noise_params_Kallinger, x);
	//std::cout << noise_params_harvey.transpose() << std::endl;
	//std::cout << "models_database.cpp: Check here that noise_params_harvey is with correct units..." << std::endl;
	//exit(EXIT_SUCCESS);

	// --- Perform rescaling of frequencies  ---
	d0l << delta0l_percent*Dnu/100., delta0l_percent*Dnu/100., delta0l_percent*Dnu/100.; // small separation l=1,2,3
	pos=where_dbl(ref_star.mode_params.col(0), 0, 1e-3); // pick l=0
	if(pos[0] != -1){
		f_ref.fl0.resize(pos.size());
		for(int n=0; n<pos.size();n++){
			f_ref.fl0[n]=ref_star.mode_params(pos[n],1);
		}
	} else{
		std::cerr << "Error in generate_cfg_from_synthese_file_Wscaled_aj : You must have at least l=0 to perform a rescale" << std::endl;
		exit(EXIT_FAILURE);
	}
	if(pos[0] != -1){
		pos=where_dbl(ref_star.mode_params.col(0), 1, 1e-3); // pick l=1
		f_ref.fl1.resize(pos.size());
		for(int n=0; n<pos.size();n++){
			f_ref.fl1[n]=ref_star.mode_params(pos[n],1);
		}
	}
	pos=where_dbl(ref_star.mode_params.col(0), 2, 1e-3); // pick l=2
	if(pos[0] != -1){
		f_ref.fl2.resize(pos.size());
		for(int n=0; n<pos.size();n++){
			f_ref.fl2[n]=ref_star.mode_params(pos[n],1);
		}
	}
	pos=where_dbl(ref_star.mode_params.col(0), 3, 1e-3); // pick l=3
	if(pos[0] != -1){
		f_ref.fl3.resize(pos.size());
		for(int n=0; n<pos.size();n++){
			f_ref.fl3[n]=ref_star.mode_params(pos[n],1);
		}
	}

	f_rescaled=rescale_freqs(Dnu, epsilon, f_ref, d0l);
	if (f_rescaled.error_status == true){
		std::cerr << "Error while rescaling: There is likely an issue in frequency tagging." << std::endl;
		std::cerr << "                       Debug in generate_cfg_from_synthese_file_Wscaled_aj required " << std::endl;
		exit(EXIT_FAILURE);
	}
	f_rescaled_lin.resize(f_ref.fl0.size()+ f_ref.fl1.size()+f_ref.fl2.size()+ f_ref.fl3.size());
	xmin=-nu_spread;
	xmax=nu_spread;
	for(int n=0; n<f_ref.fl0.size(); n++){
		mode_params(n,0)=0;
		f_rescaled_lin[n]=f_rescaled.fl0[n];
		// Adding to f_rescaled_lin[n] a uniform random quantity that is bounded by xmi and xmax
		if (nu_spread > 0)
		{
				f_rescaled_lin[n]=f_rescaled_lin[n] + (xmax-xmin)*distrib(gen);
		}
	}
	for(int n=0; n<f_ref.fl1.size(); n++){
		mode_params(n+f_ref.fl0.size(),0)=1;
		f_rescaled_lin[n+f_ref.fl0.size()]=f_rescaled.fl1[n];
		if (nu_spread > 0)
		{
				f_rescaled_lin[n+f_ref.fl0.size()]=f_rescaled_lin[n+f_ref.fl0.size()] + (xmax-xmin)*distrib(gen);
		}
	}
	for(int n=0; n<f_ref.fl2.size(); n++){
		mode_params(n+f_ref.fl0.size()+f_ref.fl1.size(),0)=2;
		f_rescaled_lin[n+f_ref.fl0.size()+f_ref.fl1.size()]=f_rescaled.fl2[n];
		if (nu_spread > 0)
		{
				f_rescaled_lin[n+f_ref.fl0.size()+f_ref.fl1.size()]=f_rescaled_lin[n+f_ref.fl0.size()+f_ref.fl1.size()] + (xmax-xmin)*distrib(gen);
		}
	}
	for(int n=0; n<f_ref.fl3.size(); n++){
		mode_params(n+f_ref.fl0.size()+f_ref.fl1.size()+f_ref.fl2.size(),0)=3;
		f_rescaled_lin[n+f_ref.fl0.size()+f_ref.fl1.size()+f_ref.fl2.size()]=f_rescaled.fl3[n];
		if (nu_spread > 0)
		{
				f_rescaled_lin[n+f_ref.fl0.size()+f_ref.fl1.size()+f_ref.fl2.size()]=f_rescaled_lin[n+f_ref.fl0.size()+f_ref.fl1.size()+f_ref.fl2.size()] + (xmax-xmin)*distrib(gen);
		}
	}

	// Compute the local noise for the new star
	local_noise_new.resize(f_rescaled_lin.size());
	local_noise_new.setZero();
	local_noise_new=harvey_like(noise_params_harvey, f_rescaled_lin, local_noise_new, 3); // Iterate on Noise_l0 to update it by putting the noise profile with one harvey profile
	/* --- This is to test with the equivalent function in Python that use directly a1,a2, etc... (and not the translation to a Harvey-like)
	std::cout << "noise_params_harvey = " << noise_params_harvey.transpose() << std::endl;
	for (int i=0;i<local_noise.size();i++){
		std::cout << local_noise[i] << std::setw(20) << local_noise_new[i] << std::endl;
	}
	VectorXd y(x.size());
	y.setZero();
	y=harvey_like(noise_params_harvey, x, y, 3); // Iterate on Noise_l0 to update it by putting the noise profile with one harvey profile
	std::ofstream debugFile("debug.txt");
	// Check if the debug file was successfully opened
	if (debugFile.is_open()) {
		for (int i = 0; i < x.size(); i++) {
			debugFile << x(i) << "\t" << y(i) << "\n";
		}
		debugFile.close();
	} else {
		std::cout << "Error: Unable to open the debug file." << std::endl;
	}
	std::cout << "Check that the local_noise_new makes sense when star_ref == star_sim" << std::endl;
	std::cout << "noise_params_harvey    = " << noise_params_harvey.transpose() << std::endl;
	std::cout << "noise_params_kallinger = " << noise_params_Kallinger.transpose() << std::endl;
	std::cout << "   - numax = " << numax_star << std::endl;
	exit(EXIT_SUCCESS);
	*/
	// ---- Rescaling Height and Width profiles ----
	HNRref=ref_star.mode_params.col(2);
	HNRref=HNRref.cwiseProduct(local_noise.cwiseInverse());
	pos=where_dbl(ref_star.mode_params.col(0), 0, 1e-3);

	HNRmaxref=0;
	for (int n=0; n<pos.size();n++){
		if (HNRmaxref < HNRref[pos[n]]){
			HNRmaxref=HNRref[pos[n]];
		}
	}

	Height_factor=HNR/HNRmaxref;  // compute the coeficient required to get the proper max(HNR)
	pos=where_dbl(HNRref, HNRmaxref, 0.001);
	if (pos[0] >= 0){
		Gamma_coef=Gamma_at_numax/ref_star.mode_params(pos[0], 3); // Correction coeficient to apply on Gamma(nu) in order to ensure that we have Gamma(nu=numax) = Gamma_at_numax
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
	if (local_noise_new.sum()!= 0){
		//N0=1; // Imposing the the Noise background is 1
		//std::cout <<  "Using N0=" << N0 << " (white noise)" << std::endl;
		std::cout << "Using the Noise profile of the New star, using provided noise parameters and Karoff et al. 2010 prescription" << std::endl;
		std::cout << "HNR of all modes:" << std::endl;
		std::cout << "     "  << std::setw(15) << "degree" << std::setw(15) << "freq_template" << std::setw(15) <<  "freq_rescaled" << std::setw(15) << "HNR" << std::setw(15) <<  "Height" << std::setw(15) << "local_noise" << std::endl;
		for(i =0; i<ref_star.mode_params.rows(); i++){
			h_star[i]=Height_factor * HNRref[i] * local_noise_new[i];
			// Adding a unifomly random value, with maximum range 2*H_spread
			if (H_spread > 0)
			{
					xmin=h_star[i]*(1. - H_spread/100.);
					xmax=h_star[i]*(1. + H_spread/100.);
					h_star[i]=xmin + (xmax-xmin)*distrib(gen);
			}
			std::cout << "     " << std::setw(15) << ref_star.mode_params(i,0) << std::setw(15) << ref_star.mode_params(i,1)  << std::setw(15)  << f_rescaled_lin[i] << std::setw(15) << Height_factor * HNRref[i]  << std::setw(15) << h_star[i] << std::setw(15)  << local_noise_new[i] << std::endl;
		}
	} else{
		std::cerr << "Error: local_noise_new has no valid value. Debug required." << std::endl;
		std::cerr << "         The program will stop now" << std::endl;
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
	
	// Filling the output Matrix. 
	//   Note on assumptions here: We assume that the template has monotonically increasing frequency, for each l. 
	//                             It is also assumed that block of l are in the strict order l=0,1,2,3

	//mode_params.col(0)=ref_star.mode_params.col(0); // List of els
	//mode_params.col(1)=ref_star.mode_params.col(1); // List of frequencies
	mode_params.col(1)=f_rescaled_lin; // Frequencies in a flat array
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

	noise_params(0,0)=noise_params_harvey[0];
	noise_params(0,1)=noise_params_harvey[1];
	noise_params(0,2)=noise_params_harvey[2];
	noise_params(1,0)=noise_params_harvey[3];
	noise_params(1,1)=noise_params_harvey[4];
	noise_params(1,2)=noise_params_harvey[5]; 
	noise_params(2, 0)=noise_params_harvey[6];
	noise_params(2, 1)=noise_params_harvey[7];
	noise_params(2, 2)=noise_params_harvey[8];
	noise_params(3, 0)=noise_params_harvey[9]; // White noise
	noise_params(3, 1)=-2;
	noise_params(3, 2)=-2;
	// A FUNCTION THAT WRITES THE Noise
	write_star_noise_params(noise_params, file_out_noise);
	}



// ------ Common ------

// Compute the effect of centrifugal distorsion on mode frequencies
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
    eta0=3.*M_PI/(rho*G); 
    return eta0;
}
