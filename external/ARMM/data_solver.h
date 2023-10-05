/*
 * data.h
 *
 * Header file that contains all kind of class/structures
 * used to process and/or encapsulate data
 * 
 *  Created on: 22 Feb 2016
 *      Author: obenomar
 */

#pragma once
#include <Eigen/Dense>
#include <string>
#include <vector>
#include "gnuplot-iostream.h"

using Eigen::MatrixXd;
using Eigen::VectorXi;
using Eigen::VectorXd;

// ----------------------------------------
// ----- For mixed modes calculation ------
// ----------------------------------------
struct Data_2vectXd{
	VectorXd vecXd1;
	VectorXd vecXd2;		
};

// For simulations only, we use a template to derive height and widths of modes
struct template_file{
	std::string ID_ref;
	double numax_ref;
	double Dnu_ref; 
	double epsilon_ref;
	MatrixXd data_ref;
};

struct Data_coresolver{
	VectorXd nu_m, ysol, nu,pnu, gnu; //
};

struct Data_eigensols{
	VectorXd nu_p, nu_g, nu_m, dnup, dPg; // The two last ones are derivatives of nu_p ~ Dnu and 1/nu_g ~ DPl
};

struct Data_eigensols_all{
	VectorXd nu_l0;
};

struct Data_rot2zone{
	long double rot_core, rot_env;
};

struct Envelope_lat_dif_rot{ // Added on 13 Sep 2023
	long double a3_l2=0; 
	long double a3_l3=0;
	long double a5_l3=0;
};

struct Envelope_asphericity{
	long double a2_l1; // Added on 13 Sep 2023
	long double a2_l2; // Added on 13 Sep 2023
	long double a2_l3; // Added on 13 Sep 2023
	long double a4_l2; // Added on 13 Sep 2023
	long double a4_l3; // Added on 13 Sep 2023
	long double a6_l3; // Added on 13 Sep 2023
};
struct Cfg_synthetic_star{
	long double Teff_star; 
	long double numax_star;
	long double Dnu_star;
	long double epsilon_star;
	long double delta0l_percent_star;
	long double beta_p_star;
	long double alpha_p_star;
	long double nmax_star;
	long double DPl_star;
	long double alpha_g_star;
	long double q_star;
	long double fmin; 
	long double fmax;
	long double maxHNR_l0;
	VectorXd noise_params_harvey_like;
	long double Gamma_max_l0;
	long double rot_env_input;
	long double rot_ratio_input; 
	long double rot_core_input;
	Envelope_lat_dif_rot env_lat_dif_rot; // Added on 13 Sep 2023
	Envelope_asphericity env_aspher;
	std::string output_file_rot;
	VectorXd Vl;
	long double H0_spread;
	std::string filetemplate;
	long double sigma_p;
	long double sigma_m;
	long double Hfactor;
	long double Wfactor;
	long double inclination;
	MatrixXd nu_nl; // Frequencies of the modes, are here if provided by a template (e.g a theoretical model) and handled by the MCMC model
	VectorXi Nf_el; // Gives the number of modes 
	bool use_nu_nl=false; // If set to true, use the nu_nl frequencies instead of computing them from the asymptotic. These must be set
};



struct Params_synthetic_star{
	VectorXd nu_l0;
	VectorXd nu_p_l1;
	VectorXd nu_g_l1; 
	VectorXd nu_m_l1;
	VectorXd nu_l2;
	VectorXd nu_l3;
	VectorXd width_l0;
	VectorXd width_l1;
	VectorXd width_l2; 
	VectorXd width_l3;
	VectorXd height_l0;
	VectorXd height_l1;
	VectorXd height_l2;
	VectorXd height_l3;
	VectorXd a1_l1;
	VectorXd a1_l2; 
	VectorXd a1_l3;
	VectorXd a2_l1; // Added on 13 Sep 2023
	VectorXd a2_l2; // Added on 13 Sep 2023
	VectorXd a2_l3; // Added on 13 Sep 2023
	VectorXd a3_l2; // Added on 13 Sep 2023
	VectorXd a3_l3; // Added on 13 Sep 2023
	VectorXd a4_l2; // Added on 13 Sep 2023
	VectorXd a4_l3; // Added on 13 Sep 2023
	VectorXd a5_l3; // Added on 13 Sep 2023
	VectorXd a6_l3; // Added on 13 Sep 2023
};