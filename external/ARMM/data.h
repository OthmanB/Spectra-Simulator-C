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

//using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXi;
using Eigen::VectorXd;

// The general structure that contain the input data (spectrum, lightcurve, ...) that has to be analysed
struct Data{
		VectorXd x; // In the case of a 1D fit (e.g. fit of Mass, [Fe/H], ...), this variable will be ignored by the model function.
		VectorXd xrange; // Contains the min and max of x... in practice, it is used to limit the data range if requested by the cfg file (e.g. freq_range option in the .MCMC)
		VectorXd y;
		VectorXd sigma_y;
		long Nx; // Ny is not checked but should be as long as x
		std::string xlabel; // label for x-axis. In ASCII file, the axis labels are identified by the symbol ! at the begining of the line
		std::string ylabel;
		std::string xunit;
		std::string yunit;
		std::vector<std::string> header; // Any kind of information that is useful. e.g. if it is a test, model info + model inputs. Any header (In ASCII, marked by #), will be put there.
};

struct Input_Data{
	std::string model_fullname; // The fullname of the model that is going to be processed
	std::vector<std::string> inputs_names;
	std::vector<std::string> priors_names;
	VectorXi priors_names_switch; // Used instead of priors_names in loops. This allows the use of the switch(x) - case statements instead of if/else.
	VectorXd inputs;
	VectorXi relax;
	MatrixXd priors;
	VectorXi plength;
	VectorXd extra_priors; // Contains extra parameters that could be used for priors
};

// A Generic structure that helps to encapsulate a Matrix of information along with some metadata
struct Data_Nd{
	MatrixXd data; // 2D array in which each columns are expected to contain the values of a given parameter
	std::vector<std::string> header;
	std::vector<std::string> labels;
	std::vector<std::string> units;
};

// Header of the parameters
// Useful to read the stand alone ASCII header written when dealing with BINARY OUTPUTS.
struct Params_hdr{
	std::vector<std::string> header; // Any comment
	int Nsamples; // Total number of samples (no necessarily the actual number of written samples)
	int Nchains;
 	int Nvars;
	int Ncons;
	//std::vector<bool> relax;
	VectorXi relax;
	VectorXi plength;
	std::vector<std::string> constant_names; // The name of the constants
	VectorXd constant_values;
	std::vector<std::string> variable_names; // The name of the variables
};

struct MCMC_files{
	std::string ID;
	double Dnu;
	double numax;
	double C_l;
	VectorXi els;
	VectorXd freq_range;
	std::vector<std::string> param_type;
	//std::vector<double> freqs_ref;conservativeResize
	VectorXd freqs_ref;
	std::vector<bool> relax_freq, relax_gamma, relax_H;

	//std::vector<double> hyper_priors;
	VectorXd hyper_priors;
	MatrixXd eigen_params;
	VectorXd noise_params;
	MatrixXd noise_s2;

	std::vector<std::string> common_names;
	std::vector<std::string> common_names_priors;
	MatrixXd modes_common;
};

// Structure that keep information of the derivatives
struct Deriv_out{
	VectorXd xderiv;
	VectorXd deriv;
	VectorXd error;
};

struct gnuplt_Data {
/*
 * This is an encapsulator for data when ploting with gnuplot-iostream.h
*/
    double x;  // x axis value
    double y1;             // y axis series 1
    double y2;             // y axis series 2
    double y3;             // y axis series 3
};

typedef std::vector<gnuplt_Data> gnuplt_Dataset;

namespace gnuplotio {
    template<>
    struct TextSender<gnuplt_Data> {
        static void send(std::ostream &stream, const gnuplt_Data &v) {
           // TextSender<std::string>::send(stream, v.x);
            stream << " ";
            TextSender<double>::send(stream, v.x);
            stream << " ";
            TextSender<double>::send(stream, v.y1);
            stream << " ";
            TextSender<double>::send(stream, v.y2);
            stream << " ";
            TextSender<float>::send(stream, v.y3);
        }
    };
}

struct Data_Basic{
	std::vector<std::string> strarr; // Any comment
	VectorXi vecXi; // Case number
};


// ----------------------------------------
// ----- For mixed modes calculation ------
// ----------------------------------------
/*struct Data_2vectXd{
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
	std::string output_file_rot;
	VectorXd Vl;
	long double H0_spread;
	std::string filetemplate;
	long double sigma_p;
	long double sigma_m;
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
};
*/