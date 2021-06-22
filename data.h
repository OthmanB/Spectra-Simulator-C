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
using Eigen::VectorXd;


struct Data_Nd{
	MatrixXd data; // 2D array in which each columns are expected to contain the values of a given parameter
	std::vector<std::string> header;
	std::vector<std::string> labels;
	std::vector<std::string> units;
};

struct Model_data{
    VectorXd params; // 1D array with the data for one single model from the read text file.
    std::vector<std::string> labels_params; // Name of the parameters within params
    MatrixXd freqs; // Table of frequencies with (l,n)
    std::vector<std::string> labels_freqs;
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

struct Star_params{
/* 
    A data structure used by the models that need to read the parameters of a reference 
    star from a .in file to operate
*/
    VectorXd spec_params; 
    MatrixXd mode_params, noise_params;
    std::string identifier;
    double cadence;
    double Tobs;
};

struct Config_Data{
/* 
 * Encapsulate the parameters used to generate the spectra
*/

	bool erase_old_files; // if set to 1 any older of same name will be erased
    bool doplots; // if set to 0, no plots
    bool write_inmodel; // if set to 0, output ascii file does not contain a column with the input model
    bool limit_data_range;
    bool do_modelfiles;
    std::string modefile_modelname; // model name writen within the model file
	std::string model_name; // model Name within this program
	std::string forest_type;
    std::string extra_params; // Any kind of extra parameters that would be encoded in a string (e.g. filenames separated by " ")
    std::vector<std::string> template_files; // List of template files used to determine Height and widths. If none, set to None
	double Tobs, Cadence, Nspectra, Nrealisation;
	std::vector<std::string> labels; 
	std::vector<double> val_min, val_max, step, forest_params;

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
};


// Structure that keep information of the derivatives
struct Deriv_out{
    VectorXd xderiv;
    VectorXd deriv;
    VectorXd error;
};


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
    double fmin;
    double fmax;
};