/*
 * artificial_spectrum.cpp
 *
 * Contains all kind of methods that is producing
 * a single artificial spectrum
 * 
 *  Created on: 05 May 2016
 *      Author: obenomar
 */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <Eigen/Dense>
//# include <random>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include "artificial_spectrum.h"
#include "noise_models.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;


void artificial_spectrum_act_asym(const double Tobs, const double Cadence, const double Nspectra, const long Nrealisation, const std::string dir_core, const std::string identifier, const bool doplots, const bool write_inmodel){

	// General variables
	Data_Nd data_modes, data_noise;
	std::string file_in_modes=dir_core + "Configurations/tmp/modes_tmp.cfg";
	std::string file_in_noise=dir_core + "Configurations/tmp/noise_tmp.cfg";
	std::string dir_common_template = dir_core + "Configurations/common_modelfile/";
	std::string fileout_spectrum;
	std::string fileout_params;
	std::string fileout_plot;

	const std::string delimiter=" ";
	const bool verbose_data=0; // SET TO 1 IF YOU WANT MORE INFO ABOUT WHAT IS READ
	int Ndata;
	double df, Delta, scoef1, scoef2;
	VectorXd freq;

	// Initialize the random generators 
 	boost::mt19937 *rng = new boost::mt19937();
  	rng->seed(time(NULL));

 	boost::normal_distribution<> distribution(0, 1);
  	boost::variate_generator< boost::mt19937, boost::normal_distribution<> > dist(*rng, distribution);
	

	// Variable of the models
	int l;
	double fc_l, H_l, gamma_l, f_s, eta, a3, b, alfa, beta_asym, angle_l, H, tau, p;
	double F1, F2;
	VectorXd ratios, s_mode, spec_modes, spec_noise, s_noise, ones, input_spec_model, spec_reg, stmp;

	data_modes=read_data_ascii_Ncols(file_in_modes, delimiter, verbose_data);
	data_noise=read_data_ascii_Ncols(file_in_noise, delimiter, verbose_data);

	df=1e6/(Tobs * 86400.);
	Delta=1e6/Cadence/2;
	Ndata=Delta/df;
	freq.setLinSpaced(Ndata, 0, Delta);

	// Build the sum of modes...
	std::cout << "    - Generating model of p modes..." << std::endl;
	spec_modes.setZero(Ndata);
	scoef1=0;
	for(int i=0; i<data_modes.data.rows(); i++){
		l=data_modes.data(i,0);
		fc_l=data_modes.data(i,1);
		H_l=data_modes.data(i,2);
		gamma_l=data_modes.data(i,3);
		f_s=data_modes.data(i,4);
		eta=data_modes.data(i,5);
		a3=data_modes.data(i,6);
		b=data_modes.data(i,7);
		alfa=data_modes.data(i,8);
		beta_asym=data_modes.data(i,9);
		angle_l=data_modes.data(i,10);
		ratios=amplitude_ratio(l, angle_l);
		s_mode=build_l_mode_act_simu(freq, H_l, fc_l, f_s, eta, a3, b, alfa, beta_asym, gamma_l, l, ratios);
		
		spec_modes=spec_modes + s_mode;

		scoef1=scoef1 + gamma_l;
	}
	scoef1=scoef1/data_modes.data.cols(); // scoef1 is based on the average mode width
	scoef2=scoef1; // scoef2 is based on the average mode width

	scoef1=1*scoef1;
	scoef2=10*scoef2;

	// Build the noise background using the Harvey like profile...
	std::cout << "    - Generating the model of noise..." << std::endl;
	spec_noise=harvey_1985(data_noise.data, freq);

	// Final spectrum and saving functions
	input_spec_model=spec_noise + spec_modes;

	// Making Nrealisation of noise
	fileout_params=dir_core + "Data/Spectra_info/" + strtrim(identifier) + ".in";
    spec_reg.resize(Ndata);
    stmp.resize(Ndata);
	//const int Nrealisation=2;
	for (int nr=0; nr<Nrealisation; nr++){
		fileout_spectrum=dir_core + "Data/Spectra_ascii/" + strtrim(identifier)  + "." + int_to_str(nr) + ".ascii";
		fileout_plot=dir_core + "Data/Spectra_plot/" + strtrim(identifier) + "." + int_to_str(nr) + ".eps" ;

    	std::cout << "    ["<< nr+1 << "/" << Nrealisation << "]" << std::endl;
    	std::cout << "        - Simulating a stochastic noise with chi(2,2p) statistics (p=" << std::floor(Nspectra) << ")..." << std::endl;
    	spec_reg.setZero();
    	for(int ns=0; ns<Nspectra; ns++){
    	    std::cout << "[" << ns << "]...";
    	    for(int i=0; i<Ndata; i++){
    	        F1=dist();
    	        F2=dist();
    	        stmp[i]=std::abs(F1*F1*input_spec_model[i]/2. + F2*F2*input_spec_model[i]/2.);
    	    }
    	    spec_reg=spec_reg + stmp/Nspectra;
    	}
    	std::cout << "Done" << std::endl;

		std::cout << "        - Saving the spectrum..." << std::endl;
		write_spectrum(freq, spec_reg, input_spec_model, fileout_spectrum, write_inmodel);

    	if(doplots == 1 ){
    	    std::cout << "        - Generating a plot of the spectrum..." << std::endl;
    	    gnuplt_model(freq, spec_reg, input_spec_model, scoef1, scoef2, fileout_plot);
    	}
    }
	std::cout << "        - Saving the input parameters..." << std::endl;
	VectorXd spec_params(2);
	spec_params << Tobs, Cadence;
	write_star_params_act_asym(spec_params, data_modes.data, data_noise.data, fileout_params, identifier);

}


void artificial_spectrum_a1a2a3asym(const double Tobs, const double Cadence, const double Nspectra, const long Nrealisation, const std::string dir_core, const std::string identifier, const bool doplots, const bool write_inmodel){

	// General variables
	Data_Nd data_modes, data_noise;
	std::string file_in_modes=dir_core + "Configurations/tmp/modes_tmp.cfg";
	std::string file_in_noise=dir_core + "Configurations/tmp/noise_tmp.cfg";
	std::string fileout_spectrum;
	std::string fileout_params;
	std::string fileout_plot;

	const std::string delimiter=" ";
	const bool verbose_data=0; // SET TO 1 IF YOU WANT MORE INFO ABOUT WHAT IS READ
	int Ndata;
	double df, Delta, scoef1, scoef2;
	VectorXd freq;

	// Initialize the random generators 
 	boost::mt19937 *rng = new boost::mt19937();
  	rng->seed(time(NULL));

 	boost::normal_distribution<> distribution(0, 1);
  	boost::variate_generator< boost::mt19937, boost::normal_distribution<> > dist(*rng, distribution);
	

	// Variable of the models
	int l;
	double fc_l, H_l, gamma_l, f_s, a2, a3, b, alfa, beta_asym, angle_l, H, tau, p;
	double F1, F2;
	VectorXd ratios, s_mode, spec_modes, spec_noise, s_noise, ones, input_spec_model, spec_reg, stmp;

	data_modes=read_data_ascii_Ncols(file_in_modes, delimiter, verbose_data);
	data_noise=read_data_ascii_Ncols(file_in_noise, delimiter, verbose_data);

	df=1e6/(Tobs * 86400.);
	Delta=1e6/Cadence/2;
	Ndata=Delta/df;
	freq.setLinSpaced(Ndata, 0, Delta);

	// Build the sum of modes...
	std::cout << "    - Generating model of p modes..." << std::endl;
	spec_modes.setZero(Ndata);
	scoef1=0;
	for(int i=0; i<data_modes.data.rows(); i++){
		l=data_modes.data(i,0);
		fc_l=data_modes.data(i,1);
		H_l=data_modes.data(i,2);
		gamma_l=data_modes.data(i,3);
		f_s=data_modes.data(i,4);
		a2=data_modes.data(i,5);
		a3=data_modes.data(i,6);
//		b=data_modes.data(i,7);
//		alfa=data_modes.data(i,8);
		beta_asym=data_modes.data(i,7);
		angle_l=data_modes.data(i,8);
		ratios=amplitude_ratio(l, angle_l);
		s_mode=build_l_mode_a1a2a3(freq, H_l, fc_l, f_s, a2, a3, beta_asym, gamma_l, l, ratios);
		
		spec_modes=spec_modes + s_mode;

		scoef1=scoef1 + gamma_l;
	}
	scoef1=scoef1/data_modes.data.cols(); // scoef1 is based on the average mode width
	scoef2=scoef1; // scoef2 is based on the average mode width

	scoef1=1*scoef1;
	scoef2=10*scoef2;

	// Build the noise background using the Harvey like profile...
	std::cout << "    - Generating the model of noise..." << std::endl;
	spec_noise=harvey_1985(data_noise.data, freq);

	// Final spectrum and saving functions
	input_spec_model=spec_noise + spec_modes;

	// Making Nrealisation of noise
	fileout_params=dir_core + "Data/Spectra_info/" + strtrim(identifier) + ".in";
    spec_reg.resize(Ndata);
    stmp.resize(Ndata);
	//const int Nrealisation=2;
	for (int nr=0; nr<Nrealisation; nr++){
		fileout_spectrum=dir_core + "Data/Spectra_ascii/" + strtrim(identifier)  + "." + int_to_str(nr) + ".ascii";
		fileout_plot=dir_core + "Data/Spectra_plot/" + strtrim(identifier) + "." + int_to_str(nr) + ".eps" ;

    	std::cout << "    ["<< nr+1 << "/" << Nrealisation << "]" << std::endl;
    	std::cout << "        - Simulating a stochastic noise with chi(2,2p) statistics (p=" << std::floor(Nspectra) << ")..." << std::endl;
    	spec_reg.setZero();
    	for(int ns=0; ns<Nspectra; ns++){
    	    std::cout << "[" << ns << "]...";
    	    for(int i=0; i<Ndata; i++){
    	        F1=dist();
    	        F2=dist();
    	        stmp[i]=std::abs(F1*F1*input_spec_model[i]/2. + F2*F2*input_spec_model[i]/2.);
    	    }
    	    spec_reg=spec_reg + stmp/Nspectra;
    	}
    	std::cout << "Done" << std::endl;

		std::cout << "        - Saving the spectrum..." << std::endl;
		write_spectrum(freq, spec_reg, input_spec_model, fileout_spectrum, write_inmodel);

    	if(doplots == 1 ){
    	    std::cout << "        - Generating a plot of the spectrum..." << std::endl;
    	    gnuplt_model(freq, spec_reg, input_spec_model, scoef1, scoef2, fileout_plot);
    	}
    }
	std::cout << "        - Saving the input parameters..." << std::endl;
	VectorXd spec_params(2);
	spec_params << Tobs, Cadence;
	write_star_params_act_asym(spec_params, data_modes.data, data_noise.data, fileout_params, identifier);

}


void artificial_spectrum_a1Alma3(const double Tobs, const double Cadence, const double Nspectra, const long Nrealisation, const std::string dir_core,
							    const std::string identifier, const bool doplots, const bool write_inmodel, const bool domodelfiles, 
							    const bool limit_data_range, const std::string modelname){

	const bool verbose_data=0; // SET TO 1 IF YOU WANT MORE INFO ABOUT WHAT IS READ
	const int common_case=0; // Case identifying this scenario
	const std::string delimiter=" ";
	const std::string file_in_modes=dir_core + "Configurations/tmp/modes_tmp.cfg";	
	const std::string file_in_noise=dir_core + "Configurations/tmp/noise_tmp.cfg";
	const std::string dir_common_template = dir_core + "Configurations/common_modelfile/";

	// General variables
	std::string fileout_spectrum;
	std::string fileout_params;
	std::string fileout_plot;	
	std::string file_out_modelfile;	
	
	int Ndata;
	double df, Delta, scoef1, scoef2;
	Data_Nd data_modes, data_noise;
	
	// Initialize the random generators 
 	boost::mt19937 *rng = new boost::mt19937();
  	rng->seed(time(NULL));

 	boost::normal_distribution<> distribution(0, 1);
  	boost::variate_generator< boost::mt19937, boost::normal_distribution<> > dist(*rng, distribution);
	
	// Variable of the models
	int l;
	double fc_l, H_l, gamma_l, a1, eta0, a3, beta_asym, angle_l, H, tau, p, epsilon_nl;
	double F1, F2;
	VectorXd ratios, thetas(2), s_mode, mode_range(2),spec_modes, spec_noise, s_noise, ones, input_spec_model, freq, spec_reg, stmp;

	mode_range << -1, -1; // Initialise the range to invalid values

	data_modes=read_data_ascii_Ncols(file_in_modes, delimiter, verbose_data);
	data_noise=read_data_ascii_Ncols(file_in_noise, delimiter, verbose_data);

	df=1e6/(Tobs * 86400.);
	Delta=1e6/Cadence/2;
	Ndata=Delta/df;
	freq.setLinSpaced(Ndata, 0, Delta);

	// Build the sum of modes...
	std::cout << "    - Generating model of p modes..." << std::endl;
	spec_modes.setZero(Ndata);
	scoef1=0;
	for(int i=0; i<data_modes.data.rows(); i++){
		l=data_modes.data(i,0);
		fc_l=data_modes.data(i,1);
		H_l=data_modes.data(i,2);
		gamma_l=data_modes.data(i,3);
		a1=data_modes.data(i,4);
		eta0=data_modes.data(i,5);
		epsilon_nl=data_modes.data(i,6);
		thetas[0]=data_modes.data(i,7);
		thetas[1]=data_modes.data(i,8);
		a3=data_modes.data(i,9);
		beta_asym=data_modes.data(i,10);
		angle_l=data_modes.data(i,11);
		ratios=amplitude_ratio(l, angle_l);
		s_mode=build_l_mode_a1etaAlma3(freq,  H_l, fc_l, a1, eta0, epsilon_nl, thetas, a3, beta_asym, gamma_l, l, ratios);
		spec_modes=spec_modes + s_mode;
		scoef1=scoef1 + gamma_l;
	}
	scoef1=scoef1/data_modes.data.cols(); // scoef1 is based on the average mode width
	scoef2=scoef1; // scoef2 is based on the average mode width

	scoef1=1*scoef1;
	scoef2=10*scoef2;

	// Build the noise background using the Harvey like profile...
	std::cout << "    - Generating the model of noise..." << std::endl;
	spec_noise=harvey_1985(data_noise.data, freq);

	// Final spectrum and saving functions
	input_spec_model=spec_noise + spec_modes;

	// Making Nrealisation of noise
	fileout_params=dir_core + "Data/Spectra_info/" + strtrim(identifier) + ".in";
    spec_reg.resize(Ndata);
    stmp.resize(Ndata);
	//const int Nrealisation=2;
	for (int nr=0; nr<Nrealisation; nr++){
		fileout_spectrum=dir_core + "Data/Spectra_ascii/" + strtrim(identifier)  + "." + int_to_str(nr) + ".data";
		fileout_plot=dir_core + "Data/Spectra_plot/" + strtrim(identifier) + "." + int_to_str(nr) + ".eps" ;
		file_out_modelfile=dir_core + "Data/Spectra_modelfile/" + strtrim(identifier) + "." + int_to_str(nr) + ".model" ;

    	std::cout << "    ["<< nr+1 << "/" << Nrealisation << "]" << std::endl;
    	std::cout << "        - Simulating a stochastic noise with chi(2,2p) statistics (p=" << std::floor(Nspectra) << ")..." << std::endl;
    	spec_reg.setZero();
    	for(int ns=0; ns<Nspectra; ns++){
    	    std::cout << "[" << ns << "]...";
    	    for(int i=0; i<Ndata; i++){
    	        F1=dist();
    	        F2=dist();
    	        stmp[i]=std::abs(F1*F1*input_spec_model[i]/2. + F2*F2*input_spec_model[i]/2.);
    	    }
    	    spec_reg=spec_reg + stmp/Nspectra;
    	}
    	std::cout << "Done" << std::endl;

    	if (domodelfiles==1){
    		std::cout << "        - Saving the model file" << std::endl;
    		mode_range=write_star_model(data_modes.data, data_noise.data, file_out_modelfile, identifier, modelname, dir_common_template);
  			if (limit_data_range == 0){ // Reset the mode_range if the full range was requested
  				mode_range[0] =-1;
  				mode_range[1] =-1;
			}
    	}
		std::cout << "        - Saving the spectrum..." << std::endl;
		write_spectrum(freq, spec_reg, input_spec_model, fileout_spectrum, write_inmodel, mode_range[0], mode_range[1]);

    	if(doplots == 1 ){
    	    std::cout << "        - Generating a plot of the spectrum..." << std::endl;
    	    gnuplt_model(freq, spec_reg, input_spec_model, scoef1, scoef2, fileout_plot);
    	}
    }
	std::cout << "        - Saving the input parameters..." << std::endl;
	VectorXd spec_params(2);
	spec_params << Tobs, Cadence;

	write_star_params_Alm(spec_params, data_modes.data, data_noise.data, fileout_params, identifier);

}

void artificial_spectrum_aj(const double Tobs, const double Cadence, const double Nspectra, const long Nrealisation, 
								 const std::string dir_core, const std::string identifier, const bool doplots, const bool write_inmodel,
								 const bool domodelfiles, const bool limit_data_range, const std::string modelname){

	const bool verbose_data=0; // SET TO 1 IF YOU WANT MORE INFO ABOUT WHAT IS READ
	const int common_case=0; // Case identifying this scenario
	const std::string delimiter=" ";
	const std::string file_in_modes=dir_core + "Configurations/tmp/modes_tmp.cfg";	
	const std::string file_in_noise=dir_core + "Configurations/tmp/noise_tmp.cfg";
	const std::string dir_common_template = dir_core + "Configurations/common_modelfile/";

	// General variables
	std::string fileout_spectrum;
	std::string fileout_params;
	std::string fileout_plot;	
	std::string file_out_modelfile;	
	
	int Ndata;
	double df, Delta, scoef1, scoef2;
	Data_Nd data_modes, data_noise;
	
	// Initialize the random generators 
 	boost::mt19937 *rng = new boost::mt19937();
  	rng->seed(time(NULL));

 	boost::normal_distribution<> distribution(0, 1);
  	boost::variate_generator< boost::mt19937, boost::normal_distribution<> > dist(*rng, distribution);
	
	// Variable of the models
	int l;
	double fc_l, H_l, gamma_l, a1, a2, a3, a4,a5,a6, beta_asym, angle_l, H, tau, p, eta0;
	double F1, F2;
	VectorXd ratios, s_mode, mode_range(2),spec_modes, spec_noise, s_noise, ones, input_spec_model, freq, spec_reg, stmp;

	mode_range << -1, -1; // Initialise the range to invalid values

	data_modes=read_data_ascii_Ncols(file_in_modes, delimiter, verbose_data);
	data_noise=read_data_ascii_Ncols(file_in_noise, delimiter, verbose_data);

	df=1e6/(Tobs * 86400.);
	Delta=1e6/Cadence/2;
	Ndata=Delta/df;
	freq.setLinSpaced(Ndata, 0, Delta);

	// Build the sum of modes...
	std::cout << "    - Generating model of p modes..." << std::endl;
	spec_modes.setZero(Ndata);
	scoef1=0;
	for(int i=0; i<data_modes.data.rows(); i++){
		l=data_modes.data(i,0);
		fc_l=data_modes.data(i,1);
		H_l=data_modes.data(i,2);
		gamma_l=data_modes.data(i,3);
		a1=data_modes.data(i,4);
		a2=data_modes.data(i,5);
		a3=data_modes.data(i,6);
		a4=data_modes.data(i,7);
		a5=data_modes.data(i,8);
		a6=data_modes.data(i,9);
		beta_asym=data_modes.data(i,10);
		angle_l=data_modes.data(i,11);
		ratios=amplitude_ratio(l, angle_l);
		eta0=0;
		s_mode=build_l_mode_aj(freq,  H_l, fc_l, a1, a2, a3, a4, a5, a6, eta0, beta_asym, gamma_l, l, ratios);
		spec_modes=spec_modes + s_mode;
		scoef1=scoef1 + gamma_l;
	}
	scoef1=scoef1/data_modes.data.cols(); // scoef1 is based on the average mode width
	scoef2=scoef1; // scoef2 is based on the average mode width

	scoef1=1*scoef1;
	scoef2=10*scoef2;

	// Build the noise background using the Harvey like profile...
	std::cout << "    - Generating the model of noise..." << std::endl;
	spec_noise=harvey_1985(data_noise.data, freq);

	// Final spectrum and saving functions
	input_spec_model=spec_noise + spec_modes;

	// Making Nrealisation of noise
	fileout_params=dir_core + "Data/Spectra_info/" + strtrim(identifier) + ".in";
    spec_reg.resize(Ndata);
    stmp.resize(Ndata);
	//const int Nrealisation=2;
	for (int nr=0; nr<Nrealisation; nr++){
		fileout_spectrum=dir_core + "Data/Spectra_ascii/" + strtrim(identifier)  + "." + int_to_str(nr) + ".data";
		fileout_plot=dir_core + "Data/Spectra_plot/" + strtrim(identifier) + "." + int_to_str(nr) + ".eps" ;
		file_out_modelfile=dir_core + "Data/Spectra_modelfile/" + strtrim(identifier) + "." + int_to_str(nr) + ".model" ;

    	std::cout << "    ["<< nr+1 << "/" << Nrealisation << "]" << std::endl;
    	std::cout << "        - Simulating a stochastic noise with chi(2,2p) statistics (p=" << std::floor(Nspectra) << ")..." << std::endl;
    	spec_reg.setZero();
    	for(int ns=0; ns<Nspectra; ns++){
    	    std::cout << "[" << ns << "]...";
    	    for(int i=0; i<Ndata; i++){
    	        F1=dist();
    	        F2=dist();
    	        stmp[i]=std::abs(F1*F1*input_spec_model[i]/2. + F2*F2*input_spec_model[i]/2.);
    	    }
    	    spec_reg=spec_reg + stmp/Nspectra;
    	}
    	std::cout << "Done" << std::endl;

    	if (domodelfiles==1){
    		std::cout << "        - Saving the model file" << std::endl;
    		mode_range=write_star_model(data_modes.data, data_noise.data, file_out_modelfile, identifier, modelname,dir_common_template);
  			if (limit_data_range == 0){ // Reset the mode_range if the full range was requested
  				mode_range[0] =-1;
  				mode_range[1] =-1;
			}
    	}
		std::cout << "        - Saving the spectrum..." << std::endl;
		write_spectrum(freq, spec_reg, input_spec_model, fileout_spectrum, write_inmodel, mode_range[0], mode_range[1]);

    	if(doplots == 1 ){
    	    std::cout << "        - Generating a plot of the spectrum..." << std::endl;
    	    gnuplt_model(freq, spec_reg, input_spec_model, scoef1, scoef2, fileout_plot);
    	}
    }
	std::cout << "        - Saving the input parameters..." << std::endl;
	VectorXd spec_params(2);
	spec_params << Tobs, Cadence;

	write_star_params_aj(spec_params, data_modes.data, data_noise.data, fileout_params, identifier);

}

/*
int main(){
	double Tobs, Cadence;
	std::string dir_core;
	std::string identifier;

	Tobs=100; // days
	Cadence=60; // seconds
	
	dir_core="/home/obenomar/Dropbox/Temporary/Cpp-Spec-Sim/";
	identifier="0000";
	artificial_spectrum_act_asym(Tobs, Cadence, dir_core, identifier);
	//std::cout << "Hello" << std::endl;
}
*/
