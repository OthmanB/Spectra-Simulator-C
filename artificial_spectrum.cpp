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
#include <Eigen/Dense>
//# include <random>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include "artificial_spectrum.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;


void artificial_spectrum_act_asym(double Tobs, double Cadence, double Nspectra, std::string dir_core, std::string identifier, bool doplots, bool write_inmodel){

	// General variables
	Data_Nd data_modes, data_noise;
	std::string file_in_modes=dir_core + "Configurations/tmp/modes_tmp.cfg";
	std::string file_in_noise=dir_core + "Configurations/tmp/noise_tmp.cfg";
	std::string fileout_spectrum=dir_core + "Data/Spectra_ascii/" + strtrim(identifier) + ".ascii";
	std::string fileout_params=dir_core + "Data/Spectra_info/" + strtrim(identifier) + ".in";
	std::string fileout_plot=dir_core + "Data/Spectra_plot/" + strtrim(identifier) + ".eps" ;

	std::string delimiter=" ";
	bool verbose_data=0; // SET TO 1 IF YOU WANT MORE INFO ABOUT WHAT IS READ
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
	//std::cout << "data_modes.data.cols()=" << data_modes.data.cols() << std::endl;
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
		//std::cout << "Here" << std::endl;
	}
	//std::cout << "end loop" << std::endl;
	scoef1=scoef1/data_modes.data.cols(); // scoef1 is based on the average mode width
	scoef2=scoef1; // scoef2 is based on the average mode width

	scoef1=1*scoef1;
	scoef2=10*scoef2;

	// Build the noise background using the Harvey like profile...
	std::cout << "    - Generating the model of noise..." << std::endl;
	ones.resize(Ndata);
	ones.setConstant(1);
	spec_noise.setZero(Ndata);
	s_noise.resize(Ndata);
	//std::cout << "data_noise.data.rows()=" << data_noise.data.rows() << std::endl;
	for(int i=0; i<data_noise.data.rows(); i++){
		H=data_noise.data(i,0);
		tau=data_noise.data(i,1);
		p=data_noise.data(i,2);
		
		/*
		std::cout << "--------" << std::endl;
		std::cout << "H=" << H << std::endl;
		std::cout << "tau=" << tau << std::endl;
		std::cout << "p=" << p << std::endl;
		std::cout << "--------" << std::endl;
		*/
		s_noise.setZero();	
		if(H>0 && (tau<0 || p<0) ){
			if(tau==-2 && p==-2){
				s_noise.setConstant(H);
				//std::cout << "Condition H ok but others at -2" << std::endl;
				spec_noise=spec_noise + s_noise;
			} else {
				std::cout << "Problem with the definition of the White noise" << std::endl;
				std::cout << "   H is defined positive but tau and/or p are not set to -2" << std::endl;
				std::cout << "   For a white noise, please set tau and p together to -2" << std::endl;
				std::cout << "The program will stop now" << std::endl;
				exit(EXIT_FAILURE);
			}
		} else {
			if(H<=0 && (tau>0 || p>0)){
				std::cout << "Problem with the definition of the Height parameter for the noise" << std::endl;
				std::cout << " if tau or p positive, H cannot be negative" << std::endl;
				std::cout << "The program will stop now" << std::endl;
				exit(EXIT_FAILURE);
			}
			if( ((tau < 0) && (p > 0)) || ((tau > 0) && (p < 0)) ){
				std::cout << "Problem with the definition of the Harvey-like profile" << std::endl;
				std::cout << "   H is defined but tau or p are set to -1" << std::endl;
				std::cout << "   For a Harvey-like profile, please set tau and p together" << std::endl;
				std::cout << "The program will stop now" << std::endl;
				exit(EXIT_FAILURE);
			}
			if( (tau >0 && p >0) ){ // If all the parameters of the Harvey noise are defined, then...
				s_noise=((1e-3)*tau*freq).array().pow(p); 
				s_noise=H*(s_noise + ones).cwiseInverse();
				//std::cout << "Condition x y z defined" << std::endl;
				spec_noise=spec_noise + s_noise;
			}			
		}
		//std::cout << "s_noise=" << s_noise << std::endl;
		//std::cout << "s_noise.size()=" << s_noise.size() << std::endl;		
	}
	
	// Final spectrum and saving functions
	input_spec_model=spec_noise + spec_modes;

    std::cout << "    - Simulating a stochastic noise with chi(2,2p) statistics (p=" << std::floor(Nspectra) << ")..." << std::endl;
    spec_reg.resize(Ndata);
    spec_reg.setZero();
    stmp.resize(Ndata);
    for(int ns=0; ns<Nspectra; ns++){
        std::cout << "[" << ns << "]...";
        for(int i=0; i<Ndata; i++){
            F1=dist();
            F2=dist();
            stmp[i]=std::abs(F1*F1*input_spec_model[i]/2. + F2*F2*input_spec_model[i]/2.);
        }
        spec_reg=spec_reg + stmp/Nspectra;
    }
    std::cout << "done" << std::endl;
/*
	for(int ii=0; ii<freq.size(); ii++){
			std::cout << freq[ii] << "   " << input_spec_model[ii]  << "  " << spec_reg[ii] << std::endl;
	}
	exit(EXIT_SUCCESS);
*/
	std::cout << "    - Saving the spectrum..." << std::endl;
	write_spectrum(freq, spec_reg, input_spec_model, fileout_spectrum, write_inmodel);

	std::cout << "    - Saving the input parameters..." << std::endl;
	VectorXd spec_params(2);
	spec_params << Tobs, Cadence;
	write_star_params_act_asym(spec_params, data_modes.data, data_noise.data, fileout_params, identifier);

    if(doplots == 1 ){
        std::cout << "    - Generating a plot of the spectrum..." << std::endl;
        gnuplt_model(freq, spec_reg, input_spec_model, scoef1, scoef2, fileout_plot);
    }
	//std::cout << "All done for the simulated star " << identifier << std::endl;

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
