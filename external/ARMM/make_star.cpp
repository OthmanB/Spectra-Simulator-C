/**
 * @file make_star.cpp
 * @brief Program to generate a simulated star spectrum.
 * 
 * This program generates a simulated star spectrum based on the provided configuration file. It uses various modules and functions to perform the simulation.
 */
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <boost/program_options.hpp>
#include "version_solver.h"
#include "data_solver.h"
#include "string_handler.h"
#include "interpol.h"
#include "noise_models.h" // get the harvey_1985 function
#include "solver_mm.h"
#include "readparams_job.h"
#include "writeparams_job.h"
#include "bump_DP.h"
#include "configure_make_star.h"

namespace po = boost::program_options;
/**
 * @brief Main function to run the make_star program.
 * 
 * @param argc The number of command-line arguments.
 * @param argv An array of command-line arguments.
 * @return int The exit status of the program.
 */
int main(int argc, char* argv[]){
	MatrixXd mode_params, noise_params(3,3);

	Cfg_synthetic_star cfg_star;
	Params_synthetic_star params_out;
    std::string cfg_file;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,H", "produce help message")
        ("config-file", po::value<std::string>(&cfg_file)->default_value("../config/make_star.cfg"), "Configuration file for a given simulated star");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << "This program allows you to run a fully simulated star spectrum" << std::endl;
        std::cout << "By default it uses the Sun Height and Width profile. And mixed mode frequencies " << std::endl;
        std::cout << "are calculated using the mixed mode solver" << std::endl;
        std::cout << "For further details on options, refers to the options explanations below " << std::endl;
        std::cout << desc << std::endl;
        return 1;
    }
    // Read the configuration file
    std::unordered_map<std::string, std::string> input_params=readParameterFile(cfg_file);

    cfg_star=configure_make_star(input_params);
    // ---- Make the star ----
	params_out=make_synthetic_asymptotic_star(cfg_star);
    // Organize the paramaters 
	mode_params=bumpoutputs_2_MatrixXd_with_aj(params_out, cfg_star.inclination); // get the output in a format that can be written with the writting function
	copy_cfg(cfg_file, input_params["file_out_modes"]); // We include the full configuration without comments inside the output file
    write_range_modes(cfg_star, params_out, input_params["file_out_modes"],true);
    write_star_l1_roots(params_out, input_params["file_out_modes"], true);
    	//el, nu, h, w, a1, a2, a3, a4, a5, a6, asym, inc
	write_star_mode_params_asympt_model(mode_params, input_params["file_out_modes"], true);

	if (cfg_star.Dnu_star <= 15){
		std::cout << "    Model with small Dnu ==> many mixed modes. This might be long to find the solutions..." << std::endl; 
	}

	double tau=cfg_star.noise_params_harvey_like[3] * pow(cfg_star.numax_star*1e-6,cfg_star.noise_params_harvey_like[4]) + cfg_star.noise_params_harvey_like[5]; // Granulation timescale (in seconds)
	double H=cfg_star.noise_params_harvey_like[0] * pow(cfg_star.numax_star*1e-6,cfg_star.noise_params_harvey_like[1]) + cfg_star.noise_params_harvey_like[2]; // Granulation Amplitude
	H=H/tau ; //This is due to the used definition for the Harvey profile (conversion from Hz to microHz)
	tau=tau/1000. ; //conversion in ksec
	noise_params(0,0)=-1;
	noise_params(0,1)=-1;
	noise_params(0,2)=-1; 
	noise_params(1,0)=H;
	noise_params(1,1)=tau;
	noise_params(1,2)=cfg_star.noise_params_harvey_like[6]; // power law:  MUST BE CLOSE TO 2
	noise_params(2, 0)=cfg_star.noise_params_harvey_like[7]; // White noise
	noise_params(2, 1)=-2;
	noise_params(2, 2)=-2;
	// A FUNCTION THAT WRITES THE Noise
	write_star_noise_params(noise_params, input_params["file_out_modes"], true);
}