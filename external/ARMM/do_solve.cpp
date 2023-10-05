/*
   A small function that compute the results of the solver as it is implemented into models.cpp for
   mixed modes models. Mostly for debug purpose
*/

#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/program_options.hpp>
#include <unordered_map>

#include "string_handler.h"
#include "noise_models.h"
#include "version_solver.h"
#include "solver_mm.h"
#include "bump_DP.h"
#include "readparams_job.h"

void showversion();
void usage();

namespace po = boost::program_options;

int main(int argc, char* argv[]){

		bool verbose;
		Data_Nd data;
		std::string filename_output, filename_params;
		std::string char0, line0;
		std::ifstream cfg_session;
		const long double sigma_p_l1=0.0;
		int l;
		double DPl, epsilon_g, delta0l, q_star, resol, fmin ,fmax;
		VectorXd fl0, h1_h0_ratio, zeta_pg;
		Data_eigensols freqs_l1;	

		po::options_description desc("Allowed options");
		desc.add_options()
				("help,H", "produce help message")
				("version,V", "show version")
				("verbose", po::value<bool>()->default_value(true),"If true, shows details on performed steps and configuration")
				("filename_in,I", po::value<std::string>(), "input filename")
				("filename_out,O", po::value<std::string>()->default_value("results.txt"), "output filename");
		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);
		if (vm.count("version")) {
			showversion();
			exit(EXIT_SUCCESS);
		}
		if (vm.count("help") || !vm.count("filename_in") || !vm.count("filename_out")) {
			usage();
			std::cout << "" << std::endl;
			std::cout << desc << std::endl;
			exit(EXIT_SUCCESS);
		}
		verbose = vm["verbose"].as<bool>();
		filename_params = vm["filename_in"].as<std::string>();
		filename_output = vm["filename_out"].as<std::string>();

		if (verbose == true){
			std::cout << "  0. Configuration: " << std::endl;
			std::cout << "      - Parameters file: " << filename_params << std::endl;
			std::cout << "      - Result file: " << filename_output << std::endl;

			std::cout << "  1. Reading the file with the parameters for the solver_mm.cpp::solve_mm_asymptotic_O2from_l0() function..." << std::endl;
		}
		std::unordered_map<std::string, std::string> parameters = readParameterFile(filename_params);
		// Assign the values to corresponding variables
		l = std::stoi(parameters["l"]);
		DPl = std::stod(parameters[	"DPl"]);
		epsilon_g = std::stod(parameters["epsilon_g"]);
		delta0l = std::stod(parameters["delta0l"]);
		q_star = std::stod(parameters["q_star"]);
		fmin = std::stod(parameters["fmin"]);
		fmax = std::stod(parameters["fmax"]);
		resol = std::stod(parameters["resol"]);
		// Parse fl0 as a comma-separated list into an Eigen::VectorXd
		fl0=str_to_Xdarr(parameters["fl0"], " \t");

		if (verbose == true){
			std::cout << "     Summary of the requested configuration: " << std::endl;
			std::cout << " step = " << resol << std::endl; 
			std::cout << " fmin = " << fmin << std::endl;
			std::cout << " fmax = " << fmax << std::endl;
			std::cout << " l = " << l << std::endl;
			std::cout << " DPl = " << DPl << std::endl;
			std::cout << " alpha_g = " << epsilon_g << std::endl;
			std::cout << " delta0l = " << delta0l << std::endl;
			std::cout << " q_star = " << q_star << std::endl;
			std::cout << " fl0 = " << fl0.transpose() << std::endl;
			std::cout << "---" << std::endl;
			std::cout << "  3. Computing l=" << l << " frequencies..." << std::endl;
		}
		freqs_l1=solve_mm_asymptotic_O2from_l0(fl0, l, delta0l, DPl, epsilon_g, q_star, sigma_p_l1, resol, true, false, fmin, fmax);

		// Generating widths profiles for l modes using the ksi function
    	if (verbose == true){
			std::cout << "  2. Computing l=" << l << " heights through the zeta function (see code for assumptions)..." << std::endl;
		}
		zeta_pg=ksi_fct2(freqs_l1.nu_m, freqs_l1.nu_p, freqs_l1.nu_g, freqs_l1.dnup, freqs_l1.dPg, q_star, "precise"); //"precise" // assume Dnu_p, DPl and q constant   	
    	h1_h0_ratio=h_l_rgb(zeta_pg); // WARNING: Valid assummption only not too evolved RGB stars (below the bump, see Kevin mail 10 August 2019)
    	
		std::ofstream outputFile(filename_output);

		if (verbose == true){
			std::cout << " 3. Writing Results on file " << filename_output << "... " << std::endl;
		}
		outputFile << " # nu_p    : Frequencies of pure p modes" << std::endl;
		outputFile << " # dnu_p   : Derivative of pure p modes dnup/dn" << std::endl;
		outputFile << " # nu_g    : Frequencies of pure g modes" << std::endl;
		outputFile << " # DPg     : Derivative of pure g modes dPg/dn" << std::endl;
		outputFile << " # nu_m    : Frequencies of mixed modes, found by the solver" << std::endl;
		outputFile << " # zeta_pg : Zeta function at the nu_m frequencies" << std::endl;
		outputFile << " # Hl/H0   : Height ratio between l=1 and l=0. It does not account of mode visibility" << std::endl;
		outputFile << " # Summary of the requested configuration: " << std::endl;
		outputFile << "! step = " << resol << std::endl; 
		outputFile << "! fmin = " << fmin << std::endl;
		outputFile << "! fmax = " << fmax << std::endl;
		outputFile << "! l = " << l << std::endl;
		outputFile << "! DPl = " << DPl << std::endl;
		outputFile << "! alpha_g = " << epsilon_g << std::endl;
		outputFile << "! delta0l = " << delta0l << std::endl;
		outputFile << "! q_star = " << q_star << std::endl;
		outputFile << "! fl0 = " << fl0.transpose() << std::endl;
		outputFile <<  "nu_p = " << freqs_l1.nu_p.transpose() << std::endl;
		outputFile <<  "dnu_p = " << freqs_l1.dnup.transpose() << std::endl;
		outputFile <<  "nu_g = " << freqs_l1.nu_g.transpose() << std::endl;
		outputFile <<  "DPg = " <<freqs_l1.dPg.transpose() << std::endl;
		outputFile <<  "nu_m = " << freqs_l1.nu_m.transpose() << std::endl;
		outputFile <<  "zeta_pg = " << zeta_pg.transpose() << std::endl;
		outputFile <<  "H1/H0 = " << h1_h0_ratio.transpose() << std::endl;
		outputFile.close();
		if (verbose == true){
			std::cout << "Finished." << std::endl;
		}
}


void showversion()
{
    std::cout << APP_NAME " ARMM solver " APP_VERSION "\n built on " __DATE__ << std::endl;

#   if defined(__clang__)
    	printf(" with clang " __clang_version__);
#   elif defined(__GNUC__)
    	printf(" with GCC");
    	printf(" %d.%d.%d", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#   elif defined(_MSC_VER)
    	printf(" with MSVC");
    	printf(" %d", MSVC_VERSION);
#   else
    printf(" unknown compiler");
#   endif

    std::cout << "\n features:";
#   if defined(__i386__) || defined(_M_IX86)
    std::cout << " i386" << std::endl;
#   elif defined(__x86_64__) || defined(_M_AMD64)
    std::cout << " x86_64" << std::endl;
#   endif
    std::cout << " Author: " << APP_COPYRIGHT << std::endl;

}

void usage() {
    std::cout << " This program allows you to use the solver as stand-alone program (in command line) " << std::endl;
    std::cout << " You need to provide at least one arguments to the program" << std::endl; 
	std::cout << " This is a configuration file with, " << std::endl;
    std::cout << "    (0) Comments on as many lines as whished ('#'), this files contains one parameter per line. Followed in this exact order by one parameter set per line: " << std::endl;
    std::cout << "    (1) l" << std::endl;
    std::cout << "    (2) DP" <<  std::endl;
    std::cout << "    (3) alpha_g" << std::endl;
    std::cout << "    (4) delta0l" << std::endl;
    std::cout << "    (5) q_star" << std::endl;
    std::cout << "    (6) fl0, a list of l=0 p mode frequencies that will be shift by l*Dnu/2 + delta0l to produce l=1 p modes" << std::endl;
    std::cout << "    (7) step / resolution of the spectrum to define the precision on the retrieve modes" <<  std::endl;
    std::cout << "    (8) fmin: Remove solutions found below fmin. Note that fl0 defines the total number of calculated frequencies. Those may/may not be below fmin" <<  std::endl;
    std::cout << "    (9) fmax: Remove solutions found above fmax. Note that fl0 defines the total number of calculated frequencies. Those may/may not be above fmax" <<  std::endl;
}
