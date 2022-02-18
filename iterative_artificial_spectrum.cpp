/*
 * iterative_artificial_spectrum.cpp
 *
 * Header file that contains all kind of methods
 * used to generate models for the pulsation/noise
 * 
 *  Created on: 05 May 2016
 *      Author: obenomar
 */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include "artificial_spectrum.h"
#include "models_database.h"
#include "models_database_grid.h"
#include "version.h"
#include "combi.h"
#include "stellar_models.h"
#include "io_star_params.h"

void showversion();
int  options(int argc, char* argv[]);

void iterative_artificial_spectrum(std::string dir_core);
//VectorXd order_input_params(VectorXd cte_params, VectorXd var_params, std::vector<std::string> cte_names, 
//		std::vector<std::string> var_names, std::vector<std::string> param_names);
void generate_random(Config_Data cfg, std::vector<std::string> param_names, std::string dir_core, std::string file_out_modes, 
		std::string file_out_noise, std::string file_out_combi, int N_model, std::string file_cfg_mm, std::string external_path,  std::string templates_dir);
void generate_grid(Config_Data cfg, bool usemodels, Data_Nd models, std::vector<std::string> param_names, std::string dir_core, std::string dir_freqs, std::string file_out_modes, 
		std::string file_out_noise, std::string file_out_combi, int N_model);

bool call_model_random(std::string model_name, VectorXd input_params, std::string file_out_modes, std::string file_out_noise, 
		 std::string file_cfg_mm, std::string dir_core, std::string identifier, Config_Data cfg, std::string external_path, std::string template_file);
bool call_model_grid(std::string model_name, VectorXd input_params, Model_data input_model, std::string file_out_modes, std::string file_out_noise, 
		std::string dir_core, std::string id_str, Config_Data cfg);

std::vector<std::string> list_dir(const std::string path, const std::string filter);


void iterative_artificial_spectrum(std::string dir_core){
/*
 * Run iteratively the program artificial_spectrum in order to generate a
 * serie of model according to a main configuration file (/Configurations/main.cfg).
*/

	bool passed;
	bool usemodels=0;
	bool verbose_data=0; // SET TO 1 IF YOU WANT MORE INFO ABOUT WHAT IS READ
	int Nmodel;
	std::vector<std::string> param_names;

	std::string delimiter=" ";
	std::string cfg_file, model_file, file_out_modes, file_out_noise, file_cfg_mm, file_out_mm, file_out_mm2, file_out_combi;
	std::string external_path, templates_dir;
	Config_Data cfg;
	Data_Nd models;
	external_path=dir_core + "external/"; 
	templates_dir=dir_core + "Configurations/templates/";
	cfg_file=dir_core + "Configurations/main.cfg";
	model_file=dir_core + "external/MESA_grid/models.params";
	file_out_modes=dir_core + "Configurations/tmp/modes_tmp.cfg";
	file_out_noise=dir_core + "Configurations/tmp/noise_tmp.cfg";
	file_cfg_mm=dir_core + "external/ARMM-solver/star_params.global"; // Used only for models with mixed modes and interfaced with the Python solver
	
	file_out_combi=dir_core + "Data/Combinations.txt";

	std::string dir_freqs=dir_core + "external/MESA_grid/frequencies/";

	std::cout << "1. Read the configuration file..." << std::endl;

	cfg=read_main_cfg(cfg_file);

	std::cout << "---------------------------------------------------------------------------------------" << std::endl;
	std::cout << "      ATTENTION: All model configuration should be contained into models_database.cpp  " << std::endl;
	std::cout << "                 All mode generators should be contained into build_lorentzian.cpp     " << std::endl;
	std::cout << "                 All model generators should be contained into artificial_spectrum.cpp " << std::endl;
	std::cout << "---------------------------------------------------------------------------------------" << std::endl;

	passed=0;
	if(cfg.model_name == "generate_cfg_asymptotic_act_asym_Hgauss"){
		Nmodel=21;
		param_names.push_back("numax"); param_names.push_back("Dnu"); param_names.push_back("epsilon"); param_names.push_back("D0"); param_names.push_back("maxH"); 
		param_names.push_back("Gamma"); param_names.push_back("lmax"); param_names.push_back("Nmax"); param_names.push_back("a1"); param_names.push_back("a3"); 
		param_names.push_back("b"); param_names.push_back("alfa"); param_names.push_back("beta"); param_names.push_back("i"); 
		param_names.push_back("Hnoise1"); param_names.push_back("tau1"); param_names.push_back("p1");
		param_names.push_back("Hnoise2"); param_names.push_back("tau2"); param_names.push_back("p2");  
		param_names.push_back("N0");
		if(param_names.size() != Nmodel){
			std::cout << "    Invalid number of parameters for model_name= 'generate_cfg_asymptotic_act_asym'" << std::endl;
			std::cout << "    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size() << std::endl;
			std::cout << "    Check your main configuration file" << std::endl;
			std::cout << "    The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}
		passed=1;
	}
	if(cfg.model_name == "generate_cfg_from_synthese_file_Wscaled_a1a2a3asymovGamma"){
		Nmodel=11;
		param_names.push_back("HNR"); 
		param_names.push_back("a1ovGamma"); 
		param_names.push_back("Gamma_at_numax"); 
		param_names.push_back("a2_0"); 
		param_names.push_back("a2_1"); 
		param_names.push_back("a2_2"); 
		param_names.push_back("a3_0"); 
		param_names.push_back("a3_1"); 
		param_names.push_back("a3_2"); 
		param_names.push_back("beta_asym");
		param_names.push_back("i");
		if(param_names.size() != Nmodel){
			std::cout << "    Invalid number of parameters for model_name= 'generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma'" << std::endl;
			std::cout << "    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size() << std::endl;
			std::cout << "    Check your main configuration file" << std::endl;
			std::cout << "    The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}
		passed=1;
	}
	if(cfg.model_name == "generate_cfg_from_synthese_file_Wscaled_Alm"){
		Nmodel=9;
		param_names.push_back("HNR"); 
		param_names.push_back("a1ovGamma"); 
		param_names.push_back("Gamma_at_numax"); 
		param_names.push_back("epsilon_nl"); 
		param_names.push_back("theta0"); 
		param_names.push_back("delta"); 
		param_names.push_back("a3"); 
		param_names.push_back("beta_asym");
		param_names.push_back("i");
		if(param_names.size() != Nmodel){
			std::cout << "    Invalid number of parameters for model_name= 'generate_cfg_from_synthese_file_Wscaled_a1a2a3asymovGamma'" << std::endl;
			std::cout << "    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size() << std::endl;
			std::cout << "    Check your main configuration file" << std::endl;
			std::cout << "    The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}
		passed=1;
	}
	if(cfg.model_name == "generate_cfg_from_synthese_file_Wscaled_aj"){
		Nmodel=10;
		param_names.push_back("HNR"); 
		param_names.push_back("a1ovGamma"); 
		param_names.push_back("Gamma_at_numax"); 
		param_names.push_back("a2"); 
		param_names.push_back("a3"); 
		param_names.push_back("a4"); 
		param_names.push_back("a5"); 
		param_names.push_back("a6"); 
		param_names.push_back("beta_asym");
		param_names.push_back("i");
		if(param_names.size() != Nmodel){
			std::cout << "    Invalid number of parameters for model_name= 'generate_cfg_from_synthese_file_Wscaled_a1a2a3asymovGamma'" << std::endl;
			std::cout << "    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size() << std::endl;
			std::cout << "    Check your main configuration file" << std::endl;
			std::cout << "    The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}
		passed=1;
	}
	if(cfg.model_name == "generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma"){
		Nmodel=6;
		param_names.push_back("HNR"); 
		param_names.push_back("a1ovGamma"); 
		param_names.push_back("Gamma_at_numax"); 
		param_names.push_back("a3"); 
		param_names.push_back("beta_asym");
		param_names.push_back("i");
		if(param_names.size() != Nmodel){
			std::cout << "    Invalid number of parameters for model_name= 'generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma'" << std::endl;
			std::cout << "    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size() << std::endl;
			std::cout << "    Check your main configuration file" << std::endl;
			std::cout << "    The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}
		passed=1;
	}
	if(cfg.model_name == "asymptotic_mm_v1"){
		//Nmodel=7;
		param_names.push_back("Teff"); 
		param_names.push_back("Dnu"); 
		param_names.push_back("epsilon"); 
		param_names.push_back("delta0l_percent"); 
		param_names.push_back("beta_p_star"); 
		param_names.push_back("nmax_spread"); 
		param_names.push_back("DP_var_percent"); 
		param_names.push_back("alpha"); 
		param_names.push_back("q");
		param_names.push_back("SNR");
		param_names.push_back("maxGamma");
		param_names.push_back("Vl1");
		param_names.push_back("Vl2");
		param_names.push_back("Vl3");
		param_names.push_back("H0_spread");
		param_names.push_back("A_Pgran");	
		param_names.push_back("B_Pgran");	
		param_names.push_back("C_Pgran");	
		param_names.push_back("A_taugran");	
		param_names.push_back("B_taugran");	
		param_names.push_back("C_taugran");	
		param_names.push_back("P");	
		param_names.push_back("N0");			
		Nmodel=param_names.size();
		if(param_names.size() != Nmodel){
			std::cout << "    Invalid number of parameters for model_name= 'asymptotic_mm_v1'" << std::endl;
			std::cout << "    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size() << std::endl;
			std::cout << "    Check your main configuration file" << std::endl;
			std::cout << "    The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}
		passed=1;
	}
	if(cfg.model_name == "asymptotic_mm_v2"){
		//Nmodel=12;
		param_names.push_back("nurot_env"); 
		param_names.push_back("nurot_ratio"); 
		param_names.push_back("Dnu"); 
		param_names.push_back("epsilon"); 
		param_names.push_back("delta0l_percent"); 
		param_names.push_back("beta_p_star"); 
		param_names.push_back("nmax_spread"); 
		param_names.push_back("DP_var_percent"); 
		param_names.push_back("alpha"); 
		param_names.push_back("q");
		param_names.push_back("SNR");
		param_names.push_back("maxGamma");
		param_names.push_back("Vl1");
		param_names.push_back("Vl2");
		param_names.push_back("Vl3");
		param_names.push_back("H0_spread");	
		param_names.push_back("A_Pgran");	
		param_names.push_back("B_Pgran");	
		param_names.push_back("C_Pgran");	
		param_names.push_back("A_taugran");	
		param_names.push_back("B_taugran");	
		param_names.push_back("C_taugran");	
		param_names.push_back("P");	
		param_names.push_back("N0");
		Nmodel=param_names.size();
		if(param_names.size() != Nmodel){
			std::cout << "    Invalid number of parameters for model_name= 'asymptotic_mm_v2'" << std::endl;
			std::cout << "    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size() << std::endl;
			std::cout << "    Check your main configuration file" << std::endl;
			std::cout << "    The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}
		passed=1;
	}
	if(cfg.model_name == "asymptotic_mm_v3"){
		//Nmodel=12;
		param_names.push_back("nurot_env"); 
		param_names.push_back("nurot_core"); 
		param_names.push_back("Dnu"); 
		param_names.push_back("epsilon"); 
		param_names.push_back("delta0l_percent"); 
		param_names.push_back("beta_p_star"); 
		param_names.push_back("nmax_spread"); 
		param_names.push_back("DP_var_percent"); 
		param_names.push_back("alpha"); 
		param_names.push_back("q");
		param_names.push_back("SNR");
		param_names.push_back("maxGamma");
		param_names.push_back("Vl1");
		param_names.push_back("Vl2");
		param_names.push_back("Vl3");
		param_names.push_back("H0_spread");	
		param_names.push_back("A_Pgran");	
		param_names.push_back("B_Pgran");	
		param_names.push_back("C_Pgran");	
		param_names.push_back("A_taugran");	
		param_names.push_back("B_taugran");	
		param_names.push_back("C_taugran");	
		param_names.push_back("P");	
		param_names.push_back("N0");
		Nmodel=param_names.size();
		if(param_names.size() != Nmodel){
			std::cout << "    Invalid number of parameters for model_name= 'asymptotic_mm_v3'" << std::endl;
			std::cout << "    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size() << std::endl;
			std::cout << "    Check your main configuration file" << std::endl;
			std::cout << "    The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}
		passed=1;
	}
	if(cfg.model_name == "asymptotic_mm_freeDp_numaxspread_curvepmodes_v1"){
		//Nmodel=12;
		param_names.push_back("Teff"); 
		param_names.push_back("Dnu"); 
		param_names.push_back("epsilon"); 
		param_names.push_back("delta0l_percent"); 
		param_names.push_back("beta_p_star"); 
		param_names.push_back("nmax_spread"); 
		param_names.push_back("DP1"); 
		param_names.push_back("alpha"); 
		param_names.push_back("q");
		param_names.push_back("SNR");
		param_names.push_back("maxGamma");
		param_names.push_back("numax_spread");
		param_names.push_back("Vl1");
		param_names.push_back("Vl2");
		param_names.push_back("Vl3");
		param_names.push_back("H0_spread");	
		param_names.push_back("A_Pgran");	
		param_names.push_back("B_Pgran");	
		param_names.push_back("C_Pgran");	
		param_names.push_back("A_taugran");	
		param_names.push_back("B_taugran");	
		param_names.push_back("C_taugran");	
		param_names.push_back("P");	
		param_names.push_back("N0");
		Nmodel=param_names.size();		
		if(param_names.size() != Nmodel){
			std::cout << "    Invalid number of parameters for model_name= 'asymptotic_mm_v1'" << std::endl;
			std::cout << "    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size() << std::endl;
			std::cout << "    Check your main configuration file" << std::endl;
			std::cout << "    The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}
		passed=1;
	}
	if(cfg.model_name == "asymptotic_mm_freeDp_numaxspread_curvepmodes_v2"){
		//Nmodel=13;
		param_names.push_back("nurot_env"); 
		param_names.push_back("nurot_ratio"); 
		param_names.push_back("Dnu"); 
		param_names.push_back("epsilon");
		param_names.push_back("delta0l_percent"); 
		param_names.push_back("beta_p_star"); 
		param_names.push_back("nmax_spread"); 
		param_names.push_back("DP1"); 
		param_names.push_back("alpha"); 
		param_names.push_back("q");
		param_names.push_back("SNR");
		param_names.push_back("maxGamma");
		param_names.push_back("numax_spread");
		param_names.push_back("Vl1");
		param_names.push_back("Vl2");
		param_names.push_back("Vl3");
		param_names.push_back("H0_spread");	
		param_names.push_back("A_Pgran");	
		param_names.push_back("B_Pgran");	
		param_names.push_back("C_Pgran");	
		param_names.push_back("A_taugran");	
		param_names.push_back("B_taugran");	
		param_names.push_back("C_taugran");	
		param_names.push_back("P");	
		param_names.push_back("N0");
		Nmodel=param_names.size();		
		if(param_names.size() != Nmodel){
			std::cout << "    Invalid number of parameters for model_name= 'asymptotic_mm_v2'" << std::endl;
			std::cout << "    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size() << std::endl;
			std::cout << "    Check your main configuration file" << std::endl;
			std::cout << "    The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}
		passed=1;
	}
	if(cfg.model_name == "asymptotic_mm_freeDp_numaxspread_curvepmodes_v3"){
		//Nmodel=13;
		param_names.push_back("nurot_env"); 
		param_names.push_back("nurot_core"); 
		param_names.push_back("Dnu"); 
		param_names.push_back("epsilon");
		param_names.push_back("delta0l_percent"); 
		param_names.push_back("beta_p_star"); 
		param_names.push_back("nmax_spread"); 
		param_names.push_back("DP1");  
		param_names.push_back("alpha"); 
		param_names.push_back("q");
		param_names.push_back("SNR");
		param_names.push_back("maxGamma");
		param_names.push_back("numax_spread");	
		param_names.push_back("Vl1");
		param_names.push_back("Vl2");
		param_names.push_back("Vl3");
		param_names.push_back("H0_spread");	
		param_names.push_back("A_Pgran");	
		param_names.push_back("B_Pgran");	
		param_names.push_back("C_Pgran");	
		param_names.push_back("A_taugran");	
		param_names.push_back("B_taugran");	
		param_names.push_back("C_taugran");	
		param_names.push_back("P");	
		param_names.push_back("N0");
		Nmodel=param_names.size();		
		if(param_names.size() != Nmodel){
			std::cout << "    Invalid number of parameters for model_name= 'asymptotic_mm_v3'" << std::endl;
			std::cout << "    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size() << std::endl;
			std::cout << "    Check your main configuration file" << std::endl;
			std::cout << "    The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}
		passed=1;
	}

	// -----------------------------------------------------------
	// ----- Models available only with the grid approach -----
	// -----------------------------------------------------------
	if(cfg.model_name == "generate_cfg_from_refstar_HWscaled"){
		//Nmodel=12;
		param_names.push_back("Model_i"); param_names.push_back("maxHNR"); param_names.push_back("Gamma_maxHNR"); param_names.push_back("a1"); param_names.push_back("i"); 
		param_names.push_back("Hnoise1"); param_names.push_back("tau1"); param_names.push_back("p1");
		param_names.push_back("Hnoise2"); param_names.push_back("tau2"); param_names.push_back("p2");  
		param_names.push_back("N0");
		Nmodel=param_names.size();
		if(param_names.size() != Nmodel){
			std::cout << "    Invalid number of parameters for model_name= 'generate_cfg_asymptotic_act_asym'" << std::endl;
			std::cout << "    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size() << std::endl;
			std::cout << "    Check your main configuration file" << std::endl;
			std::cout << "    The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}
		std::cout << " 	   - reading models from model file (requires for setting frequencies)..." << std::endl;
		models=read_data_ascii_Ncols(model_file, delimiter, verbose_data);
		usemodels=1;
		passed=1;
	}
	if(cfg.model_name == "generate_cfg_from_refstar_HWscaled_GRANscaled"){
		//Nmodel=12;
		param_names.push_back("Model_i"); param_names.push_back("maxHNR"); param_names.push_back("Gamma_maxHNR"); param_names.push_back("a1"); param_names.push_back("i"); 
		param_names.push_back("A_Pgran"); param_names.push_back("B_Pgran"); param_names.push_back("C_Pgran");
		param_names.push_back("A_taugran"); param_names.push_back("B_taugran"); param_names.push_back("C_taugran");  
		param_names.push_back("N0");
		Nmodel=param_names.size();
		if(param_names.size() != Nmodel){
			std::cout << "    Invalid number of parameters for model_name= 'generate_cfg_asymptotic_act_asym'" << std::endl;
			std::cout << "    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size() << std::endl;
			std::cout << "    Check your main configuration file" << std::endl;
			std::cout << "    The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}
		std::cout << " 	   - reading models from model file (requires for setting frequencies)..." << std::endl;
		models=read_data_ascii_Ncols(model_file, delimiter, verbose_data);
		usemodels=1;
		passed=1;
	}
 	if(passed == 0){
		std::cout << "    model_name= " << cfg.model_name << " is not a recognized keyword for models" << std::endl;
		std::cout << "    Check models_database.h to see the available model names" << std::endl;
		std::cout << "    The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

	// -----------------------------------------------------------
	// -----------------------------------------------------------
	// -----------------------------------------------------------


	std::cout << "2. Generating the models using the subroutine " << cfg.model_name << " of model_database.cpp..." << std::endl;
	if(cfg.forest_type == "random"){
		std::cout << "   Values are randomly generated into a uniform range defined in the main configuration file" << std::endl;
		generate_random(cfg, param_names, dir_core, file_out_modes, file_out_noise, file_out_combi, Nmodel, file_cfg_mm, external_path, templates_dir);

	}
	if(cfg.forest_type == "grid"){
		std::cout << "   Values are generated over a grid using all possible combinations according to inputs in the main configuration file" << std::endl;
		//generate_grid();
		generate_grid(cfg, usemodels, models, param_names, dir_core, dir_freqs, file_out_modes, file_out_noise, file_out_combi, Nmodel);
	}
	if(cfg.forest_type != "random" && cfg.forest_type != "grid"){ 
		std::cout << " Problem in the main configuration file. It is expected that the forest type parameters is either random OR grid" << std::endl;
		std::cout << " Check your main configuration file" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

}

void generate_random(Config_Data cfg, std::vector<std::string> param_names, std::string dir_core, std::string file_out_modes, 
		std::string file_out_noise, std::string file_out_combi, int N_model,  std::string file_cfg_mm, std::string external_path, std::string templates_dir){

	bool neg=0, passed=0;
	int i, ii;//, xmin_int_rgen, xmax_int_rgen;
	long lastid, id0;
	std::string id_str, str_tmp;
	std::string template_file;
	VectorXd cte_params, val_min, val_max, input_params;
	MatrixXd currentcombi, allcombi;
	std::vector<double> pos_zero, pos_one;	
	std::vector<std::string> var_names, cte_names;
	// We first check that the cfg file has a coherent setup
	ii=0;
	while(neg == 0 && (ii < N_model)){
		//std::cout << "[" << ii << "] " << cfg.val_max[ii] - cfg.val_min[ii] << std::endl;
	
		if( (cfg.val_max[ii] - cfg.val_min[ii]) < 0){ 
			neg=1;
		}
		ii=ii+1;
	}
	if(neg == 1){
		std::cout << "     Warning: val_max < val_min for some of the parameters!" << std::endl;
		std::cout << "     This is not permitted" << std::endl;
		std::cout << "     Check your main configuration file" << std::endl;
		std::cout << "     The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

	//  Define variables and constants
	pos_one=where(cfg.step, "=", 1, 0); // All positions of cfg.step that are equal to 1. Return position (last parameter is 0)
	pos_zero=where(cfg.step, "=", 0, 0); // All positions of cfg.step that are equal to 1. Return position (last parameter is 0)

	val_min.resize(pos_one.size());
	val_max.resize(pos_one.size());
	cte_params.resize(pos_zero.size());
	for(int i=0; i<pos_one.size(); i++){
		var_names.push_back(cfg.labels[pos_one[i]]);
		val_min[i]=cfg.val_min[pos_one[i]];
		val_max[i]=cfg.val_max[pos_one[i]];
	
	}
	for(int i=0; i<pos_zero.size(); i++){
		cte_names.push_back(cfg.labels[pos_zero[i]]);
		cte_params[i]=cfg.val_min[pos_zero[i]];
	}
	std::cout << "Constants: ";
	for(int i=0; i<cte_names.size(); i++){ std::cout << cte_names[i] << "  ";}
	std::cout << std::endl;

	std::cout << "Variables: ";
	for(int i=0; i<var_names.size(); i++){ std::cout << var_names[i] << "  ";}
	std::cout << std::endl;

	//  ------------ Generate the random values -----------
	// Initialize the random generators: Uniform over [0,1]. This is used for the random parameters
    boost::mt19937 generator(time(NULL));
    boost::uniform_01<boost::mt19937> dist_gr(generator);
 	
    // Generator of integers for selecting ramdomly a template file that contains a height and width profile
    switch(cfg.template_files.size()){
    	case 0: // The string is somewhat empty
    		std::cout << "Error: The template_file cannot be empty. If not used, please set it to NONE" << std::endl;
    		exit(EXIT_SUCCESS);
    	case 1:
    		template_file=templates_dir + cfg.template_files[0];
    	default:
    		if(cfg.template_files[0] == "all" || cfg.template_files[0] == "ALL" || cfg.template_files[0] == "*"){ // If the user specify that all *.template files should be used
    			cfg.template_files=list_dir(templates_dir, "template");
    		} 
    		if (cfg.template_files.size() == 0){
    			std::cout << "Could not find the template file in the expected directory" << std::endl;
    			std::cout << "Be sure to have it in Configurations/templates/ " << std::endl;
    			std::cout << "If the used mode does not require templates, specify NONE in the dedicated field " << std::endl;
    			exit(EXIT_FAILURE);
    		}
 			boost::random::uniform_int_distribution<> dist(0, cfg.template_files.size()-1);
//   			std::cout << "cfg.template_files =" << std::endl;
//    		for (int ii=0; ii< cfg.template_files.size(); ii++){
//    			std::cout << cfg.template_files[ii] << std::endl;
//    		}
    		template_file=templates_dir + cfg.template_files[dist(generator)];
    }	
    std::cout << "Selected template file: " << template_file << std::endl;	
//    exit(EXIT_SUCCESS);

 	std::cout << " -----------------------------------------------------" << std::endl;
	std::cout << " List of all combinations written iteratively into " << std::endl;
	std::cout << "       " <<  file_out_combi << std::endl; 
	std::cout << " -----------------------------------------------------" << std::endl;


	if(cfg.erase_old_files == 0){
		if(file_exists(file_out_combi) == 1){
			std::cout << "                 erase_old_files=0..." << std::endl;
			std::cout << "                 ...Older combinatory file found!" << std::endl;
			std::cout << "                 Name of the found file: " << file_out_combi << std::endl;
			std::cout << "                 Reading the combinatory file in order to determince the last value of the samples..." << std::endl;
			lastid=read_id_allcombi(file_out_combi);

			//exit(EXIT_SUCCESS);

		} else{
			lastid=0; // If no Combi file while erase_old_files is set to 1
			std::cout << "                 erase_old_files=0..."<< std::endl;
			std::cout << "                 ... but no older combi file was found" << std::endl;
			std::cout << "                 The program will therefore behave as if erase_old_files=1" << std::endl;
		}
	} else {
		lastid=0; // If no Combi file while erase_old_files is set to 1
		std::cout << "                 erase_old_files=1..."<< std::endl;
		std::cout << "                 Note that no deleting action is performed by this program" << std::endl;
		std::cout << "                 But any older file might be overwritten!" << std::endl;
		std::cout << "                 Be sure to have no important previous results in the Data directory and subdirectories!" << std::endl;
		//std::cout << "                  The program will start now" << std::endl;
		//sleep(1000*5); // give 5 seconds to kill the process
	}
	id0=lastid+1;

	//allcombi.resize(cfg.forest_params[0], pos_one.size());
	currentcombi.resize(1, pos_one.size());
	for(int c=0; c<cfg.forest_params[0]; c++){
		for(int i=0; i<pos_one.size();i++){
			//std::cout << dist_gr() << " " << std::endl;
			currentcombi(0,i)=dist_gr() * (val_max[i] - val_min[i]) + val_min[i]; // HERE allcombi HAS JUST ONE LINE
			//allcombi(c,i)=currentcombi(0,i);
		}
		id_str=write_allcombi(currentcombi, cte_params, cfg, file_out_combi, cfg.erase_old_files, c, id0, cte_names, var_names,  param_names); // update/write the combination file
		std::cout << "Combination number: " << id_str << std::endl;

		input_params=order_input_params(cte_params, currentcombi.row(0), cte_names, var_names, param_names);

		passed=call_model_random(cfg.model_name, input_params, file_out_modes, file_out_noise,  file_cfg_mm, dir_core, id_str, cfg, external_path, template_file);
		if(passed == 0){
			std::cout << "Warning: The function call_model did not generate any configuration file!" << std::endl;
			std::cout << "         It is very likely that you tried to start a model in 'random mode' while this model is not available for the random approach" << std::endl;
			std::cout << "         The program will stop now" << std::endl;
			exit(EXIT_FAILURE);
		}
		id0=id0+1;
		passed=0;

		// Erasing the temporary files to avoid those to be loaded at vitam eternam if an error in their generetation at step c=c_err happens
		str_tmp="rm " + file_out_modes;	
		const char *cmd1 = str_tmp.c_str(); 
		system(cmd1);
		str_tmp="rm " + file_out_noise;	
		const char *cmd2 = str_tmp.c_str(); 
		system(cmd2);
	}

	std::cout << "All requested tasks executed succesfully" << std::endl;
	
}

void generate_grid(Config_Data cfg, bool usemodels, Data_Nd models, std::vector<std::string> param_names, std::string dir_core, std::string dir_freqs, std::string file_out_modes, 
		std::string file_out_noise, std::string file_out_combi, int N_model){

	bool neg=0, passed=0;
	int i, ii, Nvar_params;
	long lastid, id0, Ncombi;
	std::string id_str, str_tmp;
	VectorXi pos, Nvals; // May have to use a Xd when dealing with the round off
	VectorXd tmp, cte_params, val_min, val_max, steps, input_params, delta;
	MatrixXd var_params, currentcombi, allcombi;
	std::vector<double> pos_zero, pos_one;	
	std::vector<std::string> var_names, cte_names;

	Model_data input_model;

	// We first check that the cfg file has a coherent setup
	delta.resize(N_model);
	ii=0;
	while(neg == 0 && (ii < N_model)){
		//std::cout << "[" << ii << "] " << cfg.val_max[ii] - cfg.val_min[ii] << std::endl;
		delta[ii]=cfg.val_max[ii] - cfg.val_min[ii];
		if( delta[ii] < 0){ 
			neg=1;
		}
		ii=ii+1;
	}
	if(neg == 1){
		std::cout << "     Warning: val_max < val_min for some of the parameters!" << std::endl;
		std::cout << "     This is not permitted" << std::endl;
		std::cout << "     Check your main configuration file" << std::endl;
		std::cout << "     The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

	//  Define variables and constants
	pos_one=where(cfg.step, "!=", 0, 0); // All positions of cfg.step that are not equal to 0. Return position (last parameter is 0)
	pos_zero=where(cfg.step, "<=", 0, 0); // All positions of cfg.step that are equal or less than 0. Return position (last parameter is 0)

	Nvals.resize(pos_one.size());
	for(int ii=0; ii<pos_one.size();ii++){
		//std::cout << pos_one[ii] << std::endl;
		Nvals[ii]=delta[pos_one[ii]]/cfg.step[pos_one[ii]] + 1e-15; // 1d-15 here because floor bugs sometimes due to round-off
		//std::cout << "Nvals["<< ii << "]=" << Nvals[ii] << std::endl; //Nvals=Nvals.floor();
	}
	val_min.resize(pos_one.size());
	val_max.resize(pos_one.size());
	//steps.resize(pos_one.size());
	cte_params.resize(pos_zero.size());
	var_params.resize(pos_one.size(), Nvals.maxCoeff()+1);	
	for(int i=0; i<pos_one.size(); i++){
		var_names.push_back(cfg.labels[pos_one[i]]);
		val_min[i]=cfg.val_min[pos_one[i]];
		val_max[i]=cfg.val_max[pos_one[i]];
		tmp.setLinSpaced(Nvals[i]+1, val_min[i], val_max[i]);
		for(int k=0; k<tmp.size();k++){
			var_params(i,k)=tmp[k];
		}
		Nvals[i]=Nvals[i]+1; // Number of values
	}

	for(int i=0; i<pos_zero.size(); i++){
		cte_names.push_back(cfg.labels[pos_zero[i]]);
		cte_params[i]=cfg.val_min[pos_zero[i]];
	}
	std::cout << "Constants: ";
	for(int i=0; i<cte_names.size(); i++){ std::cout << cte_names[i] << "  ";}
	std::cout << std::endl;

	std::cout << "Variables: ";
	Nvar_params=var_names.size();
	for(int i=0; i<Nvar_params; i++){ std::cout << var_names[i] << "  ";}
	std::cout << std::endl;

	//  ------------ Generate the grid -----------
	std::cout << "       - Generating all the requested combinations. This may take a while..." << std::endl;
	var_params.transposeInPlace();
	allcombi=define_all_combinations(var_params, Nvals, Nvar_params);
	Ncombi=allcombi.rows();
	std::cout << allcombi << std::endl;

	std::cout << "Ncombi = " << Ncombi << std::endl;
	std::cout << " -----------------------------------------------------" << std::endl;
	std::cout << " List of all combinations written iteratively into " << std::endl;
	std::cout << "       " <<  file_out_combi << std::endl; 
	std::cout << " -----------------------------------------------------" << std::endl;

	if(cfg.erase_old_files == 0){
		if(file_exists(file_out_combi) == 1){
			std::cout << "                 erase_old_files=0..." << std::endl;
			std::cout << "                 ...Older combinatory file found!" << std::endl;
			std::cout << "                 Name of the found file: " << file_out_combi << std::endl;
			std::cout << "                 Reading the combinatory file in order to determince the last value of the samples..." << std::endl;
			lastid=read_id_allcombi(file_out_combi);

			//exit(EXIT_SUCCESS);

		} else{
			lastid=-1; // If no Combi file while erase_old_files is set to 1
			std::cout << "                 erase_old_files=0..."<< std::endl;
			std::cout << "                 ... but no older combi file was found" << std::endl;
			std::cout << "                 The program will therefore behave as if erase_old_files=0" << std::endl;
		}
	} else {
		lastid=-1; // If no Combi file while erase_old_files is set to 1
		std::cout << "                 erase_old_files=1..."<< std::endl;
		std::cout << "                 Note that no deleting action is performed by this program" << std::endl;
		std::cout << "                 But any older file might be overwritten!" << std::endl;
		std::cout << "                 Be sure to have no important previous results in the Data directory and subdirectories!" << std::endl;
	}
	id0=lastid+1;
	//id0=1; // DEBUG ONLY

	currentcombi.resize(1, pos_one.size());
	for(int c=0; c<Ncombi; c++){
		for(int i=0; i<pos_one.size();i++){
			currentcombi.row(0)=allcombi.row(id0); // HERE currentcombi HAS JUST ONE LINE
		}
		id_str=write_allcombi(currentcombi, cte_params, cfg, file_out_combi, cfg.erase_old_files, c, id0, cte_names, var_names,  param_names); // update/write the combination file
		std::cout << "Combination number: " << id_str  << "  (last is: " << Ncombi-1 << ")" << std::endl;
	
		input_params=order_input_params(cte_params, currentcombi.row(0), cte_names, var_names, param_names);
		if(usemodels == 1){
			pos=where_strXi(param_names, "Model_i"); // Look for the position where Model_i is given (should be 0)	
			//std::cout << pos << std::endl;
			if(pos[0] !=-1){
				input_model=get_model_param(models, input_params[pos[0]], dir_freqs);
			} else{
				std::cout << "Could not find the column that contains the 'Model_i' (model index)" << std::endl;
				std::cout << "Debug required. The program will exit now" << std::endl;
				exit(EXIT_SUCCESS);
			}
		}
		passed=call_model_grid(cfg.model_name, input_params, input_model, file_out_modes, file_out_noise, dir_core, id_str, cfg);
		if(passed == 0){
			std::cout << "Warning: The function call_model_grid did not generate any configuration file!" << std::endl;
			std::cout << "         It is very likely that you tried to start a model in 'grid mode' while this model is not available for the grid approach" << std::endl;
			std::cout << "         The program will stop now" << std::endl;
			exit(EXIT_FAILURE);
		}
		id0=id0+1;
		passed=0;

		// Erasing the temporary files to avoid those to be loaded at vitam eternam if an error in their generetation at step c=c_err happens
		str_tmp="rm " + file_out_modes;	
		const char *cmd1 = str_tmp.c_str(); 
		system(cmd1);
		str_tmp="rm " + file_out_noise;	
		const char *cmd2 = str_tmp.c_str(); 
		system(cmd2);
	}

	std::cout << "All requested tasks executed succesfully" << std::endl;
//exit(EXIT_SUCCESS);	
}



// This is the routine that calls models that are valid on the random approach
bool call_model_random(std::string model_name, VectorXd input_params, std::string file_out_modes, std::string file_out_noise, 
		std::string file_cfg_mm, std::string dir_core, std::string id_str, Config_Data cfg, std::string external_path, std::string template_file){

	std::string str;
	bool passed=0, subpassed=0;
	//std::string id_str;

	if(model_name == "generate_cfg_asymptotic_act_asym_Hgauss"){
		generate_cfg_asymptotic_act_asym_Hgauss(input_params, file_out_modes, file_out_noise);
		artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel);
		passed=1;
	}
	if(model_name =="generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma"){
		generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma(input_params, file_out_modes,  file_out_noise, cfg.extra_params); // extra_params must points towards a .in file
		artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel);
		passed=1;		
	}
	if(model_name =="generate_cfg_from_synthese_file_Wscaled_a1a2a3asymovGamma"){
		generate_cfg_from_synthese_file_Wscaled_a1a2a3asymovGamma(input_params, file_out_modes,  file_out_noise, cfg.extra_params); // extra_params must points towards a .in file
		artificial_spectrum_a1a2a3asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel);
		passed=1;		
	}
	if(model_name =="generate_cfg_from_synthese_file_Wscaled_Alm"){
		generate_cfg_from_synthese_file_Wscaled_Alm(input_params, file_out_modes,  file_out_noise, cfg.extra_params); // extra_params must points towards a .in file
		artificial_spectrum_a1Alma3(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, 
									cfg.write_inmodel, cfg.do_modelfiles, cfg.limit_data_range, cfg.modefile_modelname);
		passed=1;		
	}
	if(model_name =="generate_cfg_from_synthese_file_Wscaled_aj"){
		generate_cfg_from_synthese_file_Wscaled_aj(input_params, file_out_modes,  file_out_noise, cfg.extra_params); // extra_params must points towards a .in file
		artificial_spectrum_aj(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, 
									cfg.write_inmodel, cfg.do_modelfiles, cfg.limit_data_range, cfg.modefile_modelname);
		passed=1;		
	}
	if(model_name == "asymptotic_mm_v1" || model_name == "asymptotic_mm_v2" || model_name == "asymptotic_mm_v3" ||
		model_name == "asymptotic_mm_freeDp_numaxspread_curvepmodes_v1" || model_name == "asymptotic_mm_freeDp_numaxspread_curvepmodes_v2" ||
		model_name == "asymptotic_mm_freeDp_numaxspread_curvepmodes_v3"){
		if(model_name =="asymptotic_mm_v1"){
			asymptotic_mm_v1(input_params, file_out_modes, file_out_noise,  file_cfg_mm, external_path, template_file);
			subpassed=1;
		}
		if(model_name =="asymptotic_mm_v2"){
			asymptotic_mm_v2(input_params, file_out_modes, file_out_noise,  file_cfg_mm, external_path, template_file);
			subpassed=1;
		}
		if(model_name =="asymptotic_mm_v3"){
			asymptotic_mm_v3(input_params, file_out_modes, file_out_noise,  file_cfg_mm, external_path, template_file);
			subpassed=1;
		}
		if(model_name =="asymptotic_mm_freeDp_numaxspread_curvepmodes_v1"){ 
			asymptotic_mm_freeDp_numaxspread_curvepmodes_v1(input_params, file_out_modes, file_out_noise,  file_cfg_mm, external_path, template_file);
			subpassed=1;
		}
		if(model_name =="asymptotic_mm_freeDp_numaxspread_curvepmodes_v2"){ 
			asymptotic_mm_freeDp_numaxspread_curvepmodes_v2(input_params, file_out_modes, file_out_noise,  file_cfg_mm, external_path, template_file);
			subpassed=1;
		}
		if(model_name =="asymptotic_mm_freeDp_numaxspread_curvepmodes_v3"){ 
			asymptotic_mm_freeDp_numaxspread_curvepmodes_v3(input_params, file_out_modes, file_out_noise,  file_cfg_mm, external_path, template_file);
			subpassed=1;
		}
		if(subpassed == 0){
			std::cout << "Warning: The function call_model did not generate any configuration file!" << std::endl;
			std::cout << "         Debug required in the section handling the models asymptotic_mm_vX" << std::endl;
			std::cout << "         The program will stop now" << std::endl;
			exit(EXIT_FAILURE);
		}
		
		//artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, dir_core, id_str, cfg.doplots, cfg.write_inmodel);
		artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel);
		std::cout << "   - [Warning DIRTY CODE IN LINE 281, iterative_artificial_spectrum.cpp] Hard coded path for saving python3 cfg in the spectra_info dir" << std::endl;
		std::cout << "                                                                         Stable final version should avoid this" << std::endl;
		str="cp " + file_cfg_mm + " " + dir_core + "Data/Spectra_info/" + strtrim(id_str) + ".global";
		const char *command = str.c_str(); 
		system(command);
		str="cp " + dir_core + "external/ARMM-solver/star_params.range " + dir_core + "Data/Spectra_info/" + strtrim(id_str) + ".range";
		const char *command2 = str.c_str(); 
		system(command2);
		str="cp " + dir_core + "external/ARMM-solver/star_params.rot " + dir_core + "Data/Spectra_info/" + strtrim(id_str) + ".rot";
		const char *command3 = str.c_str(); 
		system(command3);
		passed=1;
	}

	//exit(EXIT_SUCCESS);
	return passed;
}


// This is the routine that calls models that are valid on the grid approach
bool call_model_grid(std::string model_name, VectorXd input_params, Model_data input_model, std::string file_out_modes, std::string file_out_noise, 
		std::string dir_core, std::string id_str, Config_Data cfg){

	bool passed=0;
	std::string file_ref_star, modelID_str;
	std::ostringstream strs;

	if(model_name == "generate_cfg_asymptotic_act_asym_Hgauss"){
		generate_cfg_asymptotic_act_asym_Hgauss(input_params, file_out_modes, file_out_noise);
		//id_str=identifier2chain(identifier);
		modelID_str="NONE";
		artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel);
		//artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, dir_core, id_str, modelID_str, cfg.doplots, cfg.write_inmodel);
		passed=1;
	}
	if(model_name == "generate_cfg_from_refstar_HWscaled"){
		file_ref_star=dir_core + "Configurations/ref_spectra.params"; // The spectra used for the reference star
		generate_cfg_from_refstar_HWscaled(input_params, input_model, file_ref_star, file_out_modes, file_out_noise);
		strs << input_model.params[0];
		modelID_str=format_freqname(strs.str());
		artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel);
		//artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, dir_core, id_str, modelID_str, cfg.doplots, cfg.write_inmodel);
		passed=1;
	}
	if(model_name == "generate_cfg_from_refstar_HWscaled_GRANscaled"){
		file_ref_star=dir_core + "Configurations/ref_spectra.params"; // The spectra used for the reference star
		generate_cfg_from_refstar_HWscaled_GRANscaled(input_params, input_model, file_ref_star, file_out_modes, file_out_noise);
		strs << input_model.params[0];
		modelID_str=format_freqname(strs.str());
		artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel);
		//artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, dir_core, id_str, modelID_str, cfg.doplots, cfg.write_inmodel);
		passed=1;
	}

	if(model_name =="generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma"){
		generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma(input_params, file_out_modes,  file_out_noise, cfg.extra_params); // extra_params must points towards a .in file
		artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel);
		passed=1;		
	}
	if(model_name =="generate_cfg_from_synthese_file_Wscaled_a1a2a3asymovGamma"){
		generate_cfg_from_synthese_file_Wscaled_a1a2a3asymovGamma(input_params, file_out_modes,  file_out_noise, cfg.extra_params); // extra_params must points towards a .in file
		artificial_spectrum_a1a2a3asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel);
		passed=1;		
	}
	if(model_name =="generate_cfg_from_synthese_file_Wscaled_Alm"){
		generate_cfg_from_synthese_file_Wscaled_Alm(input_params, file_out_modes,  file_out_noise, cfg.extra_params); // extra_params must points towards a .in file
		artificial_spectrum_a1Alma3(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, 
									cfg.write_inmodel, cfg.do_modelfiles, cfg.limit_data_range, cfg.modefile_modelname);
		passed=1;		
	}
	if(model_name =="generate_cfg_from_synthese_file_Wscaled_aj"){
		generate_cfg_from_synthese_file_Wscaled_aj(input_params, file_out_modes,  file_out_noise, cfg.extra_params); // extra_params must points towards a .in file
		artificial_spectrum_aj(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, 
									cfg.write_inmodel, cfg.do_modelfiles, cfg.limit_data_range, cfg.modefile_modelname);
		passed=1;		
	}
	return passed;
}


// -----------------------------------------------------------------------
// ------------------------- SECONDARY METHODS ---------------------------
// -----------------------------------------------------------------------
// Retrieve a list of files within path. 
// Note that path cannot contain any system-based filter such as '*.*' etc... 
// Those have to be put in the extension string
std::vector<std::string> list_dir(const std::string path, const std::string extension){

	const std::string tmp_file="dir.list";
	const std::string cmd="ls -p " + path + "| grep -v / >> " + tmp_file; // Get a list of files and store them into the temporary dir.list file
	const std::string cmd_erase="rm " + tmp_file; // Once finished with the temporary file, we erase it

	std::vector<std::string> files, line_split;
	std::string line;
	std::ifstream file_in;

	//std::cout << "command: " << cmd << std::endl; 
	system(cmd.c_str());
	
	file_in.open(tmp_file.c_str());
   	if (file_in.is_open()) {
   		while(!file_in.eof()){
   			std::getline(file_in, line);
			line_split=strsplit(line, ".");
			if (line_split[line_split.size()-1] == extension){
				files.push_back(strtrim(line)); // remove any white space at the begining/end of the string
   			} else{
   			}
   		}
   		system(cmd_erase.c_str());
   	} else{
   		std::cout << "I/O Error for the temporary temporary file " + tmp_file + "!"  << std::endl;
   		std::cout << "Cannot proceed. Check that you have the proper I/O rights on the root program directory..." << std::endl;
   		exit(EXIT_FAILURE);
   	}
   	return files;
}

/*
long read_id_allcombi(std::string file_combi){

	std::string lastline;
	std::vector<std::string> vals_last;

	lastline=read_lastline_ascii(file_combi);
	vals_last=strsplit(strtrim(lastline), " ");

	std::cout << "lastline=" << lastline << std::endl;
	for(int i=0; i<vals_last.size(); i++){
		std::cout << "vals_last[" << i << "]=" << vals_last[i] << std::endl;	
	}
	return str_to_long(vals_last[0]);
}
*/

/*
std::string read_lastline_ascii(std::string filename){
//
// Code that jumps to last line and read it
//
//

	std::string line;
	int i;
	std::string lastline;

	std::ifstream myfile;
        myfile.open(filename.c_str());

	while (getline (myfile,line)){
  		lastline=line;
	} 
        myfile.close();
   
	return lastline;
  
}
*/

/*
std::string write_allcombi(MatrixXd allcombi, VectorXd cte_params, Config_Data cfg, std::string fileout, bool erase_old_file, long iter, long id0, 
		    std::vector<std::string> cte_names, std::vector<std::string> var_names, std::vector<std::string> param_names){
	
	int Nchars, precision;
	std::string id_str;
	VectorXd input_params, var_params;
	std::ofstream outfile;

	Nchars = 14;
	precision = 5;

	if(erase_old_file == 1 && iter == 0) {
		outfile.open(fileout.c_str()); // write a new file
	} else{
		outfile.open(fileout.c_str(), std::ios::app); // append
	}
	if(outfile.is_open()){
		//std::cout << "File opened" << std::endl;
		if(erase_old_file == 1 && iter == 0) { // Write Header only if we do not erase the old file AND this is the first execution of the function
			outfile << "model_name= " << cfg.model_name << std::endl;
			outfile << " --------------------------" << std::endl;
			outfile << "  List of all combinations " << std::endl;
			outfile << " --------------------------" << std::endl;
			outfile << "#" << std::setw(5) << "id  ";
			for(int s=0; s<param_names.size(); s++){
				outfile << std::setw(14) << param_names[s];
			}
			outfile << std::endl;
		} 

		for(int i=0; i<allcombi.rows(); i++){
			id_str=identifier2chain(i + id0); // The identifier corresponds to the index of the current process + the initial id0
			outfile << std::setw(7) << id_str;
			//std::cout << "id_str=" << id_str << std::endl;
			
			var_params=allcombi.row(i).transpose();
			input_params=order_input_params(cte_params, var_params, cte_names, var_names, param_names);
			//std::cout << "input_params=" << input_params << std::endl;
			for(int j=0; j<input_params.size(); j++){			
				outfile << std::setw(Nchars) << std::setprecision(precision) << input_params(j);	
			}
			outfile << std::endl;
		}
	outfile.close();
	}  
	else {
		std::cout << " Unable to open file " << fileout << std::endl;	
		std::cout << " Check that the full path exists" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

	return id_str;
}
*/

/*
std::string identifier2chain(long identifier){

	std::string out;

	if  (identifier < 10) { out="000000" + strtrim(lng_to_str(identifier)); }
	if ((identifier >= 10) && (identifier < 100)) { out="00000" + strtrim(lng_to_str(identifier));}
	if ((identifier >= 100) && (identifier < 1000)) { out="0000" + strtrim(lng_to_str(identifier));}
	if ((identifier >= 1000) && (identifier < 10000)) { out="000" + strtrim(lng_to_str(identifier));}
	if ((identifier >= 10000) && (identifier  < 100000)) { out= "00" + strtrim(lng_to_str(identifier)); }
	if ((identifier >= 100000) && (identifier< 1000000)) { out="0" + strtrim(lng_to_str(identifier));}
	if (identifier  >  10000000) { 
		std::cout << "Warning: This cannot handle greater number than 99999" << std::endl;
		std::cout << "Pursuing will lead to an improper naming of the models" << std::endl;
		std::cout << "Please update identifier2chain in order to handle greater numbers" << std::endl;
		std::cout << "The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
	}
	return out;
}

std::string identifier2chain(std::string identifier){

	std::string out;

	if  (str_to_dbl(identifier) < 10) { out="000000" + strtrim(identifier); }
	if ((str_to_dbl(identifier) >= 10) && (str_to_dbl(identifier)    < 100)) { out="00000" + strtrim(identifier);}
	if ((str_to_dbl(identifier) >= 100) && (str_to_dbl(identifier)   < 1000)) { out="0000" + strtrim(identifier);}
	if ((str_to_dbl(identifier) >= 1000) && (str_to_dbl(identifier)  < 10000)) { out="000" + strtrim(identifier);}
	if ((str_to_dbl(identifier)  >= 10000) && (str_to_dbl(identifier)  < 100000)) { out= "00" + strtrim(identifier); }
	if ((str_to_dbl(identifier) >= 100000) && (str_to_dbl(identifier)< 1000000)) { out="0" + strtrim(identifier);}
	if (str_to_dbl(identifier)  >  10000000) { 
		std::cout << "Warning: This cannot handle greater number than 9999999" << std::endl;
		std::cout << "Pursuing will lead to an improper naming of the models" << std::endl;
		std::cout << "Please update identifier2chain in order to handle greater numbers" << std::endl;
		std::cout << "The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
	}
	return out;
}
*/

void showversion()
{
    std::cout << APP_NAME " " APP_VERSION "\n built on " __DATE__ << std::endl;

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

int options(int argc, char* argv[]){

	std::string arg1, arg2;
	int val;
	
	val=-2;
	arg1="";
	
	//std::cout << argc << std::endl;
	if(argc == 1){
		val=-1; // Code for do nothing here (no options)
		//std::cout << "No option passed, continuing..." << std::endl;
		std::cout << " ------------" << std::endl;
		showversion();
		std::cout << " ------------" << std::endl << std::endl;
		
		
	}
	if(argc > 1){
		arg1=argv[1];
	}
	if(argc > 2){
		std::cout << "Too many arguments. Allowed number: 1. " << std::endl;
		std::cout << "Extra arguments will be ignored" << std::endl;	
		//arg2=argv[2];
	}
	if(argc == 2){
		if(arg1 == "version"){
			 val=0;
		} else{
			std::cout << "Unknown argument. Allowed arguments: "<< std::endl;
			std::cout << "    version : Returns the version of the code and exit" << std::endl;
		}
	}

	if (val == -2){ // Error code
		exit(EXIT_FAILURE);
	} 
	if (val == 0){ // Version code
		showversion(); 
		exit(EXIT_SUCCESS);
	}
	if (val >= -1 ){
		return val; // Execution code val
	} else{
		return -2; // Default value is to return an error code
	}
}


int main(int argc, char* argv[]){

	boost::filesystem::path full_path( boost::filesystem::current_path() );

// -------- Options Handling ------
	int msg_code;
	msg_code=options(argc, argv);
	if(msg_code == -2){
		std::cout << "Error detected in options. Cannot proceed. Debug required." << std::endl;
		std::cout << "The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	} 
	// -------------------------------

	iterative_artificial_spectrum(full_path.string() + "/");

}
