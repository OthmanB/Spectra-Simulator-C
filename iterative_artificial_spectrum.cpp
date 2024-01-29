/*
 * iterative_artificial_spectrum.cpp
 *
 * Header file that contains all kind of methods
 * used to generate iteratively artificial spectra
 * 
 *  Created on: 05 May 2016
 *      Author: obenomar
 */

/**
 * @file iterative_artificial_spectrum.cpp
 * @brief Header file for the artificial spectra generation
 *
 * Header file for the functions used to generate iteratively artificial spectra
 * 
 *
 * @date 05 May 2016
 * @author obenomar
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/program_options.hpp>
#include "artificial_spectrum.h"
#include "models_database.h"
#include "models_database_grid.h"
#include "version.h"
#include "combi.h"
#include "stellar_models.h"
#include "io_star_params.h"

/**
 * @brief Creates directories and subdirectories.
 *
 * This function creates the specified directories and subdirectories. It can either create the directories if they don't exist or check if they exist and return an error if they don't.
 *
 * @param output_dir The output directory path.
 * @param force_mkdir Flag indicating whether to create the directories if they don't exist (true) or check if they exist (false).
 * @return Returns true if the directories were created successfully or if they already exist, and false otherwise.
 */
bool createDirectories(const std::string& output_dir, bool force_mkdir);

/**
 * @brief Display the version information of the application.
 *
 * This function displays the version information of the application, including the application name, version number, and build date. It also displays information about the compiler used to build the application and the features supported by the platform.
 */
void showversion();

/**
 * @brief Run iteratively the program artificial_spectrum to generate a series of models.
 *
 * This function runs the program artificial_spectrum iteratively in order to generate a series of models.
 * It takes a directory path and a configuration file path as input parameters.
 * The directory path specifies the location where the generated models will be saved.
 * The configuration file path specifies the main configuration file that will be used for generating the models.
 *
 * @param dir_core The directory path where the generated models will be saved.
 * @param cfg_file The path of the main configuration file.
 * @param cfg_noise_file The full path to a noise configuration file
 * @param data_path The root path on which the data are going to be written
 * @date 29 January 2024 (Update)
 */
void iterative_artificial_spectrum(const std::string dir_core, const std::string cfg_file, const std::string cfg_noise_file, const std::string data_path);

/**
 * @brief Function to check validity of the configuration parameters
 * 	
 * 	The function check the validity of parameters in the cfg structure. If it is not as expected, it forces exit of the program
 * 
 * @param cfg Configuration structure of type Config_Data
 * @param N_model Number of model parameters
 * @date 6 December 2023
 */
void check_params(Config_Data cfg, const int N_model);

/**
 * @brief Aggregates the main configuration data with the noise configuration data.
 * 
 * @param cfg_main The main configuration data.
 * @param cfg_noise The noise configuration data.
 * @return The aggregated configuration data.
 * 
 * @details This function aggregates the main configuration data with the noise configuration data. It copies all the content of the cfg_main structure into the cfg structure using the update_cfg function with priority 1. Then, it updates the cfg structure by adding the cfg_noise at the end.
 * 
 * If the forest_type in cfg_main is "random", it iterates through the noise configuration data and adds the corresponding variables to the cfg structure. It also sets the step value based on the distribution type.
 * 
 * If the forest_type in cfg_main is "grid", it iterates through the noise configuration data and adds the corresponding variables to the cfg structure. It sets the step value to the kerror_grid value for the "Uniform" distribution type. If the distribution type is not "Uniform", it throws an error and exits.
 * 
 * If the forest_type in cfg_main is neither "random" nor "grid", it throws an error and exits.
 * 
 * The function returns the aggregated cfg structure.
 * @date 6 December 2023
 */
Config_Data agregate_maincfg_noisecfg(Config_Data cfg_main, Config_Noise cfg_noise);

/**
 * @brief Updates the target configuration data with the source configuration data based on priority.
 *
 * @param cfg_target The target configuration data.
 * @param cfg_source The source configuration data.
 * @param priority4common The priority for common variables.
 *                        -1: Check and enforce the equality of common variables between cfg_target and cfg_source.
 *                         0: No priority. No action is taken.
 *                         1: Priority on cfg_source. Common variables are updated from cfg_source.
 *
 * @return The updated configuration data.
 *
 * @details This function updates the target configuration data with the source configuration data based on priority.
 *          If priority4common is -1, it checks and enforces the equality of common variables between cfg_target and cfg_source.
 *          If any common variable is not equal, it exits with an error.
 * 			If priority4common is 0, it uses the cfg_target values as they are and increase the length of the variables 
 * 			"labels", "distrib", "val_min", "val_max", "step" using the information in cfg_source
 *          If priority4common is 1, it overwrites the common variables in cfg_target with the values from cfg_source.
 *          It also updates the "labels", "distrib", "val_min", "val_max", "step"  variables in cfg_target with the 
 * 			corresponding values from cfg_source.
 *          The function returns the updated cfg_target.
 * @date 6 December 2023
 */
Config_Data update_cfg(Config_Data cfg_target, const Config_Data cfg_source, const int priority4common);

/**
* @brief Generate random models based on a configuration file.
*
* This function generates random models based on a given configuration file. It takes several input parameters including the configuration file, directory path, and file paths for output files. The function generates a series of random values for the parameters specified in the configuration file and uses them to generate the models. It also selects a template file containing height and width profiles for the models.
*
* @param cfg The configuration data.
* @param param_names The names of the parameters.
* @param dir_core The directory path where the generated models will be saved.
* @param file_out_modes The file path for the output modes.
* @param file_out_noise The file path for the output noise.
* @param file_out_combi The file path for the output combinations.
* @param N_model The number of models to generate.
* @param file_cfg_mm The file path for the main configuration file.
* @param external_path The external path.
* @param templates_dir The directory path for the template files. 
*/
void generate_random(Config_Data cfg, std::vector<std::string> param_names, std::string dir_core, std::string file_out_modes, 
		std::string file_out_noise, std::string file_out_combi, int N_model, std::string file_cfg_mm, std::string external_path,  std::string templates_dir, std::string data_path);

/**
 * @brief Generate a grid of combinations based on a configuration file.
 *
 * This function generates a grid of combinations for the variables of the model, based on a given configuration file. It takes several input parameters including the configuration data, whether to use models, the models data, parameter names, directory paths, and file paths for output files. The function generates all possible combinations of values for the parameters specified in the configuration file and uses them to generate the grid. It also calls the model grid function to generate the models for each combination.
 *
 * @param cfg The configuration data.
 * @param usemodels Whether to use models.
 * @param models The models data.
 * @param param_names The names of the parameters.
 * @param dir_core The directory path where the generated models will be saved.
 * @param dir_freqs The directory path for the frequency files.
 * @param file_out_modes The file path for the output modes.
 * @param file_out_noise The file path for the output noise.
 * @param file_out_combi The file path for the output combinations.
 * @param N_model The number of models.
 */
void generate_grid(Config_Data cfg, bool usemodels, Data_Nd models, std::vector<std::string> param_names, std::string dir_core, std::string dir_freqs, std::string file_out_modes, 
		std::string file_out_noise, std::string file_out_combi, int N_model, std::string data_path);


/**
 * @brief Call a model based on the given model name and input parameters.
 *
 * This function calls a specific model based on the given model name and input parameters. Compared to call_model_grid(), only the models compatible with variables generated randomly are listed. It also generates output files for the modes and noise. The function checks the model name and performs the corresponding actions for each model.
 *
 * @param model_name The name of the model to call.
 * @param input_params The input parameters for the model.
 * @param file_out_modes The file path for the output modes.
 * @param file_out_noise The file path for the output noise.
 * @param file_cfg_mm The file path for the main configuration file.
 * @param dir_core The directory path where the generated models will be saved.
 * @param identifier The ID of the star for the current combination.
 * @param cfg The configuration data.
 * @param external_path The external path.
 * @param template_file The template file path.
 * @return True if the model was called successfully, false otherwise.
 */
bool call_model_random(std::string model_name, VectorXd input_params, std::string file_out_modes, std::string file_out_noise, 
		 std::string file_cfg_mm, std::string dir_core, std::string identifier, Config_Data cfg, std::string external_path, std::string template_file, std::string data_path);

/**
 * @brief Call a model based on the given model name and input parameters.
 *
 * This function calls a specific model based on the given model name and input parameters. It also generates output files for the modes and noise. The function checks the model name and performs the corresponding actions for each model.
 *
 * @param model_name The name of the model to call.
 * @param input_params The input parameters for the model.
 * @param input_model The input model data.
 * @param file_out_modes The file path for the output modes.
 * @param file_out_noise The file path for the output noise.
 * @param dir_core The directory path where the generated models will be saved.
 * @param id_str The ID string for the current combination.
 * @param cfg The configuration data.
 * @return True if the model was called successfully, false otherwise.
 */
bool call_model_grid(std::string model_name, VectorXd input_params, Model_data input_model, std::string file_out_modes, std::string file_out_noise, 
		std::string dir_core, std::string id_str, Config_Data cfg, std::string data_path);

/**
 * @brief Retrieve a list of files within a specified path.
 *
 * This function retrieves a list of files within the specified path. The path should not contain any system-based filters such as '*.*'. These filters should be included in the extension string parameter.
 *
 * @param path The path from which to retrieve the list of files.
 * @param filter Term allowing to return syntax-specific name, eg. kplr*
 * @return A vector of strings containing the names of the files in the specified path.
 */
std::vector<std::string> list_dir(const std::string path, const std::string filter);


void iterative_artificial_spectrum(const std::string dir_core, const std::string cfg_file, const std::string cfg_noise_file, const std::string data_out_path){
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
	std::string model_file, file_out_modes, file_out_noise, file_cfg_mm, file_out_mm, file_out_mm2, file_out_combi;
	std::string data_path, external_path, templates_dir;
	Config_Data cfg;
	Config_Noise cfg_noise;
	Data_Nd models;
	external_path=dir_core + "external/"; 
	templates_dir=dir_core + "Configurations/templates/";
	model_file=dir_core + "external/MESA_grid/models.params"; // This is for models of MS stars from a grid
	file_out_modes=dir_core + "Configurations/tmp/modes_tmp.cfg";
	file_out_noise=dir_core + "Configurations/tmp/noise_tmp.cfg";
	file_cfg_mm=dir_core + "Configurations/MixedModes_models/star_params.theoretical"; // This is for models of Evolved stars. Used only when the user has already a set of frequencies and want to use them to compute a model
	
	if(data_out_path == ""){
		data_path=dir_core + "/Data/";
	} else{
		data_path=data_out_path;
	}
	file_out_combi=data_path + "Combinations.txt";

	std::string dir_freqs=dir_core + "external/MESA_grid/frequencies/";

	std::cout << "1. Read the configuration file..." << std::endl;

	cfg=read_main_cfg(cfg_file);
	cfg_noise=readNoiseConfigFile(cfg_noise_file); // Used to handle models with a separate noise, essentially the Kallinger+2014 noise 

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
			std::cout << "    Invalid number of parameters for model_name= 'generate_cfg_from_synthese_file_Wscaled_Alm'" << std::endl;
			std::cout << "    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size() << std::endl;
			std::cout << "    Check your main configuration file" << std::endl;
			std::cout << "    The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}
		passed=1;
	}
	if(cfg.model_name == "generate_cfg_from_synthese_file_Wscaled_aj"){
		Nmodel=17;
		param_names.push_back("Dnu");
		param_names.push_back("epsilon"); 
		param_names.push_back("delta0l_percent"); 
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
		param_names.push_back("H_spread");
		param_names.push_back("nu_spread");
		param_names.push_back("Gamma_spread");
		param_names.push_back("do_flat_noise");
		if(param_names.size() != Nmodel){
			std::cout << "    Invalid number of parameters for model_name= 'generate_cfg_from_synthese_file_Wscaled_aj'" << std::endl;
			std::cout << "    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size() << std::endl;
			std::cout << "    Check your main configuration file" << std::endl;
			std::cout << "    The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}
		passed=1;
	}
	if(cfg.model_name == "generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled"){
		Nmodel=24;
		param_names.push_back("Dnu");
		param_names.push_back("epsilon"); 
		param_names.push_back("delta0l_percent"); 
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
		param_names.push_back("A_Pgran"); // Coefficients for scaling the noise with numax. See Karoff 2010 or Kallinger 2014.
		param_names.push_back("B_Pgran");	
		param_names.push_back("C_Pgran");	
		param_names.push_back("A_taugran");	
		param_names.push_back("B_taugran");	
		param_names.push_back("C_taugran");	
		param_names.push_back("P");	
		param_names.push_back("N0");
		param_names.push_back("numax_spread"); // Add a spread to numax to avoid to have a noise that stricly follow the Karoff et al. 2010 relation
		param_names.push_back("H_spread");
		param_names.push_back("nu_spread");
		if(param_names.size() != Nmodel){
			std::cout << "    Invalid number of parameters for model_name= 'generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled'" << std::endl;
			std::cout << "    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size() << std::endl;
			std::cout << "    Check your main configuration file" << std::endl;
			std::cout << "    The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}
		passed=1;
	}
	if(cfg.model_name == "generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled_Kallinger2014"){
		// Warning: This model uses the file noise_Kallinger2014.cfg to set the noise parameters
		const int Nmodel_modes=16;
		const int Nmodel_noise=14; // 6 Dec 2023: Many of these parameters are in fact generated according to a Gaussian 
		// Agregate the old configuration with the noise configuration
		cfg=agregate_maincfg_noisecfg(cfg, cfg_noise);	
		param_names.push_back("Dnu");
		param_names.push_back("epsilon"); 
		param_names.push_back("delta0l_percent"); 
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
		param_names.push_back("numax_spread"); // Add a spread to numax to avoid to have a noise that stricly follow the Karoff et al. 2010 relation
		param_names.push_back("H_spread");
		param_names.push_back("nu_spread");
		//k_Agran         s_Agran         k_taugran       s_taugran       c0              ka              ks              k1              s1              c1              k2              s2              c2              N0
		param_names.push_back("k_Agran");
		param_names.push_back("s_Agran");
		param_names.push_back("k_taugran");
		param_names.push_back("s_taugran");
		param_names.push_back("c0");
		param_names.push_back("ka");
		param_names.push_back("ks");
		param_names.push_back("k1");
		param_names.push_back("s1");
		param_names.push_back("c1");
		param_names.push_back("k2");
		param_names.push_back("s2");
		param_names.push_back("c2");
		param_names.push_back("N0");
		if(param_names.size() != Nmodel_modes+Nmodel_noise){
			std::cout << "    Invalid number of parameters for model_name= 'generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled_Kallinger2014'" << std::endl;
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
		param_names.push_back("Hfactor");
		param_names.push_back("Wfactor");			
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
		param_names.push_back("Hfactor");
		param_names.push_back("Wfactor");			
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
		param_names.push_back("Hfactor");
		param_names.push_back("Wfactor");			
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
		param_names.push_back("Hfactor");
		param_names.push_back("Wfactor");			
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
		param_names.push_back("a2_l1_core"); 
		param_names.push_back("a2_l1_env"); 
		param_names.push_back("a2_l2_env"); 
		param_names.push_back("a2_l3_env"); 
		param_names.push_back("a3_l2_env"); 
		param_names.push_back("a3_l3_env"); 
		param_names.push_back("a4_l2_env"); 
		param_names.push_back("a4_l3_env"); 
		param_names.push_back("a5_l3_env"); 
		param_names.push_back("a6_l3_env"); 
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
		param_names.push_back("Hfactor");
		param_names.push_back("Wfactor");			
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
		param_names.push_back("nurot_env"); 
		param_names.push_back("nurot_core"); 
		param_names.push_back("a2_l1_core"); 
		param_names.push_back("a2_l1_env"); 
		param_names.push_back("a2_l2_env"); 
		param_names.push_back("a2_l3_env"); 
		param_names.push_back("a3_l2_env"); 
		param_names.push_back("a3_l3_env"); 
		param_names.push_back("a4_l2_env"); 
		param_names.push_back("a4_l3_env"); 
		param_names.push_back("a5_l3_env"); 
		param_names.push_back("a6_l3_env");
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
		param_names.push_back("Hfactor");
		param_names.push_back("Wfactor");			
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

	if(cfg.model_name == "asymptotic_mm_freeDp_numaxspread_curvepmodes_v3_GRANscaled_Kallinger2014"){
		// Warning: This model uses the file noise_Kallinger2014.cfg to set the noise parameters
		const int Nmodel_modes=29;
		const int Nmodel_noise=14; // 6 Dec 2023: Many of these parameters are in fact generated according to a Gaussian 
		// Agregate the old configuration with the noise configuration
		cfg=agregate_maincfg_noisecfg(cfg, cfg_noise);	
		param_names.push_back("nurot_env"); 
		param_names.push_back("nurot_core"); 
		param_names.push_back("a2_l1_core"); 
		param_names.push_back("a2_l1_env"); 
		param_names.push_back("a2_l2_env"); 
		param_names.push_back("a2_l3_env"); 
		param_names.push_back("a3_l2_env"); 
		param_names.push_back("a3_l3_env"); 
		param_names.push_back("a4_l2_env"); 
		param_names.push_back("a4_l3_env"); 
		param_names.push_back("a5_l3_env"); 
		param_names.push_back("a6_l3_env");
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
		//k_Agran         s_Agran         k_taugran       s_taugran       c0              ka              ks              k1              s1              c1              k2              s2              c2              N0
		param_names.push_back("k_Agran");
		param_names.push_back("s_Agran");
		param_names.push_back("k_taugran");
		param_names.push_back("s_taugran");
		param_names.push_back("c0");
		param_names.push_back("ka");
		param_names.push_back("ks");
		param_names.push_back("k1");
		param_names.push_back("s1");
		param_names.push_back("c1");
		param_names.push_back("k2");
		param_names.push_back("s2");
		param_names.push_back("c2");
		param_names.push_back("N0");
		param_names.push_back("Hfactor");
		param_names.push_back("Wfactor");	
		if(param_names.size() != Nmodel_modes+Nmodel_noise){
			std::cout << "    Invalid number of parameters for model_name= 'generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled_Kallinger2014'" << std::endl;
			std::cout << "    Expecting " << Nmodel_modes + Nmodel_noise << " parameters, but found " << cfg.val_min.size() << std::endl;
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
		generate_random(cfg, param_names, dir_core, file_out_modes, file_out_noise, file_out_combi, Nmodel, file_cfg_mm, external_path, templates_dir, data_path);

	}
	if(cfg.forest_type == "grid"){
		std::cout << "   Values are generated over a grid using all possible combinations according to inputs in the main configuration file" << std::endl;
		generate_grid(cfg, usemodels, models, param_names, dir_core, dir_freqs, file_out_modes, file_out_noise, file_out_combi, Nmodel,data_path);
	}
	if(cfg.forest_type != "random" && cfg.forest_type != "grid"){ 
		std::cout << " Problem in the main configuration file. It is expected that the forest type parameters is either random OR grid" << std::endl;
		std::cout << " Check your main configuration file" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

}

void check_params(Config_Data cfg, const int N_model){
	bool neg=0;
	int ii=0;
	while(neg == 0 && (ii < N_model)){
		if( (cfg.distrib[ii] == "Uniform") && (cfg.val_max[ii] - cfg.val_min[ii]) < 0){ 
			neg=1;
			std::cout << "List of the parameters and values :" << std::endl;
			for (int i=0;i<cfg.labels.size();i++){
				std::cout << "[" << i << "] " << std::setw(20) << cfg.labels[i] << std::setw(12)<< cfg.val_min[i] << std::setw(12) << cfg.val_max[i] << std::endl;
			}
			std::cout << " ------- " << std::endl;
			std::cout << "       ii = " << ii << std::endl;
			std::cout << "       name:    " << cfg.labels[ii] << std::endl;
			std::cout << "       val_max =" << cfg.val_max[ii] << std::endl;
			std::cout << "       val_min =" << cfg.val_min[ii] << std::endl;
		}
		ii=ii+1;
	}
	if(neg == 1){
		std::cout << "     Warning: val_max < val_min for some of the parameters while a uniform distribution is requested!" << std::endl;
		std::cout << "     This is not permitted" << std::endl;
		std::cout << "     Check your main configuration file" << std::endl;
		std::cout << "     The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}
}

Config_Data update_cfg(Config_Data cfg_target, const Config_Data cfg_source, const int priority4common){
	bool pass=false;
	if (priority4common <-1 || priority4common >1){
		std::cerr << "Error: Invalid pririty4common value. Set it to -1, 0 or 1." << std::endl;
		exit(EXIT_FAILURE);		
	}
	// Case where there is no priority... then we check that common variables in cfg_target and cfg_source are the
	// same. If it is not the case, we force the exit with error
	if (priority4common == -1) {
		// First check the single elements
		if (cfg_target.erase_old_files != cfg_source.erase_old_files ||
			cfg_target.doplots != cfg_source.doplots ||
			cfg_target.write_inmodel != cfg_source.write_inmodel ||
			cfg_target.limit_data_range != cfg_source.limit_data_range ||
			cfg_target.do_modelfiles != cfg_source.do_modelfiles ||
			cfg_target.modefile_modelname != cfg_source.modefile_modelname ||
			cfg_target.model_name != cfg_source.model_name ||
			cfg_target.extra_params != cfg_source.extra_params ||
			cfg_target.Tobs != cfg_source.Tobs ||
			cfg_target.Cadence != cfg_source.Cadence ||
			cfg_target.Nspectra != cfg_source.Nspectra ||
			cfg_target.Nrealisation != cfg_source.Nrealisation ||
			cfg_target.forest_type != cfg_source.forest_type)
			{
				std::cerr << "Error: Common variables in cfg_target and cfg_source are not the same." << std::endl;
				exit(EXIT_FAILURE);
			}
	    // Check each element of the template_files vector
    	for (size_t i = 0; i < cfg_target.template_files.size(); i++) {
			if (cfg_target.template_files[i] != cfg_source.template_files[i]) {
				std::cerr << "Error: template_files in cfg_target and cfg_source are not the same." << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		// Check each element of the forest_params vector
    	for (size_t i = 0; i < cfg_target.forest_params.size(); i++) {
			if (cfg_target.forest_params[i] != cfg_source.forest_params[i]) {
				std::cerr << "Error: forest_params in cfg_target and cfg_source are not the same." << std::endl;
				exit(EXIT_FAILURE);
			}
		}
	}
	// Case where the priority is on cfg_target. Then, common variables are just those of cfg_target
	// --> nothing to do
	// Case where the priority is on cfg_source. Then common variables are those from cfg_source
	if (priority4common == 1){
		cfg_target.erase_old_files=cfg_source.erase_old_files;
		cfg_target.doplots=cfg_source.doplots;
		cfg_target.write_inmodel=cfg_source.write_inmodel;
		cfg_target.limit_data_range=cfg_source.limit_data_range;
		cfg_target.do_modelfiles=cfg_source.do_modelfiles;
		cfg_target.modefile_modelname=cfg_source.modefile_modelname;
		cfg_target.model_name=cfg_source.model_name;
		cfg_target.forest_type=cfg_source.forest_type;
		cfg_target.extra_params=cfg_source.extra_params;
		cfg_target.template_files=cfg_source.template_files;
		cfg_target.Tobs=cfg_source.Tobs;
		cfg_target.Cadence=cfg_source.Cadence;
		cfg_target.Nspectra=cfg_source.Nspectra;
		cfg_target.Nrealisation=cfg_source.Nrealisation;
		cfg_target.forest_params=cfg_source.forest_params;
		pass=true;
	}
	// Update the arrays of variables
	for (size_t i = 0; i < cfg_source.labels.size(); i++) {
		cfg_target.labels.push_back(cfg_source.labels[i]);
		cfg_target.distrib.push_back(cfg_source.distrib[i]);
		cfg_target.val_min.push_back(cfg_source.val_min[i]);
		cfg_target.val_max.push_back(cfg_source.val_max[i]);
		cfg_target.step.push_back(cfg_source.step[i]);
	}
	return cfg_target;
}


Config_Data agregate_maincfg_noisecfg(Config_Data cfg_main, Config_Noise cfg_noise){
	bool pass=false;
	Config_Data cfg;
	// Copy all the content of the cfg_main structure into the cfg structure
	cfg=update_cfg(cfg, cfg_main, 1); 
	// Update the cfg structure by adding the cfg_noise at the end
	if (cfg.forest_type == "random"){
		for (int i=0; i<cfg_noise.name_random.size();i++){
			cfg.labels.push_back(cfg_noise.name_random[i]);
			cfg.distrib.push_back(cfg_noise.distrib_random[i]);
			cfg.val_min.push_back(cfg_noise.x1_random[i]);
			if (cfg_noise.distrib_random[i] == "Gaussian"){
				cfg.val_max.push_back(cfg_noise.kerror_random[i]*cfg_noise.x2_random[i]);
				cfg.step.push_back(1);
			} else if (cfg_noise.distrib_random[i] == "Uniform"){
				cfg.val_max.push_back(cfg_noise.x2_random[i]);
				cfg.step.push_back(1);
			} else if (cfg_noise.distrib_random[i] == "Fix"){
				cfg.val_max.push_back(cfg_noise.x2_random[i]);
				cfg.step.push_back(0);
			} else{
				std::cerr << "Error in iterative_artificial_spectrum::agregate_maincfg_noisecfg():" << std::endl;
				std::cerr << "    In random mode and when agregation of a main.cfg with a noise.cfg is requested, only Gaussian, Fix and Uniform is a valid entry for the distribution parameter" << std::endl;
				std::cerr << "    Check the content of the used noise configuration file." << std::endl;
				std::cerr << "    The program will exit now" << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		pass=true;
	}
	if (cfg.forest_type == "grid"){
		for (int i=0; i<cfg_noise.name_random.size();i++){
			cfg.labels.push_back(cfg_noise.name_grid[i]);
			cfg.distrib.push_back(cfg_noise.distrib_grid[i]);
			cfg.val_min.push_back(cfg_noise.x1_grid[i]);
			cfg.val_max.push_back(cfg_noise.x2_grid[i]);
			if (cfg_noise.distrib_grid[i] == "Fix"){
				cfg.step.push_back(0);
			} else if (cfg_noise.distrib_grid[i] != "Uniform"){
					std::cerr << "Error in iterative_artificial_spectrum::agregate_maincfg_noisecfg():" << std::endl;
					std::cerr << "    In grid mode and when agregation of a main.cfg with a noise.cfg is requested, only Fix and Uniform is a valid entry for the distribution parameter" << std::endl;
					std::cerr << "    Check the content of the used noise configuration file." << std::endl;
					std::cerr << "    The program will exit now" << std::endl;
					exit(EXIT_FAILURE);
			} else{ // In the uniform case, the kerror_grid is used to set the step
				cfg.step.push_back(cfg_noise.kerror_grid[i]);
				cfg.step.push_back(1);
			}
		}
		pass=true;
	}
	if (pass == false){
		std::cerr << "Error in iterative_artificial_spectrum::agregate_maincfg_noisecfg():" << std::endl;
		std::cerr << "      Wrong argument found for forest_type while agregating noise and main cfg files" << std::endl;
		exit(EXIT_FAILURE);
	}
	return cfg;
}

void generate_random(Config_Data cfg, std::vector<std::string> param_names, std::string dir_core, std::string file_out_modes, 
		std::string file_out_noise, std::string file_out_combi, int N_model,  std::string file_cfg_mm, std::string external_path, std::string templates_dir, std::string data_path){

	bool neg=0, passed=0;
	int i;
	long lastid, id0;
	std::string id_str, str_tmp;
	std::string template_file;
	VectorXd cte_params, val_min, val_max, input_params;
	MatrixXd currentcombi, allcombi;
	std::vector<double> pos_zero, pos_one;	
	std::vector<std::string> var_names, cte_names, distrib;
	// We first check that the cfg file has a coherent setup
	check_params(cfg, N_model);

	//  Define variables and constants
	pos_one=where(cfg.step, "=", 1, 0); // All positions of cfg.step that are equal to 1. Return position (last parameter is 0)
	pos_zero=where(cfg.step, "=", 0, 0); // All positions of cfg.step that are equal to 1. Return position (last parameter is 0)

	val_min.resize(pos_one.size());
	val_max.resize(pos_one.size());
	distrib.resize(pos_one.size()); // added on 6 Dec 2023
	cte_params.resize(pos_zero.size());
	for(int i=0; i<pos_one.size(); i++){
		var_names.push_back(cfg.labels[pos_one[i]]);
		val_min[i]=cfg.val_min[pos_one[i]];
		val_max[i]=cfg.val_max[pos_one[i]];
		distrib[i]=cfg.distrib[pos_one[i]]; // added on 6 Dec 2023
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
 	boost::normal_distribution<> dist_gr_gauss(1,1); // 6Dec2023: Gaussian random number of mean=1 and std=1 
	boost::variate_generator<boost::mt19937&, boost::normal_distribution<>> gaussian(generator, dist_gr_gauss);

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
    		template_file=templates_dir + cfg.template_files[dist(generator)];
    }	
    std::cout << "Selected template file: " << template_file << std::endl;	
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
	}
	id0=lastid+1;

	currentcombi.resize(1, pos_one.size());
	for(int c=0; c<cfg.forest_params[0]; c++){
		for(int i=0; i<pos_one.size();i++){
			if (distrib[i] == "Uniform"){
				currentcombi(0,i)=dist_gr() * (val_max[i] - val_min[i]) + val_min[i]; // HERE allcombi HAS JUST ONE LINE
			} else if (distrib[i] == "Gaussian"){
				currentcombi(0,i)=gaussian()*val_max[i] + val_min[i]; // In this context, val_min == mean, val_max == stddev
			} else{
				std::cerr << "Error while attempting to generate random numbers in iterative_artifical_spectrum::generate_random()" << std::endl;
				std::cerr << "     You are allowed to only have Uniform or Gaussian distributions. Found distribution: " << distrib[i] << std::endl;
			}
		}
		id_str=write_allcombi(currentcombi, cte_params, cfg, file_out_combi, cfg.erase_old_files, c, id0, cte_names, var_names,  param_names); // update/write the combination file
		std::cout << "Combination number: " << id_str << std::endl;

		input_params=order_input_params(cte_params, currentcombi.row(0), cte_names, var_names, param_names);

		passed=call_model_random(cfg.model_name, input_params, file_out_modes, file_out_noise,  file_cfg_mm, dir_core, id_str, cfg, external_path, template_file, data_path);
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
		std::string file_out_noise, std::string file_out_combi, int N_model, std::string data_path){

	bool passed=0;
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
	check_params(cfg, N_model);

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
		passed=call_model_grid(cfg.model_name, input_params, input_model, file_out_modes, file_out_noise, dir_core, id_str, cfg, data_path);
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
		std::string file_cfg_mm, std::string dir_core, std::string id_str, Config_Data cfg, std::string external_path, 
		std::string template_file, std::string data_path){

	std::string str;
	bool passed=0, subpassed=0;
	//std::string id_str;

	if(model_name == "generate_cfg_asymptotic_act_asym_Hgauss"){
		generate_cfg_asymptotic_act_asym_Hgauss(input_params, file_out_modes, file_out_noise);
		artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel, data_path);
		passed=1;
	}
	if(model_name =="generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma"){
		generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma(input_params, file_out_modes,  file_out_noise, cfg.extra_params); // extra_params must points towards a .in file
		artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel, data_path);
		passed=1;		
	}
	if(model_name =="generate_cfg_from_synthese_file_Wscaled_a1a2a3asymovGamma"){
		generate_cfg_from_synthese_file_Wscaled_a1a2a3asymovGamma(input_params, file_out_modes,  file_out_noise, cfg.extra_params); // extra_params must points towards a .in file
		artificial_spectrum_a1a2a3asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel, data_path);
		passed=1;		
	}
	if(model_name =="generate_cfg_from_synthese_file_Wscaled_Alm"){
		generate_cfg_from_synthese_file_Wscaled_Alm(input_params, file_out_modes,  file_out_noise, cfg.extra_params); // extra_params must points towards a .in file
		artificial_spectrum_a1Alma3(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, 
									cfg.write_inmodel, cfg.do_modelfiles, cfg.limit_data_range, cfg.modefile_modelname, data_path);
		passed=1;		
	}
	if(model_name =="generate_cfg_from_synthese_file_Wscaled_aj"){
		generate_cfg_from_synthese_file_Wscaled_aj(input_params, file_out_modes,  file_out_noise, cfg.extra_params); // extra_params must points towards a .in file
		artificial_spectrum_aj(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, 
									cfg.write_inmodel, cfg.do_modelfiles, cfg.limit_data_range, cfg.modefile_modelname, data_path, "harvey_like");
		passed=1;		
	}
	if(model_name =="generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled"){
		generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled(input_params, file_out_modes,  file_out_noise, cfg.extra_params); // extra_params must points towards a .in file
		artificial_spectrum_aj(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, 
									cfg.write_inmodel, cfg.do_modelfiles, cfg.limit_data_range, cfg.modefile_modelname, data_path);
		passed=1;		
	}
	if(model_name =="generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled_Kallinger2014"){
		size_t old_size=input_params.size();
		input_params.conservativeResize(old_size+2);
		input_params[old_size]=cfg.Tobs;
		input_params[old_size+1]=cfg.Cadence;
		generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled_Kallinger2014(input_params, file_out_modes,  file_out_noise, cfg.extra_params); // extra_params must points towards a .in file
		artificial_spectrum_aj(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, 
									cfg.write_inmodel, cfg.do_modelfiles, cfg.limit_data_range, cfg.modefile_modelname, data_path, "harvey_like");
		passed=1;		
	}
	if(model_name == "asymptotic_mm_v1" || model_name == "asymptotic_mm_v2" || model_name == "asymptotic_mm_v3" ||
		model_name == "asymptotic_mm_freeDp_numaxspread_curvepmodes_v1" || model_name == "asymptotic_mm_freeDp_numaxspread_curvepmodes_v2" ||
		model_name == "asymptotic_mm_freeDp_numaxspread_curvepmodes_v3" || model_name == "asymptotic_mm_freeDp_numaxspread_curvepmodes_v3_GRANscaled_Kallinger2014"){
		if(model_name =="asymptotic_mm_v1"){
			asymptotic_mm_v1(input_params, file_out_modes, file_out_noise,  file_cfg_mm, external_path, template_file);
			artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel, data_path);
			subpassed=1;
		}
		if(model_name =="asymptotic_mm_v2"){
			asymptotic_mm_v2(input_params, file_out_modes, file_out_noise,  file_cfg_mm, external_path, template_file);
			artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel, data_path);
			subpassed=1;
		}
		if(model_name =="asymptotic_mm_v3"){
			asymptotic_mm_v3(input_params, file_out_modes, file_out_noise,  file_cfg_mm, external_path, template_file);
			artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel, data_path);
			subpassed=1;
		}
		if(model_name =="asymptotic_mm_freeDp_numaxspread_curvepmodes_v1"){ 
			file_cfg_mm=""; // By-pass the definition to avoid issues
			asymptotic_mm_freeDp_numaxspread_curvepmodes_v1(input_params, file_out_modes, file_out_noise,  file_cfg_mm, external_path, template_file);
			artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel, data_path);
			subpassed=1;
		}
		if(model_name =="asymptotic_mm_freeDp_numaxspread_curvepmodes_v2"){ 
			file_cfg_mm=""; // By-pass the definiton to avoid issues
			asymptotic_mm_freeDp_numaxspread_curvepmodes_v2(input_params, file_out_modes, file_out_noise,  file_cfg_mm, external_path, template_file);
			artificial_spectrum_aj(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, 
									cfg.write_inmodel, cfg.do_modelfiles, cfg.limit_data_range, cfg.modefile_modelname, data_path);
			subpassed=1;
		}
		if(model_name =="asymptotic_mm_freeDp_numaxspread_curvepmodes_v3"){ 
			asymptotic_mm_freeDp_numaxspread_curvepmodes_v3(input_params, file_out_modes, file_out_noise,  file_cfg_mm, external_path, template_file);
			artificial_spectrum_aj(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, 
									cfg.write_inmodel, cfg.do_modelfiles, cfg.limit_data_range, cfg.modefile_modelname, data_path);
			subpassed=1;
		}
		if(model_name =="asymptotic_mm_freeDp_numaxspread_curvepmodes_v3_GRANscaled_Kallinger2014"){ 
			size_t old_size=input_params.size();
			input_params.conservativeResize(old_size+2);
			input_params[old_size]=cfg.Tobs;
			input_params[old_size+1]=cfg.Cadence;
			asymptotic_mm_freeDp_numaxspread_curvepmodes_v3_GRANscaled_Kallinger2014(input_params, file_out_modes, file_out_noise,  file_cfg_mm, external_path, template_file);
			artificial_spectrum_aj(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, 
									cfg.write_inmodel, cfg.do_modelfiles, cfg.limit_data_range, cfg.modefile_modelname, data_path, "harvey_like");
			subpassed=1;
		}
		if(subpassed == 0){
			std::cout << "Warning: The function call_model did not generate any configuration file!" << std::endl;
			std::cout << "         Debug required in the section handling the models asymptotic_mm_vX" << std::endl;
			std::cout << "         The program will stop now" << std::endl;
			exit(EXIT_FAILURE);
		}
		
		if(file_cfg_mm != ""){
			str="cp " + file_cfg_mm + " " + data_path + "/Spectra_info/" + strtrim(id_str) + ".global";
			const char *command = str.c_str(); 
			system(command);
		}
		passed=1;
	}

	return passed;
}


// This is the routine that calls models that are valid on the grid approach
bool call_model_grid(std::string model_name, VectorXd input_params, Model_data input_model, std::string file_out_modes, std::string file_out_noise, 
		std::string dir_core, std::string id_str, Config_Data cfg, std::string data_path){

	bool passed=0;
	std::string file_ref_star, modelID_str;
	std::ostringstream strs;

	if(model_name == "generate_cfg_asymptotic_act_asym_Hgauss"){
		generate_cfg_asymptotic_act_asym_Hgauss(input_params, file_out_modes, file_out_noise);
		//id_str=identifier2chain(identifier);
		modelID_str="NONE";
		artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel, data_path);
		//artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, dir_core, id_str, modelID_str, cfg.doplots, cfg.write_inmodel);
		passed=1;
	}
	if(model_name == "generate_cfg_from_refstar_HWscaled"){
		file_ref_star=dir_core + "Configurations/ref_spectra.params"; // The spectra used for the reference star
		generate_cfg_from_refstar_HWscaled(input_params, input_model, file_ref_star, file_out_modes, file_out_noise);
		strs << input_model.params[0];
		modelID_str=format_freqname(strs.str());
		artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel, data_path);
		//artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, dir_core, id_str, modelID_str, cfg.doplots, cfg.write_inmodel);
		passed=1;
	}
	if(model_name == "generate_cfg_from_refstar_HWscaled_GRANscaled"){
		file_ref_star=dir_core + "Configurations/ref_spectra.params"; // The spectra used for the reference star
		generate_cfg_from_refstar_HWscaled_GRANscaled(input_params, input_model, file_ref_star, file_out_modes, file_out_noise);
		strs << input_model.params[0];
		modelID_str=format_freqname(strs.str());
		artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel, data_path);
		//artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, dir_core, id_str, modelID_str, cfg.doplots, cfg.write_inmodel);
		passed=1;
	}

	if(model_name =="generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma"){
		generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma(input_params, file_out_modes,  file_out_noise, cfg.extra_params); // extra_params must points towards a .in file
		artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel, data_path);
		passed=1;		
	}
	if(model_name =="generate_cfg_from_synthese_file_Wscaled_a1a2a3asymovGamma"){
		generate_cfg_from_synthese_file_Wscaled_a1a2a3asymovGamma(input_params, file_out_modes,  file_out_noise, cfg.extra_params); // extra_params must points towards a .in file
		artificial_spectrum_a1a2a3asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel, data_path);
		passed=1;		
	}
	if(model_name =="generate_cfg_from_synthese_file_Wscaled_Alm"){
		generate_cfg_from_synthese_file_Wscaled_Alm(input_params, file_out_modes,  file_out_noise, cfg.extra_params); // extra_params must points towards a .in file
		artificial_spectrum_a1Alma3(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, 
									cfg.write_inmodel, cfg.do_modelfiles, cfg.limit_data_range, cfg.modefile_modelname, data_path);
		passed=1;		
	}
	if(model_name =="generate_cfg_from_synthese_file_Wscaled_aj"){
		generate_cfg_from_synthese_file_Wscaled_aj(input_params, file_out_modes,  file_out_noise, cfg.extra_params); // extra_params must points towards a .in file
		artificial_spectrum_aj(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, 
									cfg.write_inmodel, cfg.do_modelfiles, cfg.limit_data_range, cfg.modefile_modelname, data_path);
		passed=1;		
	}
	if(model_name =="generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled"){
		generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled(input_params, file_out_modes,  file_out_noise, cfg.extra_params); // extra_params must points towards a .in file
		artificial_spectrum_aj(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, 
									cfg.write_inmodel, cfg.do_modelfiles, cfg.limit_data_range, cfg.modefile_modelname, data_path);
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
#	elif (defined(__arm64__) && defined(__APPLE__)) || defined(__aarch64__)
		std::cout << " arm64 / Apple" << std::endl;
#   elif 
		std::cout << " Unknown" << std::endl;
#   endif
    std::cout << " Author: " << APP_COPYRIGHT << std::endl;
}


bool createDirectories(const std::string& output_dir, bool force_mkdir) {
    boost::filesystem::path out_dir(output_dir);
    if (force_mkdir == true) {
        try {
            boost::filesystem::create_directory(out_dir);
        } catch (const std::exception& e) {
            std::cerr << "Error creating directory: " << e.what() << std::endl;
            return false;
        }
    } else {
        if (!boost::filesystem::exists(out_dir)) {
            std::cerr << "Error: Output directory does not exist." << std::endl;
            return false;
        }
    }
    std::vector<std::string> subdirectories = {"Spectra_ascii", "Spectra_info", "Spectra_modelfile", "Spectra_plot"};
    for (const auto& subdirectory : subdirectories) {
        boost::filesystem::path subdirectory_path = out_dir / subdirectory;
        if (force_mkdir) {
            try {
                boost::filesystem::create_directory(subdirectory_path);
            } catch (const std::exception& e) {
                std::cerr << "Error creating subdirectory: " << e.what() << std::endl;
                return false;
            }
        } else {
            if (!boost::filesystem::exists(subdirectory_path)) {
                std::cerr << "Error: Subdirectory '" << subdirectory << "' does not exist." << std::endl;
                return false;
            }
        }
    }
    return true;
}


/**
 * @brief The main entry point of the program.
 *
 * This function is the main entry point of the program. It takes command line arguments and performs the necessary operations based on the provided options.
 *
 * @param argc The number of command line arguments.
 * @param argv An array of strings containing the command line arguments.
 * @return An integer representing the exit status of the program.
 */
int main(int argc, char* argv[]){

	boost::filesystem::path full_path( boost::filesystem::current_path() );

	// -------- Options Handling ------
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
		("version,v", "show program version")
        ("main_file,f", boost::program_options::value<std::string>()->default_value("main.cfg"), "Filename for the main configuration file. If not set, use the default filename.")
		("noise_file,n", boost::program_options::value<std::string>()->default_value("noise_Kallinger2014.cfg"), "Filename for the noise configuration file. If not set, use the default filename. Note that this is only for models with Kallinger+2014 noise at the moment.")
		("main_dir,g", boost::program_options::value<std::string>()->default_value("Configurations/"), "Full path for the main configuration file. If not set, use the default sub-directory 'Configurations/.")
		("out_dir,o", boost::program_options::value<boost::filesystem::path>()->default_value("Data/"), "Full path for the main configuration file. If not set, use the default sub-directory 'Data/.")
		("force-create-output-dir", boost::program_options::value<bool>()->default_value(0), "If set to 1=true, it will create the output directory defined by output_dir and all the required subdirectory. If set to 0=false (default), it will not create the directories, but only check if they exist and stop the program if they do not");	
	boost::program_options::variables_map vm;
	try {
        boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
        boost::program_options::notify(vm);
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return -1;
    }
	if (vm.count("help")) {
		std::cout << desc << std::endl;
		return 1;
	}
	if (vm.count("version")) {
        showversion();
		return 1;
    }
	// -------------------------------

	std::string cfg_file, cfg_noise_file;
	std::string main_f = vm["main_file"].as<std::string>();
	std::string noise_f = vm["noise_file"].as<std::string>();
	std::string main_dir = vm["main_dir"].as<std::string>();
	boost::filesystem::path out_dir = vm["out_dir"].as<boost::filesystem::path>();
	bool force_mkdir = vm["force-create-output-dir"].as<bool>();
	if(main_dir == "Configurations/"){
		cfg_file=full_path.string() + "/" + main_dir + main_f;
		cfg_noise_file=full_path.string() + "/" + main_dir + noise_f;
	} else{
		cfg_file=main_dir + main_f;
		cfg_noise_file=main_dir + noise_f;
	}
    if (out_dir.is_relative()) {
        boost::filesystem::path current_path = boost::filesystem::current_path();
        boost::filesystem::path full_path = current_path / out_dir;
        out_dir = full_path;
    }
    if (force_mkdir == true){
		std::cout << "           --force-create-output-dir = 1, Create the directories at the destination whenever possible: " << std::endl;
		std::cout << "          " << out_dir.string() << "..." << std::endl;
		std::cout << std::endl;
	} else{
		std::cout << "            --force-create-output-dir = 0, Checking if the directories at " << out_dir.string() << "  do already exist. The process will fail they don't... " << std::endl;
		std::cout << std::endl;
	}
	if (createDirectories(out_dir.string(), force_mkdir)) {
		std::cout << "Directories created successfully." << std::endl;
	} else {
		std::cerr << "Error creating directories." << std::endl;
		exit(EXIT_FAILURE);
	}
	iterative_artificial_spectrum(full_path.string() + "/", cfg_file, cfg_noise_file, out_dir.string());

}
