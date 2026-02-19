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
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cctype>
#include <random>
#include <cmath>
#include <filesystem>
#include <Eigen/Dense>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "artificial_spectrum.h"
#include "models_database.h"
#include "models_database_grid.h"
#include "version.h"
#include "combi.h"
#include "stellar_models.h"
#include "io_star_params.h"
#include "rng.h"
#include "logging.h"

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

static void cfg_validation_fail(const std::string& cfg_file, const Config_Data& cfg, const std::string& stage, const std::string& msg){
	LOG_ERROR("Configuration validation failed");
	LOG_ERROR("  stage     : " << stage);
	LOG_ERROR("  cfg_file  : " << cfg_file);
	LOG_ERROR("  model_name: " << cfg.model_name);
	LOG_ERROR("  details   : " << msg);
	exit(EXIT_FAILURE);
}

static void validate_cfg_vector_sizes(const Config_Data& cfg, const std::string& cfg_file, const std::string& stage){
	const size_t n = cfg.labels.size();
	if(n == 0){
		cfg_validation_fail(cfg_file, cfg, stage, "cfg.labels is empty");
	}
	if(cfg.val_min.size() != n || cfg.val_max.size() != n || cfg.step.size() != n || cfg.distrib.size() != n){
		std::ostringstream oss;
		oss << "Vector sizes mismatch: labels=" << n
			<< " val_min=" << cfg.val_min.size()
			<< " val_max=" << cfg.val_max.size()
			<< " step=" << cfg.step.size()
			<< " distrib=" << cfg.distrib.size();
		cfg_validation_fail(cfg_file, cfg, stage, oss.str());
	}
}

static void validate_cfg_forest_params(const Config_Data& cfg, const std::string& cfg_file, const std::string& stage){
	if(cfg.forest_params.size() < 1){
		cfg_validation_fail(cfg_file, cfg, stage, "cfg.forest_params is empty (expected at least one value)");
	}
	if(cfg.forest_type == "random"){
		if(cfg.forest_params[0] <= 0){
			cfg_validation_fail(cfg_file, cfg, stage, "random forest requires forest_params[0] > 0 (number of samples)");
		}
	}
}

static void validate_cfg_step_semantics(const Config_Data& cfg, const std::string& cfg_file, const std::string& stage){
	if(cfg.forest_type == "random"){
		for(size_t i=0; i<cfg.step.size(); i++){
			const double s = cfg.step[i];
			const bool is0 = std::abs(s - 0.0) < 1e-12;
			const bool is1 = std::abs(s - 1.0) < 1e-12;
			if(!is0 && !is1){
				std::ostringstream oss;
				oss << "random forest requires step values in {0,1}. Found step=" << s
					<< " for label='" << cfg.labels[i] << "' (index " << i << ")";
				cfg_validation_fail(cfg_file, cfg, stage, oss.str());
			}
		}
	}
}

static void validate_cfg_labels_unique(const Config_Data& cfg, const std::string& cfg_file, const std::string& stage){
	std::unordered_map<std::string, int> counts;
	for(const auto& s : cfg.labels){
		counts[s] += 1;
	}
	std::vector<std::string> dups;
	dups.reserve(cfg.labels.size());
	for(const auto& kv : counts){
		if(kv.second > 1){
			dups.push_back(kv.first);
		}
	}
	if(!dups.empty()){
		std::ostringstream oss;
		oss << "Duplicate labels are not allowed. Duplicates:";
		for(const auto& d : dups){
			oss << " " << d;
		}
		cfg_validation_fail(cfg_file, cfg, stage, oss.str());
	}
}

static void validate_cfg_for_model(const Config_Data& cfg, const std::string& cfg_file, const std::string& stage, const std::vector<std::string>& param_names){
	// param_names is the canonical ordered list expected by the code for the selected model.
	if(param_names.size() != cfg.labels.size()){
		std::unordered_set<std::string> provided(cfg.labels.begin(), cfg.labels.end());
		std::unordered_set<std::string> expected(param_names.begin(), param_names.end());

		std::vector<std::string> missing;
		std::vector<std::string> extra;
		missing.reserve(param_names.size());
		extra.reserve(cfg.labels.size());

		for(const auto& p : param_names){
			if(provided.find(p) == provided.end()){
				missing.push_back(p);
			}
		}
		for(const auto& l : cfg.labels){
			if(expected.find(l) == expected.end()){
				extra.push_back(l);
			}
		}

		std::ostringstream oss;
		oss << "Parameter list mismatch for model. expected=" << param_names.size() << " provided=" << cfg.labels.size();
		if(!missing.empty()){
			oss << "\n  Missing labels:";
			for(const auto& m : missing){ oss << " " << m; }
		}
		if(!extra.empty()){
			oss << "\n  Extra labels:";
			for(const auto& e : extra){ oss << " " << e; }
		}
		oss << "\n  Expected (ordered):";
		for(const auto& p : param_names){ oss << " " << p; }
		oss << "\n  Provided:";
		for(const auto& l : cfg.labels){ oss << " " << l; }
		cfg_validation_fail(cfg_file, cfg, stage, oss.str());
	}

	// Size matches; now verify set equality and duplicates.
	validate_cfg_labels_unique(cfg, cfg_file, stage);
	std::unordered_set<std::string> provided(cfg.labels.begin(), cfg.labels.end());
	std::unordered_set<std::string> expected(param_names.begin(), param_names.end());

	std::vector<std::string> missing;
	std::vector<std::string> extra;
	missing.reserve(param_names.size());
	extra.reserve(cfg.labels.size());

	for(const auto& p : param_names){
		if(provided.find(p) == provided.end()){
			missing.push_back(p);
		}
	}
	for(const auto& l : cfg.labels){
		if(expected.find(l) == expected.end()){
			extra.push_back(l);
		}
	}
	if(!missing.empty() || !extra.empty()){
		std::ostringstream oss;
		oss << "Parameter names do not match the selected model";
		if(!missing.empty()){
			oss << "\n  Missing labels:";
			for(const auto& m : missing){ oss << " " << m; }
		}
		if(!extra.empty()){
			oss << "\n  Extra labels:";
			for(const auto& e : extra){ oss << " " << e; }
		}
		oss << "\n  Expected (ordered):";
		for(const auto& p : param_names){ oss << " " << p; }
		oss << "\n  Provided:";
		for(const auto& l : cfg.labels){ oss << " " << l; }
		cfg_validation_fail(cfg_file, cfg, stage, oss.str());
	}
}

static void warn_negative_delta0l_percent(const Config_Data& cfg, const std::string& cfg_file, const std::string& stage){
	for(size_t i=0; i<cfg.labels.size(); i++){
		if(cfg.labels[i] == "delta0l_percent"){
			if(cfg.val_min[i] < 0 || cfg.val_max[i] < 0){
				LOG_WARN("Warning: delta0l_percent is negative in cfg");
				LOG_WARN("  stage     : " << stage);
				LOG_WARN("  cfg_file  : " << cfg_file);
				LOG_WARN("  model_name: " << cfg.model_name);
				LOG_WARN("  val_min   : " << cfg.val_min[i]);
				LOG_WARN("  val_max   : " << cfg.val_max[i]);
				LOG_WARN("  Note      : Convention expects a positive magnitude; values are negated internally.");
			}
			break;
		}
	}
}

static std::string to_lower_copy(const std::string& input){
	std::string out=input;
	std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
	return out;
}

static bool is_none_token(const std::string& value){
	const std::string v=to_lower_copy(strtrim(value));
	return (v == "none");
}

static bool is_all_token(const std::string& value){
	const std::string v=to_lower_copy(strtrim(value));
	return (v == "all" || v == "*");
}

static bool validate_template_file(const std::string& file_path, std::string* reason){
	std::ifstream in(file_path.c_str());
	if(!in.is_open()){
		if(reason){ *reason="cannot open file"; }
		return false;
	}

	bool has_numax=false;
	bool has_dnu=false;
	bool has_epsilon=false;
	bool has_data=false;

	std::string line;
	// Skip initial comment and empty lines
	while(std::getline(in, line)){
		line=strtrim(line);
		if(line.empty()){
			continue;
		}
		if(line[0] == '#'){
			continue;
		}
		break;
	}

	if(in.eof() && strtrim(line).empty()){
		if(reason){ *reason="file contains no content"; }
		return false;
	}

	// Read key/value lines until the data header (a line starting with '#')
	while(true){
		line=strtrim(line);
		if(line.empty()){
			if(!std::getline(in, line)){
				break;
			}
			continue;
		}
		if(line[0] == '#'){
			break;
		}
		const size_t pos=line.find('=');
		if(pos == std::string::npos){
			if(reason){ *reason="expected key=value entries before data table"; }
			return false;
		}
		const std::string key=strtrim(line.substr(0, pos));
		if(key == "numax_ref"){
			has_numax=true;
		} else if(key == "Dnu_ref"){
			has_dnu=true;
		} else if(key == "epsilon_ref"){
			has_epsilon=true;
		}
		if(!std::getline(in, line)){
			line="";
			break;
		}
	}
	if(!has_numax || !has_dnu || !has_epsilon){
		if(reason){
			std::ostringstream oss;
			oss << "missing required keyword(s):";
			if(!has_numax){ oss << " numax_ref"; }
			if(!has_dnu){ oss << " Dnu_ref"; }
			if(!has_epsilon){ oss << " epsilon_ref"; }
			*reason=oss.str();
		}
		return false;
	}

	// If the current line is a data line (not a comment), validate it too.
	if(!line.empty() && line[0] != '#'){
		std::vector<std::string> cols=strsplit(line, " \t");
		std::vector<std::string> filtered;
		for(size_t i=0; i<cols.size(); i++){
			const std::string c=strtrim(cols[i]);
			if(c != ""){ filtered.push_back(c); }
		}
		if(filtered.size() != 3){
			if(reason){ *reason="data row does not have exactly 3 columns"; }
			return false;
		}
		for(size_t i=0; i<filtered.size(); i++){
			std::istringstream iss(filtered[i]);
			double v=0;
			if(!(iss >> v)){
				if(reason){ *reason="non-numeric value in data row"; }
				return false;
			}
		}
		has_data=true;
	}

	while(std::getline(in, line)){
		line=strtrim(line);
		if(line.empty()){
			continue;
		}
		if(line[0] == '#'){
			continue;
		}
		std::vector<std::string> cols=strsplit(line, " \t");
		std::vector<std::string> filtered;
		for(size_t i=0; i<cols.size(); i++){
			const std::string c=strtrim(cols[i]);
			if(c != ""){ filtered.push_back(c); }
		}
		if(filtered.size() != 3){
			if(reason){ *reason="data row does not have exactly 3 columns"; }
			return false;
		}
		for(size_t i=0; i<filtered.size(); i++){
			std::istringstream iss(filtered[i]);
			double v=0;
			if(!(iss >> v)){
				if(reason){ *reason="non-numeric value in data row"; }
				return false;
			}
		}
		has_data=true;
	}
	if(!has_data){
		if(reason){ *reason="no data rows found"; }
		return false;
	}
	return true;
}


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

	// Ensure the temporary configuration directory exists.
	// Several generators write into Configurations/tmp/*.cfg and will fail if the directory is missing.
	{
		std::filesystem::path tmp_dir(dir_core + "Configurations/tmp");
		if (!std::filesystem::exists(tmp_dir)) {
			try {
				std::filesystem::create_directories(tmp_dir);
			} catch (const std::exception& e) {
				LOG_ERROR("Error creating temporary directory: " << tmp_dir.string());
				LOG_ERROR("Reason: " << e.what());
				exit(EXIT_FAILURE);
			}
		}
	}
	
	if(data_out_path == ""){
		data_path=dir_core + "/Data/";
	} else{
		data_path=data_out_path;
	}
	file_out_combi=data_path + "/Combinations.txt";

	std::string dir_freqs=dir_core + "external/MESA_grid/frequencies/";

	LOG_INFO("1. Read the configuration file...");

	cfg=read_main_cfg(cfg_file);
	validate_cfg_vector_sizes(cfg, cfg_file, "read_main_cfg");
	validate_cfg_forest_params(cfg, cfg_file, "read_main_cfg");
	validate_cfg_step_semantics(cfg, cfg_file, "read_main_cfg");
	validate_cfg_labels_unique(cfg, cfg_file, "read_main_cfg");
	cfg_noise=readNoiseConfigFile(cfg_noise_file); // Used to handle models with a separate noise, essentially the Kallinger+2014 noise 

	LOG_INFO("---------------------------------------------------------------------------------------");
	LOG_INFO("      ATTENTION: All model configuration should be contained into models_database.cpp  ");
	LOG_INFO("                 All mode generators should be contained into build_lorentzian.cpp     ");
	LOG_INFO("                 All model generators should be contained into artificial_spectrum.cpp ");
	LOG_INFO("---------------------------------------------------------------------------------------");

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
			LOG_ERROR("    Invalid number of parameters for model_name= 'generate_cfg_asymptotic_act_asym'");
			LOG_ERROR("    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size());
			LOG_ERROR("    Check your main configuration file");
			LOG_ERROR("    The program will exit now");
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
			LOG_ERROR("    Invalid number of parameters for model_name= 'generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma'");
			LOG_ERROR("    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size());
			LOG_ERROR("    Check your main configuration file");
			LOG_ERROR("    The program will exit now");
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
			LOG_ERROR("    Invalid number of parameters for model_name= 'generate_cfg_from_synthese_file_Wscaled_Alm'");
			LOG_ERROR("    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size());
			LOG_ERROR("    Check your main configuration file");
			LOG_ERROR("    The program will exit now");
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
			LOG_ERROR("    Invalid number of parameters for model_name= 'generate_cfg_from_synthese_file_Wscaled_aj'");
			LOG_ERROR("    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size());
			LOG_ERROR("    Check your main configuration file");
			LOG_ERROR("    The program will exit now");
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
			LOG_ERROR("    Invalid number of parameters for model_name= 'generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled'");
			LOG_ERROR("    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size());
			LOG_ERROR("    Check your main configuration file");
			LOG_ERROR("    The program will exit now");
			exit(EXIT_FAILURE);
		}
		passed=1;
	}
	if(cfg.model_name == "generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled_Kallinger2014"){
		// Warning: This model uses the file noise_Kallinger2014.cfg to set the noise parameters
		const int Nmodel_modes=16;
		const int Nmodel_noise=14; // 6 Dec 2023: Many of these parameters are in fact generated according to a Gaussian 
		Nmodel=Nmodel_modes+Nmodel_noise;
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
		if(param_names.size() != Nmodel){
			LOG_ERROR("    Invalid number of parameters for model_name= 'generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled_Kallinger2014'");
			LOG_ERROR("    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size());
			LOG_ERROR("    Check your main configuration file");
			LOG_ERROR("    The program will exit now");
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
			LOG_ERROR("    Invalid number of parameters for model_name= 'generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma'");
			LOG_ERROR("    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size());
			LOG_ERROR("    Check your main configuration file");
			LOG_ERROR("    The program will exit now");
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
			LOG_ERROR("    Invalid number of parameters for model_name= 'asymptotic_mm_v1'");
			LOG_ERROR("    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size());
			LOG_ERROR("    Check your main configuration file");
			LOG_ERROR("    The program will exit now");
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
			LOG_ERROR("    Invalid number of parameters for model_name= 'asymptotic_mm_v2'");
			LOG_ERROR("    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size());
			LOG_ERROR("    Check your main configuration file");
			LOG_ERROR("    The program will exit now");
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
			LOG_ERROR("    Invalid number of parameters for model_name= 'asymptotic_mm_v3'");
			LOG_ERROR("    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size());
			LOG_ERROR("    Check your main configuration file");
			LOG_ERROR("    The program will exit now");
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
			LOG_ERROR("    Invalid number of parameters for model_name= 'asymptotic_mm_v1'");
			LOG_ERROR("    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size());
			LOG_ERROR("    Check your main configuration file");
			LOG_ERROR("    The program will exit now");
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
			LOG_ERROR("    Invalid number of parameters for model_name= 'asymptotic_mm_v2'");
			LOG_ERROR("    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size());
			LOG_ERROR("    Check your main configuration file");
			LOG_ERROR("    The program will exit now");
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
			LOG_ERROR("    Invalid number of parameters for model_name= 'asymptotic_mm_v3'");
			LOG_ERROR("    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size());
			LOG_ERROR("    Check your main configuration file");
			LOG_ERROR("    The program will exit now");
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
			LOG_ERROR("    Invalid number of parameters for model_name= 'generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled_Kallinger2014'");
			LOG_ERROR("    Expecting " << Nmodel_modes + Nmodel_noise << " parameters, but found " << cfg.val_min.size());
			LOG_ERROR("    Check your main configuration file");
			LOG_ERROR("    The program will exit now");
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
			LOG_ERROR("    Invalid number of parameters for model_name= 'generate_cfg_asymptotic_act_asym'");
			LOG_ERROR("    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size());
			LOG_ERROR("    Check your main configuration file");
			LOG_ERROR("    The program will exit now");
			exit(EXIT_FAILURE);
		}
		LOG_INFO(" 	   - reading models from model file (requires for setting frequencies)...");
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
			LOG_ERROR("    Invalid number of parameters for model_name= 'generate_cfg_asymptotic_act_asym'");
			LOG_ERROR("    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size());
			LOG_ERROR("    Check your main configuration file");
			LOG_ERROR("    The program will exit now");
			exit(EXIT_FAILURE);
		}
		LOG_INFO(" 	   - reading models from model file (requires for setting frequencies)...");
		models=read_data_ascii_Ncols(model_file, delimiter, verbose_data);
		usemodels=1;
		passed=1;
	}
 	if(passed == 0){
		LOG_INFO("    model_name= " << cfg.model_name << " is not a recognized keyword for models");
		LOG_INFO("    Check models_database.h to see the available model names");
		LOG_ERROR("    The program will exit now");
		exit(EXIT_FAILURE);
	}

	// -----------------------------------------------------------
	// -----------------------------------------------------------
	// -----------------------------------------------------------


	// Validate cfg after any model-specific augmentation (e.g. noise cfg aggregation).
	validate_cfg_vector_sizes(cfg, cfg_file, "post_model_selection");
	validate_cfg_forest_params(cfg, cfg_file, "post_model_selection");
	validate_cfg_step_semantics(cfg, cfg_file, "post_model_selection");
	validate_cfg_for_model(cfg, cfg_file, "post_model_selection", param_names);
	warn_negative_delta0l_percent(cfg, cfg_file, "post_model_selection");
	Nmodel=static_cast<int>(param_names.size());

	LOG_INFO("2. Generating the models using the subroutine " << cfg.model_name << " of model_database.cpp...");
	if(cfg.forest_type == "random"){
		LOG_INFO("   Values are randomly generated into a uniform range defined in the main configuration file");
		generate_random(cfg, param_names, dir_core, file_out_modes, file_out_noise, file_out_combi, Nmodel, file_cfg_mm, external_path, templates_dir, data_path);

	}
	if(cfg.forest_type == "grid"){
		LOG_INFO("   Values are generated over a grid using all possible combinations according to inputs in the main configuration file");
		generate_grid(cfg, usemodels, models, param_names, dir_core, dir_freqs, file_out_modes, file_out_noise, file_out_combi, Nmodel,data_path);
	}
	if(cfg.forest_type != "random" && cfg.forest_type != "grid"){ 
		LOG_INFO(" Problem in the main configuration file. It is expected that the forest type parameters is either random OR grid");
		LOG_ERROR(" Check your main configuration file");
		LOG_ERROR(" The program will exit now");
		exit(EXIT_FAILURE);
	}

}

void check_params(Config_Data cfg, const int N_model){
	bool neg=0;
	int ii=0;
	while(neg == 0 && (ii < N_model)){
		if( (cfg.distrib[ii] == "Uniform") && (cfg.val_max[ii] - cfg.val_min[ii]) < 0){ 
			neg=1;
			LOG_INFO("List of the parameters and values :");
			for (int i=0;i<cfg.labels.size();i++){
				LOG_INFO("[" << i << "] " << std::setw(20) << cfg.labels[i] << std::setw(12)<< cfg.val_min[i] << std::setw(12) << cfg.val_max[i]);
			}
			LOG_INFO(" ------- ");
			LOG_INFO("       ii = " << ii);
			LOG_INFO("       name:    " << cfg.labels[ii]);
			LOG_INFO("       val_max =" << cfg.val_max[ii]);
			LOG_INFO("       val_min =" << cfg.val_min[ii]);
		}
		ii=ii+1;
	}
	if(neg == 1){
		LOG_WARN("     Warning: val_max < val_min for some of the parameters while a uniform distribution is requested!");
		LOG_INFO("     This is not permitted");
		LOG_ERROR("     Check your main configuration file");
		LOG_ERROR("     The program will exit now");
		exit(EXIT_FAILURE);
	}
}

Config_Data update_cfg(Config_Data cfg_target, const Config_Data cfg_source, const int priority4common){
	bool pass=false;
	if (priority4common <-1 || priority4common >1){
		LOG_ERROR("Error: Invalid pririty4common value. Set it to -1, 0 or 1.");
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
				LOG_ERROR("Error: Common variables in cfg_target and cfg_source are not the same.");
				exit(EXIT_FAILURE);
			}
	    // Check each element of the template_files vector
    	for (size_t i = 0; i < cfg_target.template_files.size(); i++) {
			if (cfg_target.template_files[i] != cfg_source.template_files[i]) {
				LOG_ERROR("Error: template_files in cfg_target and cfg_source are not the same.");
				exit(EXIT_FAILURE);
			}
		}
		// Check each element of the forest_params vector
    	for (size_t i = 0; i < cfg_target.forest_params.size(); i++) {
			if (cfg_target.forest_params[i] != cfg_source.forest_params[i]) {
				LOG_ERROR("Error: forest_params in cfg_target and cfg_source are not the same.");
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
				LOG_ERROR("Error in iterative_artificial_spectrum::agregate_maincfg_noisecfg():");
				LOG_ERROR("    In random mode and when agregation of a main.cfg with a noise.cfg is requested, only Gaussian, Fix and Uniform is a valid entry for the distribution parameter");
				LOG_ERROR("    Check the content of the used noise configuration file.");
				LOG_ERROR("    The program will exit now");
				exit(EXIT_FAILURE);
			}
		}
		pass=true;
	}
	if (cfg.forest_type == "grid"){
		for (int i=0; i<cfg_noise.name_grid.size();i++){
			cfg.labels.push_back(cfg_noise.name_grid[i]);
			cfg.distrib.push_back(cfg_noise.distrib_grid[i]);
			cfg.val_min.push_back(cfg_noise.x1_grid[i]);
			cfg.val_max.push_back(cfg_noise.x2_grid[i]);
			if (cfg_noise.distrib_grid[i] == "Fix"){
				cfg.step.push_back(0);
			} else if (cfg_noise.distrib_grid[i] != "Uniform"){
					LOG_ERROR("Error in iterative_artificial_spectrum::agregate_maincfg_noisecfg():");
					LOG_ERROR("    In grid mode and when agregation of a main.cfg with a noise.cfg is requested, only Fix and Uniform is a valid entry for the distribution parameter");
					LOG_ERROR("    Check the content of the used noise configuration file.");
					LOG_ERROR("    The program will exit now");
					exit(EXIT_FAILURE);
			} else{ // In the uniform case, the kerror_grid is used to set the step
				cfg.step.push_back(cfg_noise.kerror_grid[i]);
			}
		}
		pass=true;
	}
	if (pass == false){
		LOG_ERROR("Error in iterative_artificial_spectrum::agregate_maincfg_noisecfg():");
		LOG_ERROR("      Wrong argument found for forest_type while agregating noise and main cfg files");
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
	{
		std::ostringstream oss;
		oss << "Constants: ";
		for(int i=0; i<cte_names.size(); i++){ oss << cte_names[i] << "  ";}
		LOG_INFO(oss.str());
	}

	{
		std::ostringstream oss;
		oss << "Variables: ";
		for(int i=0; i<var_names.size(); i++){ oss << var_names[i] << "  ";}
		LOG_INFO(oss.str());
	}

	//  ------------ Generate the random values -----------
	// Initialize the random generators: Uniform over [0,1]. This is used for the random parameters
	std::mt19937& generator = global_rng();
	std::uniform_real_distribution<double> dist_gr(0.0, 1.0);
	std::normal_distribution<double> dist_gr_gauss(1.0, 1.0); // 6Dec2023: Gaussian random number of mean=1 and std=1 

	// Generator of integers for selecting randomly a template file that contains a height and width profile
	std::vector<std::string> template_candidates=cfg.template_files;
	if(template_candidates.size() == 0){
		LOG_ERROR("Error: The template_file cannot be empty. If not used, please set it to NONE");
		exit(EXIT_SUCCESS);
	}
	if(template_candidates.size() == 1 && is_none_token(template_candidates[0])){
		template_file="";
	} else{
		if(is_all_token(template_candidates[0])){ // If the user specifies that all *.template files should be used
			template_candidates=list_dir(templates_dir, "template");
		}
		if (template_candidates.size() == 0){
			LOG_INFO("Could not find the template file in the expected directory");
			LOG_INFO("Be sure to have it in Configurations/templates/ ");
			LOG_INFO("If the used mode does not require templates, specify NONE in the dedicated field ");
			exit(EXIT_FAILURE);
		}
		std::vector<std::string> valid_templates;
		for(size_t i=0; i<template_candidates.size(); i++){
			const std::string candidate=strtrim(template_candidates[i]);
			if(candidate == "" || is_none_token(candidate)){
				continue;
			}
			const std::string full_path=templates_dir + candidate;
			std::string reason;
			if(validate_template_file(full_path, &reason)){
				valid_templates.push_back(candidate);
			} else{
				LOG_ERROR("Warning: Skipping invalid template file: " << full_path);
				LOG_ERROR("         Reason: " << reason);
			}
		}
		if(valid_templates.size() == 0){
			LOG_ERROR("Error: No valid template files found after validation");
			LOG_ERROR("       Check templates in: " << templates_dir);
			LOG_ERROR("       If the used mode does not require templates, specify NONE in the dedicated field");
			exit(EXIT_FAILURE);
		}
		std::uniform_int_distribution<int> dist(0, static_cast<int>(valid_templates.size()-1));
		template_file=templates_dir + valid_templates[dist(generator)];
	}
	if(template_file == ""){
		LOG_INFO("Selected template file: NONE");
	} else{
		LOG_INFO("Selected template file: " << template_file);
	}
 	LOG_INFO(" -----------------------------------------------------");
	LOG_INFO(" List of all combinations written iteratively into ");
	LOG_INFO("       " <<  file_out_combi);
	LOG_INFO(" -----------------------------------------------------");

	if(cfg.erase_old_files == 0){
		if(file_exists(file_out_combi) == 1){
			LOG_INFO("                 erase_old_files=0...");
			LOG_INFO("                 ...Older combinatory file found!");
			LOG_INFO("                 Name of the found file: " << file_out_combi);
			LOG_INFO("                 Reading the combinatory file in order to determince the last value of the samples...");
			lastid=read_id_allcombi(file_out_combi);
		} else{
			lastid=0; // If no Combi file while erase_old_files is set to 1
			LOG_INFO("                 erase_old_files=0...");
			LOG_INFO("                 ... but no older combi file was found");
			LOG_INFO("                 The program will therefore behave as if erase_old_files=1");
		}
	} else {
		lastid=0; // If no Combi file while erase_old_files is set to 1
		LOG_INFO("                 erase_old_files=1...");
		LOG_INFO("                 Note that no deleting action is performed by this program");
		LOG_INFO("                 But any older file might be overwritten!");
		LOG_INFO("                 Be sure to have no important previous results in the Data directory and subdirectories!");
	}
	id0=lastid+1;

	currentcombi.resize(1, pos_one.size());
	for(int c=0; c<cfg.forest_params[0]; c++){
		for(int i=0; i<pos_one.size();i++){
			if (distrib[i] == "Uniform"){
				currentcombi(0,i)=dist_gr(generator) * (val_max[i] - val_min[i]) + val_min[i]; // HERE allcombi HAS JUST ONE LINE
			} else if (distrib[i] == "Gaussian"){
				currentcombi(0,i)=dist_gr_gauss(generator)*val_max[i] + val_min[i]; // In this context, val_min == mean, val_max == stddev
			} else{
				LOG_ERROR("Error while attempting to generate random numbers in iterative_artifical_spectrum::generate_random()");
				LOG_ERROR("     You are allowed to only have Uniform or Gaussian distributions. Found distribution: " << distrib[i]);
			}
		}
		id_str=write_allcombi(currentcombi, cte_params, cfg, file_out_combi, cfg.erase_old_files, c, id0, cte_names, var_names,  param_names); // update/write the combination file
		LOG_INFO("Combination number: " << id_str);

		input_params=order_input_params(cte_params, currentcombi.row(0), cte_names, var_names, param_names);

		passed=call_model_random(cfg.model_name, input_params, file_out_modes, file_out_noise,  file_cfg_mm, dir_core, id_str, cfg, external_path, template_file, data_path);
		if(passed == 0){
			LOG_WARN("Warning: The function call_model did not generate any configuration file!");
			LOG_INFO("         It is very likely that you tried to start a model in 'random mode' while this model is not available for the random approach");
			LOG_INFO("         The program will stop now");
			exit(EXIT_FAILURE);
		}
		id0=id0+1;
		passed=0;

	// Erasing the temporary files to avoid those to be loaded at vitam eternam if an error in their generetation at step c=c_err happens
	std::error_code ec;
	std::filesystem::remove(file_out_modes, ec);
	std::filesystem::remove(file_out_noise, ec);
	}

	LOG_INFO("All requested tasks executed succesfully");
	
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
		//LOG_INFO(pos_one[ii]);
		Nvals[ii]=delta[pos_one[ii]]/cfg.step[pos_one[ii]] + 1e-15; // 1d-15 here because floor bugs sometimes due to round-off
		//LOG_INFO("Nvals["<< ii << "]=" << Nvals[ii] ; //Nvals=Nvals.floor());
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
	{
		std::ostringstream oss;
		oss << "Constants: ";
		for(int i=0; i<cte_names.size(); i++){ oss << cte_names[i] << "  ";}
		LOG_INFO(oss.str());
	}

	{
		std::ostringstream oss;
		oss << "Variables: ";
		Nvar_params=var_names.size();
		for(int i=0; i<Nvar_params; i++){ oss << var_names[i] << "  ";}
		LOG_INFO(oss.str());
	}

	//  ------------ Generate the grid -----------
	LOG_INFO("       - Generating all the requested combinations. This may take a while...");
	var_params.transposeInPlace();
	allcombi=define_all_combinations(var_params, Nvals, Nvar_params);
	Ncombi=allcombi.rows();
	LOG_INFO(allcombi);

	LOG_INFO("Ncombi = " << Ncombi);
	LOG_INFO(" -----------------------------------------------------");
	LOG_INFO(" List of all combinations written iteratively into ");
	LOG_INFO("       " <<  file_out_combi);
	LOG_INFO(" -----------------------------------------------------");

	if(cfg.erase_old_files == 0){
		if(file_exists(file_out_combi) == 1){
			LOG_INFO("                 erase_old_files=0...");
			LOG_INFO("                 ...Older combinatory file found!");
			LOG_INFO("                 Name of the found file: " << file_out_combi);
			LOG_INFO("                 Reading the combinatory file in order to determince the last value of the samples...");
			lastid=read_id_allcombi(file_out_combi);

			//exit(EXIT_SUCCESS);

		} else{
			lastid=-1; // If no Combi file while erase_old_files is set to 1
			LOG_INFO("                 erase_old_files=0...");
			LOG_INFO("                 ... but no older combi file was found");
			LOG_INFO("                 The program will therefore behave as if erase_old_files=0");
		}
	} else {
		lastid=-1; // If no Combi file while erase_old_files is set to 1
		LOG_INFO("                 erase_old_files=1...");
		LOG_INFO("                 Note that no deleting action is performed by this program");
		LOG_INFO("                 But any older file might be overwritten!");
		LOG_INFO("                 Be sure to have no important previous results in the Data directory and subdirectories!");
	}
	id0=lastid+1;
	//id0=1; // DEBUG ONLY

	currentcombi.resize(1, pos_one.size());
	for(int c=0; c<Ncombi; c++){
		for(int i=0; i<pos_one.size();i++){
			currentcombi.row(0)=allcombi.row(id0); // HERE currentcombi HAS JUST ONE LINE
		}
		id_str=write_allcombi(currentcombi, cte_params, cfg, file_out_combi, cfg.erase_old_files, c, id0, cte_names, var_names,  param_names); // update/write the combination file
		LOG_INFO("Combination number: " << id_str  << "  (last is: " << Ncombi-1 << ")");
	
		input_params=order_input_params(cte_params, currentcombi.row(0), cte_names, var_names, param_names);
		if(usemodels == 1){
			pos=where_strXi(param_names, "Model_i"); // Look for the position where Model_i is given (should be 0)	
			//LOG_INFO(pos);
			if(pos[0] !=-1){
				input_model=get_model_param(models, input_params[pos[0]], dir_freqs);
			} else{
				LOG_INFO("Could not find the column that contains the 'Model_i' (model index)");
				LOG_ERROR("Debug required. The program will exit now");
				exit(EXIT_SUCCESS);
			}
		}
		passed=call_model_grid(cfg.model_name, input_params, input_model, file_out_modes, file_out_noise, dir_core, id_str, cfg, data_path);
		if(passed == 0){
			LOG_WARN("Warning: The function call_model_grid did not generate any configuration file!");
			LOG_INFO("         It is very likely that you tried to start a model in 'grid mode' while this model is not available for the grid approach");
			LOG_INFO("         The program will stop now");
			exit(EXIT_FAILURE);
		}
		id0=id0+1;
		passed=0;

		// Erasing the temporary files to avoid those to be loaded at vitam eternam if an error in their generetation at step c=c_err happens
		std::error_code ec;
		std::filesystem::remove(file_out_modes, ec);
		std::filesystem::remove(file_out_noise, ec);
	}

	LOG_INFO("All requested tasks executed succesfully");
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
			bool failed=asymptotic_mm_v1(input_params, file_out_modes, file_out_noise,  file_cfg_mm, external_path, template_file);
			if(failed == false){
				artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel, data_path);
			}
			subpassed=1;
		}
		if(model_name =="asymptotic_mm_v2"){
			bool failed=asymptotic_mm_v2(input_params, file_out_modes, file_out_noise,  file_cfg_mm, external_path, template_file);
			if(failed == false){
				artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel, data_path);
			}
			subpassed=1;
		}
		if(model_name =="asymptotic_mm_v3"){
			bool failed=asymptotic_mm_v3(input_params, file_out_modes, file_out_noise,  file_cfg_mm, external_path, template_file);
			if(failed == false){
				artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel, data_path);
			}
			subpassed=1;
		}
		if(model_name =="asymptotic_mm_freeDp_numaxspread_curvepmodes_v1"){ 
			file_cfg_mm=""; // By-pass the definition to avoid issues
			bool failed=asymptotic_mm_freeDp_numaxspread_curvepmodes_v1(input_params, file_out_modes, file_out_noise,  file_cfg_mm, external_path, template_file);
			if(failed == false){
				artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, cfg.write_inmodel, data_path);
			}
			subpassed=1;
		}
		if(model_name =="asymptotic_mm_freeDp_numaxspread_curvepmodes_v2"){ 
			file_cfg_mm=""; // By-pass the definiton to avoid issues
			bool failed=asymptotic_mm_freeDp_numaxspread_curvepmodes_v2(input_params, file_out_modes, file_out_noise,  file_cfg_mm, external_path, template_file);
			if(failed == false){
				artificial_spectrum_aj(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, 
									cfg.write_inmodel, cfg.do_modelfiles, cfg.limit_data_range, cfg.modefile_modelname, data_path);
			}
			subpassed=1;
		}
		if(model_name =="asymptotic_mm_freeDp_numaxspread_curvepmodes_v3"){ 
			bool failed=asymptotic_mm_freeDp_numaxspread_curvepmodes_v3(input_params, file_out_modes, file_out_noise,  file_cfg_mm, external_path, template_file);
			if(failed == false){
				artificial_spectrum_aj(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, 
									cfg.write_inmodel, cfg.do_modelfiles, cfg.limit_data_range, cfg.modefile_modelname, data_path);
			}
			subpassed=1;
		}
		if(model_name =="asymptotic_mm_freeDp_numaxspread_curvepmodes_v3_GRANscaled_Kallinger2014"){ 
			size_t old_size=input_params.size();
			input_params.conservativeResize(old_size+2);
			input_params[old_size]=cfg.Tobs;
			input_params[old_size+1]=cfg.Cadence;
			bool failed=asymptotic_mm_freeDp_numaxspread_curvepmodes_v3_GRANscaled_Kallinger2014(input_params, file_out_modes, file_out_noise,  file_cfg_mm, external_path, template_file);
			if(failed == false){
				artificial_spectrum_aj(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, dir_core, id_str, cfg.doplots, 
									cfg.write_inmodel, cfg.do_modelfiles, cfg.limit_data_range, cfg.modefile_modelname, data_path, "harvey_like");
			}
			subpassed=1;
		}
		if(subpassed == 0){
			LOG_WARN("Warning: The function call_model did not generate any configuration file!");
			LOG_INFO("         Debug required in the section handling the models asymptotic_mm_vX");
			LOG_INFO("         The program will stop now");
			exit(EXIT_FAILURE);
		}
		
		if(file_cfg_mm != ""){
			std::filesystem::path src(file_cfg_mm);
			std::filesystem::path dst = std::filesystem::path(data_path) / "Spectra_info" / (strtrim(id_str) + ".global");
			std::error_code ec;
			std::filesystem::copy_file(src, dst, std::filesystem::copy_options::overwrite_existing, ec);
			if(ec){
				LOG_ERROR("Warning: failed to copy mixed-mode cfg file: " << src.string());
				LOG_ERROR("         to: " << dst.string());
				LOG_ERROR("         reason: " << ec.message());
			}
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

	std::vector<std::string> files;
	std::filesystem::path dir_path(path);
	if(!std::filesystem::exists(dir_path)){
		LOG_INFO("I/O Error: directory does not exist: " << path);
		exit(EXIT_FAILURE);
	}
	if(!std::filesystem::is_directory(dir_path)){
		LOG_INFO("I/O Error: not a directory: " << path);
		exit(EXIT_FAILURE);
	}

	const std::string ext_lower=to_lower_copy(extension);
	for(const auto& entry : std::filesystem::directory_iterator(dir_path)){
		if(!entry.is_regular_file()){
			continue;
		}
		const std::filesystem::path p=entry.path();
		if(!p.has_extension()){
			continue;
		}
		std::string ext=p.extension().string();
		if(!ext.empty() && ext[0] == '.'){
			ext=ext.substr(1);
		}
		if(to_lower_copy(ext) == ext_lower){
			files.push_back(p.filename().string());
		}
	}
	std::sort(files.begin(), files.end());
    return files;
}

void showversion()
{
	std::ostringstream oss;
	oss << APP_NAME << " " << APP_VERSION << "\n built on " << __DATE__;

#   if defined(__clang__)
	oss << " with clang " << __clang_version__;
#   elif defined(__GNUC__)
	oss << " with GCC";
	oss << " " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__;
#   elif defined(_MSC_VER)
	oss << " with MSVC";
	oss << " " << MSVC_VERSION;
#   else
	oss << " unknown compiler";
#   endif

	oss << "\n features:";
#   if defined(__i386__) || defined(_M_IX86)
	oss << " i386";
#   elif defined(__x86_64__) || defined(_M_AMD64)
	oss << " x86_64";
#	elif (defined(__arm64__) && defined(__APPLE__)) || defined(__aarch64__)
	oss << " arm64 / Apple";
#   else
	oss << " Unknown";
#   endif

	oss << "\n Author: " << APP_COPYRIGHT;
	LOG_INFO(oss.str());
}


bool createDirectories(const std::string& output_dir, bool force_mkdir) {
    boost::filesystem::path out_dir(output_dir);
    if (force_mkdir == true) {
        try {
            boost::filesystem::create_directory(out_dir);
        } catch (const std::exception& e) {
            LOG_ERROR("Error creating directory: " << e.what());
            return false;
        }
    } else {
        if (!boost::filesystem::exists(out_dir)) {
            LOG_ERROR("Error: Output directory does not exist.");
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
                LOG_ERROR("Error creating subdirectory: " << e.what());
                return false;
            }
        } else {
            if (!boost::filesystem::exists(subdirectory_path)) {
                LOG_ERROR("Error: Subdirectory '" << subdirectory << "' does not exist.");
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
		("out_dir,o", boost::program_options::value<boost::filesystem::path>()->default_value("Data/"), "Full path or relative path for the outputs. If not set, use the default sub-directory 'Data/.")
		("force-create-output-dir", boost::program_options::value<bool>()->default_value(0), "If set to 1=true, it will create the output directory defined by output_dir and all the required subdirectory. If set to 0=false (default), it will not create the directories, but only check if they exist and stop the program if they do not")
		("log-level", boost::program_options::value<std::string>()->default_value("info"), "Log level: debug, info, warn, error")
		("seed", boost::program_options::value<long long>()->default_value(-1), "Seed for RNG (>=0). If set, results are deterministic across runs.");	
	boost::program_options::variables_map vm;
	try {
        boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
        boost::program_options::notify(vm);
	} catch (const std::exception& ex) {
		std::cerr << ex.what() << std::endl;
		return -1;
	}

	LogLevel log_level = LogLevel::info;
	std::string log_level_str = vm["log-level"].as<std::string>();
	if (!try_parse_log_level(log_level_str, &log_level)) {
		std::cerr << "Invalid --log-level value: " << log_level_str << std::endl;
		std::cerr << "Valid values: debug, info, warn, error" << std::endl;
		return -1;
	}
	init_logging(log_level);
	LOG_INFO("Log level set to " << log_level_str);

	if (vm.count("help")) {
		set_log_level(LogLevel::info);
		LOG_INFO(desc);
		return 1;
	}
	if (vm.count("version")) {
		set_log_level(LogLevel::info);
		showversion();
		return 1;
	}
	long long seed_value = vm["seed"].as<long long>();
	if(seed_value >= 0){
		set_global_seed(static_cast<uint64_t>(seed_value));
		LOG_INFO("Using RNG seed: " << seed_value);
	}
	// -------------------------------

	std::string cfg_file, cfg_noise_file;
	std::string main_f = vm["main_file"].as<std::string>();
	std::string noise_f = vm["noise_file"].as<std::string>();
	std::string main_dir = vm["main_dir"].as<std::string>();
	boost::filesystem::path out_dir = vm["out_dir"].as<boost::filesystem::path>();
	bool force_mkdir = vm["force-create-output-dir"].as<bool>();
	boost::filesystem::path main_dir_path(main_dir);
	boost::filesystem::path main_cfg_path(main_f);
	if(main_cfg_path.is_absolute() || main_cfg_path.has_parent_path()){
		cfg_file = main_cfg_path.string();
	}else{
		cfg_file = (main_dir_path / main_cfg_path).string();
	}

	boost::filesystem::path noise_cfg_path(noise_f);
	if(noise_cfg_path.is_absolute() || noise_cfg_path.has_parent_path()){
		cfg_noise_file = noise_cfg_path.string();
	}else{
		boost::filesystem::path noise_in_main_dir = main_dir_path / noise_cfg_path;
		boost::filesystem::path noise_in_default_dir = full_path / "Configurations" / noise_cfg_path;
		if(boost::filesystem::exists(noise_in_main_dir)){
			cfg_noise_file = noise_in_main_dir.string();
		}else if(boost::filesystem::exists(noise_in_default_dir)){
			cfg_noise_file = noise_in_default_dir.string();
		}else{
			cfg_noise_file = noise_in_main_dir.string();
		}
	}
    if (out_dir.is_relative()) {
        boost::filesystem::path current_path = boost::filesystem::current_path();
        boost::filesystem::path full_path = current_path / out_dir;
        out_dir = full_path;
    }
    if (force_mkdir == true){
		LOG_INFO("           --force-create-output-dir = 1, Create the directories at the destination whenever possible: ");
		LOG_INFO("          " << out_dir.string() << "...");
		LOG_INFO("");
	} else{
		LOG_INFO("            --force-create-output-dir = 0, Checking if the directories at " << out_dir.string() << "  do already exist. The process will fail they don't... ");
		LOG_INFO("");
	}
	if (createDirectories(out_dir.string(), force_mkdir)) {
		LOG_INFO("Directories created successfully.");
	} else {
		LOG_ERROR("Error creating directories.");
		exit(EXIT_FAILURE);
	}
	iterative_artificial_spectrum(full_path.string() + "/", cfg_file, cfg_noise_file, out_dir.string());

}
