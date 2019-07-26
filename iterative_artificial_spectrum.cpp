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
#include <Eigen/Dense>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include "artificial_spectrum.h"
#include "models_database.h"

void iterative_artificial_spectrum(std::string dir_core);
std::vector<double> where(std::vector<double> vec, std::string condition, double value, bool return_values);
VectorXi where_str(std::vector<std::string> vec, std::string value);
VectorXd order_input_params(VectorXd cte_params, VectorXd var_params, std::vector<std::string> cte_names, 
		std::vector<std::string> var_names, std::vector<std::string> param_names);
void generate_random(Config_Data cfg, std::vector<std::string> param_names, std::string dir_core, std::string file_out_modes, 
		std::string file_out_noise, std::string file_out_combi, int N_model, std::string file_cfg_mm, std::string external_path);

void generate_grid();
std::string write_allcombi(MatrixXd allcombi, VectorXd cte_params, Config_Data cfg, std::string fileout, bool erase_old_file, long iter, long id0, 
		    std::vector<std::string> cte_names, std::vector<std::string> var_names, std::vector<std::string> param_names);
bool call_model(std::string model_name, VectorXd input_params, std::string file_out_modes, std::string file_out_noise, 
		 std::string file_cfg_mm, std::string dir_core, std::string identifier, Config_Data cfg, std::string external_path);
std::string identifier2chain(std::string identifier);
std::string identifier2chain(long identifier);
long read_id_allcombi(std::string file_combi);
std::string read_lastline_ascii(std::string filename);

void iterative_artificial_spectrum(std::string dir_core){
/*
 * Run iteratively the program artificial_spectrum in order to generate a
 * serie of model according to a main configuration file (/Configurations/main.cfg).
*/

	bool passed;
	int Nmodel;
	std::vector<std::string> param_names;

	std::string cfg_file, file_out_modes, file_out_noise, file_cfg_mm, file_out_mm, file_out_combi;
	std::string external_path;
	Config_Data cfg;

	external_path=dir_core + "external/";
	cfg_file=dir_core + "Configurations/main.cfg";
	file_out_modes=dir_core + "Configurations/tmp/modes_tmp.cfg";
	file_out_noise=dir_core + "Configurations/tmp/noise_tmp.cfg";
	file_cfg_mm=dir_core + "external/ARMM-solver/star_params.global"; // Used only for models with mixed modes and interfaced with the Python solver
	file_out_combi=dir_core + "Data/Combinations.txt";

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
	if(cfg.model_name == "asymptotic_mm_v1"){
		Nmodel=6;
		param_names.push_back("Dnu"); 
		param_names.push_back("epsilon"); 
		param_names.push_back("alpha"); 
		param_names.push_back("q");
		param_names.push_back("SNR");
		param_names.push_back("maxGamma");
		if(param_names.size() != Nmodel){
			std::cout << "    Invalid number of parameters for model_name= 'asymptotic_mm_v1'" << std::endl;
			std::cout << "    Expecting " << Nmodel << " parameters, but found " << cfg.val_min.size() << std::endl;
			std::cout << "    Check your main configuration file" << std::endl;
			std::cout << "    The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}
		passed=1;
	}

	if(passed == 0){
		std::cout << "    model_name= " << cfg.model_name << " is not a recognized keyword for models" << std::endl;
		std::cout << "    Check models_database.h to see the available model names" << std::endl;
		std::cout << "    The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

	std::cout << "2. Generating the models using the subroutine " << cfg.model_name << " of model_database.cpp..." << std::endl;
	if(cfg.forest_type == "random"){
		std::cout << "   Values are randomly generated into a uniform range defined in the main configuration file" << std::endl;
		generate_random(cfg, param_names, dir_core, file_out_modes, file_out_noise, file_out_combi, Nmodel, file_cfg_mm, external_path);

	}
	if(cfg.forest_type == "grid"){
		std::cout << "   Values are generated over a grid using all possible combinations according to inputs in the main configuration file" << std::endl;
		generate_grid();
	}
	if(cfg.forest_type != "random" && cfg.forest_type != "grid"){ 
		std::cout << " Problem in the main configuration file. It is expected that the forest type parameters is either random OR grid" << std::endl;
		std::cout << " Check your main configuration file" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

	

}


void generate_random(Config_Data cfg, std::vector<std::string> param_names, std::string dir_core, std::string file_out_modes, 
		std::string file_out_noise, std::string file_out_combi, int N_model,  std::string file_cfg_mm, std::string external_path){

	bool neg=0, passed=0;
	int i, ii;
	long lastid, id0;
	std::string id_str;
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

	// Initialize the random generators: Uniform over [0,1]
    boost::mt19937 generator(time(NULL));
    boost::uniform_01<boost::mt19937> dist_gr(generator);
 
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
			std::cout << "                 The program will therefore behave as if erase_old_files=0" << std::endl;
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

		passed=call_model(cfg.model_name, input_params, file_out_modes, file_out_noise,  file_cfg_mm, dir_core, id_str, cfg, external_path);
		if(passed == 0){
			std::cout << "Warning: The function call_model did not generate any configuration file!" << std::endl;
			std::cout << "         Debug required" << std::endl;
			std::cout << "         The program will stop now" << std::endl;
			exit(EXIT_FAILURE);
		}
		id0=id0+1;
		passed=0;
	}

	std::cout << "All requested tasks executed succesfully" << std::endl;
	
}


void generate_grid(){

	std::cout << "The program does not handle yet the grid case" << std::endl;
	std::cout << "...Need to be completed" << std::endl;
	std::cout << "Use the random mode instead" << std::endl;
	std::cout << "The program will stop now" << std::endl;
	exit(EXIT_SUCCESS);
}

bool call_model(std::string model_name, VectorXd input_params, std::string file_out_modes, std::string file_out_noise, 
		std::string file_cfg_mm, std::string dir_core, std::string id_str, Config_Data cfg, std::string external_path){

	std::string str;
	bool passed=0;
	//std::string id_str;

	if(model_name == "generate_cfg_asymptotic_act_asym_Hgauss"){
		generate_cfg_asymptotic_act_asym_Hgauss(input_params, file_out_modes, file_out_noise);
		//id_str=identifier2chain(identifier);
		artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, dir_core, id_str, cfg.doplots, cfg.write_inmodel);
		passed=1;
	}
	if(model_name == "asymptotic_mm_v1"){
		asymptotic_mm_v1(input_params, file_out_modes, file_out_noise,  file_cfg_mm, external_path);
		artificial_spectrum_act_asym(cfg.Tobs, cfg.Cadence, cfg.Nspectra, dir_core, id_str, cfg.doplots, cfg.write_inmodel);
		std::cout << "   - [Warning DIRTY CODE IN LINE 281, iterative_artificial_spectrum.cpp] Hard coded path for saving python3 cfg in the spectra_info dir" << std::endl;
		std::cout << "                                                                         Stable final version should avoid this" << std::endl;
		str="cp " + file_cfg_mm + " " + dir_core + "Data/Spectra_info/" + strtrim(id_str) + ".python.global";
		const char *command = str.c_str(); 
		system(command);
		passed=1;
	}

	return passed;
}


// -----------------------------------------------------------------------
// ------------------------- SECONDARY METHODS ---------------------------
// -----------------------------------------------------------------------

VectorXi where_str(std::vector<std::string> vec, std::string value){
/*
 * Gives the indexes of values of an array that match the value
 *
*/
   int cpt;
   VectorXi index(vec.size());

	cpt=0;
	for(int i=0; i<vec.size(); i++){
		if(vec[i] == value){
			index[cpt]=i;
			cpt=cpt+1;
		}		
	}
	index.conservativeResize(cpt);
/*
	std::cout << "in where_str()" << std::endl;
	std::cout << index << std::endl;
	exit(EXIT_SUCCESS);
*/
	return index;
 
}

std::vector<double> where(std::vector<double> vec, std::string condition, double value, bool return_values){
/*
 * If return_values == 0:
 * 	Gives the indexes of values of an array that fullfil the condition
 * If return_values == 1:
 * 	Gives the values of an array that fullfil the condition
 *
 * The condition can be the operator "=", "!=", ">", "<", ">=" or "<="
*/
	std::vector<double> index, values;
	
   if(return_values == 0){
	if(condition == "="){
		for(int i=0; i<vec.size(); i++){
			if(vec[i] == value){
				index.push_back(i);
				//std::cout << "vec[" << i << "]= " << vec[i] << std::endl;
			}		
		}
	}

	if(condition == "!="){
		for(int i=0; i<vec.size(); i++){
			if(vec[i] != value){
				index.push_back(i);
			}		
		}
	}
	if(condition == ">"){
		for(int i=0; i<vec.size(); i++){
			if(vec[i] > value){
				index.push_back(i);
			}		
		}
	}
	if(condition == "<"){
		for(int i=0; i<vec.size(); i++){
			if(vec[i] < value){
				index.push_back(i);
			}		
		}
	}
	if(condition == ">="){
		for(int i=0; i<vec.size(); i++){
			if(vec[i] >= value){
				index.push_back(i);
			}		
		}
	}
	if(condition == "<="){
		for(int i=0; i<vec.size(); i++){
			if(vec[i] <= value){
				index.push_back(i);
			}		
		}
	}

	return index;
   } else {
	if(condition == "="){
		for(int i=0; i<vec.size(); i++){
			if(vec[i] == value){
				values.push_back(vec[i]);
			}		
		}
	}
	if(condition == "!="){
		for(int i=0; i<vec.size(); i++){
			if(vec[i] != value){
				values.push_back(vec[i]);
			}		
		}
	}
	if(condition == ">"){
		for(int i=0; i<vec.size(); i++){
			if(vec[i] > value){
				values.push_back(vec[i]);
			}		
		}
	}
	if(condition == "<"){
		for(int i=0; i<vec.size(); i++){
			if(vec[i] < value){
				values.push_back(vec[i]);
			}		
		}
	}
	if(condition == ">="){
		for(int i=0; i<vec.size(); i++){
			if(vec[i] >= value){
				values.push_back(vec[i]);
			}		
		}
	}
	if(condition == "<="){
		for(int i=0; i<vec.size(); i++){
			if(vec[i] <= value){
				values.push_back(vec[i]);
			}		
		}
	}

	return values;
   }

}

long read_id_allcombi(std::string file_combi){

 /* bool keyword_found=0, group_found=0;
    int cpt;
    std::string line0, char0;
    std::vector<std::string> word;
 
    std::ifstream cfg_file_session;
    cfg_file_session.open(cfg_file.c_str());
    if (cfg_file_session.is_open()) {
		cpt=0;
		std::getline(cfg_file_session, line0);
		char0=strtrim(line0.substr(0, 1));
		while(char0 != "#" && file_sesson.eof()){ // Skip lines until we reach the "#" symbol
			std::getline(cfg_file_session, line0);
			char0=strtrim(line0.substr(0, 1));
		}

    }
*/
	std::string lastline;
	std::vector<std::string> vals_last;

	lastline=read_lastline_ascii(file_combi);
	vals_last=strsplit(strtrim(lastline), " ");

	std::cout << "lastline=" << lastline << std::endl;
	for(int i=0; i<vals_last.size(); i++){
		std::cout << "vals_last[" << i << "]=" << vals_last[i] << std::endl;	
	}
	return str_to_lng(vals_last[0]);
}

std::string read_lastline_ascii(std::string filename){
/*
 * Code that jumps to last line and read it
 *
*/

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

std::string write_allcombi(MatrixXd allcombi, VectorXd cte_params, Config_Data cfg, std::string fileout, bool erase_old_file, long iter, long id0, 
		    std::vector<std::string> cte_names, std::vector<std::string> var_names, std::vector<std::string> param_names){
	
	int Nchars, precision;
	std::string id_str;
	VectorXd input_params, var_params;
	std::ofstream outfile;

	Nchars = 11;
	precision = 5;

	if(erase_old_file == 1 && iter == 0) {
		outfile.open(fileout.c_str()); // write a new file
	} else{
		outfile.open(fileout.c_str(), std::ios::app); // append
	}
	/*
	if (erase_old_file == 1) {
		outfile.open(fileout.c_str()); // write a new file
	} else {
		if(iter != 0){
			outfile.open(fileout.c_str(), std::ios::app); // append
		} else {
			outfile.open(fileout.c_str()); // write a new file
		}
	}
	*/
	if(outfile.is_open()){
		//std::cout << "File opened" << std::endl;
		if(erase_old_file == 1 && iter == 0) { // Write Header only if we do not erase the old file AND this is the first execution of the function
			outfile << "model_name= " << cfg.model_name << std::endl;
			outfile << " --------------------------" << std::endl;
			outfile << "  List of all combinations " << std::endl;
			outfile << " --------------------------" << std::endl;
			outfile << "#" << std::setw(5) << "id  ";
			for(int s=0; s<param_names.size(); s++){
				outfile << std::setw(11) << param_names[s];
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


VectorXd order_input_params(VectorXd cte_params, VectorXd var_params, std::vector<std::string> cte_names, 
		std::vector<std::string> var_names, std::vector<std::string> param_names){


	VectorXd input(cte_params.size() + var_params.size());
	VectorXd input2(cte_params.size() + var_params.size());
	std::vector<std::string> names;
	VectorXi order(cte_params.size() + var_params.size()), ind;
	
	//std::cout << "  in order_input_params..." << std::endl;
	input.segment(0, cte_params.size())=cte_params;
	input.segment(cte_params.size(), var_params.size())=var_params;

	for(int i=0; i<cte_names.size(); i++){
		names.push_back(cte_names[i]);
	}
	for(int i=0; i<var_names.size(); i++){
		names.push_back(var_names[i]);
	}
        //std::cout << "names.size()=" << names.size() << std::endl;
	for(int i=0; i<names.size(); i++){
		//std::cout << "param_names[i]=" << param_names[i] << std::endl;

		ind=where_str(names, strtrim(param_names[i])); // Assumes that only one value matches the criteria
		//std::cout << "ind =" << ind << std::endl;
		if(ind.size() == 1){
			order[i]=ind[0];
			input2[i]=input[order[i]];
	
		} else {
			if(ind.size() == 0){
				std::cout << "Some values of cfg.labels could not be matched with param_names" << std::endl;
				std::cout << "This is likely due to a mispelling" << std::endl;
				std::cout << "Debug is required. The program will stop now" << std::endl;
			} else {
				std::cout << "Problem when organizing the parameters and varialbe in the correct order" << std::endl;
				std::cout << "Multiple identical parameter names found" << std::endl;
				std::cout << "This is prohibited" << std::endl;
				std::cout << "Debug required: Keywords from the main.cfg must match keywords hardcoded in the program"  << std::endl;
				std::cout << "The program will stop now" << std::endl;
			}
			exit(EXIT_FAILURE);
		}
	}

	return input2;
}

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



int main(){

	boost::filesystem::path full_path( boost::filesystem::current_path() );

	iterative_artificial_spectrum(full_path.string() + "/");

}
