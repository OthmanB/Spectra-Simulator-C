/*
 * combi.cpp
 *
 *  Function that handle combination generation
 *  and read/write into files
 * 
 *  Created on: 10 Oct 2017
 *      Author: obenomar
 */
/**
 * @file combi.cpp
 * @brief Function that handles combination generation and read/write into files
 * 
 * This file contains functions for generating combinations and reading/writing them into files.
 * 
 * @date 10 Oct 2017
 * @author obenomar
 */ 

# include <iostream>
# include <iomanip>
#include <fstream>
# include <Eigen/Dense>
# include <vector>
#include <string>
# include "data.h"
# include "ioproc.h"
# include "format.h"

using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;

/*
* Function that generates all possible combinations recursively
*/
void generate_combinations(const MatrixXd& Values, const VectorXi& Nvalues, const int Nparams, MatrixXd& allcombi, long& z, VectorXd& current_combi, int current_param) {
    if (current_param == Nparams) {
        allcombi.row(z) = current_combi;
        z++;
        return;
    }
    
    for (int i = 0; i < Nvalues[current_param]; i++) {
        current_combi(current_param) = Values(i, current_param);
        generate_combinations(Values, Nvalues, Nparams, allcombi, z, current_combi, current_param + 1);
    }
}

/*
* Function that generates all possible combinations
*/
MatrixXd define_all_combinations(const MatrixXd& Values, const VectorXi& Nvalues, const int Nparams) {
    long Ncombi = Nvalues.prod();
    MatrixXd allcombi(Ncombi, Nparams);
    long z = 0;
    VectorXd current_combi(Nparams);
  
    generate_combinations(Values, Nvalues, Nparams, allcombi, z, current_combi, 0);
    
    return allcombi;
}


long read_id_allcombi(std::string file_combi){

	std::string lastline;
	std::vector<std::string> vals_last;

	lastline=read_lastline_ascii(file_combi);
	vals_last=strsplit(strtrim(lastline), " ");

	std::cout << "lastline=" << lastline << std::endl;
	//for(int i=0; i<vals_last.size(); i++){
	//	std::cout << "vals_last[" << i << "]=" << vals_last[i] << std::endl;	
	//}
	return str_to_lng(vals_last[0]);
}

std::string write_allcombi(MatrixXd& allcombi, VectorXd& cte_params, Config_Data cfg, std::string fileout, bool erase_old_file, long iter, long id0, 
		    std::vector<std::string> cte_names, std::vector<std::string> var_names, std::vector<std::string> param_names){
	
	int Nchars, precision;
	std::string id_str;
	VectorXd input_params, var_params;
	std::ofstream outfile;

	Nchars = 17;
	precision = 5;

	if(erase_old_file == 1 && iter == 0) {
		outfile.open(fileout.c_str()); // write a new file
	} else{
		outfile.open(fileout.c_str(), std::ios::app); // append
	}
	if(outfile.is_open()){
		if(erase_old_file == 1 && iter == 0) { // Write Header only if we do not erase the old file AND this is the first execution of the function
			outfile << "model_name= " << cfg.model_name << std::endl;
			outfile << " --------------------------" << std::endl;
			outfile << "  List of all combinations " << std::endl;
			outfile << " --------------------------" << std::endl;
			outfile << "#" << std::setw(7) << "id  ";
			for(int s=0; s<param_names.size(); s++){
				outfile << std::setw(Nchars) << param_names[s];
			}
			outfile << std::endl;
		} 
		for(int i=0; i<allcombi.rows(); i++){
			id_str=identifier2chain(i + id0); // The identifier corresponds to the index of the current process + the initial id0
			outfile << std::setw(7) << id_str;			
			var_params=allcombi.row(i).transpose();
			input_params=order_input_params(cte_params, var_params, cte_names, var_names, param_names);
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


