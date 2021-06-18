/*
 * main.cpp
 *
 * This procedure:
 *    (1) Read a MESA summary file as defined by Kuldeep
 *    (2) Goes through a directory and check if freq files exist for each Model ID of the summary file
 *    (3) Report missing freq files in an output file
 * 
 *  Created on: 16 Oct 2017
 *      Author: obenomar
 */
# include <iostream>
# include <iomanip>
#include <fstream>
//# include <Eigen/Dense>
# include <vector>
#include <string>
# include "data.h"
# include "format.h"
# include "io_star_params.h"

bool check_freqfile(const std::string file_in_name, const bool verbose_data);

int main(){

	const bool verbose_data=0;
	bool status;
	long j;
	std::ostringstream strs;
	std::string dir_core, dir_freqs, delimiter=" ";
	std::string model_file, model_name, freq_file, file_out;
	Data_Nd models;

	dir_core="/home/obenomar/Dropbox/Temporary/Spectra-Simulator-Cpp_Wgrid/external/MESA_grid/";
	dir_freqs=dir_core + "frequencies/";
	//model_file=dir_core + "models_sample.params";
	model_file=dir_core + "models.params";

	file_out= "report.txt";

	models=read_data_ascii_Ncols(model_file, delimiter, verbose_data);
	
   	std::ofstream file_out_session;

        file_out_session.open(file_out.c_str());
        if (file_out_session.is_open()){
		j=0;
		for(long i=0; i<models.data.rows(); i++){
			// 1. Determine the model name and frequency filename
			//std::cout << std::fixed << std::setprecision(2) << models.data(i, 0) << std::endl;
	
			strs.str(std::string());
			strs << std::fixed << std::setprecision(2) << models.data(i,0);
			model_name=format_freqname(strs.str());
			freq_file=dir_freqs + "freq." + model_name;
	
			// 2. Read the file with frequencies
			status=check_freqfile(freq_file, verbose_data);
		
			// 3. Report the missing files if necessary
			//std::cout << "Stellar Model ID: " << model_name << "   " << status << std::endl;
			if(status == 0){
				file_out_session << j << "   " << model_name << std::endl;
				std::cout << j << "   " << model_name << std::endl;// Open a file and write on it the missing model name
				//std::cout << "Stellar Model ID: " << model_name << "   " << status << std::endl;
				j=j+1;
			}
		}
		file_out_session.close();
	} else{
		std::cout << " Could not open the output file " << std::endl;
		std::cout << " Filename: " << file_out << std::endl;
		std::cout << " Check that the output directory is valid " << std::endl;
		std::cout << " The program will exit now" << std::endl;
	}

}

bool check_freqfile(const std::string file_in_name, const bool verbose_data){
/*
 * This function attempt to open a file. If successful, it returns 1. Otherwise it returns 0
*/
   
    bool status;

    std::ifstream file_in;
    file_in.open(file_in_name.c_str());
    if (file_in.is_open()) {
	status=1;
	file_in.close();
     } else {
	status=0;
     }
     
     if (verbose_data == 1){ 
	std::cout << file_in_name << "    " << status << std::endl; 
     }
return status;
}


