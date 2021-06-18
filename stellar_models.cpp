/*
 * stellar_models.cpp
 *
 *  High level functions that handle read/write
 *  of files from stellar models (ADIPLS, MESA)
 * 
 *  Created on: 10 Oct 2017
 *      Author: obenomar
 */

#include <math.h>
#include <Eigen/Dense>
# include <iostream>
# include <iomanip>
# include <vector>
# include "stellar_models.h"
# include "data.h"
# include "ioproc.h" // contains the string handlers
# include "format.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;


/* 
* Function that extract the relevant information from the table of models
* given as input of the function. 
* Condition: 
*	(1) Data are arranged by rows so that the function extracts a given 
*	    requested row.
*	(2) The first columns contains the model index so that the program can find
*	    the relevant position for the model
*	(3) The frequency files should have the following format: freq.[Model_index]
*           Example: freq.000101.24
* Inputs: 
*	(1) table: The input table that contains information from the models (Model_index, M, R, Fe/H, numax, Dnu, ...)
*	(2) table_index: column index that is requested to extract
*	(3) dir_freqs: director that contains all of the frequencies files for each model
* Output: A structure of type Model_data (see data.h for further information)
*/
Model_data get_model_param(const Data_Nd table, const long table_index, const std::string dir_freqs){

	bool verbose_data;
	std::string model_name, freq_file;
	std::ostringstream strs;
	Model_data model_param;	
	Data_Nd model_freqs;
	
	model_param.params.resize(table.data.cols());	

	model_param.params=table.data.row(table_index);
	for(int i=0; i<table.labels.size();i++){
		model_param.labels_params.push_back(table.labels[i]);
	}
	
	// ------- Dealing with frequencies -------
	// 1. Determine the model name and frequency filename
	strs << std::fixed << std::setprecision(2) << model_param.params[0];
	model_name=format_freqname(strs.str());
	freq_file=dir_freqs + "freq." + model_name;
	
	//std::cout << std::fixed << model_param.params[0] << std::endl;
	std::cout << "Stellar Model ID: " << model_name << std::endl;

	// 2. Read the file with frequencies
	verbose_data=0;
	model_freqs=read_freq(freq_file, verbose_data);

	// 3. export table of frequencies into model_param (Model_data  structure)
	model_param.freqs.resize(model_freqs.data.rows(), model_freqs.data.cols());
	model_param.freqs=model_freqs.data;
	for(int i=0; i<model_freqs.labels.size();i++){
		model_param.labels_freqs.push_back(model_freqs.labels[i]);
	}

return model_param;
}

Data_Nd read_freq(const std::string file_in_name, const bool verbose_data){
/*
 * This function is specifically to read frequencies outputs from adipls.
 * It is an adaptation of read_data_ascii_Ncols(const std::string file_in_name, const std::string delimiter, const bool verbose_data)
 * There was a need for an adaptation because the header/comments in freq files are on a side column, with frequencies (not on top)
 * Conditions for use: 
 *        - Read 3 first columns of a matrix-format file (extra columns are ignored)
 *        - The separator must be an empty space (at least)
 *        - Can handle up to 10000 lines (ie frequencies)
*/
    
    const std::string delimiter=" ";
    const int Ncols=3;

    int cpt, Nrows;
    std::string line, line0, subline0; //token
    std::vector<std::string> data_str;
    long double tmp_val;
    MatrixXd data;
    Data_Nd all_data_out; // The structure that encapsulate all the data, the header, the labels and units
    int data_Maxsize=10000;
 
    std::ifstream file_in;
    if (verbose_data == 1) {
 	   std::cout << "  Assumptions for the data file: " << std::endl;
	   std::cout << "       - No top header or extra lines" << std::endl;
 	   std::cout << "       - Read 3 first columns of a matrix-format file (extra columns are ignored)" <<std::endl;
	   std::cout << "       - The separator must be an empty space (at least)" <<std::endl;
 	   std::cout << "       - Maximum number of lines for the data: " << data_Maxsize << std::endl;
	}
    file_in.open(file_in_name.c_str());
    if (file_in.is_open()) {
	if (verbose_data == 1) {std::cout << "Data File opened... processing lines" << std::endl;}

	   //  Read the data...
	   if (verbose_data == 1) {std::cout <<  "   [4] Now processing the data..." << std::endl;}

	   Nrows=0;
	   std::getline(file_in, line0);
	   while(!file_in.eof()){
		//data_str=strsplit(strtrim(line0), " \t");
		data_str=strsplit(strtrim(line0), " ");  
		if (Nrows == 0) {
			data.resize(data_Maxsize, Ncols);
			data.setConstant(-2);
		}
		for(int i=0; i<Ncols;i++){
			if ( ! (std::istringstream(data_str[i]) >> tmp_val) ){tmp_val = nan("");} // If the number can be converted, then tmp_val=value. Otherwise tmp_val = NaN
			data(Nrows, i)=tmp_val;		
		}
		if (verbose_data == 1) {std::cout << data.row(Nrows) << std::endl;} // Show all entries only if requested
		std::getline(file_in, line0);
		Nrows=Nrows+1;
	    }
	file_in.close();
	data.conservativeResize(Nrows, data.cols());
	if (verbose_data == 1) { std::cout << "         - Number of lines found: " << Nrows << std::endl; }
     } else {
	std::cout << "Could not open the data file!" << std::endl;
	std::cout << "file: " << file_in_name << std::endl;
	std::cout << "use the full path of the file in order to read it without issue " << std::endl;
	exit(EXIT_FAILURE);
     }

     all_data_out.data=data;
     all_data_out.header.push_back("Frequencies computed using ADIPLS outputs. Details given in" + file_in_name);
     all_data_out.labels.push_back("degree");
     all_data_out.labels.push_back("radial order");
     all_data_out.labels.push_back("frequency");
     //all_data_out.labels.push_back("Inertia (?)");

     all_data_out.units.push_back("");
     all_data_out.units.push_back("");
     all_data_out.units.push_back("muHz");
     //all_data_out.units.push_back("");

return all_data_out;
}


