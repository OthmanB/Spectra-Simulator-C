/*
 * ioproc.cpp
 *
 * Various functions handling strings and sorting
 * 
 *  Created on: 10 Oct 2017
 *      Author: obenomar
 */

#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include "ioproc.h"
#ifndef _WIN32
#include <unistd.h>
#endif

using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;

std::vector<int> where_str(const std::vector<std::string> vec, const std::string value){
/*
 * Gives the indexes of values of an array that match the value
 *
*/
   int cpt;
  
   std::vector<int> index;
   index.resize(vec.size());
	
	cpt=0;
	for(int i=0; i<vec.size(); i++){
		if(vec[i] == value){
			index[cpt]=i;
			cpt=cpt+1;
		}		
	}
	index.resize(cpt);

	return index;
 
}

/*
CANNOT BE USED BECAUSE 'functions that differ only in their return type cannot be overloaded'
VectorXi where_str(const std::vector<std::string> vec, const std::string value){
//
// Gives the indexes of values of an array that match the value
//
//
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
	return index;
}
*/

VectorXi where_strXi(const std::vector<std::string> vec, const std::string value){
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
	return index;
}

std::vector<double> where(const VectorXd& vec, const std::string condition, const double value, bool return_values){
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

std::vector<int> where_index(const VectorXd& vec, const std::string condition, const double value){
/*
 * If return_values == 0:
 * 	Gives the indexes of values of an array that fullfil the condition
 * If return_values == 1:
 * 	Gives the values of an array that fullfil the condition
 *
 * The condition can be the operator "=", "!=", ">", "<", ">=" or "<="
*/
	std::vector<int> index;
	
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
}

std::vector<int> where_index(const std::vector<double> vec, const std::string condition, const double value){
/*
 * If return_values == 0:
 * 	Gives the indexes of values of an array that fullfil the condition
 * If return_values == 1:
 * 	Gives the values of an array that fullfil the condition
 *
 * The condition can be the operator "=", "!=", ">", "<", ">=" or "<="
*/
	std::vector<int> index;
	
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
}


std::vector<double> where(const std::vector<double> vec, const std::string condition, const double value, const bool return_values){
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

VectorXi where_int(const VectorXi& vec, const int value){
/*
 * Gives the indexes of values of an array that match the value.
 *
*/
   int cpt;
   VectorXi index_out;
  
   std::vector<int> index;
   index.resize(vec.size());
	
	cpt=0;
	for(int i=0; i<vec.size(); i++){
		if(vec[i] == value){
			index[cpt]=i;
			cpt=cpt+1;
		}		
	}
	if(cpt >=1){
		index_out.resize(cpt);
		for(int i=0; i<cpt; i++){
			index_out[i]=index[i];
		}
	} else{
		index_out.resize(1);
		index_out[0]=-1;
	}
	return index_out;
}

std::vector<int> where_int(const std::vector<int> vec, const int value){
/*
 * Gives the indexes of values of an array that match the value.
 *
*/
   int cpt;
  
   std::vector<int> index;
   index.resize(vec.size());
	
	cpt=0;
	for(int i=0; i<vec.size(); i++){
		if(vec[i] == value){
			index[cpt]=i;
			cpt=cpt+1;
		}		
	}
	if(cpt >=1){
		index.resize(cpt);
	} else{
		index.resize(1);
		index[0]=-1;
	}

	return index;
}

std::vector<int> where_dbl(const std::vector<double> vec, const double value, const double tolerance){
/*
 * Gives the indexes of values of an array that match the value.
 * A tolerance parameter allows you to control how close the match
 * is considered as acceptable. The tolerance is in the same unit
 * as the value
 *
*/
   int cpt;
  
   std::vector<int> index;
   index.resize(vec.size());
	
	cpt=0;
	for(int i=0; i<vec.size(); i++){
		if(vec[i] > value - tolerance && vec[i] < value + tolerance){
			index[cpt]=i;
			cpt=cpt+1;
		}		
	}
	if(cpt >=1){
		index.resize(cpt);
	} else{
		index.resize(1);
		index[0]=-1;
	}

	return index;
}

VectorXi where_dbl(const VectorXd& vec, double value, const double tolerance){
/*
 * Gives the indexes of values of an array that match the value.
 * A tolerance parameter allows you to control how close the match
 * is considered as acceptable. The tolerance is in the same unit
 * as the value
 *
*/
   int cpt;
   VectorXi index_out;
  
   //std::vector<int> index;
   index_out.resize(vec.size());
	
	cpt=0;
	for(int i=0; i<vec.size(); i++){
		if(vec[i] > value - tolerance && vec[i] < value + tolerance){
			index_out[cpt]=i;
			cpt=cpt+1;
		}		
	}
	if(cpt >=1){
		index_out.conservativeResize(cpt);
		//for(int i=0; i<cpt; i++){
		//	index_out[i]=index[i];
		//}
	} else{
		index_out.resize(1);
		index_out[0]=-1;
	}
	return index_out;
}

VectorXi where_dbl(const VectorXd& vec, double value, const double tolerance, const int imin_search, const int imax_search){
/*
 * Gives the indexes of values of an array that match the value.
 * A tolerance parameter allows you to control how close the match
 * is considered as acceptable. The tolerance is in the same unit
 * as the value
 *
*/
   int cpt;
   VectorXi index_out;
  
   //std::vector<int> index;
   index_out.resize(vec.size());
	
	cpt=0;
	for(int i=imin_search; i<=imax_search; i++){
		if(vec[i] > value - tolerance && vec[i] < value + tolerance){
			index_out[cpt]=i;
			cpt=cpt+1;
		}		
	}
	if(cpt >=1){
		index_out.conservativeResize(cpt);
		//for(int i=0; i<cpt; i++){
		//	index_out[i]=index[i];
		//}
	} else{
		index_out.resize(1);
		index_out[0]=-1;
	}
	return index_out;
}

VectorXi where_in_range(const VectorXd& vec, const double value_min, const double value_max, const bool strict){
/*
 * Gives the indexes of values of an array within a specified range 
 *
*/
   int cpt;
   VectorXi index(vec.size());
  
   //std::vector<int> index;
   //index.resize(vec.size());
	
	cpt=0;
	for(int i=0; i<vec.size(); i++){
		if (strict == 1){
			if(vec[i] > value_min && vec[i] < value_max){
				index[cpt]=i;
				cpt=cpt+1;
			}
		} else{
			if(vec[i] >= value_min && vec[i] <= value_max){
				index[cpt]=i;
				cpt=cpt+1;
			}		
		}		
	}
	if(cpt >=1){
		index.conservativeResize(cpt);
	} else{
		index.resize(1);
		index[0]=-1;
	}
	return index;
}

std::vector<int> where_in_range(const std::vector<double> vec, const double value_min, const double value_max, const bool strict){
/*
 * Gives the indexes of values of an array within a range.
 * If strict = 1 values have to be strictly within the range. 
 *
*/
   int cpt;
  
   std::vector<int> index;
   index.resize(vec.size());
	
	cpt=0;
	for(int i=0; i<vec.size(); i++){
		if (strict == 1){
			if(vec[i] > value_min && vec[i] < value_max){
				index[cpt]=i;
				cpt=cpt+1;
			}
		} else{
			if(vec[i] >= value_min && vec[i] <= value_max){
				index[cpt]=i;
				cpt=cpt+1;
			}
		}	
	}
	if(cpt >=1){
		index.resize(cpt);
	} else{
		index.resize(1);
		index[0]=-1;
	}

	return index;
}


VectorXi str_to_Xiarr(const std::string str, const std::string delimiters){

	int int_v;
	std::vector<std::string> words;
	std::vector<int> arri;
	VectorXi Xiarr;

	words=strsplit(str, delimiters);
	for(int i=0; i<words.size(); i++){
		if(strtrim(words[i]) != ""){
			std::stringstream(words[i]) >> int_v;
			arri.push_back(int_v);
		}
	}
	Xiarr.resize(arri.size());
	for(int i=0; i<arri.size(); i++){
		Xiarr[i]=arri[i];
	}

return Xiarr;
}


// Version of str_to_Xdarr imported from Diagnotics.cpp
VectorXd str_to_Xdarr(const std::string str, const std::string delimiters){

	std::vector<std::string> words;
	Eigen::VectorXd Xdarr;

	words=strsplit(str, delimiters);
	Xdarr.resize(words.size());
	for(int i=0; i<words.size(); i++){
		if(strtrim(words[i]) != ""){ // Added on 1 May 2018 ... check the impact
			std::stringstream(words[i]) >> Xdarr[i];
		}
	}

return Xdarr;
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

/* 
 * Function that sort a vector of variable and constant in the same
 * order as a vector of variable names and constant names
 *
*/
VectorXd order_input_params(const VectorXd& cte_params, const VectorXd& var_params, const std::vector<std::string> cte_names, 
		const std::vector<std::string> var_names, const std::vector<std::string> param_names){


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

		ind=where_strXi(names, strtrim(param_names[i])); // Assumes that only one value matches the criteria
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

int str_to_int(const std::string str){

	int int_out;
	std::stringstream(strtrim(str)) >> int_out;
return int_out;
}


bool str_to_bool(const std::string str){

	bool bool_out;

	std::stringstream(strtrim(str)) >> bool_out;
return bool_out;
}


long str_to_long(const std::string str){
	return str_to_lng(str);
}

long str_to_lng(const std::string str){

	long lng_out;

	std::stringstream(strtrim(str)) >> lng_out;
return lng_out;
}

double str_to_dbl(const std::string str){

	//std::stringstream d0;
	double dbl_out;

	//d0 <<  std::setprecision(9) << str; strtrim(str);
	//d0 >> dbl_out;
	std::stringstream(strtrim(str)) >> dbl_out;
	//std::cout << "!!!! " << dbl_out;
return dbl_out;
}

std::vector<int> str_to_arrint(const std::string str, const std::string delimiters){
	/* 
	 * Used to extract arrays of values from a string. 
	 * The terminator indicates the symbol that indicate the end of a line (beyond it is comments)
	 * The delimiter indicates how the values are separated (e.g. with a "," or with " ")
	*/
	int int_v;
	std::string str0;
	std::vector<std::string> vals_strs;
	std::vector<int> int_out;

	str0=strtrim(str); // remove blanks and convert to a stream

	vals_strs=strsplit(str0, delimiters); // get all values in a string

	for (int i=0; i<vals_strs.size(); i++){
		if(strtrim(vals_strs[i]) != ""){
				std::stringstream(vals_strs[i]) >> int_v;
				int_out.push_back(int_v);
		}
	}


	return int_out;
}

std::vector<double> str_to_arrdbl(const std::string str, const std::string delimiters){
	/* 
	 * Used to extract arrays of values from a string. 
	 * The terminator indicates the symbol that indicate the end of a line (beyond it is comments)
	 * The delimiter indicates how the values are separated (e.g. with a "," or with " ")
	*/
	std::string str0;
	std::vector<std::string> vals_strs;
	std::vector<double> dbl_out;

	str0=strtrim(str); // remove blanks and convert to a stream

	vals_strs=strsplit(str0, delimiters); // get all values in a string
	dbl_out.resize(vals_strs.size());
	for (int i=0; i<vals_strs.size(); i++){
		std::stringstream(vals_strs[i]) >> dbl_out[i];
	}

	return dbl_out;
}

std::vector<double> str_to_dblarr(const std::string str, const std::string delimiters){

	std::vector<std::string> words;
	std::vector<double> dblarr;

	words=strsplit(str, delimiters);
	dblarr.resize(words.size());
	for(int i=0; i<words.size(); i++){
		std::stringstream(strtrim(words[i])) >> dblarr[i];
	}

return dblarr;
}


std::vector<double> arrstr_to_arrdbl(const std::vector<std::string> vals_strs){
	/* 
	 * Used to extract arrays of values from a string. 
	 * The terminator indicates the symbol that indicate the end of a line (beyond it is comments)
	 * The delimiter indicates how the values are separated (e.g. with a "," or with " ")
	*/

	std::vector<double> dbl_out;

	dbl_out.resize(vals_strs.size());
	for (int i=0; i<vals_strs.size(); i++){
		std::stringstream(vals_strs[i]) >> dbl_out[i];
	}

	return dbl_out;
}

VectorXd arrstr_to_Xdarrdbl(const std::vector<std::string> vals_strs){
	/* 
	 * Used to extract VectorXd from an array of strings (vector<string>). 
	*/

	VectorXd dbl_out;

	dbl_out.resize(vals_strs.size());
	for (int i=0; i<vals_strs.size(); i++){
		std::stringstream(vals_strs[i]) >> dbl_out[i];
	}

	return dbl_out;
}


VectorXi arrstr_to_Xiarr(const std::vector<std::string> vals_strs){
	/* 
	 * Used to extract VectorXd from an array of strings (vector<string>). 
	*/

	VectorXi int_out;

	int_out.resize(vals_strs.size());
	for (int i=0; i<vals_strs.size(); i++){
		std::stringstream(vals_strs[i]) >> int_out[i];
	}

	return int_out;
}


bool file_exists(const std::string& name) {
    return ( access( name.c_str(), F_OK ) != -1 );
}


std::string strtrim(const std::string& str){
/*
 * Small function that remove white space at the end and at
 * the begining of a string. 
 * The original program was taken from http://stackoverflow.com/questions/1798112/removing-leading-and-trailing-spaces-from-a-string
 * Modified in order to not use the C++11 standard:
 *  	- const auto --> const string
 *      - optional argument for the separator is now hardcoded
*/
    std::string whitespace = " \t";
    size_t strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    size_t strEnd = str.find_last_not_of(whitespace);
    size_t strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

/*
// WARNING ON 17/06/2021: THIS FUNCTION WAS MOVED HERE FROM STRING_HANDLER.CPP
//                        WHEN TRYING TO FUSION IT WITH IOPROC.CPP
//                        IT SLIGHTLY DIFFERS FROM THE ABOVE STRSPLIT() AND THE EFFECT IS YET TO ASSES
std::vector<std::string> strsplit(const std::string str, const std::string delimiters){
//
// Take a string and split it each time one of the listed delimiters is detected
//
	std::string str0=strtrim(str);
	size_t pos=0;
	std::vector<std::string> str_splitted;

	while ((pos = str0.find_first_of(delimiters)) != std::string::npos) {
		    
		str_splitted.push_back(str0.substr(0, pos)); // get the substring
		str0.erase(0, pos + delimiters.length());
		str0=strtrim(str0); // remove any extra white space at the begining before iterating
	}
	str_splitted.push_back(str0); // do not forget to add the end of the string

return str_splitted;
}
*/

std::vector<std::string> strsplit(const std::string str, const std::string delimiters){
//
// Take a string and split it each time one of the listed delimiters is detected
//
	std::string str0=strtrim(str);
	size_t pos=0;
	std::vector<std::string> str_splitted;

	while ((pos = str0.find_first_of(delimiters)) != std::string::npos) {
		if (pos !=0){
			str_splitted.push_back(strtrim(str0.substr(0, pos))); // get the substring
		}
			str0.erase(0, pos +1);
			str0=strtrim(str0); // remove any extra white space at the begining before iterating
	}
	if (pos !=0){
		str_splitted.push_back(str0); // do not forget to add the end of the string
	}
return str_splitted;
}

std::vector<double> vectXd_to_vec(const VectorXd& vecXd_in){

   std::vector<double> vec;
   vec.resize(vecXd_in.size());
   VectorXd::Map(&vec[0], vecXd_in.size()) = vecXd_in;
   return vec;
}

std::string int_to_str(const int value){
	
	std::ostringstream convert;   // stream used for the conversion
	convert << value;      // insert the textual representation of 'Number' in the characters in the stream

	return convert.str(); // set 'Result' to the contents of the stream
}


std::string dbl_to_str(const double ind){

    std::stringstream ss;

    ss.str(std::string());
    ss << ind;

    return ss.str();
}

std::string lng_to_str(const long ind){

    std::stringstream ss;

    ss.str(std::string());
    ss << ind;

    return ss.str();
}


std::vector<double> filter_range(const std::vector<double> param_in, const std::vector<double> f_in, const double fmin, const double fmax){
/* This function takes a vector of parameters (param_in) and returns only its values at
   The indexes that matches the range condition fmin < f_in < fmax
   Inputs:
   	param_in: input vector to be filtered
   	f_in: vector that contains the values that define the filtering condition
   	fmin: minimum valid value of f_in (lower bound for the rejection)
   	fmax: maximum valid value of f_in (upper bound for the rejection)
*/
	std::vector<int> pos_OK;
	std::vector<double> param_out;
	
	pos_OK=where_in_range(f_in, fmin, fmax, 1);
		
	for (int i=0; i<pos_OK.size(); i++){
		if (pos_OK[i] != -1){
			param_out.push_back(param_in[pos_OK[i]]);
		}
	}
	return param_out;
}

std::vector<bool> filter_range(const std::vector<bool> param_in, const std::vector<double> f_in, const double fmin, const double fmax){
/* This function takes a vector of parameters (param_in) and returns only its values at
   The indexes that matches the range condition fmin < f_in < fmax
   Inputs:
   	param_in: input vector to be filtered
   	f_in: vector that contains the values that define the filtering condition
   	fmin: minimum valid value of f_in (lower bound for the rejection)
   	fmax: maximum valid value of f_in (upper bound for the rejection)
*/
	std::vector<int> pos_OK;
	std::vector<bool> param_out;
	
	pos_OK=where_in_range(f_in, fmin, fmax, 1);
		
	for (int i=0; i<pos_OK.size(); i++){
		if (pos_OK[i] != -1){
			param_out.push_back(param_in[pos_OK[i]]);
		}
	}
	return param_out;
}
