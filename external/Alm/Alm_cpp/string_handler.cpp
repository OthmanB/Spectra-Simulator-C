/*
 * string_handler.cpp
 *
 * Contains all kind of functions
 * used to arrange/handle the strings
 * 
 *  Created on: 20 Jun 2016
 *      Author: obenomar
 */

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "string_handler.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

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
	 * Used to extract vector<double> from an array of strings (vector<string>). 
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


std::vector<double> str_to_arrdbl(const std::string str, const std::string delimiters){
	/* 
	 * Used to extract arrays of values from a string. 
	 * The terminator indicates the symbol that indicate the end of a line (beyond it is comments)
	 * The delimiter indicates how the values are separated (e.g. with a "," or with " ")
	*/
	double dbl_v;
	std::string str0;
	std::vector<std::string> vals_strs;
	std::vector<double> dbl_out;

	str0=strtrim(str); // remove blanks and convert to a stream

	vals_strs=strsplit(str0, delimiters); // get all values in a string

	for (int i=0; i<vals_strs.size(); i++){
		if(strtrim(vals_strs[i]) !=""){
			std::stringstream(vals_strs[i]) >> dbl_v;
			dbl_out.push_back(dbl_v);
		}
	}

	return dbl_out;
}

std::string dbl_to_str(const double ind){

    std::stringstream ss;

    ss.str(std::string());
    ss << ind;

    return ss.str();
}

std::string int_to_str(const int value){
	
	std::ostringstream convert;   // stream used for the conversion
	convert << value;      // insert the textual representation of 'Number' in the characters in the stream

	return convert.str(); // set 'Result' to the contents of the stream
}


double str_to_dbl(const std::string str){

	double dbl_out;

	std::stringstream(strtrim(str)) >> dbl_out;
return dbl_out;
}

int str_to_int(const std::string str){

	int int_out;

	std::stringstream(strtrim(str)) >> int_out;
return int_out;
}

long str_to_long(const std::string str){

	long long_out;

	std::stringstream(strtrim(str)) >> long_out;
return long_out;
}

bool str_to_bool(const std::string str){

	bool bool_out;

	std::stringstream(strtrim(str)) >> bool_out;
return bool_out;
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
