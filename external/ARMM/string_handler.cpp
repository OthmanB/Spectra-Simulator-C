/**
 * @file string_handler.cpp
 * @brief Header file containing utility functions for string manipulation.
 * @date 20 Jun 2016
 * @author obenomar
 * This file contains functions for trimming strings, splitting strings, and finding indexes of values in arrays.
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

/**
 * @brief Removes leading and trailing white space from a string.
 * @param str The input string.
 * @return The trimmed string.
 *
 * This function removes any leading and trailing white space characters from the input string.
 */
std::string strtrim(const std::string& str){
    std::string whitespace = " \t";
    size_t strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    size_t strEnd = str.find_last_not_of(whitespace);
    size_t strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}


/**
 * @brief Splits a string into substrings based on specified delimiters.
 * @param str The input string.
 * @param delimiters The delimiters used for splitting.
 * @return A vector of substrings.
 *
 * This function splits the input string into substrings whenever one of the specified delimiters is found.
 */
std::vector<std::string> strsplit(const std::string str, const std::string delimiters){
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

/**
 * @brief Finds the indexes of values in an array that match a specified value.
 * @param vec The input array.
 * @param value The value to match.
 * @return A vector of indexes.
 *
 * This function returns a vector containing the indexes of values in the input array that match the specified value.
 */
std::vector<int> where_str(const std::vector<std::string> vec, const std::string value){
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

/**
 * @brief Finds the indexes of values in an array that match a specified value within a tolerance.
 * @param vec The input array.
 * @param value The value to match.
 * @param tolerance The tolerance for matching.
 * @return A vector of indexes.
 *
 * This function returns a vector containing the indexes of values in the input array that match the specified value within the specified tolerance.
 */
std::vector<int> where_dbl(const std::vector<double> vec, const double value, const double tolerance){
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

/**
 * @brief Gives the indexes of values of an array that match the value within a tolerance.
 * 
 * @param vec The input array.
 * @param value The value to match.
 * @param tolerance The tolerance for matching.
 * @return VectorXi The vector of indexes.
 * 
 * This function returns a VectorXi containing the indexes of values in the input array that match the specified value within the specified tolerance.
 */
VectorXi where_dbl(const VectorXd& vec, double value, const double tolerance){
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

/**
 * @brief Gives the indexes of values of an array that match the value within a tolerance, within a specified range.
 * 
 * @param vec The input array.
 * @param value The value to match.
 * @param tolerance The tolerance for matching.
 * @param imin_search The minimum index to search.
 * @param imax_search The maximum index to search.
 * @return VectorXi The vector of indexes.
 * 
 * This function returns a VectorXi containing the indexes of values in the input array that match the specified value within the specified tolerance, within the specified range.
 */
VectorXi where_dbl(const VectorXd& vec, double value, const double tolerance, const int imin_search, const int imax_search){

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

/**
 * @brief Gives the indexes of values of an array within a specified range.
 * 
 * @param vec The input array.
 * @param value_min The minimum value of the range.
 * @param value_max The maximum value of the range.
 * @param strict A flag indicating whether the range is strict or inclusive.
 * @return VectorXi The vector of indexes.
 * 
 * This function returns a VectorXi containing the indexes of values in the input array that fall within the specified range. The strict flag determines whether the range is strict (exclusive) or inclusive.
 */
VectorXi where_in_range(const VectorXd& vec, const double value_min, const double value_max, const bool strict){

   int cpt;
   VectorXi index(vec.size());
  
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

/**
 * @brief Gives the indexes of values of an array within a specified range.
 * 
 * @param vec The input array.
 * @param value_min The minimum value of the range.
 * @param value_max The maximum value of the range.
 * @param strict A flag indicating whether the range is strict or inclusive.
 * @return std::vector<int> The vector of indexes.
 * 
 * This function returns a std::vector<int> containing the indexes of values in the input array that fall within the specified range. The strict flag determines whether the range is strict (exclusive) or inclusive.
 */
std::vector<int> where_in_range(const std::vector<double> vec, const double value_min, const double value_max, const bool strict){
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

/**
 * @brief Gives the indexes of values of an array that match the value.
 *
 * @param vec The input array.
 * @param value The value to match.
 * @return VectorXi The vector of indexes.
 *
 * This function returns a VectorXi containing the indexes of values in the input array that match the specified value.
 */
VectorXi where_int(const VectorXi& vec, const int value){
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

/**
 * @brief Gives the indexes of values of an array that match the value.
 *
 * @param vec The input array.
 * @param value The value to match.
 * @return std::vector<int> The vector of indexes.
 *
 * This function returns a std::vector<int> containing the indexes of values in the input array that match the specified value.
 */
std::vector<int> where_int(const std::vector<int> vec, const int value){
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


/**
 * @brief Converts a string to a VectorXi array.
 *
 * @param str The input string.
 * @param delimiters The delimiters used for splitting.
 * @return VectorXi The converted VectorXi array.
 *
 * This function converts a string to a VectorXi array by splitting the string using the specified delimiters and converting each substring to an integer.
 */
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


/**
 * @brief Converts a string to a VectorXd array.
 *
 * @param str The input string.
 * @param delimiters The delimiters used for splitting.
 * @return VectorXd The converted VectorXd array.
 *
 * This function converts a string to a VectorXd array by splitting the string using the specified delimiters and converting each substring to a double.
 */
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

/**
 * @brief Converts a string to a std::vector<double> array.
 *
 * @param str The input string.
 * @param delimiters The delimiters used for splitting.
 * @return std::vector<double> The converted std::vector<double> array.
 *
 * This function converts a string to a std::vector<double> array by splitting the string using the specified delimiters and converting each substring to a double.
 */
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

/**
 * @brief Converts an array of strings to an array of doubles.
 *
 * @param vals_strs The array of strings.
 * @return std::vector<double> The converted array of doubles.
 *
 * This function converts an array of strings to an array of doubles by converting each string to a double.
 */
std::vector<double> arrstr_to_arrdbl(const std::vector<std::string> vals_strs){

	std::vector<double> dbl_out;

	dbl_out.resize(vals_strs.size());
	for (int i=0; i<vals_strs.size(); i++){
		std::stringstream(vals_strs[i]) >> dbl_out[i];
	}

	return dbl_out;
}

/**
 * @brief Converts an array of strings to a VectorXd array of doubles.
 *
 * @param vals_strs The array of strings.
 * @return VectorXd The converted VectorXd array of doubles.
 *
 * This function converts an array of strings to a VectorXd array of doubles by converting each string to a double.
 */
VectorXd arrstr_to_Xdarrdbl(const std::vector<std::string> vals_strs){

	VectorXd dbl_out;

	dbl_out.resize(vals_strs.size());
	for (int i=0; i<vals_strs.size(); i++){
		std::stringstream(vals_strs[i]) >> dbl_out[i];
	}

	return dbl_out;
}

/**
 * @brief Extracts a VectorXi from an array of strings.
 *
 * @param vals_strs The array of strings.
 * @return VectorXi The extracted VectorXi.
 *
 * This function extracts a VectorXi from an array of strings by converting each string to an integer.
 */
VectorXi arrstr_to_Xiarr(const std::vector<std::string> vals_strs){
	VectorXi int_out;

	int_out.resize(vals_strs.size());
	for (int i=0; i<vals_strs.size(); i++){
		std::stringstream(vals_strs[i]) >> int_out[i];
	}

	return int_out;
}

/**
 * @brief Extracts an array of doubles from a string.
 *
 * @param str The input string.
 * @param delimiters The delimiters used for splitting.
 * @return std::vector<double> The extracted array of doubles.
 *
 * This function extracts an array of doubles from a string by splitting the string using the specified delimiters and converting each substring to a double.
 * The terminator indicates the symbol that indicates the end of a line (beyond it is comments).
 * The delimiter indicates how the values are separated (e.g. with a "," or with " ").
 */
std::vector<double> str_to_arrdbl(const std::string str, const std::string delimiters){
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

/**
 * @brief Converts a double to a string.
 *
 * @param ind The input double.
 * @return std::string The converted string.
 *
 * This function converts a double to a string using a stringstream.
 */
std::string dbl_to_str(const double ind){

    std::stringstream ss;

    ss.str(std::string());
    ss << ind;

    return ss.str();
}

/**
 * @brief Converts an integer to a string.
 *
 * @param value The input integer.
 * @return std::string The converted string.
 *
 * This function converts an integer to a string using a stringstream.
 */
std::string int_to_str(const int value){
	
	std::ostringstream convert;   // stream used for the conversion
	convert << value;      // insert the textual representation of 'Number' in the characters in the stream

	return convert.str(); // set 'Result' to the contents of the stream
}

/**
 * @brief Converts a string to a double.
 *
 * @param str The input string.
 * @return double The converted double.
 *
 * This function converts a string to a double using a stringstream.
 */
double str_to_dbl(const std::string str){

	double dbl_out;

	std::stringstream(strtrim(str)) >> dbl_out;
return dbl_out;
}

/**
 * @brief Converts a string to an integer.
 *
 * @param str The input string.
 * @return int The converted integer.
 *
 * This function converts a string to an integer using a stringstream.
 */
int str_to_int(const std::string str){

	int int_out;

	std::stringstream(strtrim(str)) >> int_out;
return int_out;
}

/**
 * @brief Converts a string to a long.
 *
 * @param str The input string.
 * @return long The converted long.
 *
 * This function converts a string to a long using a stringstream.
 */
long str_to_long(const std::string str){

	long long_out;

	std::stringstream(strtrim(str)) >> long_out;
return long_out;
}

/**
 * @brief Converts a string to a boolean.
 *
 * @param str The input string.
 * @return bool The converted boolean.
 *
 * This function converts a string to a boolean using a stringstream.
 */
bool str_to_bool(const std::string str){

	bool bool_out;

	std::stringstream(strtrim(str)) >> bool_out;
return bool_out;
}

/**
 * @brief Extracts an array of integers from a string.
 *
 * This function is used to extract arrays of values from a string. The terminator indicates the symbol that indicates the end of a line (beyond it is comments). The delimiter indicates how the values are separated (e.g. with a "," or with " ").
 *
 * @param str The input string.
 * @param delimiters The delimiters used for splitting.
 * @return std::vector<int> The extracted array of integers.
 */
std::vector<int> str_to_arrint(const std::string str, const std::string delimiters){

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

/**
 * @brief Filters a vector of parameters based on a range condition.
 *
 * This function takes a vector of parameters and returns only its values at the indexes that match the range condition fmin < f_in < fmax.
 *
 * @param param_in The input vector to be filtered.
 * @param f_in The vector that contains the values that define the filtering condition.
 * @param fmin The minimum valid value of f_in (lower bound for the rejection).
 * @param fmax The maximum valid value of f_in (upper bound for the rejection).
 * @return std::vector<double> The filtered vector of parameters.
 */
std::vector<double> filter_range(const std::vector<double> param_in, const std::vector<double> f_in, const double fmin, const double fmax){

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

/**
 * @brief Filters a vector of parameters based on a range condition.
 *
 * This function takes a vector of parameters and returns only its values at the indexes that match the range condition fmin < f_in < fmax.
 *
 * @param param_in The input vector to be filtered.
 * @param f_in The vector that contains the values that define the filtering condition.
 * @param fmin The minimum valid value of f_in (lower bound for the rejection).
 * @param fmax The maximum valid value of f_in (upper bound for the rejection).
 * @return std::vector<bool> The filtered vector of parameters.
 */
std::vector<bool> filter_range(const std::vector<bool> param_in, const std::vector<double> f_in, const double fmin, const double fmax){

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

