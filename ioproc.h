/*
 * ioproc.h
 *
 * Various functions handling strings and sorting
 * 
 *  Created on: 10 Oct 2017
 *      Author: obenomar
 */


/**
 * @file ioproc.h
 * @brief Header file for string handling and conversion to numbers
 *
 * Various functions handling strings and sorting
 * 
 *
 * @date 10 Oct 2017
 * @author obenomar
 */

// --- MERGED WITH STRING_HANDLER.CPP ON 17/06/2021 ---

#pragma once

#include <math.h>
#include <Eigen/Dense>
# include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>

using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;

/**
 * @brief Gives the indexes of values of an array that match the value
 * 
 * @param vec The input vector
 * @param value The value to match
 * @return std::vector<int> The indexes of values that match the value
 */
std::vector<int> where_str(const std::vector<std::string> vec, const std::string value);

/**
 * @brief Gives the indexes of values of an array that match the value
 * 
 * @param vec The input vector
 * @param value The value to match
 * @return Eigen::VectorXi The indexes of values that match the value
 */
VectorXi where_strXi(const std::vector<std::string> vec, const std::string value);

/**
 * @brief Gives the indexes or values of an array that fulfill a condition
 *
 * If return_values is false, this function returns the indexes of values in the input vector that fulfill the condition.
 * If return_values is true, this function returns the values in the input vector that fulfill the condition.
 *
 * The condition can be one of the following operators: "=", "!=", ">", "<", ">=", or "<=".
 *
 * @param vec The input vector
 * @param condition The condition to be checked
 * @param value The value to compare with
 * @param return_values Flag indicating whether to return the indexes or values (default: false)
 * @return std::vector<double> The indexes or values that fulfill the condition
 */
std::vector<double> where(const std::vector<double> vec, const std::string condition, const double value, const bool return_values);


/**
 * @brief Gives the indexes or values of an array that fulfill a condition
 *
 * If return_values is 0, this function returns the indexes of values in the input vector that fulfill the condition.
 * If return_values is 1, this function returns the values in the input vector that fulfill the condition.
 *
 * The condition can be one of the following operators: "=", "!=", ">", "<", ">=", or "<=".
 *
 * @param vec The input vector
 * @param condition The condition to be checked
 * @param value The value to compare with
 * @return std::vector<int> The indexes or values that fulfill the condition
 */
std::vector<int> where_index(const std::vector<double> vec, const std::string condition, const double value);

/**
 * @brief Gives the indexes or values of an array that fulfill a condition
 *
 * If return_values is 0, this function returns the indexes of values in the input vector that fulfill the condition.
 * If return_values is 1, this function returns the values in the input vector that fulfill the condition.
 *
 * The condition can be one of the following operators: "=", "!=", ">", "<", ">=", or "<=".
 *
 * @param vec The input Eigen vector
 * @param condition The condition to be checked
 * @param value The value to compare with
 * @return std::vector<int> The indexes or values that fulfill the condition
 */
std::vector<int> where_index(const VectorXd& vec, const std::string condition, const double value);


/**
 * @brief Gives the indexes of values of an array that match the value.
 *
 * This function returns the indexes of values in the input vector that match the specified value.
 *
 * @param vec The input vector
 * @param value The value to match
 * @return std::vector<int> The indexes of values that match the value
 */
std::vector<int> where_int(const std::vector<int> vec, const int value);

/**
 * @brief Gives the indexes of values of an array that match the value.
 *
 * @param vec The input vector
 * @param value The value to match
 * @return VectorXi The indexes of values that match the value
 */
VectorXi where_int(const VectorXi& vec, const int value);

/**
 * @brief Gives the indexes of values of an array that match the value with a tolerance.
 *
 * This function returns the indexes of values in the input vector that are within a specified tolerance of the given value.
 *
 * @param vec The input vector
 * @param value The value to match
 * @param tolerance The tolerance within which the values are considered a match
 * @return std::vector<int> The indexes of values that match the value within the tolerance
 */
std::vector<int> where_dbl(const std::vector<double> vec, const double value, const double tolerance);

/**
 * @brief Gives the indexes of values of an array that match the value within a specified tolerance.
 *
 * This function returns the indexes of values in the input vector that are within a specified tolerance of the given value.
 * The tolerance parameter allows you to control how close the match is considered as acceptable. The tolerance is in the same unit as the value.
 *
 * @param vec The input vector
 * @param value The value to match
 * @param tolerance The tolerance within which the values are considered a match
 * @return VectorXi The indexes of values that match the value within the tolerance
 */
VectorXi where_dbl(const VectorXd& vec, const double value, const double tolerance);

/**
 * @brief Gives the indexes of values of an array that match the value within a specified tolerance.
 *
 * This function returns the indexes of values in the input vector that are within a specified tolerance of the given value.
 * The tolerance parameter allows you to control how close the match is considered as acceptable. The tolerance is in the same unit as the value.
 *
 * @param vec The input vector
 * @param value The value to match
 * @param tolerance The tolerance within which the values are considered a match
 * @param imin_search The minimum index to search within the vector
 * @param imax_search The maximum index to search within the vector
 * @return VectorXi The indexes of values that match the value within the tolerance
 */
VectorXi where_dbl(const VectorXd& vec, double value, const double tolerance, const int imin_search, const int imax_search);

/**
 * @brief Reads the last line of a text file in ASCII format.
 *
 * This function opens a text file specified by the filename and reads the last line of the file.
 *
 * @param filename The name of the file to read.
 * @return The last line of the file as a string.
 */
std::string read_lastline_ascii(const std::string filename);

/**
 * @brief Sorts a vector of variable and constant values in the same order as a vector of variable and constant names.
 *
 * This function takes in a vector of constant values (cte_params) and a vector of variable values (var_params), along with vectors of constant names (cte_names), variable names (var_names), and parameter names (param_names). It sorts the values in the same order as the names and returns the sorted values in a new vector (input2).
 *
 * @param cte_params The vector of constant values
 * @param var_params The vector of variable values
 * @param cte_names The vector of constant names
 * @param var_names The vector of variable names
 * @param param_names The vector of parameter names
 * @return The sorted vector of values (input2)
 */
VectorXd order_input_params(const VectorXd& cte_params, const VectorXd& var_params, const std::vector<std::string> cte_names, 
		const std::vector<std::string> var_names, const std::vector<std::string> param_names);


/**
 * @brief Gives the indexes of values of an array within a range.
 *
 * This function takes in a vector of double values (vec), along with a minimum value (value_min) and a maximum value (value_max). It returns a vector of integers (index) containing the indexes of the values in the input vector that fall within the specified range. The strict parameter determines whether the values have to be strictly within the range or can be equal to the minimum or maximum values.
 *
 * @param vec The vector of double values
 * @param value_min The minimum value of the range
 * @param value_max The maximum value of the range
 * @param strict Flag indicating whether the values have to be strictly within the range
 * @return The vector of indexes of values within the range (index)
 */
std::vector<int> where_in_range(const std::vector<double> vec, const double value_min, const double value_max, const bool strict);

/**
 * @brief Gives the indexes of values of an array within a specified range.
 *
 * This function takes in a vector of double values (vec), along with a minimum value (value_min) and a maximum value (value_max). It returns a VectorXi containing the indexes of the values in the input vector that fall within the specified range. The strict parameter determines whether the values have to be strictly within the range or can be equal to the minimum or maximum values.
 *
 * @param vec The vector of double values
 * @param value_min The minimum value of the range
 * @param value_max The maximum value of the range
 * @param strict Flag indicating whether the values have to be strictly within the range
 * @return The VectorXi containing the indexes of values within the range
 */
VectorXi where_in_range(const VectorXd& vec, const double value_min, const double value_max, const bool strict);


/**
 * @brief Converts a string to a VectorXi array.
 *
 * This function takes in a string (str) and a string of delimiters (delimiters). It splits the string into words using the delimiters and converts each word to an integer. The resulting integers are stored in a std::vector<int> (arri) and then copied to a VectorXi (Xiarr).
 *
 * @param str The input string
 * @param delimiters The string of delimiters used to split the input string
 * @return The VectorXi array containing the converted integers
 */
VectorXi str_to_Xiarr(const std::string str, const std::string delimiters);

/**
*
* @brief Converts a string to a VectorXd array.
* This function takes in a string (str) and a string of delimiters (delimiters). It splits the string into words using the delimiters and converts each word to a double. The resulting doubles are stored in a std::vector (words) and then copied to a VectorXd (Xdarr).
* @param str The input string
* @param delimiters The string of delimiters used to split the input string
* @return The VectorXd array containing the converted doubles 
*/
VectorXd str_to_Xdarr(const std::string str, const std::string delimiters);


/**
 * @brief Converts a string to an integer.
 * 
 * This function takes in a string and converts it to an integer using std::stringstream.
 * The resulting integer is returned.
 * 
 * @param str The input string to be converted
 * @return The converted integer
 */
int str_to_int(const std::string str);

/**
* @brief Converts a string to a boolean value.
* This function takes in a string and converts it to a boolean value using std::stringstream.
* The resulting boolean value is returned.
* @param str The input string to be converted.
* @return The converted boolean value. 
*/
bool str_to_bool(const std::string str);

/**
 * @brief Converts a string to a long integer.
 *
 * This function takes in a string and converts it to a long integer using std::stringstream.
 * The resulting long integer is returned.
 *
 * @param str The input string to be converted.
 * @return The converted long integer.
 */
long str_to_lng(const std::string str);

/**
 * @brief Converts a string to a long integer.
 *
 * This function takes in a string and converts it to a long integer using stdstringstream.
 * The resulting long integer is returned.
 *
 * @param str The input string to be converted.
 * @return The converted long integer.
 */
long str_to_long(const std::string str);

/**
 * @brief Converts a string to a double.
 *
 * This function takes in a string and converts it to a double using std::stringstream.
 * The resulting double is returned.
 *
 * @param str The input string to be converted.
 * @return The converted double.
 */
double str_to_dbl(const std::string str);

/**
 * @brief Converts a string to a double array.
 *
 * This function takes in a string and a set of delimiters, and converts the string into a double array.
 * The string is split into individual words using the specified delimiters, and each word is converted to a double using std::stringstream.
 * The resulting double array is returned.
 *
 * @param str The input string to be converted.
 * @param delimiters The delimiters used to split the string into words.
 * @return The converted double array.
 */
std::vector<double> str_to_dblarr(const std::string str, const std::string delimiters);

/**
 * @brief Converts an integer to a string.
 *
 * This function takes an integer value and converts it to a string using std::ostringstream.
 * The integer value is inserted into the stream and the resulting string is returned.
 *
 * @param value The integer value to be converted.
 * @return The converted string.
 */
std::string int_to_str(const int value);

/**
 * @brief Converts a double to a string.
 *
 * This function takes a double value and converts it to a string using std::stringstream.
 * The double value is inserted into the stream and the resulting string is returned.
 *
 * @param ind The double value to be converted.
 * @return The converted string.
 */
std::string dbl_to_str(const double ind);

/**
 * @brief Converts a long integer to a string.
 *
 * This function takes a long integer value and converts it to a string using std::stringstream.
 * The long integer value is inserted into the stream and the resulting string is returned.
 *
 * @param ind The long integer value to be converted.
 * @return The converted string.
 */
std::string lng_to_str(const long ind);

/**
 * @brief Converts a string to an array of double values.
 *
 * This function takes a string and extracts an array of double values from it. The string is trimmed to remove any leading or trailing whitespace. The values in the string are separated by the specified delimiters. Each value is converted to a double using std::stringstream and stored in an array. The resulting array of double values is returned.
 *
 * @param str The string to extract the double values from.
 * @param delimiters The delimiters used to separate the double values in the string.
 * @return The array of double values extracted from the string.
 */
std::vector<double> str_to_arrdbl(const std::string str, const std::string delimiters);

/**
* @brief Converts a string to an array of integer values.
* This function takes a string and extracts an array of integer values from it. The string is trimmed to remove any leading or trailing whitespace. The values in the string are separated by the specified delimiters. Each value is converted to an integer using stdstringstream and stored in an array. The resulting array of integer values is returned.
* @param str The string to extract the integer values from.
* @param delimiters The delimiters used to separate the integer values in the string.
* @return The array of integer values extracted from the string. 
*/ 
std::vector<int> str_to_arrint(const std::string str, const std::string delimiters);

/**
 * @brief Converts an array of strings to a VectorXd array of double values.
 *
 * This function takes an array of strings and extracts a VectorXd array of double values from it. Each string value in the array is converted to a double using std::stringstream and stored in the corresponding element of the VectorXd array. The resulting VectorXd array of double values is returned.
 *
 * @param vals_strs The array of strings to extract the double values from.
 * @return The VectorXd array of double values extracted from the array of strings.
 */
VectorXd arrstr_to_Xdarrdbl(const std::vector<std::string> vals_strs);

/**
 * @brief Converts an array of strings to a VectorXi array of integer values.
 *
 * This function takes an array of strings and extracts a VectorXi array of integer values from it. Each string value in the array is converted to an integer using std::stringstream and stored in the corresponding element of the VectorXi array. The resulting VectorXi array of integer values is returned.
 *
 * @param vals_strs The array of strings to extract the integer values from.
 * @return The VectorXi array of integer values extracted from the array of strings.
 */
VectorXi arrstr_to_Xiarr(const std::vector<std::string> vals_strs);


/**
 * @brief Converts an array of strings to an array of double values.
 *
 * This function takes an array of strings and extracts an array of double values from it. Each string value in the array is converted to a double using std::stringstream and stored in the corresponding element of the double array. The resulting array of double values is returned.
 *
 * @param vals_strs The array of strings to extract the double values from.
 * @return The array of double values extracted from the array of strings.
 */
std::vector<double> arrstr_to_arrdbl(const std::vector<std::string> vals_strs);

/**
 * @brief Converts an Eigen VectorXd to a std::vector<double>.
 *
 * This function takes an Eigen VectorXd as input and converts it to a std::vector<double>.
 *
 * @param vecXd_in The input Eigen VectorXd to be converted.
 * @return The converted std::vector<double>.
 */
std::vector<double> vectXd_to_vec(const VectorXd& vecXd_in);

/**
 * @brief Checks if a file exists.
 *
 * This function checks if a file with the given name exists in the file system.
 *
 * @param name The name of the file to check.
 * @return True if the file exists, false otherwise.
 */
bool file_exists(const std::string& name);

/**
 * @brief Removes leading and trailing white spaces from a string.
 *
 * This function removes any leading and trailing white spaces from the given string. It uses the find_first_not_of and find_last_not_of functions to find the first and last non-white space characters in the string. The resulting substring, without the leading and trailing white spaces, is returned.
 *
 * @param str The string to remove leading and trailing white spaces from.
 * @return The modified string with leading and trailing white spaces removed.
 */
std::string strtrim(const std::string& str);


/**
 * @brief Splits a string into substrings using specified delimiters.
 *
 * This function takes a string and splits it into substrings each time one of the specified delimiters is detected.
 * The function first trims any leading and trailing white spaces from the input string.
 * It then iterates through the string, finding the positions of the delimiters using the find_first_of function.
 * For each delimiter found, the function extracts the substring before the delimiter, trims any leading and trailing white spaces from it, and adds it to the resulting vector of substrings.
 * After all delimiters have been processed, the remaining part of the string is also added to the vector of substrings.
 *
 * @param str The string to be split.
 * @param delimiters A string containing the delimiters to split the string by.
 * @return A vector of substrings obtained by splitting the input string.
 */
std::vector<std::string> strsplit(const std::string str, const std::string delimiters);

/**
 * @brief Filters a vector of parameters based on a range condition.
 *
 * This function takes a vector of parameters and returns only the values that satisfy the range condition fmin < f_in < fmax.
 * It first calls the where_in_range function to obtain a vector of indexes that match the range condition.
 * Then, it iterates through the indexes and adds the corresponding parameter values to the output vector.
 *
 * @param param_in The input vector of parameters to be filtered.
 * @param f_in The vector of values that define the filtering condition.
 * @param fmin The minimum valid value of f_in (lower bound for the rejection).
 * @param fmax The maximum valid value of f_in (upper bound for the rejection).
 * @return A vector of filtered parameters that satisfy the range condition.
 */
std::vector<double> filter_range(const std::vector<double> param_in, const std::vector<double> f_in, const double fmin, const double fmax);

/**
 * @brief Filters a vector of parameters based on a range condition.
 *
 * This function takes a vector of parameters and returns only the values that satisfy the range condition fmin < f_in < fmax.
 * It first calls the where_in_range function to obtain a vector of indexes that match the range condition.
 * Then, it iterates through the indexes and adds the corresponding parameter values to the output vector.
 *
 * @param param_in The input vector of parameters to be filtered.
 * @param f_in The vector of values that define the filtering condition.
 * @param fmin The minimum valid value of f_in (lower bound for the rejection).
 * @param fmax The maximum valid value of f_in (upper bound for the rejection).
 * @return A vector of filtered parameters that satisfy the range condition.
 */
std::vector<bool> filter_range(const std::vector<bool> param_in, const std::vector<double> f_in, const double fmin, const double fmax);


