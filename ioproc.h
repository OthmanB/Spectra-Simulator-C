/*
 * ioproc.h
 *
 * Various functions handling strings and sorting
 * 
 *  Created on: 10 Oct 2017
 *      Author: obenomar
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

std::vector<int> where_str(const std::vector<std::string> vec, const std::string value);
//VectorXi where_str(const std::vector<std::string> vec, const std::string value); // Cannot be used here because 'functions that differ only in their return type cannot be overloaded'
VectorXi where_strXi(const std::vector<std::string> vec, const std::string value);
std::vector<double> where(const std::vector<double> vec, const std::string condition, const double value, const bool return_values);
std::vector<double> where(const VectorXd& vec, const std::string condition, const double value, const bool return_values);
std::vector<int> where_index(const std::vector<double> vec, const std::string condition, const double value);
std::vector<int> where_index(const VectorXd& vec, const std::string condition, const double value);


std::string read_lastline_ascii(const std::string filename);
VectorXd order_input_params(const VectorXd& cte_params, const VectorXd& var_params, const std::vector<std::string> cte_names, 
		const std::vector<std::string> var_names, const std::vector<std::string> param_names);
std::vector<int> where_int(const std::vector<int> vec, const int value);
VectorXi where_int(const VectorXi& vec, const int value);
VectorXi where_dbl(const VectorXd& vec, const double value, const double tolerance);
VectorXi where_dbl(const VectorXd& vec, double value, const double tolerance, const int imin_search, const int imax_search);
std::vector<int> where_dbl(const std::vector<double> vec, const double value, const double tolerance);
std::vector<int> where_in_range(const std::vector<double> vec, const double value_min, const double value_max, const bool strict);
VectorXi where_in_range(const VectorXd& vec, const double value_min, const double value_max, const bool strict);

VectorXi str_to_Xiarr(const std::string str, const std::string delimiters);
VectorXd str_to_Xdarr(const std::string str, const std::string delimiters);

int str_to_int(const std::string str);
bool str_to_bool(const std::string str);
long str_to_lng(const std::string str);
long str_to_long(const std::string str);
double str_to_dbl(const std::string str);
std::vector<double> str_to_dblarr(const std::string str, const std::string delimiters);
std::string int_to_str(const int value);
std::string dbl_to_str(const double ind);
std::string lng_to_str(const long ind);
std::vector<double> str_to_arrdbl(const std::string str, const std::string delimiters);
std::vector<int> str_to_arrint(const std::string str, const std::string delimiters);

VectorXd arrstr_to_Xdarrdbl(const std::vector<std::string> vals_strs);
VectorXi arrstr_to_Xiarr(const std::vector<std::string> vals_strs);
std::vector<double> arrstr_to_arrdbl(const std::vector<std::string> vals_strs);

std::vector<double> vectXd_to_vec(const VectorXd& vecXd_in);
bool file_exists(const std::string& name);
std::string strtrim(const std::string& str);
std::vector<std::string> strsplit(const std::string str, const std::string delimiters);

std::vector<double> filter_range(const std::vector<double> param_in, const std::vector<double> f_in, const double fmin, const double fmax);
std::vector<bool> filter_range(const std::vector<bool> param_in, const std::vector<double> f_in, const double fmin, const double fmax);


//std::vector<std::string> strsplit(const std::string str, const std::string delimiters); // Alias of strsplit2... used in diagnostics.cpp. THERE IS A WARNING ON IT PUT ON 17/06/2021

