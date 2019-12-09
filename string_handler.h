/*
 * string_handler.h
 *
 * Contains all kind of header for the functions
 * used to arrange/handle the strings
 * 
 *  Created on: 20 Jun 2016
 *      Author: obenomar
 */
#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

std::string strtrim(const std::string& str);
std::vector<int> where_str(std::vector<std::string> vec, std::string value);
std::vector<int> where_int(std::vector<int> vec, int value);
VectorXi where_int(VectorXi vec, int value);
VectorXi where_dbl(VectorXd vec, double value, double tolerance);
std::vector<int> where_dbl(std::vector<double> vec, double value, double tolerance);
std::vector<int> where_in_range(std::vector<double> vec, double value_min, double value_max, bool strict);
VectorXi where_in_range(VectorXd vec, double value_min, double value_max, bool strict);
VectorXi str_to_Xiarr(const std::string str, const std::string delimiters);
VectorXd str_to_Xdarr(const std::string str, const std::string delimiters);
std::vector<double> arrstr_to_arrdbl(const std::vector<std::string> vals_strs);
VectorXd arrstr_to_Xdarrdbl(const std::vector<std::string> vals_strs);
VectorXi arrstr_to_Xiarr(const std::vector<std::string> vals_strs);
std::vector<double> str_to_arrdbl(const std::string str, const std::string delimiters);
double str_to_dbl(const std::string str);
int str_to_int(const std::string str);
long str_to_long(const std::string str);
bool str_to_bool(const std::string str);
std::vector<int> str_to_arrint(const std::string str, const std::string delimiters);
// -- Transfered from various other functions on 1 May 2018 --
std::string dbl_to_str(const double ind);
std::vector<std::string> strsplit(const std::string str, const std::string delimiters); // Alias of strsplit2... used in diagnostics.cpp
std::string int_to_str(const int value);
std::vector<double> str_to_dblarr(const std::string str, const std::string delimiters);
std::vector<double> filter_range(const std::vector<double> param_in, const std::vector<double> f_in, const double fmin, const double fmax);
std::vector<bool> filter_range(const std::vector<bool> param_in, const std::vector<double> f_in, const double fmin, const double fmax);
