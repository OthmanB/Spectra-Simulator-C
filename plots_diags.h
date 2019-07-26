#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include "gnuplot-iostream.h"
#include "data.h"

//#include <gsl/gsl_histogram.h>
//#include <gsl/gsl_rng.h>


using Eigen::MatrixXd;
using Eigen::VectorXd;

void gnuplt_modelv0(VectorXd x, VectorXd y, VectorXd model, double data_scoef1, double data_scoef2, std::string file_model); // old type of function which is quite slow
void gnuplt_model(VectorXd x, VectorXd y, VectorXd model, double scoef1, double scoef2, std::string file_model);
std::vector<double> vectXd_to_vec(VectorXd vecXd_in);
std::string dbl_to_str(const double ind);
std::string lng_to_str(const long ind);


