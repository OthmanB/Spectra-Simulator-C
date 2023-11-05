#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <vector>
using Eigen::VectorXd;


VectorXd linspace(const long double start_in, const long double end_in, const long num_in);
std::vector<double> linspace_vec(const long double start, const long double end, const long num);
