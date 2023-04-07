#pragma once 
#include <Eigen/Dense>
#include <vector>

using Eigen::VectorXd;
using Eigen::MatrixXd;

struct Data_asympt_p{
	bool error_status; // Used to detect if there was a problem when writting the data
	VectorXd n;
	double Dnu;
	double epsilon;
	double d0l;
	VectorXd O2_term;	
};

struct Freq_modes{
	bool error_status; // Used to detect if there was a problem when writting the data
    VectorXd fl0;
    VectorXd fl1;
    VectorXd fl2;
    VectorXd fl3;
	Data_asympt_p p_asymptotic;
};

struct Data_file{
	bool error_status; // Used to detect if there was a problem when writting the data
	std::vector<std::string> header;
	bool do_rescale;
	double Dnu_target;
	double epsilon_target;
	VectorXd d0l_target;
	MatrixXd data;
};
