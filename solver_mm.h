#pragma once
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>

#include "version_solver.h"
#include "data.h"
//#include "string_handler.h" // Replaced by ioproc.h
#include "ioproc.h"
#include "interpol.h"
#include "derivatives_handler.h"

using Eigen::MatrixXd;
using Eigen::VectorXi;
using Eigen::VectorXd;

VectorXd linspace(const long double start_in, const long double end_in, const long num_in);
VectorXi sign_change(const VectorXd& x, bool return_indices=true);
VectorXd pnu_fct(const VectorXd& nu, const long double nu_p);
long double pnu_fct(const long double nu, const long double nu_p);
VectorXd gnu_fct(const VectorXd& nu, const long double nu_g, const long double Dnu_p, const long double DPl, const long double q);
long double gnu_fct(const long double nu, const long double nu_g, const long double Dnu_p, const long double DPl, const long double q);
long double asympt_nu_p(const long double Dnu_p, const int np, const long double epsilon, const int l, 
	const long double delta0l, const long double alpha, const long double nmax, long double r=0);
long double asympt_nu_p_from_l0(const VectorXd& nu_l0, const long double Dnu_p, const int np, const long double epsilon, const int l, 
	const long double delta0l, long double r=0);
VectorXd asympt_nu_p_from_l0_Xd(const VectorXd& nu_l0, const long double Dnu_p, const int l, const long double delta0l, long double fmin=-1, long double fmax=-1);
long double asympt_nu_g(const long double DPl, const int ng, const long double alpha, long double r=0);
Data_coresolver solver_mm(const long double nu_p, const long double nu_g, const long double Dnu_p, const long double DPl, const long double q, 
	const long double numin, const long double numax, const long double resol, const bool returns_axis=false, const bool verbose=false, const long double factor=0.05);
void test_sg_solver_mm();
Data_eigensols solve_mm_asymptotic_O2p(const long double Dnu_p, const long double epsilon, const int el, const long double delta0l, const long double alpha_p, 
	const long double nmax, const long double DPl, const long double alpha, const long double q, const long double sigma_p, 
	const long double fmin, const long double fmax, const long double resol, bool returns_pg_freqs=true, bool verbose=false);
Data_eigensols solve_mm_asymptotic_O2from_l0(const VectorXd& nu_l0_in, const int el, const long double delta0l, 
    const long double DPl, const long double alpha, const long double q, const long double sigma_p, 
	const long double resol, bool returns_pg_freqs=true, bool verbose=false, const long double freq_min=0, const long double freq_max=1e6);
void test_rgb_solver_mm();
void test_asymptotic_rgb_O2();
Data_eigensols test_asymptotic_sg_O2();
Data_eigensols test_asymptotic_sg_O2from_l0();
