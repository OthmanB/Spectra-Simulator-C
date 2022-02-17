#pragma once
#include <math.h>
#include <Eigen/Dense>
#include "function_rot.h"
# include <iostream>
# include <iomanip>
#include <string>

using Eigen::VectorXd;
using Eigen::VectorXi;

double Qlm(const int l, const int m);
VectorXi set_imin_imax(const VectorXd& x, const int l, const double fc_l, const double gamma_l, const double f_s, const double c, const double step);


VectorXd build_l_mode_a1l_etaa3(const VectorXd& x_l, const double H_l,  const double fc_l,  const double f_s1,  const double f_s2,  const double eta,  const double a3,  const double asym,  const double gamma_l, const int l,  const VectorXd& V);
VectorXd build_l_mode_a1etaa3(const VectorXd& x_l,  const double H_l,  const double fc_l,  const double f_s,  const double eta, const double a3,  const double asym,  const double gamma_l, const int l, const VectorXd& V);
VectorXd build_l_mode_a1acta3(const VectorXd& x_l,  const double H_l,  const double fc_l,  const double f_s,  const double eta, const double a3,  const double b,  const double alpha,  const double asym, double gamma_l, const int l,  const VectorXd& V);
VectorXd build_l_mode_a1etaa3_v2(const VectorXd& x_l, const VectorXd& H_lm,  const double fc_l,  const double f_s,  const double eta,  const double a3,  const double asym,  const double gamma_l, const int l);
VectorXd build_l_mode_a1l_etaa3_v2(const VectorXd& x_l, const VectorXd& H_lm,  const double fc_l,  const double f_s1,  const double f_s2,  const double eta,  const double a3,  const double asym, double gamma_l, const int l);
VectorXd build_l_mode_a1l_a2a3(const VectorXd& x_l, const double H_l,  const double fc_l,  const double f_s1,  const double f_s2,  const double a2,  const double a3,  const double asym,  const double gamma_l, const int l,  const VectorXd& V);
VectorXd build_l_mode_a1a2a3(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s, const double a2, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V);

VectorXd optimum_lorentzian_calc_a1l_etaa3(const VectorXd& x, const VectorXd& y, const  double H_l, const  double fc_l,  const double f_s1,  const double f_s2,  const double eta,  const double a3,  const double asym,  const double gamma_l, const int l,  const VectorXd& V,  const double step, const double c);
VectorXd optimum_lorentzian_calc_a1etaa3(const VectorXd& x, const VectorXd& y,  const double H_l,  const double fc_l,  const double f_s,  const double eta,  const double a3,  const double asym,  const double gamma_l, const int l,  const VectorXd& V,  const double step, const double c);
VectorXd optimum_lorentzian_calc_a1acta3(const VectorXd& x, const VectorXd& y,  const double H_l,  const double fc_l,  const double f_s,  const double eta,  const double a3, 
		 const double b,  const double alpha,  const double asym,  const double gamma_l, const int l,  const VectorXd& V,  const double step, const double c);

VectorXd optimum_lorentzian_calc_a1l_etaa3_v2(const VectorXd& x,  const VectorXd& y,  const VectorXd& H_lm,  const double fc_l,  const double f_s1,  const double f_s2,  const double eta,  const double a3,  const double asym,  const double gamma_l, const int l, double step, const double c);
VectorXd optimum_lorentzian_calc_a1etaa3_v2(const VectorXd& x,  const VectorXd& y,  const VectorXd& H_lm,  const double fc_l,  const double f_s,  const double eta,  const double a3,  const double asym,  const double gamma_l, const int l,  const double step, const double c);

VectorXd optimum_lorentzian_calc_a1l_a2a3(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, const double f_s1, const double f_s2, const double a2, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c);
VectorXd optimum_lorentzian_calc_a1a2a3(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, const double f_s, const double a2, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c);

VectorXd build_l_mode_a1etaAlma3(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s, 
    const double eta0, const double epsilon_nl, const VectorXd& thetas, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V);
VectorXd optimum_lorentzian_calc_a1etaAlma3(const VectorXd& x, const VectorXd& y,  const double H_l,  const double fc_l,  const double f_s, 
		const double eta0, const double epsilon_nl, const VectorXd& thetas, const double a3,  const double asym,  const double gamma_l, const int l,  const VectorXd& V,  const double step, const double c);

VectorXd build_l_mode_act_simu(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s, const double eta, const double a3, 
        const double b, const double alpha, const double asym, const double gamma_l, const int l, const VectorXd& V); // Only here for compatibiility reasons

// Generalized form of a-coefficients
VectorXd build_l_mode_aj(const VectorXd& x_l, const double H_l, const double fc_l, 
        const double a1, const double a2, const double a3, const double a4, const double a5, const double a6, 
        const double eta, const double asym, const double gamma_l, const int l, const VectorXd& V);
VectorXd optimum_lorentzian_calc_aj(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, 
        const double a1, const double a2, const double a3, const double a4, const double a5, const double a6, 
        const double eta0, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c);
