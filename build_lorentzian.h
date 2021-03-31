#include <math.h>
#include <Eigen/Dense>
#include "function_rot.h"

using Eigen::VectorXd;

VectorXd build_l_mode_a1etaa3_simu(const VectorXd x_l, double H_l, double fc_l, double f_s, double eta, double a3, double gamma_l, const int l, VectorXd V);
VectorXd build_l_mode_act_simu(const VectorXd x_l, double H_l, double fc_l, double f_s, double eta, double a3, double b, double alpha, double asym, double gamma_l, const int l, VectorXd V);

VectorXd build_l_mode_a1a2a3(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s, const double a2, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V);
