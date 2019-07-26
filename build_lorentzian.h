#include <math.h>
#include <Eigen/Dense>
#include "function_rot.h"

using Eigen::VectorXd;

VectorXd build_l_mode_a1etaa3_simu(const VectorXd x_l, double H_l, double fc_l, double f_s, double eta, double a3, double gamma_l, const int l, VectorXd V);
VectorXd build_l_mode_act_simu(const VectorXd x_l, double H_l, double fc_l, double f_s, double eta, double a3, double b, double alpha, double asym, double gamma_l, const int l, VectorXd V);
