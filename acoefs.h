/*
Date: 16 Nov 2021
Functions that handle the acoefficients
It can create nu(n,l,m) from acoeffs for l<=3, j={1,2,3,4,5,6}
Or it can decompose into aj coefficients
Adapted from acoefs.py from the acoefs_check project (see github)
*/
#pragma once 
#include <math.h>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

long double Hslm_Ritzoller1991(const int s,const int l, const int m);
long double Pslm(const int s,const int l,const int m);
VectorXd Tnlm(VectorXd& nu_nlm, const int l);
VectorXd Snlm(VectorXd& nu_nlm, const int l);
VectorXd nunlm_from_acoefs(const long double nunl0, const int l, 
	const long double a1, const long double a2, const long double a3, const long double a4, const long double a5, const long double a6);
VectorXd eval_acoefs(const int l, VectorXd& nu_nls);
