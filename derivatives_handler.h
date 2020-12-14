//#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include "data.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;


// ---------------------------------------------------------
// -------------- **** 1st Derivative ***** ----------------
// ---------------------------------------------------------

Deriv_out Fstder_backward_1err_reggrid(const VectorXd& y);
Deriv_out Fstder_backward_1err_reggrid(const VectorXd& y, const VectorXd& x);
Deriv_out Fstder_forward_1err_reggrid(const VectorXd& y);
Deriv_out Fstder_forward_1err_reggrid(const VectorXd& y, const VectorXd& x);
Deriv_out Fstder_centered_2err_reggrid(const VectorXd& y);
Deriv_out Fstder_centered_2err_reggrid(const VectorXd& y, const VectorXd& x);
Deriv_out Fstder_centered_4err_reggrid(const VectorXd y);
Deriv_out Fstder_centered_4err_reggrid(const VectorXd& y, const VectorXd& x);
Deriv_out Scndder_backward_1err_reggrid(const VectorXd& y);
Deriv_out Scndder_backward_1err_reggrid(const VectorXd& y, const VectorXd& x);
Deriv_out Scndder_forward_1err_reggrid(const VectorXd& y);
Deriv_out Scndder_forward_1err_reggrid(const VectorXd& y, const VectorXd& x);
Deriv_out Scndder_centered_2err_reggrid(const VectorXd& y);
Deriv_out Scndder_centered_2err_reggrid(const VectorXd& y, const VectorXd& x);

// Main Functions
Deriv_out Frstder_adaptive_reggrid(const VectorXd& y);
Deriv_out Frstder_adaptive_reggrid(const VectorXd& y, const VectorXd& x);
Deriv_out Scndder_adaptive_reggrid(const VectorXd& y);
Deriv_out Scndder_adaptive_reggrid(const VectorXd& y, const VectorXd& x);
