/* 
 *	Simple linear fit algorithm 
 * 
 *  Created on: 11 Oct 2017
 *      Author: obenomar
*/

# include <Eigen/Dense>
# include <iostream>
# include <iomanip>
# include "linfit.h"

using Eigen::VectorXd;

// Supposedly more optimised. Tested and compared with the old version in unit_test/linfit
// Gains are of x2
VectorXd linfit(const VectorXd& x, const VectorXd& y) {
    if (x.size() != y.size()) {
        std::cout << "x and y do not have the same size!" << std::endl;
        std::cout << "Cannot proceed. The program will exit now." << std::endl;
        exit(EXIT_SUCCESS);
    }
    double sx = x.sum();
    double sy = y.sum();
    double n = x.size();
    double mean_x = sx / n;
    
    VectorXd t = x.array() - mean_x;
    VectorXd tmp = t.array() * y.array();
    double out0 = tmp.sum() / (t.array() * t.array()).sum();
    double out1 = (sy - sx * out0) / n;
    VectorXd out(2);
    out << out0, out1;
    return out;
}



