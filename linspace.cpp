/**
 * @file linspace.cpp
 * @brief Header file for the linspace function
 *
 * Header file for the function that generate equally spaced eigen vector of values
 * 
 *
 * @date 20 Apr 2016
 * @author obenomar
 */

#include <Eigen/Dense>
#include <iostream>

using Eigen::VectorXd;

// My linspace function
VectorXd linspace(const long double start_in, const long double end_in, const long num_in)
{
	if (num_in == 0) {
		std::cout << " num_in in linspace is 0. Cannot create a linspace vector. Returning -1." << std::endl;
		VectorXd linspaced(1);
		linspaced[0]=-1;
		return linspaced;
	}
	VectorXd linspaced(num_in);

	const long double delta = (end_in - start_in) / (num_in - 1);
	for(long i=0 ; i< num_in ; i++){
		linspaced[i]=start_in + delta*i;
	}

	return linspaced;
}
