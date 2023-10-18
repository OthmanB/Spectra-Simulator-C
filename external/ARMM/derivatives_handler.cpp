/**
 * @file derivatives_handler.cpp
 * @brief A set of function to compute derivation
 * 
 * This is a set of function that can compute the first and second derivatives for any discrete function
 */
#include <iostream>
#include <Eigen/Dense>
#include "derivatives_handler.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;


// ---------------------------------------------------------
// -------------- **** 1st Derivative ***** ----------------
// ---------------------------------------------------------


/**
 * @brief Computes the first order backward derivative with 1st order precision.
 * 
 * @param y A vector containing 2 points: y at x-1 (y[0]) and y at x (y[1]).
 * @return Deriv_out A structure containing the computed derivative.
 * 
 * @note The vector y must have a size of at least 2.
 *       If the size is less than 2, the program will exit with an error message.
 */
Deriv_out Fstder_backward_1err_reggrid(const VectorXd& y){
	Deriv_out der;
	
	der.deriv.resize(y.size()-1);
	
	if(y.size() < 2){
		std::cout << "Vectors must be at least of size 2 to compute the first order backward derivative at first order precision" << std::endl;
		std::cout << "The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
	} else{
		for(int i=0; i<y.size()-1; i++){
			der.deriv[i]=y[i+1] - y[i];
		}	
	}
return der;	
}


/**
 * @brief Computes the first order backward derivative with 1st order precision.
 *
 * @param y A vector containing 2 points: y at x-1 (y[0]) and y at x (y[1]).
 * @param x A vector containing the corresponding x values for the points in y.
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vectors y and x must have a size of at least 2.
 * If the size is less than 2, the program will exit with an error message.
 */
Deriv_out Fstder_backward_1err_reggrid(const VectorXd& y, const VectorXd& x){
	Deriv_out der;
	
	der.deriv.resize(y.size()-1);
	if(y.size() < 2){
		std::cout << "Vectors must be at least of size 2 to compute the first order backward derivative at first order precision" << std::endl;
		std::cout << "The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
	} else{
		for(int i=0; i<y.size()-1; i++){
			der.deriv[i]=(y[i+1] - y[i])/(x[i+1] - x[i]);
		}	
	}
return der;	
}

/**
 * @brief Computes the first order forward derivative with 1st order precision.
 *
 * @param y A vector containing 2 points: y at x (y[0]) and y at x+1 (y[1]).
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vector y must have a size of at least 2.
 * If the size is less than 2, the program will exit with an error message.
 */
Deriv_out Fstder_forward_1err_reggrid(const VectorXd& y){
	Deriv_out der;
	
	der.deriv.resize(y.size()-1);

	if(y.size() < 2){
		std::cout << "Vectors must be at least of size 2 to compute the 1st order forward derivative at first order precision" << std::endl;
		std::cout << "The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
	} else{
		for(int i=0; i<y.size()-1; i++){
			der.deriv[i]=y[i+1] - y[i];
		}	
	}
return der;	
}

/**
 * @brief Computes the first order forward derivative with 1st order precision.
 *
 * @param y A vector containing 2 points: y at x (y[0]) and y at x+1 (y[1]).
 * @param x A vector containing the corresponding x values for the points in y.
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vectors y and x must have a size of at least 2.
 * If the size is less than 2, the program will exit with an error message.
 */
Deriv_out Fstder_forward_1err_reggrid(const VectorXd& y, const VectorXd& x){
	Deriv_out der;
	
	der.deriv.resize(y.size()-1);
	
	if(y.size() < 2){
		std::cout << "Vectors must be at least of size 2 to compute the 1st order forward derivative at first order precision" << std::endl;
		std::cout << "The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
	} else{
		for(int i=0; i<y.size()-1; i++){
			der.deriv[i]=(y[i+1] - y[i])/(x[i+1] - x[i]);
		}	
	}
return der;	
}

/**
 * @brief Computes the first order centered derivative with 2nd order precision.
 *
 * @param y A vector containing 3 points: y at x-1 (y[0]), y at x (y[1]), and y at x+1 (y[2]).
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vector y must have a size of at least 3.
 * If the size is less than 3, the program will exit with an error message.
 */
Deriv_out Fstder_centered_2err_reggrid(const VectorXd& y){
	Deriv_out der;
	
	der.deriv.resize(y.size()-2);
	
	if(y.size() < 3){
		std::cout << "Vectors must be at least of size 3 to compute the 1st order centered derivative at second order precision" << std::endl;
		std::cout << "The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
	} else{
		for(int i=0; i<y.size()-2; i++){
			der.deriv[i]=(y[i+2] - y[i])/2.;
		}	
	}
return der;	
}


/**
 * @brief Computes the first order centered derivative with 2nd order precision.
 *
 * @param y A vector containing 3 points: the position y at x-1 (y[0]), y at x (y[1]), and y at x+1 (y[2]).
 * @param x A vector containing the corresponding x values for the points in y.
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vectors y and x must have a size of at least 3.
 * If the size is less than 3, the program will exit with an error message.
 */
Deriv_out Fstder_centered_2err_reggrid(const VectorXd& y, const VectorXd& x){
	Deriv_out der;
	
	der.deriv.resize(y.size()-2);
	
	if(y.size() < 3){
		std::cout << "Vectors must be at least of size 3 to compute the 1st order centered derivative at second order precision" << std::endl;
		std::cout << "The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
	} else{
		for(int i=0; i<y.size()-2; i++){
			der.deriv[i]=(y[i+2] - y[i])/(x[i+2] - x[i]);
		}	
	}
return der;	
}

/**
 * @brief Computes the first order centered derivative with 4th order precision.
 *
 * @param y A vector containing 5 points: the position y at x-2 (y[0]) to y at x+2 (y[4]).
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vector y must have a size of at least 5.
 * If the size is less than 5, the program will exit with an error message.
 */
Deriv_out Fstder_centered_4err_reggrid(const VectorXd& y){
	Deriv_out der;
	
	der.deriv.resize(y.size()-4);
	
	if(y.size() < 5){
		std::cout << "Vectors must be at least of size 5 to compute the 1st order centered derivative at fourth order precision" << std::endl;
		std::cout << "The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
	} else{
		for(int i=0; i<y.size()-4; i++){
			der.deriv[i]=(-y[i+4] + 8.*y[i+3] - 8.*y[i+1] + y[i])/12.;
		}	
	}
return der;	
}

/**
 * @brief Computes the first order centered derivative with 4th order precision.
 *
 * @param y A vector containing 5 points: the position y at x-2 (y[0]) to y at x+2 (y[4]).
 * @param x A vector containing the corresponding x values for the points in y.
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vectors y and x must have a size of at least 5.
 * If the size is less than 5, the program will exit with an error message.
 */
Deriv_out Fstder_centered_4err_reggrid(const VectorXd& y, const VectorXd& x){
	double dx;
	Deriv_out der;
	
	der.deriv.resize(y.size()-4);
	
	if(y.size() < 5){
		std::cout << "Vectors must be at least of size 5 to compute the 1st order centered derivative at fourth order precision" << std::endl;
		std::cout << "The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
	} else{
		for(int i=0; i<y.size()-4; i++){
			dx=(x[i+4] - x[i])/4;
			der.deriv[i]=(-y[i+4] + 8.*y[i+3] - 8.*y[i+1] + y[i])/(12.*dx);
		}	
	}
return der;	
}
// ---------------------------------------------------------
// -------------- **** 2nd Derivative ***** ----------------
// ---------------------------------------------------------

/**
 * @brief Computes the second order backward derivative with 1st order precision.
 *
 * @param y A vector containing 3 points: the position y at x-2 (y[0]) to y at x (y[2]).
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vector y must have a size of at least 3.
 * If the size is less than 3, the program will exit with an error message.
 */
Deriv_out Scndder_backward_1err_reggrid(const VectorXd& y){
	Deriv_out der;
	
	der.deriv.resize(y.size()-2);
	//der.error.resize(y.size()-2);
	
	if(y.size() < 3){
		std::cout << "Vectors must be at least of size 3 to compute the second order backward derivative at first order precision" << std::endl;
		std::cout << "The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
	} else{
		for(int i=0; i<y.size()-2; i++){
			der.deriv[i]=y[i+2] - 2.*y[i+1] + y[i];
			//der.error[i]=sqrt(pow(err_y[i+2],2) + 4.*pow(err_y[i+1],2) + pow(err_y[i],2));
		}	
	}
return der;	
}

/**
 * @brief Computes the second order backward derivative with 1st order precision.
 *
 * @param y A vector containing 3 points: the position y at x-2 (y[0]) to y at x (y[2]).
 * @param x A vector containing the corresponding x values for the points in y.
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vectors y and x must have a size of at least 3.
 * If the size is less than 3, the program will exit with an error message.
 */
Deriv_out Scndder_backward_1err_reggrid(const VectorXd& y, const VectorXd& x){
	double dx;
	Deriv_out der;
	
	der.deriv.resize(y.size()-2);
	
	if(y.size() < 3){
		std::cout << "Vectors must be at least of size 3 to compute the second order backward derivative at first order precision" << std::endl;
		std::cout << "The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
	} else{
		for(int i=0; i<y.size()-2; i++){
			dx=(x[i+2] - x[i])/2.;
			der.deriv[i]=(y[i+2] - 2.*y[i+1] + y[i])/pow(dx,2.);
		}	
	}
return der;	
}

/**
 * @brief Computes  the second order forward derivative with 1st order precision.
 *
 * @param y A vector containing 3 points: the position y at x (y[0]) to y at x+2 (y[2]).
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vector y must have a size of at least 3.
 * If the size is less than 3, the program will exit with an error message.
 */
Deriv_out Scndder_forward_1err_reggrid(const VectorXd& y){
	Deriv_out der;
	
	der.deriv.resize(y.size()-2);
	
	if(y.size() < 3){
		std::cout << "Vectors must be at least of size 3 to compute the second order forward derivative at first order precision" << std::endl;
		std::cout << "The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
	} else{
		for(int i=0; i<y.size()-2; i++){
			der.deriv[i]=y[i+2] - 2.*y[i+1] + y[i];
		}	
	}
return der;	
}

/**
 * @brief Computes the second order forward derivative with 1st order precision.
 *
 * @param y A vector containing 3 points: the position y at x (y[0]) to y at x+2 (y[2]).
 * @param x A vector containing the corresponding x values for the points in y.
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vectors y and x must have a size of at least 3.
 * If the size is less than 3, the program will exit with an error message.
 */
Deriv_out Scndder_forward_1err_reggrid(const VectorXd& y, const VectorXd& x){

	double dx;
	Deriv_out der;
	
	der.deriv.resize(y.size()-2);
	//der.error.resize(y.size()-2);
	
	if(y.size() < 3){
		std::cout << "Vectors must be at least of size 3 to compute the second order forward derivative at first order precision" << std::endl;
		std::cout << "The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
	} else{
		for(int i=0; i<y.size()-2; i++){
			dx=(x[i+2] - x[i])/2.;
			der.deriv[i]=(y[i+2] - 2.*y[i+1] + y[i])/pow(dx,2.);
			//der.error[i]=sqrt(pow(err_y[i+2],2) + 4.*pow(err_y[i+1],2) + pow(err_y[i],2));
		}	
	}
return der;	
}

/**
 * @brief Computes the second order centered derivative with second order precision.
 *
 * @param y A vector containing 3 points: the position y at x-1 (y[0]) to y at x+1 (y[2]).
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vector y must have a size of at least 3.
 * If the size is less than 3, the program will exit with an error message.
 */
Deriv_out Scndder_centered_2err_reggrid(const VectorXd& y){
	Deriv_out der;
	
	der.deriv.resize(y.size()-2);
	//der.error.resize(y.size()-2);
	
	if(y.size() < 3){
		std::cout << "Vectors must be at least of size 3 to compute the second order centered derivative at second order precision" << std::endl;
		std::cout << "The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
	} else{
		for(int i=0; i<y.size()-2; i++){
			der.deriv[i]=y[i+2] - 2.*y[i+1] + y[i];
			//der.error[i]=sqrt(pow(err_y[i+2],2) + 4.*pow(err_y[i+1],2) + pow(err_y[i],2));
		}	
	}
return der;	
}


/**
 * @brief Computes the second order centered derivative with second order precision.
 *
 * @param y A vector containing 3 points: the position y at x-1 (y[0]) to y at x+1 (y[2]).
 * @param x A vector containing the corresponding x values for the points in y.
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vectors y and x must have a size of at least 3.
 * If the size is less than 3, the program will exit with an error message.
 */
Deriv_out Scndder_centered_2err_reggrid(const VectorXd& y, const VectorXd& x){
	double dx;
	Deriv_out der;
	
	der.deriv.resize(y.size()-2);
	
	if(y.size() < 3){
		std::cout << "Vectors must be at least of size 3 to compute the second order centered derivative at second order precision" << std::endl;
		std::cout << "The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
	} else{
		for(int i=0; i<y.size()-2; i++){
			dx=(x[i+2] - x[i])/2.;
			der.deriv[i]=(y[i+2] - 2.*y[i+1] + y[i])/pow(dx,2.);
		}	
	}
return der;	
}

// ------------- Adaptive derivative -----------

/**
 * @brief Computes the first derivative of a list of values with adaptive precision.
 * 
 * The computation is adaptive in the sense that the precision order of the derivative is automatically adjusted depending on the region (lower edge, center, upper edge) in which the derivative is computed.
 * 
 * @param y A vector containing the list of values.
 * @return Deriv_out A structure containing the computed derivative.
 * 
 * @note The vector y must have a size greater than or equal to 2.
 * If the size is less than 2, the behavior is undefined.
 */
Deriv_out Frstder_adaptive_reggrid(const VectorXd& y){
	Deriv_out r_edge0, r_edge1, r_mid;
	Deriv_out der;
	
	der.deriv.resize(y.size());
	//der.error.resize(y.size());

	r_edge0=Fstder_forward_1err_reggrid(y.segment(0, 2));
	r_edge1=Fstder_backward_1err_reggrid(y.tail(2));
	
	// Put values on the edges
	der.deriv.segment(0,1)=r_edge0.deriv;
	der.deriv.tail(1)=r_edge1.deriv;
	
	// Put values on the center
	r_mid=Fstder_centered_2err_reggrid(y);
	der.deriv.segment(1, r_mid.deriv.size())=r_mid.deriv;

return der;
}

/**
 * @brief Computes the first derivative of a list of values with adaptive precision.
 *
 * The computation is adaptive in the sense that the precision order of the derivative is automatically adjusted depending on the region (lower edge, center, upper edge) in which the derivative is computed.
 *
 * @param y A vector containing the list of values.
 * @param x A vector containing the corresponding x values for the points in y.
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vectors y and x must have the same size and have a size greater than or equal to 2.
 * If the size is less than 2 or the sizes of y and x are not equal, the behavior is undefined.
 */
Deriv_out Frstder_adaptive_reggrid(const VectorXd& y, const VectorXd& x){
	Deriv_out r_edge0, r_edge1, r_mid;
	Deriv_out der;
	
	der.deriv.resize(y.size());

	r_edge0=Fstder_forward_1err_reggrid(y.segment(0, 2), x.segment(0, 2));
	r_edge1=Fstder_backward_1err_reggrid(y.tail(2), x.tail(2));
	
	// Put values on the edges
	der.deriv.segment(0,1)=r_edge0.deriv;
	der.deriv.tail(1)=r_edge1.deriv;
	
	// Put values on the center
	r_mid=Fstder_centered_2err_reggrid(y,x);
	der.deriv.segment(1, r_mid.deriv.size())=r_mid.deriv;

return der;
}

/**
 * @brief Computes the second derivative of a list of values with adaptive precision.
 *
 * The computation is adaptive in the sense that the precision order of the derivative is automatically adjusted depending on the region (lower edge, center, upper edge) in which the derivative is computed.
 *
 * @param y A vector containing the list of values.
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vector y must have a size greater than or equal to 3.
 * If the size is less than 3, the behavior is undefined.
 */
Deriv_out Scndder_adaptive_reggrid(const VectorXd& y){
	Deriv_out r_edge0, r_edge1, r_mid;
	Deriv_out der;
	
	der.deriv.resize(y.size());

	r_edge0=Scndder_forward_1err_reggrid(y.segment(0, 3));
	r_edge1=Scndder_backward_1err_reggrid(y.tail(3));
	
	// Put values on the edges
	der.deriv.segment(0,1)=r_edge0.deriv;
	der.deriv.tail(1)=r_edge1.deriv;
	
	// Put values on the center
	r_mid=Scndder_centered_2err_reggrid(y);
	der.deriv.segment(1, r_mid.deriv.size())=r_mid.deriv;

return der;
}


/**
 * @brief Computes the second derivative of a list of values with adaptive precision.
 *
 * The computation is adaptive in the sense that the precision order of the derivative is automatically adjusted depending on the region (lower edge, center, upper edge) in which the derivative is computed.
 *
 * @param y A vector containing the list of values.
 * @param x A vector containing the corresponding x values for the points in y.
 * @return Deriv_out A structure containing the computed derivative.
 *
 * @note The vectors y and x must have the same size and have a size greater than or equal to 3.
 * If the size is less than 3 or the sizes of y and x are not equal, the behavior is undefined.
 */
Deriv_out Scndder_adaptive_reggrid(const VectorXd& y, const VectorXd& x){
	Deriv_out r_edge0, r_edge1, r_mid;
	Deriv_out der;
	
	der.deriv.resize(y.size());
	//der.xderiv=x;

	r_edge0=Scndder_forward_1err_reggrid(y.segment(0, 3), x.segment(0, 3));
	r_edge1=Scndder_backward_1err_reggrid(y.tail(3), x.tail(3));
	
	// Put values on the edges
	der.deriv.segment(0,1)=r_edge0.deriv;
	der.deriv.tail(1)=r_edge1.deriv;
	
	// Put values on the center
	r_mid=Scndder_centered_2err_reggrid(y,x);
	der.deriv.segment(1, r_mid.deriv.size())=r_mid.deriv;

return der;
}

// Test Function
/*
#define BOOST_TEST_MODULE MyTest
#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include <Eigen/Dense>

// Include the header file for the functions being tested
#include "my_functions.h"

BOOST_AUTO_TEST_CASE(test_derivatives)
{
    VectorXd beta, dx, x;
    Deriv_out dout, dout2;

    dx.resize(100);
    x.resize(100);
    beta.resize(100);

    for (int i = 0; i < beta.size(); i++)
    {
        dx[i] = 1;
        if (i != 0)
            x[i] = x[i - 1] + dx[i];
        else
            x[i] = 0;

        beta[i] = x[i] * x[i] + 10;
    }

    dout = Frstder_adaptive_reggrid(beta);
    dout2 = Scndder_adaptive_reggrid(beta);

    std::cout << " ----- Without using a X-vector ----" << std::endl;
    std::cout << "First Derivatives=" << std::endl;
    for (int i = 0; i < beta.size(); i++)
        std::cout << "x[" << i << "]=" << x[i] << "   beta[" << i << "]=" << beta[i] << "   deriv1[" << i << "]=" << dout.deriv[i] << std::endl;

    std::cout << "Second Derivative=" << std::endl;
    for (int i = 0; i < beta.size(); i++)
        std::cout << "x[" << i << "]=" << x[i] << "   beta[" << i << "]=" << beta[i] << "   deriv2[" << i << "]=" << dout2.deriv[i] << std::endl;

    for (long i = 0; i < beta.size(); i++)
    {
        dx[i] = 0.34435; //- 0.001*i; // 0.1*i;
        if (i != 0)
            x[i] = x[i - 1] + dx[i];
        else
            x[i] = 0;

        beta[i] = x[i] * x[i] + 10;
    }

    dout = Frstder_adaptive_reggrid(beta, x);
    dout2 = Scndder_adaptive_reggrid(beta, x);
    std::cout << " ----- Using X-vector (regular spacing) ----" << std::endl;
    std::cout << "First Derivatives=" << std::endl;
    for (int i = 0; i < beta.size(); i++)
        std::cout << "x[" << i << "]=" << x[i] << "   beta[" << i << "]=" << beta[i] << "   deriv1[" << i << "]=" << dout.deriv[i] << std::endl;

    std::cout << "Second Derivative=" << std::endl;
    for (int i = 0; i < beta.size(); i++)
        std::cout << "x[" << i << "]=" << x[i] << "   beta[" << i << "]=" << beta[i] << "   deriv2[" << i << "]=" << dout2.deriv[i] << std::endl;

    // Add your assertions here
    BOOST_CHECK_EQUAL(dout.deriv.size(), beta.size());
    BOOST_CHECK_EQUAL(dout2.deriv.size(), beta.size());
}
*/