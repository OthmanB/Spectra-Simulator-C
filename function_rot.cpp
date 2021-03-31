/*
 * function_rot.cpp
 *
 *  Created on: 11 Feb 2016
 *      Author: obenomar
 */
#include <math.h>
#include <Eigen/Dense>
//#include "function_rot.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

double combi(int n, int r);
double dmm(int l, int m1, int m2, double beta);
int factorial(int n);
MatrixXd function_rot( int l, double beta);
VectorXd amplitude_ratio(int l, double beta);

VectorXd amplitude_ratio(int l, double beta){
/* Main function that calculates the mode visibilities
   requires: 
	- l: degree of the mode
	- beta: stellar inclination in degrees
   returns:
	- A vector of size 2*l+1 which contains the 
	  visibilities of the m componnents
   dependences:
	- math.h
	- the eigen library. This program was written with Eigen 3.1.
*/
	const int dim=2*l+1;
	const double PI = 3.141592653589793238462643;
	double angle=PI*beta/180.;
        double norm;
        MatrixXd zz;
        VectorXd vv;
	VectorXd pv= VectorXd::Constant(dim, 2);

        
        zz=function_rot(l, angle);
        vv=zz.col(l);
        
        for(int i=0; i<dim; i++){vv(i)=vv(i)*vv(i);}
	
return vv;        
}

MatrixXd function_rot(int l, double beta){
	
	const long dim=2*l+1;
	MatrixXd mat(dim,dim);

	for(int i=0; i<=l; i++){
		for(int j=-i; j<= i; j++){
			mat(i+l,j+l)=dmm(l, i, j, beta);
		}
	}

	for(int i=-l; i<=0 ; i++){
		for(int j=i; j<= -i; j++){
			mat(i+l,j+l)=mat(-i+l,-j+l)*pow(-1, i-j);
		}
	}

	for(int j=0; j<=l; j++){
		for(int i=-j; i<= j; i++){
			mat(i+l,j+l)=dmm(l, j, i, -beta);
		}
	}

	for(int j=-l; j<=0; j++){
		for(int i=j; i<= -j; i++){
			mat(i+l,j+l)=mat(-i+l,-j+l)*pow(-1, i-j);
		}
	}

	return mat;
}

double dmm(int l, int m1, int m2, double beta){
	double sum=0;
	double var=0;
	for(long s=0; s<= l-m1; s++ ){
		var=combi(l+m2, l-m1-s)*combi(l-m2, s)*pow(-1, l-m1-s);
		var=var*pow(cos(beta/2.),2*s+m1+m2)*pow(sin(beta/2.), 2*l-2*s-m1-m2);
		sum=sum+var;
	}
	sum=sum*sqrt(factorial(l+m1)*factorial(l-m1));
	sum=sum/sqrt(factorial(l+m2)*factorial(l-m2));

	return sum;
}

double combi(int n, int r){
	return factorial(n)/factorial(n-r)/factorial(r);
}

int factorial(int n) {
    long factorial=1;

    for (long i = 1; i <= n; i++) {
        factorial = factorial*i; 
    }
    return factorial;
}

