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

VectorXd linfit(const VectorXd& x, const VectorXd& y){

   double sx, sy, st2;
   VectorXd v(x.size()), t(x.size()), tmp(x.size()), out(2);
   
   if (y.size() != x.size()){
	std::cout << " x and y do not have same size!" << std::endl;
	std::cout << " Cannot proceed. The program will exit now" << std::endl;
	exit(EXIT_SUCCESS);
   }

   sx = x.sum();
   sy = y.sum();

 
   v.setConstant(sx/x.size());
   //std::cout << "Here3" << std::endl;

   t = x - v;
   //std::cout << "t=" << t << std::endl;
   
   for(int i=0;i<x.size();i++){
   	tmp[i]=t[i]*y[i];
   }
   out[0] = tmp.sum(); // a coefficient

   //std::cout << "out[0]=" << out[0] << std::endl;

   for(int i=0; i<x.size();i++){ tmp[i]=t[i]*t[i];}
   //std::cout << t << std::endl;
   st2=tmp.sum();
   //std::cout << "Here4" << std::endl;

   // parameter estimates
   out[0] = out[0] / st2;
   out[1] = (sy - sx * out[0]) / x.size(); // b coefficient

   //std::cout << out << std::endl;
  return out;
}


