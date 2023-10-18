/**
 * @file interpol.cpp
 * @brief A collection of functions for interpolation and resampling of arrays.
 */
#include <Eigen/Dense>
#include <cmath>
#include "interpol.h"
#include <iostream>
#include <iomanip>

using Eigen::VectorXd;

/**
 * @brief Linearly interpolates a curve (x, y) at the position x_int.
 * 
 * @param x The x-coordinates of the curve.
 * @param y The y-coordinates of the curve.
 * @param x_int The x-coordinate at which to interpolate.
 * @return double The interpolated y-coordinate.
 */
double lin_interpol(const VectorXd& x, const VectorXd& y, const double x_int){

	long i=0, Nx=x.size();
	double a=0,b=0; 
	VectorXd xtmp, ytmp;
	
	// ---- case of an actual interpolation -----
	if(x_int >= x.head(1)(0) && x_int <= x.tail(1)(0)){ 
		while((x_int < x[i] || x_int > x[i+1])){
			//std::cout << "x_int < x[i] || x_int > x[i+1]) && i<Nx-1 =" << ((x_int < x[i] || x_int > x[i+1]) && i<Nx-1) << std::endl;
			i=i+1;
		}
		if(i==0 && (x_int < x[i] || x_int > x[i+1])){i=i+1;} // case where we never passed by the loop because x_int < x[0] || x_int > x[1] = True	        
		a=(y[i+1] - y[i])/(x[i+1] - x[i]); // slope
		b=y[i] - a*x[i]; // ordinate at origin
	}
	// ---- case of an extrapolation toward the lower edge ----
	if(x_int  < x.head(1)(0)){
		a=(y[1] - y[0])/(x[1] - x[0]); // slope
		b=y[0] - a*x[0]; // ordinate at origin 
	}
	// ---- case of an extrapolation toward the upper edge ----
	if(x_int  > x.tail(1)(0)){
		ytmp=y.tail(2);
		xtmp=x.tail(2);
		a=(ytmp[1] - ytmp[0])/(xtmp[1] - xtmp[0]); // slope
		b=ytmp[0] - a*xtmp[0]; // ordinate at origin 
    }
return a*x_int+b;
}

/**
 * @brief Quadratically interpolates the array a and resamples it to size m.
 * 
 * @param a The array to be interpolated.
 * @param m The desired size of the resampled array.
 * @return VectorXd The resampled array.
 */
VectorXd quad_interpol( const VectorXd& a, const int m ){


    long double step = double( a.size() - 1 ) / (m - 1);
    VectorXd b(m);

    for( int j = 0; j < m; j ++ ){
        b[j] = interp2( j*step, a);
    }
    return b;
}

// -------------- Extra routines -------------------

/**
 * @brief Linearly interpolates a value x in an array a.
 * 
 * @param x The value to be interpolated.
 * @param a The array to interpolate.
 * @return long double The interpolated value.
 */
long double interp1( const long double x,  const VectorXd& a)
{
    int n=a.size();

    if( x <= 0 )  return a[0];
    if( x >= n - 1 )  return a[n-1];
    int j = int(x);
    return a[j] + (x - j) * (a[j+1] - a[j]);
}

/**
 * @brief Linearly interpolates an array a and resamples it to size m.
 * 
 * @param a The array to be interpolated.
 * @param m The desired size of the resampled array.
 * @return VectorXd The resampled array.
 */
VectorXd  inter1parray( const VectorXd& a, const int m )
{
    long double step = double( a.size() - 1 ) / (m - 1);
    VectorXd b(m);

    for( int j = 0; j < m; j ++ ){
        b[j] = interp1( j*step, a);
    }
    return b;
}

/**
 * @brief Calculates the value of a parabola passing through three points.
 * 
 * @param x The x-coordinate at which to evaluate the parabola.
 * @param f_1 The y-coordinate of the point with x-coordinate -1.
 * @param f0 The y-coordinate of the point with x-coordinate 0.
 * @param f1 The y-coordinate of the point with x-coordinate 1.
 * @return long double The value of the parabola at x.
 */
long double parabola( const long double x, const long double f_1, const long double f0, const long double f1 ){
    if( x <= -1 )  return f_1; 
    if( x >= 1 )  return f1; 
    long double l = f0 - x * (f_1 - f0);
    long double r = f0 + x * (f1 - f0);
    return (l + r + x * (r - l)) / 2;
}

/**
 * @brief Quadratically interpolates a value x in an array a.
 * 
 * @param x The value to be interpolated.
 * @param a The array to interpolate.
 * @return long double The interpolated value.
 */
long double interp2( const long double x, const VectorXd& a){
    if( x <= .5  ||  x >= a.size() - 1.5 )
        return interp1( x, a);
    int j = int( x + .5 );
    long double t = 2 * (x - j);  // -1 .. 1
    return parabola( t, (a[j-1] + a[j]) / 2, a[j], (a[j] + a[j+1]) / 2 );
}
