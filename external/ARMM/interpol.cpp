#include <Eigen/Dense>
#include <cmath>
#include "interpol.h"
# include <iostream>
# include <iomanip>

using Eigen::VectorXd;

const long double interp1( const long double x,  const VectorXd& a);
long double parabola( const long double x, const long double f_1, const long double f0, const long double f1 );
long double interp2( const long double x, const VectorXd& a);

double lin_interpol(const VectorXd& x, const VectorXd& y, const double x_int){
	/* Very simple function that linearly interpolate at the position x_int
	   a curve (x,y). Here we assume that x.size() = y.size()
	*/

	long i=0, Nx=x.size();
	double a=0,b=0; 
	VectorXd xtmp, ytmp;
	
	// ---- case of an actual interpolation -----
	//std::cout << "x_int=" << x_int << std::endl;

	if(x_int >= x.head(1)(0) && x_int <= x.tail(1)(0)){ 
		//while((x_int < x[i] || x_int > x[i+1]) && i<Nx-1){
		//while((x_int < x[i] || x_int > x[i+1]) && i<Nx-1){
		while((x_int < x[i] || x_int > x[i+1])){
			//std::cout << "x_int < x[i] || x_int > x[i+1]) && i<Nx-1 =" << ((x_int < x[i] || x_int > x[i+1]) && i<Nx-1) << std::endl;
			i=i+1;
		}
		if(i==0 && (x_int < x[i] || x_int > x[i+1])){i=i+1;} // case where we never passed by the loop because x_int < x[0] || x_int > x[1] = True	        
        //a=(y[i] - y[i-1])/(x[i] - x[i-1]); // slope
		//b=y[i-1] - a*x[i-1]; // ordinate at origin
		a=(y[i+1] - y[i])/(x[i+1] - x[i]); // slope
		b=y[i] - a*x[i]; // ordinate at origin
		//std::cout << "[" << i << "/" << x.size()-1 << "]  x_int=" << x_int << "   a=" << a << "   b=" << b << std::endl;
	}
	// ---- case of an extrapolation toward the lower edge ----
	if(x_int  < x.head(1)(0)){
		//std::cout << "Warning: we have to perform an extrapolation toward the lower edge" << std::endl;
		a=(y[1] - y[0])/(x[1] - x[0]); // slope
		b=y[0] - a*x[0]; // ordinate at origin 
	}
	// ---- case of an extrapolation toward the upper edge ----
	if(x_int  > x.tail(1)(0)){
		//std::cout << "Warning: we have to perform an extrapolation toward the upper edge" << std::endl;
		ytmp=y.tail(2);
		xtmp=x.tail(2);
		a=(ytmp[1] - ytmp[0])/(xtmp[1] - xtmp[0]); // slope
		b=ytmp[0] - a*xtmp[0]; // ordinate at origin 
    }
return a*x_int+b;
}


VectorXd quad_interpol( const VectorXd& a, const int m ){
   /* Function that quadratically interpolate the array a. 
    * The new size of the array is m, such that this 
    * function actually resample the array
   */

    long double step = double( a.size() - 1 ) / (m - 1);
    VectorXd b(m);

    for( int j = 0; j < m; j ++ ){
        b[j] = interp2( j*step, a);
    }
    return b;
}

// -------------- Extra routines -------------------

  // linear interpolate x in an array
const long double interp1( const long double x,  const VectorXd& a)
{
    int n=a.size();

    if( x <= 0 )  return a[0];
    if( x >= n - 1 )  return a[n-1];
    int j = int(x);
    return a[j] + (x - j) * (a[j+1] - a[j]);
}

    // linear interpolate array a[] -> array b[]
VectorXd  inter1parray( const VectorXd& a, const int m )
{
    long double step = double( a.size() - 1 ) / (m - 1);
    VectorXd b(m);

    for( int j = 0; j < m; j ++ ){
        b[j] = interp1( j*step, a);
    }
    return b;
}
    // parabola through 3 points, -1 < x < 1
long double parabola( const long double x, const long double f_1, const long double f0, const long double f1 ){
    if( x <= -1 )  return f_1; 
    if( x >= 1 )  return f1; 
    long double l = f0 - x * (f_1 - f0);
    long double r = f0 + x * (f1 - f0);
    return (l + r + x * (r - l)) / 2;
}

    // quadratic interpolate x in an array
long double interp2( const long double x, const VectorXd& a){
    if( x <= .5  ||  x >= a.size() - 1.5 )
        return interp1( x, a);
    int j = int( x + .5 );
    long double t = 2 * (x - j);  // -1 .. 1
    return parabola( t, (a[j-1] + a[j]) / 2, a[j], (a[j] + a[j+1]) / 2 );
}

/*
int main(){
        // a.out [n m] --
    int n = 6, m = 60;
    Eigen::VectorXd xa(n), a(n), xb, b;

    for (int i=0; i<n; i++){
	a[i]=i*i;
	xa[i]=0.1*i;
    }    
    xb.resize(m);
    for (int i=0; i<m; i++){
	xb[i]=0.01*i;
    }
    xb=quad_interpol(xa, m);  // Resample the array xa ==> Transform it to an array of size m
    b=quad_interpol(a, m);  // Resample the array xa ==> Transform it to an array of size m



    std::cout << " xa       a" << std::endl;
    for(int i=0; i<xa.size(); i++){
	std::cout << xa[i] << "    "  << a[i] <<  std::endl;
    }
    std::cout << " xb       b" << std::endl;
    for(int i=0; i<xb.size(); i++){
	std::cout << xb[i] << "    "  << b[i] <<  std::endl;
    }
 }
*/
