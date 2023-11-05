# integrate Alm~F(theta, phi).|Ylm|^2
Small program for integrating Glm(theta, phi), term that can be used to define activity effect on mode splittings. It uses modified functions taken from https://github.com/CD3/libIntegrate to perform integration in C++ using Gaussian Legendre Quadrature

This code can calculate Alm=Integral(Ylm^2*Filter). With Filter ~ Gate function, the solution are as defined in Gizon 2002, AN, 323, 251–253
But can be used to have more flexible activity zone hypothesis than that paper (Filter : Gauss, Gate or Triangle).

The given code requires libraries such as Boost, Iostream, and Complex. 

If the user provides the required input parameters for the program to execute, it executes the Alm function that calculates Alm using Integral Computation or Interpolation of grids. It can also choose between different shapes for the active zones (gate, triangle, gauss) depending upon the filter type selected by the user. 

Finally, the calculated values are printed with the respective values of l, m, and the coefficients.

The code comes with unit and property-based tests in the test directory.
The code also comes with a suite of tests and visualisation tools in python in ../python_rendering.

# How to use
Use cmake to compile into a build directory:

       mkdir build
       
       cd build
       
       cmake ..
       
       make

This will generate two executables in the ../bin directory.

The Alm executable computes the integral with parameters as per defined into the main() of main.cpp 

The Filter_Alm returns a file with the shape of the filter, for the user-defined parameters (as defined by the get_filter.cpp program). This is more a test program to evaluate if the filtering is properly made (see python routines in ../python_rendering).

The make_Alm_grids allows you to generate grids in the format required to perform interpolations with https://github.com/OthmanB/2D_interp. This grid is also a required component for the TAMCMC code : https://github.com/OthmanB/TAMCMC-C

Additional functions in python are used to cross-verify that everything behave properly. 

Othman Benomar
