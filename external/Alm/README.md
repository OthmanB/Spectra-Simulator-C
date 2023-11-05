# Integrate Alm~F(theta, phi).|Ylm|^2 and suite of tools related to Alm
This is a project that compute the activity effect on a-coefficients, assuming a pure geometrical relationship between the activity
and the stellar shape. The perturbation on the even a-coefficients is determined following a prescription that is first described in Gizon 2002, AN, 323, 251–253, but generalized to three cases for the active region shape: A gate function, a gaussian and a triangular function. Note that the triangular function is in fact the most realistic to reproduce solar-spots.

# What is here?

First, the directory Alm_cpp contains the source codes in C++ that compute the Alm term using a modified functions taken from https://github.com/CD3/libIntegrate to perform integration in C++ using Gaussian Legendre Quadrature. If you want to determine the perturbation, you may get to this directory and read further the README.md there in order to know how to use it.

Secondly, the python_rendering directory contains some python implementation of the the same integral computation using the scipy integrate, double quadratic function. An additional set of functions are to make grids of Alm, to translate into a-coefficients and to see their dependence with the active zone size and extension. Finally, there is also functions overlaps the solar butterfly diagram with the filter function of your choice (e.g triangle).

Here are more details on each files:

       - show_filter.py: A small function that calls the C++ program Filter_Alm in order to visualise the filter function that is implemented in C++. Mostly for debug purpose.
       
       - test_convergence.py: A function execute both the python and C++ version of the Alm computation and compares them in term of precision. Mostly for debug purpose.
       
       - make_Alm_grid.py: Tools to make grids of Alm for l=1,2,3. This calls the C++ function to get faster
       
       - acoefs.py: small routines to perform the conversion between frequencies and a-coefficients and to compute the symmetric/asymmetric splittings
       
       - activity.py: Implements the python version of the C++ code. Note that this is ~10 - 100 times slower than the C++ code, so it is to be used for tests purpose only
       
       - eval_aj_range.py: functions that can generate the min/max range of some a-coefficients, provided a grid of Alm terms and for a given set of stellar parameters
       
       - show_aj_fct_theta_delta.py: Use a grid of Alm in order to show what is the dependence of aj in function of theta0 and delta
       
       - show_butterflydiagram.py: Used to show the butterfly diagram of the Sun with the filter function of your choice superimposed.


Othman Benomar
