# integrate Glm~F(theta, phi).|Ylm|^2
Small program for integrating Glm(theta, phi), term that can be used to define activity effect on mode splittings. It uses modified functions taken from https://github.com/CD3/libIntegrate to perform integration in C++ using Gaussian Legendre Quadrature

This code can calculate Glm as defined in Gizon 2002, AN, 323, 251–253
But can be used to have more flexible activity zone hypothesis than that paper.

# How to use
Use cmake to compile into a build directory:
       'mkdir build'
       'cd build'
       'cmake ..'
       'make'

This will generate an executable with parameters as per defined into the main() of ylm.cpp


Othman Benomar
