# integrate Glm~|Ylm|^2
Small program for integrating |Ylm|^2 in C++ using Gaussian Legendre Quadrature


This code is for calculating Glm as defined in Gizon 2002, AN, 323, 251–253
It is implemented with explicit limitation to l<=3

# How to use
Use cmake to compile into a build directory:
       'mkdir build'
       'cd build'
       'cmake ..'
       'make'

This will generate an executable with parameters as per defined into the main() of ylm.cpp
Note that this main is just for tests. A proper usage would require to include ylm.cpp into a larger projet,
after commenting the main().

Othman Benomar
