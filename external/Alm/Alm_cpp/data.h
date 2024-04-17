#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <gsl/gsl_interp2d.h>

// Define a struct to hold grid data
struct GridData_Alm {
    int l;
    int Ntheta;
    int Ndelta;
    double resol_theta;
    double resol_delta;
    std::vector<double> theta;
    std::vector<double> delta;
    std::vector<std::vector<std::vector<double>>> Alm_grid;
};

// Define a struct to hold grid data for interpolation
struct GridData {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<std::vector<double>> z;
    int n_rows;
    int n_cols;
};

// Faster retrieval, but fix structure limited to lmax=3
struct GridData_Alm_fast{
    bool error;
    GridData A10; // stores l=1, m=0
    GridData A11; // stores l=1, m=+/-1
    GridData A20; // stores l=2, m=0
    GridData A21; // stores l=2, m=+/-1
    GridData A22; // stores l=2, m=+/-2
    GridData A30; // stores l=3, m=0
    GridData A31; // stores l=3, m=+/-1
    GridData A32; // stores l=3, m=+/-2
    GridData A33; // stores l=3, m=+/-3    
};

 // Flat grid for direct inputs to gsl_interp2d_eval_e
 // This avoids to flatten the array at every iteration in an
 // iterative call of interpolate
struct GridData4gsl {
    double* x;
    double* y;
    double* z;
    int nx;
    int ny;
};

// structure containing the gsl initialised grid functions
// This is similar to GridData_Alm_fast but with the gsl function instead of the matrix
struct gsl_funcs{
    bool valid;
    gsl_interp2d* interp_A10;
    gsl_interp2d* interp_A11;
    gsl_interp2d* interp_A20;
    gsl_interp2d* interp_A21;
    gsl_interp2d* interp_A22;
    gsl_interp2d* interp_A30;
    gsl_interp2d* interp_A31;
    gsl_interp2d* interp_A32;
    gsl_interp2d* interp_A33;
    GridData4gsl flat_grid_A10;
    GridData4gsl flat_grid_A11;
    GridData4gsl flat_grid_A20;
    GridData4gsl flat_grid_A21;
    GridData4gsl flat_grid_A22;
    GridData4gsl flat_grid_A30;
    GridData4gsl flat_grid_A31;
    GridData4gsl flat_grid_A32;
    GridData4gsl flat_grid_A33;
};