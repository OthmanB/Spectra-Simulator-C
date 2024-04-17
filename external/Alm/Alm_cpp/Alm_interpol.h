#pragma once 
#include <iostream>
#include <iomanip>
#include "bilinear_interpol.h"

// Load all the grids required for l=1,2,3. To be used jointly with Alm_interp4iter()
GridData_Alm_fast loadAllData(const std::string grid_dir, const std::string ftype);

// Use interplation over precomputed grids to compute Alm
double Alm_interp(const int l, const int m, const long double theta0, const long double delta, const std::string ftype, const std::string grid_dir);
double Alm_interp_iter(const int l, const int m, const long double theta0, const long double delta, const std::string ftype, GridData_Alm_fast grids);
double Alm_interp_iter_preinitialised(
                        const int l, const int m, const long double theta0, const long double delta, 
                        const std::string ftype, gsl_funcs interp_funcs);