#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <cmath>
#include <random>
#include "activity.h"
#include "Alm_interpol.h"
#include "bilinear_interpol.h"

//#include <boost/program_options.hpp>

int random_int(int min, int max) {
  std::random_device rd;  // obtain a random number from hardware
  std::mt19937 gen(rd()); // seed the generator
  
  // create a distribution
  std::uniform_int_distribution<> distr(min, max);
  
  // generate and return a random integer within the range [min, max]
  return distr(gen);
}

double random_double(double min, double max) {
  std::random_device rd;  // obtain a random number from hardware
  std::mt19937 gen(rd()); // seed the generator

  // create a distribution
  std::uniform_real_distribution<> distr(min, max);

  // generate and return a random double within the range [min, max]\n  
  return distr(gen);
}

std::string Alm_test(const int l, const double theta0, const double delta, const std::string ftype, const std::string grid_dir, const GridData_Alm_fast grids, const gsl_funcs interp_funcs,  const uint8_t use_grids) {
	const long double delta_limit=0.001, theta0_limit=M_PI;
	double r, Alm_norm;
    std::string str;
    if (delta >=delta_limit && theta0>=0 && theta0< theta0_limit){
		if (ftype == "gauss"){
			Alm_norm=gauss_filter_cte(theta0, delta);
		} else{
			Alm_norm=1;
		}
		for (int m=-l;m<=l; m++){
			if (use_grids == 0){
				r=Alm(l, m, theta0, delta, ftype); // full computation of the integral
			    str = std::to_string(l) + "   " + std::to_string(m) + "  " + std::to_string(r/Alm_norm) + "\n";
			} else{
                if (use_grids == 1){
                    if (grid_dir.empty() == false){
                        r=Alm_interp(l, m, theta0, delta, ftype, grid_dir); // interpolation on a provided on-disk grid
                    } else{
                        std::cerr << "Error: An Empty grid with use_grids = 1 is not allowed." << std::endl;
                        std::cerr << "       Either provide use_grids = 0 (no grid) or  = 2 (with a preloaded grid)" << std::endl;
                        return "0";
                    }
                }
                if (use_grids == 2){ 
                    try{
                        r=Alm_interp_iter(l, m, theta0, delta, ftype, grids); // interpolation on a provided preload grid
                    }
                    catch (exception& e) {
		                std::cerr << "Error: " << e.what() << "\n";
                        return "0";
	                }
                }
                if (use_grids == 3){ 
                    try{
                        //std::cout << "interp_funcs ( " << l << " , " << m << " )" << std::endl;
                        r=Alm_interp_iter_preinitialised(l, m, theta0, delta, ftype, interp_funcs); // interpolation on a provided preload grid
                    }
                    catch (exception& e) {
		                std::cerr << "Error: " << e.what() << "\n";
                        return "0";
	                }
                }                
                if (use_grids > 3){
                        std::cerr << "Error: use_grids must be either 0 (no grid), 1 (on-disk grid) or 2 (with a preloaded grid)" << std::endl;
                        return "0";
                    }                    
                str = std::to_string(l) + "   " + std::to_string(m) + "  " + std::to_string(r) + "\n"; // Here r is already normalised in the table
            }
		}
		return str;
	}
	if (delta <delta_limit && theta0>=0  && theta0< theta0_limit){
		for (int m=-l;m<=l; m++){
			str = std::to_string(l) + "   " + std::to_string(m) + "  \n";
		}
		return str;
	}
	if (theta0<0 || theta0 > theta0_limit){
		for (int m=-l;m<=l; m++){
			str = std::to_string(l) + "   " + std::to_string(m) + "  -9999\n";
		}
		return str;
	}	
    return "";   
}

int core_test(const std::string& ftype, const std::string& grid_dir, const int Niter, const int lmax) {
    const double theta_min = 0.0;
    const double theta_max = M_PI/2.0;  // fixed: added ".0" to specify that the division should be done in double
    const double delta_min = 0.0;
    const double delta_max = M_PI/4.0;  // fixed: added ".0" to specify that the division should be done in double
    const int Nshow=std::ceil(Niter/10);
    uint8_t use_grids;
    int l;
    double theta0, delta;
    auto start = std::chrono::high_resolution_clock::now();
    GridData_Alm_fast grids=loadAllData(grid_dir, ftype);
    std::string out;
    // Pre-initialisation of the grid into gsl : Flattening + gsl init
    gsl_funcs funcs_data;
    funcs_data.flat_grid_A10=flatten_grid(grids.A10);
    funcs_data.flat_grid_A11=flatten_grid(grids.A11);
    funcs_data.flat_grid_A20=flatten_grid(grids.A20);
    funcs_data.flat_grid_A21=flatten_grid(grids.A21);
    funcs_data.flat_grid_A22=flatten_grid(grids.A22);
    funcs_data.flat_grid_A30=flatten_grid(grids.A30);
    funcs_data.flat_grid_A31=flatten_grid(grids.A31);
    funcs_data.flat_grid_A32=flatten_grid(grids.A32);
    funcs_data.flat_grid_A33=flatten_grid(grids.A33);
    funcs_data.interp_A10=init_2dgrid(funcs_data.flat_grid_A10);
    funcs_data.interp_A11=init_2dgrid(funcs_data.flat_grid_A11);
    funcs_data.interp_A20=init_2dgrid(funcs_data.flat_grid_A20);
    funcs_data.interp_A21=init_2dgrid(funcs_data.flat_grid_A21);
    funcs_data.interp_A22=init_2dgrid(funcs_data.flat_grid_A22);
    funcs_data.interp_A30=init_2dgrid(funcs_data.flat_grid_A30);
    funcs_data.interp_A31=init_2dgrid(funcs_data.flat_grid_A31);
    funcs_data.interp_A32=init_2dgrid(funcs_data.flat_grid_A32);
    funcs_data.interp_A33=init_2dgrid(funcs_data.flat_grid_A33);
    
    start = std::chrono::high_resolution_clock::now();
    std::cout << " - Direct method... ";
    use_grids=0;
    for (int i = 0; i < Niter; i++) {
        l = random_int(1, lmax);
        theta0 = random_double(theta_min, theta_max);
        delta = random_double(delta_min, delta_max);
        out=Alm_test(l, theta0, delta, ftype, "", grids, funcs_data, use_grids);
        // Show iteration count every Nshow
        if (i % Nshow == 0) {
            std::cout << i << "...";
            std::cout.flush();
        }
    }
    std::cout << std::endl;
    auto stop = std::chrono::high_resolution_clock::now();
    auto durationA = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

    std::cout << " - Interpolation method (grids read on disk at each iteration)... ";
    start = std::chrono::high_resolution_clock::now();
    use_grids=1;
    for (int i = 0; i < Niter; i++) {
        l = random_int(1, lmax);
        theta0 = random_double(theta_min, theta_max);
        delta = random_double(delta_min, delta_max);
        out=Alm_test(l, theta0, delta, ftype, grid_dir, grids,funcs_data, use_grids);
        // Show iteration count every Nshow
        if (i % Nshow == 0) {
            std::cout << i << "...";
            std::cout.flush();
        }
    }
    std::cout << std::endl;
    stop = std::chrono::high_resolution_clock::now();
    auto durationB = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

    std::cout << " - Fast Interpolation method (grids stored in memory at each iteration)... ";
    start = std::chrono::high_resolution_clock::now();
    use_grids=2;
    for (int i = 0; i < Niter; i++) {
        l = random_int(1, lmax);
        theta0 = random_double(theta_min, theta_max);
        delta = random_double(delta_min, delta_max);
        out=Alm_test(l, theta0, delta, ftype, grid_dir, grids, funcs_data, use_grids);
        // Show iteration count every Nshow
        if (i % Nshow == 0) {
            std::cout << i << "...";
            std::cout.flush();
        }
    }
    std::cout << std::endl;
    stop = std::chrono::high_resolution_clock::now();
    auto durationC = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

    std::cout << " - Fastest Interpolation method (grids stored in memory + preflattened + preloaded in GSL)... " << std::endl;
    start = std::chrono::high_resolution_clock::now();
    use_grids=3;
    for (int i = 0; i < Niter; i++) {
        l = random_int(1, lmax);
        theta0 = random_double(theta_min, theta_max);
        delta = random_double(delta_min, delta_max);
        out=Alm_test(l, theta0, delta, ftype, grid_dir, grids, funcs_data, use_grids);
        // Show iteration count every Nshow
        if (i % Nshow == 0) {
            std::cout << i << "...";
            std::cout.flush();
        }
    }
    std::cout << std::endl;
    stop = std::chrono::high_resolution_clock::now();
    auto durationD = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

    // Output the results
    std::cout << "      Duration Alm                              " << durationA << std::endl;
    std::cout << "      Duration Alm_interp                       " << durationB << std::endl;
    std::cout << "      Duration Alm_interp_iter                  " << durationC << std::endl;
    std::cout << "      Duration Alm_interp_iter_preinitialised   " << durationD << std::endl;
    std::cout << "      Alm_interp                     " << std::fixed << std::setprecision(2) << (static_cast<double>(durationA)/durationB) << " faster than computing Alm directly " << std::endl;  // fixed: used std::fixed and std::setprecision to output the percentage with exactly 2 decimal places, and casted durationB to double to avoid integer division
    std::cout << "      Alm_interp_iter                " << std::fixed << std::setprecision(2) << (static_cast<double>(durationA)/durationC) << " faster than computing Alm directly " << std::endl;  // fixed: used std::fixed and std::setprecision to output the percentage with exactly 2 decimal places, and casted durationB to double to avoid integer division
    std::cout << "      Alm_interp_iter_preinitialised " << std::fixed << std::setprecision(2) << (static_cast<double>(durationA)/durationD) << " faster than computing Alm directly " << std::endl;  // fixed: used std::fixed and std::setprecision to output the percentage with exactly 2 decimal places, and casted durationB to double to avoid integer division
    return 0;
}


int main(){
    const int Niter = 1000;
    const int lmax = 2;
   int error;
   std::string ftype;
    std::string grid_dir;
    std::string grid_resol;
    grid_dir = "../../../data/Alm_grids_CPP/test_grids/";  // 2 degree resolution grid. Errors at 9e-3 max
    grid_resol="5 degrees";
    ftype="gate";
    std::cout << " Requested Niter = " << Niter << std::endl;
    std::cout << "           lmax  = " << lmax  << std::endl;
    std::cout << "----------------------- Testing ------------------------" << std::endl;
    std::cout << "--------         Activity shape  " << ftype << "----------------- " << std::endl;
    std::cout << "--------         Grid Resolution " << grid_resol << "------------- " << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    error=core_test(ftype, grid_dir, Niter, lmax);
    //
    ftype="triangle";
    std::cout << " Requested Niter = " << Niter << std::endl;
    std::cout << "           lmax  = " << lmax  << std::endl;
    std::cout << "----------------------- Testing ------------------------" << std::endl;
    std::cout << "--------         Activity shape  " << ftype << "----------------- " << std::endl;
    std::cout << "--------         Grid Resolution " << grid_resol << "------------- " << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    error=core_test(ftype, grid_dir, Niter, lmax);
    //
    grid_dir = "../../../data/Alm_grids_CPP/2deg_grids/";  // 2 degree resolution grid. Errors at 9e-3 max
    grid_resol="2 degrees";
    ftype="gate";
    std::cout << " Requested Niter = " << Niter << std::endl;
    std::cout << "           lmax  = " << lmax  << std::endl;
    std::cout << "----------------------- Testing ------------------------" << std::endl;
    std::cout << "--------         Activity shape  " << ftype << "----------------- " << std::endl;
    std::cout << "--------         Grid Resolution " << grid_resol << "------------- " << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    error=core_test(ftype, grid_dir, Niter, lmax);
    //
    ftype="triangle";
    std::cout << " Requested Niter = " << Niter << std::endl;
    std::cout << "           lmax  = " << lmax  << std::endl;
    std::cout << "----------------------- Testing ------------------------" << std::endl;
    std::cout << "--------         Activity shape  " << ftype << "----------------- " << std::endl;
    std::cout << "--------         Grid Resolution " << grid_resol << "------------- " << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    error=core_test(ftype, grid_dir, Niter, lmax);

    grid_dir = "../../../data/Alm_grids_CPP/1deg_grids/";  // 2 degree resolution grid. Errors at 9e-3 max
    grid_resol="1 degrees";
    ftype="gate";
    std::cout << " Requested Niter = " << Niter << std::endl;
    std::cout << "           lmax  = " << lmax  << std::endl;
    std::cout << "----------------------- Testing ------------------------" << std::endl;
    std::cout << "--------         Activity shape  " << ftype << "----------------- " << std::endl;
    std::cout << "--------         Grid Resolution " << grid_resol << "------------- " << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    error=core_test(ftype, grid_dir, Niter, lmax);
    //
    ftype="triangle";
    std::cout << " Requested Niter = " << Niter << std::endl;
    std::cout << "           lmax  = " << lmax  << std::endl;
    std::cout << "----------------------- Testing ------------------------" << std::endl;
    std::cout << "--------         Activity shape  " << ftype << "----------------- " << std::endl;
    std::cout << "--------         Grid Resolution " << grid_resol << "------------- " << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    error=core_test(ftype, grid_dir, Niter, lmax);

    return error;
}
