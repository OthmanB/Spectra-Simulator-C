#include <iostream>
#include <iomanip>
#include <filesystem>
#include "bilinear_interpol.h"
#include <algorithm>
#include <vector>

// Function that load all grids for l=1 to lmax in order 
// to be used by Alm_interp_fast
// This program expect that the grids exist
GridData_Alm_fast loadAllData(const std::string grid_dir, const std::string ftype){
    GridData_Alm_fast grids;
	const int lmax=3;
    int lm;
    grids.error=false;
    try{
        if (!std::filesystem::is_directory(grid_dir + "/" + ftype)) {
              throw std::runtime_error("Invalid grid directory or file type");
        }
		for(int l=1; l<=lmax;l++){
            for(int m=0;m<=l;m++){
                // The expected file name is following the one at creation by GridMaker: 
                // (l=1, m=-1) --> A10.gz, (l=1, m=0) --> A11.gz etc...
                // NOTE: WE HERE ONLY READ THE m>0 ... thus A11.gz, A12.gz, A22.gz, A23 etc...
                // The root directory of the grid is expected to contain subdirectories 
                // of the same name as the allowed ftype
                std::string file_gz = grid_dir + "/" + ftype + "/A" + std::to_string(l) + std::to_string(m+l) + ".gz";
                // check if file exists
                if (!std::filesystem::exists(file_gz)) {
                    throw std::runtime_error("Grid file for (l=" + std::to_string(l) + ", m=" + std::to_string(m) + " does not exist");
                }                
                //std::cout << "       File : " << file_gz << std::endl;
                // Load the correct grid
                GridData data = loadGridData(file_gz);
                lm = 10.*l + 1.*m; // m=-1 or +1 have the same results, hence the loop only on 0<m<l
                //std::cout << "l= " << l << "   lm = " << lm << std::endl;
                switch(lm){
                    // Interpolation on a case by case basis
                    case 10:
                        grids.A10=data;
                        break;
                    case 11:
                        grids.A11=data;
                        break;
                    case 20:
                        grids.A20=data;
                        break;
                    case 21:
                        grids.A21=data;
                        break;
                    case 22:
                        grids.A22=data;
                        break;
                    case 30:
                        grids.A30=data;
                        break;
                    case 31:
                        grids.A31=data;
                        break;
                    case 32:
                        grids.A32=data;
                        break;
                    case 33:
                        grids.A33=data;
                        break;
                    default:
                        throw std::invalid_argument("Invalid lm combination");
                } 
            }
        }       
    }
	catch (exception& e) {
		cerr << "Error: " << e.what() << "\n";
        grids.error=true;
	}
	catch (...) {
		cerr << "Exception of unknown type!\n";
        grids.error=true;
	}
    return grids;
}

// Use interpolation over precomputed grids to compute Alm
double Alm_interp(const int l, const int m, 
                const long double theta0, const long double delta, 
                const std::string ftype, const std::string grid_dir){
    if (l <= 0){
        std::cerr << "Error : The function Alm_interp() does not allow l<=0. " << std::endl;
        return -9998;
    }	
    try{
		// The expected file name is following the one at creation by GridMaker: 
		// (l=1, m=-1) --> A10.gz, (l=1, m=0) --> A11.gz etc...
		// The root directory of the grid is expected to contain subdirectories 
		// of the same name as the allowed ftype
		std::string file_gz = grid_dir + "/" + ftype + "/A" + std::to_string(l) + std::to_string(m+l) + ".gz";
		//std::cout << "       File : " << file_gz << std::endl;
        // Load the correct grid
		GridData data = loadGridData(file_gz);
		//show_mat(data);
		// Interpolation
         double c = interpolate(data, theta0, delta);
		//std::cout << "Interpolated value at (" << theta0 << ", " << delta << ") = " << c << '\n';
		return c;
	}
	catch (exception& e) {
		cerr << "error: " << e.what() << "\n";
		return -8888;
	}
	catch (...) {
		cerr << "Exception of unknown type!\n";
		return -9999;
	}
	return -9999;
}

// Use interpolation over precomputed grids to compute Alm
// Instead of getting a grid_dir as input, it takes a GridData_Alm_fast structure that
// should contain all of the table A10, A11 etc... in a preset way. This allows much faster calls 
// to the interpolation compared to Alm_interp as the function does not need to unpack the gz file
// It will be particularly usefull if Alm must be calculated in a loop
// However, it is limited to lmax=3 
// boost::math::bicubic_b_spline<double>& interp, 
double Alm_interp_iter(
                        const int l, const int m, const long double theta0, const long double delta, 
                        const std::string ftype, GridData_Alm_fast grids){
	const int lmax=3;
    int lm;
    double c;
    if (l <= 0){
        std::cerr << "Error : The function Alm_interp() does not allow l<=0. " << std::endl;
        return -9998;
    }	
    if (l > lmax){
        std::cerr << "Error : The function Alm_interp4iter() is not suitable for lmax>3. " << std::endl;
        std::cerr << "        Please use Alm_interp() instead" << std::endl;
        return -9998;
    }
	try{
        lm = 10*l + std::abs(m); // m=-1 or +1 have the same results, hence the abs(m)
        switch(lm){
            // Interpolation on a case by case basis
            case 10:
                c = interpolate(grids.A10, theta0, delta);
                break;
            case 11:
                c = interpolate(grids.A11, theta0, delta);
                break;
            case 20:
                c = interpolate(grids.A20, theta0, delta);
                break;
            case 21:
                c = interpolate(grids.A21, theta0, delta);
                break;
            case 22:
                c = interpolate(grids.A22, theta0, delta);
                break;
            case 30:
                c = interpolate(grids.A30, theta0, delta);
                break;
            case 31:
                c = interpolate(grids.A31, theta0, delta);
                break;
            case 32:
                c = interpolate(grids.A32, theta0, delta);
                break;
            case 33:
                c = interpolate(grids.A33, theta0, delta);
                break;
            default:
                throw std::invalid_argument("Invalid lm combination");
        }
		//std::cout << "Interpolated value at (" << theta0 << ", " << delta << ") = " << c << '\n';
		return c;
	}
	catch (exception& e) {
		cerr << "Error: " << e.what() << "\n";
		return -8888;
	}
	catch (...) {
		cerr << "Exception of unknown type!\n";
		return -9999;
	}
	return -9999;
}

// The fastest possible: Does not take the grid but the gsl object that contain the initialised grid
double Alm_interp_iter_preinitialised(
                        const int l, const int m, const long double theta0, const long double delta, 
                        const std::string ftype, gsl_funcs interp_funcs){
	const int lmax=3;
    int lm;
    double c;
    if (l <= 0){
        std::cerr << "Error : The function Alm_interp() does not allow l<=0. " << std::endl;
        return -9998;
    }	
    if (l > lmax){
        std::cerr << "Error : The function Alm_interp4iter() is not suitable for lmax>3. " << std::endl;
        std::cerr << "        Please use Alm_interp() instead" << std::endl;
        return -9998;
    }
	try{
        lm = 10*l + std::abs(m); // m=-1 or +1 have the same results, hence the abs(m)
        switch(lm){
            // Interpolation on a case by case basis
            case 10:
                        /*std::cout << "    flat_grid_A10.x = ";
                        std::cout.flush();
                        for(int kl=0; kl<interp_funcs.flat_grid_A10.nx;kl++){
                            std::cout  << interp_funcs.flat_grid_A10.x[kl] << "  ";
                            std::cout.flush();
                        }
                        std::cout << std::endl << "    flat_grid_A10.y = ";
                        std::cout.flush();
                        for(int kl=0; kl<interp_funcs.flat_grid_A10.ny;kl++){
                            std::cout  << interp_funcs.flat_grid_A10.y[kl] << "  ";
                            std::cout.flush();
                        }
                        std::cout << std::endl;*/
                c = interpolate_core(interp_funcs.interp_A10, interp_funcs.flat_grid_A10, theta0, delta);
                //std::cout << " ===----> " << c << std::endl;
                break;
            case 11:
                        /*std::cout << "    flat_grid_A11.x = ";
                        std::cout.flush();
                        for(int kl=0; kl<interp_funcs.flat_grid_A11.nx;kl++){
                            std::cout  << interp_funcs.flat_grid_A11.x[kl] << "  ";
                            std::cout.flush();
                        }
                        std::cout << std::endl << "    flat_grid_A11.y = ";
                        std::cout.flush();
                        for(int kl=0; kl<interp_funcs.flat_grid_A11.ny;kl++){
                            std::cout  << interp_funcs.flat_grid_A11.y[kl] << "  ";
                            std::cout.flush();
                        }
                        std::cout << std::endl;*/
                c = interpolate_core(interp_funcs.interp_A11, interp_funcs.flat_grid_A11, theta0, delta);
                //std::cout << " ===----> " << c << std::endl;
                break;
            case 20:
                        /*std::cout << "    flat_grid_A20.x = ";
                        std::cout.flush();
                        for(int kl=0; kl<interp_funcs.flat_grid_A20.nx;kl++){
                            std::cout  << interp_funcs.flat_grid_A20.x[kl] << "  ";
                            std::cout.flush();
                        }
                        std::cout << std::endl << "    flat_grid_A20.y = ";
                        std::cout.flush();
                        for(int kl=0; kl<interp_funcs.flat_grid_A20.ny;kl++){
                            std::cout  << interp_funcs.flat_grid_A20.y[kl] << "  ";
                            std::cout.flush();
                        }
                        std::cout << std::endl;*/
                c = interpolate_core(interp_funcs.interp_A20, interp_funcs.flat_grid_A20, theta0, delta);
                //std::cout << " ===----> " << c << std::endl;
                break;
            case 21:
                        /*std::cout << "    flat_grid_A21.x = ";
                        std::cout.flush();
                        for(int kl=0; kl<interp_funcs.flat_grid_A21.nx;kl++){
                            std::cout  << interp_funcs.flat_grid_A21.x[kl] << "  ";
                            std::cout.flush();
                        }
                        std::cout << std::endl << "    flat_grid_A21.y = ";
                        std::cout.flush();
                        for(int kl=0; kl<interp_funcs.flat_grid_A21.ny;kl++){
                            std::cout  << interp_funcs.flat_grid_A21.y[kl] << "  ";
                            std::cout.flush();
                        }
                        std::cout << std::endl;*/
                c = interpolate_core(interp_funcs.interp_A21, interp_funcs.flat_grid_A21, theta0, delta);
                //std::cout << " ===----> " << c << std::endl;
                break;
            case 22:
                        /*std::cout << std::endl << "    flat_grid_A22.x = ";
                        std::cout.flush();
                        for(int kl=0; kl<interp_funcs.flat_grid_A22.nx;kl++){
                            std::cout  << interp_funcs.flat_grid_A22.x[kl] << "  ";
                            std::cout.flush();
                        }
                        std::cout << std::endl << "    flat_grid_A22.y = ";
                        std::cout.flush();
                        for(int kl=0; kl<interp_funcs.flat_grid_A22.ny;kl++){
                            std::cout  << interp_funcs.flat_grid_A22.y[kl] << "  ";
                            std::cout.flush();
                        }
                        std::cout << std::endl;*/
                c = interpolate_core(interp_funcs.interp_A22, interp_funcs.flat_grid_A22, theta0, delta);
                //std::cout << " ===----> " << c << std::endl;
                break;
            case 30:
                        /*std::cout << std::endl << "    flat_grid_A30.x = ";
                        std::cout.flush();
                        for(int kl=0; kl<interp_funcs.flat_grid_A30.nx;kl++){
                            std::cout  << interp_funcs.flat_grid_A30.x[kl] << "  ";
                            std::cout.flush();
                        }
                        std::cout <<  std::endl << "    flat_grid_A30.y = ";
                        std::cout.flush();
                        for(int kl=0; kl<interp_funcs.flat_grid_A30.ny;kl++){
                            std::cout  << interp_funcs.flat_grid_A30.y[kl] << "  ";
                            std::cout.flush();
                        }
                        std::cout << std::endl;*/
                c = interpolate_core(interp_funcs.interp_A30, interp_funcs.flat_grid_A30, theta0, delta);
                //std::cout << " ===----> " << c << std::endl;
                break;
            case 31:
                        /*std::cout << std::endl << "    flat_grid_A31.y = ";
                        std::cout.flush();
                        for(int kl=0; kl<interp_funcs.flat_grid_A31.nx;kl++){
                            std::cout  << interp_funcs.flat_grid_A31.x[kl] << "  ";
                            std::cout.flush();
                        }
                        std::cout <<  std::endl << "    flat_grid_A31.y = ";
                        std::cout.flush();
                        for(int kl=0; kl<interp_funcs.flat_grid_A31.ny;kl++){
                            std::cout  << interp_funcs.flat_grid_A31.y[kl] << "  ";
                            std::cout.flush();
                        }
                        std::cout << std::endl;*/
                c = interpolate_core(interp_funcs.interp_A31, interp_funcs.flat_grid_A31, theta0, delta);
                //std::cout << " ===----> " << c << std::endl;
                break;
            case 32:
                c = interpolate_core(interp_funcs.interp_A32, interp_funcs.flat_grid_A32, theta0, delta);
                break;
            case 33:
                c = interpolate_core(interp_funcs.interp_A33, interp_funcs.flat_grid_A33, theta0, delta);
                break;
            default:
                std::cerr << "Invalid lm combination" << std::endl;
                return -9998;
        }
		//std::cout << "Interpolated value at (" << theta0 << ", " << delta << ") = " << c << '\n';
		return c;
	}
	catch (exception& e) {
		cerr << "Error: " << e.what() << "\n";
		return -8888;
	}
	catch (...) {
		cerr << "Exception of unknown type!\n";
		return -9999;
	}
	return -9999;
}


/*
// A quick test
int main(int argc, char** argv){
    int l=3;
    int m=-1;
    double theta0=0.999816;
    double delta=0.233657;
    std::string ftype="gate";
    std::string grid_dir="/Users/obenomar/tmp/Alm/data/Alm_grids_CPP/2deg_grids/";
    std::cout << Alm_interp(l, m, theta0, delta, ftype, grid_dir) << std::endl;
}
*/
