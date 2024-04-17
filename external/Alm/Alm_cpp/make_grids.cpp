// Program to make a grid the is compatible with the 2d interpoler program
// see https://github.com/OthmanB/2D_interp
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <cmath>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include "linspace.h"
#include "data.h"
#include "activity.h"
#include "gzip_compress.h"

bool saveAlm(const GridData_Alm& grid, const std::string& outdir, const std::string& prefix, const std::string& suffix) {
    std::string file_txt="filetmp.txt";
    std::string file_gz;
    if (grid.Alm_grid.empty()) {
        std::cerr << "Error: Grid is empty" << std::endl;
        return 0;
    }   
    for (unsigned int m = 0; m < grid.Alm_grid.size(); ++m) {
        file_gz = outdir + "/" + prefix + std::to_string(m) + suffix;
        bool err1=writeToFile(grid.theta, grid.delta, grid.Alm_grid[m], file_txt);
        bool err2=compress_gzip(file_txt, file_gz);
        bool err3=eraseFile(file_txt);
        if(err1 || err2 ||  err3){
            std::cerr << "Error: Check for errors in saveAlm()" << std::endl;
            return 0;
        }
    }
    return 1;
}

GridData_Alm make_Alm_grid(const int l, const double theta_min, const double theta_max, 
        const double delta_min, const double delta_max, 
        const double resol, const std::string ftype, bool verbose){
    //
    const long double delta_limit=0.001, theta0_limit=M_PI;
    if (theta_min < 0 || theta_max > M_PI || theta_max <= theta_min){
        std::cerr << "Error: theta_min must be > 0 and theta_max must < Pi with also theta_max>theta_min" << std::endl;
        std::cerr << "       theta_min = " << theta_min << std::endl; 
        std::cerr << "       theta_max = " << theta_max << std::endl; 
        exit(EXIT_FAILURE);
    }
    if (delta_min < 0 || delta_max > M_PI || delta_max <= delta_min){
        std::cerr << "Error: delta_min must be > 0 and delta_max must < Pi with also delta_max>delta_min" << std::endl;
        std::cerr << "       delta_min = " << delta_min << std::endl; 
        std::cerr << "       delta_max = " << delta_max << std::endl; 
        exit(EXIT_FAILURE);
    }
    if (ftype != "gauss" && ftype != "triangle" && ftype != "gate"){
        std::cerr << "Error: ftype must be either 'gate', 'triangle' or 'gauss' " << std::endl;
        std::cerr << "       ftype = " << ftype << std::endl;
        exit(EXIT_FAILURE);        
    }
    if (resol <= 0){
        std::cerr << "Error: resol must be > 0" << std::endl;
        std::cerr << "       resol = " << resol << std::endl;
        exit(EXIT_FAILURE);
    }
    // Initialisations
    int i,j,tot;
    double Alm_norm, r;
    int Ntheta = ceil((theta_max - theta_min) / resol); // convert float to integer using ceil function
    int Ndelta = ceil((delta_max - delta_min) / resol); // convert float to integer using ceil function
    std::vector<double> theta=linspace_vec(theta_min, theta_max, Ntheta);
    std::vector<double> delta=linspace_vec(delta_min, delta_max, Ndelta);
    std::vector<std::vector<std::vector<double>>> Alm_all(2*l + 1);
    for (auto& ai : Alm_all) {
        ai.assign(Ndelta, std::vector<double>(Ntheta, 0.0));
    }
    if(verbose == true){
        std::cout << "Number of data point on the theta axis:" << Ntheta << std::endl;
        std::cout << "Number of data point on the delta axis:" <<  Ndelta << std::endl;
        std::cout << "Resolution: " <<  resol << std::endl;
        std::cout << "Theta range: " << "[" << theta_min << " , " << theta_max << ']'  << std::endl;
        std::cout << "Delta range: " << "[" << delta_min << " , " <<  delta_max << ']' << std::endl;
        std::cout << "ftype     = " << ftype << std::endl;
        std::cout << "l         = " << l << std::endl;
    }
    i=0; // Index in theta
	j=0; // Index in delta
	for (auto theta0 : theta){
        if(verbose == true){
            std::cout << "theta0 = " << theta0 << "   index:" << i+1 << "/" << Ntheta << "  (timestamp: " << std::time(NULL) << ")" << std::endl;
            std::cout << "          index : [ 1 , " << Ndelta << " ]" << std::endl;
        }
        for (auto delta0 : delta){
            if (delta0 >=delta_limit && theta0>=0 && theta0< theta0_limit){
                for (int m=-l;m<=l; m++){
                    r=Alm(l, m, theta0, delta0, ftype);
                    if (ftype == "gauss"){
                        Alm_norm= gauss_filter_cte(theta0, delta0);
                    }
                    else{
                        Alm_norm=1;
                    }
                    r=r/Alm_norm;
                    Alm_all[m+l][j][i] = r;
                }
            }
            if (delta0 <delta_limit && theta0>=0  && theta0< theta0_limit){
                for (int m=-l;m<=l; m++){
                    Alm_all[m+l][j][i] = 0;
                }
            }
            j++;
        }
        j=0;
        i++;
    }
    i=0;	
	const double resol_theta=theta[1]-theta[0]; // Use the two first element of the regular grid to evaluate the effective resolution of the grid as it may differ (larger) than resol
	const double resol_delta=delta[1]-delta[0];

    GridData_Alm grid;
    grid.theta=theta;
    grid.delta=delta;
    grid.Alm_grid=Alm_all;
    grid.resol_theta=resol_theta;
    grid.resol_delta=resol_delta;
    grid.Ntheta=Ntheta;
    grid.Ndelta=Ndelta;
    grid.l=l;
	return grid;
}

