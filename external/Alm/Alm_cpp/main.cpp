#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include "GaussLegendre2D.hpp"
#include <Eigen/Dense>
#include "activity.h"
//#include "linspace.h"
#include "version.h"
#include "bilinear_interpol.h"
#include "Alm_interpol.h"

namespace po = boost::program_options;

void show_version(){
    std::cout << " " << APP_NAME << " " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << "\n built on " __DATE__ << std::endl;
#   if defined(__clang__)
    	printf(" with clang " __clang_version__);
#   elif defined(__GNUC__)
    	printf(" with GCC");
    	printf(" %d.%d.%d", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#   elif defined(_MSC_VER)
    	printf(" with MSVC");
    	printf(" %d", MSVC_VERSION);
#   else
    printf(" unknown compiler");
#   endif

    std::cout << "\n architecture:";
#   if defined(__i386__) || defined(_M_IX86)
    	std::cout << " i386" << std::endl;
#   elif defined(__x86_64__) || defined(_M_AMD64)
    	std::cout << " x86_64" << std::endl;
#	elif (defined(__arm64__) && defined(__APPLE__)) || defined(__aarch64__)
		std::cout << " arm64 / Apple" << std::endl;
#   elif 
		std::cout << " Unknown" << std::endl;
#   endif
    std::cout << " Author: " << APP_COPYRIGHT << std::endl;
}

int main(int argc, char* argv[]){

	const long double delta_limit=0.001, theta0_limit=M_PI;
	double r, Alm_norm;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
		("version,v", "show program version")
        ("degree,d", po::value<int>(), "mode degree")
        ("theta0", po::value<double>(), "active region colatitude theta0 in rad")
        ("delta", po::value<double>(), "active region width in rad")
        ("filter_type,f", po::value<std::string>(), "filter type")
		("grid_dir, g", po::value<std::string>()->default_value(""), "If set, the computation is made by interpolation of grids mades with the GridMaker. \nOtherwise, perform a full integral computation");
	po::variables_map vm;
	//po::store(po::parse_command_line(argc, argv, desc), vm);
	//po::notify(vm);
	try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return -1;
    }
	if (vm.count("help")) {
		std::cout << desc << std::endl;
		return 1;
	}
	if (vm.count("version")) {
        show_version();
		return 1;
    }
	if (!(vm.count("degree") && vm.count("theta0") && vm.count("delta") && vm.count("filter_type"))) {
		std::cerr << "Missing required arguments!\n" << desc << "\n";
		return -2;
	}
	int l = vm["degree"].as<int>();
	double theta0 = vm["theta0"].as<double>();
	double delta = vm["delta"].as<double>();
	std::string ftype = vm["filter_type"].as<std::string>();
	std::string grid_dir = vm["grid_dir"].as<std::string>();

	if (ftype == "gauss" && grid_dir.empty() == false){
		std::cout << "Warning: The gauss approach with grid showed to be very imprecise even with a 1 deg grid (up to 10 percent error )" << std::endl;
		std::cout << "         The reason being unknown, please do not use the grid approach with a gauss. " <<std::endl;
		std::cout << "         Use instead a direct calculation. If you require speed, prefer the triangle approximation" << std::endl;
		std::cout << "         It is a very good approximation for the Sun. And in fact better than the gaussian spot zone model" << std::endl;
		return 0;
	}
	std::cout << "#---------------" << std::endl;
	std::cout << "#Configuration: " << std::endl;
	std::cout << "#l=" << l << std::endl;
	std::cout << "#theta0 =" << theta0 << std::endl;
	std::cout << "#delta = " << delta << std::endl;
	std::cout << "#ftype = " << ftype << std::endl;
	std::cout << "#--------------"  << std::endl;
	std::cout << "#l     m      Alm" << std::endl;
	if (delta >=delta_limit && theta0>=0 && theta0< theta0_limit){
		if (ftype == "gauss"){
			Alm_norm=gauss_filter_cte(theta0, delta);
		} else{
			Alm_norm=1;
		}
		//std::cout << "Alm_norm =" << Alm_norm << std::endl;
		for (int m=-l;m<=l; m++){
			if (grid_dir.empty()){
				r=Alm(l, m, theta0, delta, ftype); // full computation of the integral
				std::cout << l << "   " << m << "  " << r/Alm_norm << std::endl;
			} else{
				r=Alm_interp(l, m, theta0, delta, ftype, grid_dir); // interpolation on a provided gri
				std::cout << l << "   " << m << "  " << r << std::endl; // Here r is already normalised in the table
				//if(r <= -8888){
					//std::cerr << "Error : The interpolation of Alm returned the error code " << r << std::endl;
				//}
			}
		}	
		return 1;
	}
	if (delta <delta_limit && theta0>=0  && theta0< theta0_limit){
		for (int m=-l;m<=l; m++){
			std::cout << l << "   " << m << "  " << 0 << std::endl;
		}
		return 1;
	}
	if (theta0<0 || theta0 > theta0_limit){
		for (int m=-l;m<=l; m++){
			std::cout << l << "   " << m << "  " << -9999 << std::endl;
		}
		return 0;
	}	
	return 0; // error if we passed by none of the if statements
}
