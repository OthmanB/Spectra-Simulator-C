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
#include "make_grids.h"
//#include "linspace.h"
#include "data.h"
//#include "activity.h"
#include <sys/stat.h>
#include <sys/types.h>

bool do_grids(int argc, char** argv) {
    // Declare the supported options.
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("outdir,d", boost::program_options::value<std::string>()->default_value("./grids/"), "file type (default is ./grids/)")
        ("resol,r", boost::program_options::value<double>()->default_value(M_PI/180.), "resolution of grid (default is Pi/180= 1 Deg)")
        ("lmax", boost::program_options::value<int>()->default_value(3), "maximum l value (default is l=3)")
        ("ftype,f", boost::program_options::value<std::string>()->default_value(""), "Filter type: Use 'gate', 'triangle' or 'gauss'")
        ("verbose,v", boost::program_options::value<bool>()->default_value(true), "show information on the created grid (true) or not (default is true)");
            
    boost::program_options::variables_map vm;
    try {
        boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
        boost::program_options::notify(vm);
        if (vm.count("help")) {
            std::cout << desc << "\n";
            return 0;
        }

        if (vm["ftype"].as<std::string>() != "gate" && vm["ftype"].as<std::string>() != "triangle" && vm["ftype"].as<std::string>() != "gauss") {
            std::cerr << "Error: Invalid filter type specified. ";
            std::cerr << desc << "\n";
            return 1;
        }
    } catch (boost::program_options::error& e) {
        std::cerr << "Error: " << e.what() << "\n";
        std::cerr << desc << "\n";
        return 1;
    }         
        const std::string outdir = vm["outdir"].as<std::string>();
        const double resol = vm["resol"].as<double>();
        const int lmax = vm["lmax"].as<int>();
        const std::string filter_type = vm["ftype"].as<std::string>();
        const bool verbose = vm["verbose"].as<bool>();
        const double theta_min = 0;
        const double theta_max = M_PI/2;
        const double delta_min = 0;
        const double delta_max = M_PI/4 + resol*2; // Added on 27/04/2023: As a margin of security on the edge, we add 2*resol. This because we do some rounding when prioritising the resolution
        
        // Create directory if it does not already exist
        struct stat info;
        if (stat(outdir.c_str(), &info) != 0) {
            if(mkdir(outdir.c_str(), S_IRWXU) !=0){
                std::cerr << "Error: Could not create output directory. ";
                std::cerr << "\n";
                return 1;                
            } //| S_IRWXG | S_IROTH | S_IXOTH);
        } else{
            std::cout << "Warning: The directory already exist. Any previous grid may be overwritten" << std::endl;
        }
        // Generate the grid
        for(int l=1; l<=lmax;l++){
            GridData_Alm grid= make_Alm_grid(l, theta_min, theta_max, delta_min, delta_max, resol, filter_type, verbose);
            // Save it in multiple files within the outdir directory and return an error if it failed
            if (saveAlm(grid, outdir, "A"+std::to_string(l), ".gz")) {
                std::cout << "All subgrids for l= " << l << " saved successfully in the  designated directory: " << outdir << "\n";
            } else {
                std::cerr << "Error: could not save subgrid to " << outdir << "\n";
                return 1;
            }
        }
            std::cout << " All grids from l=1 to lmax=" << lmax << " Done." << std::endl;
    return 0;
}