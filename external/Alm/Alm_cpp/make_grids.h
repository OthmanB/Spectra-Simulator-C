#pragma once 

#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <cmath>
#include <string>
#include <vector>
#include "linspace.h"
#include "data.h"
#include "activity.h"
#include <sys/stat.h>
#include <sys/types.h>

bool saveAlm(const GridData_Alm& grid, const std::string& outdir, const std::string& prefix, const std::string& suffix);
GridData_Alm make_Alm_grid(const int l, const double theta_min, const double theta_max, 
        const double delta_min, const double delta_max, 
        const double resol, const std::string ftype, bool verbose = true);
