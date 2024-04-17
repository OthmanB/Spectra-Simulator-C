#pragma once 

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <gsl/gsl_interp2d.h>
#include "string_handler.h"
#include "data.h"

using namespace std;

GridData loadGridData(const string& filename);
double interpolate(const GridData& data, double x, double y);
gsl_interp2d* init_2dgrid(const GridData4gsl& data_flatten);
double interpolate_core(gsl_interp2d* interp, const GridData4gsl& data_flatten, double x, double y);
GridData4gsl flatten_grid(const GridData& data);