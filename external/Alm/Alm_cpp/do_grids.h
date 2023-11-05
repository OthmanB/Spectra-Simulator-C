// Program to make a grid the is compatible with the 2d interpoler program
// see https://github.com/OthmanB/2D_interp
#pragma once 
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

bool do_grids(int argc, char** argv);