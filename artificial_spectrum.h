#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/Dense>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include "io_star_params.h"
#include "build_lorentzian.h"
#include "plots_diags.h"
#include <string>
void artificial_spectrum_act_asym(double Tobs, double Cadence, double Nspectra, std::string dir_core, std::string identifier, bool doplots, bool write_inmodel);
