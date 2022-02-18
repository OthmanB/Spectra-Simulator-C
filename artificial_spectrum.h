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
//void artificial_spectrum_act_asym(double Tobs, double Cadence, double Nspectra, std::string dir_core, std::string identifier, bool doplots, bool write_inmodel);
void artificial_spectrum_act_asym(const double Tobs, const double Cadence, const double Nspectra, const long Nrealisation, const std::string dir_core, const std::string identifier, const bool doplots, const bool write_inmodel);
void artificial_spectrum_a1a2a3asym(const double Tobs, const double Cadence, const double Nspectra, const long Nrealisation, const std::string dir_core, const std::string identifier, const bool doplots, const bool write_inmodel);
void artificial_spectrum_a1Alma3(const double Tobs, const double Cadence, const double Nspectra, const long Nrealisation, 
								 const std::string dir_core, const std::string identifier, const bool doplots, const bool write_inmodel,
								 const bool domodelfiles, const bool limit_data_range, const std::string modelname);
void artificial_spectrum_aj(const double Tobs, const double Cadence, const double Nspectra, const long Nrealisation, 
								 const std::string dir_core, const std::string identifier, const bool doplots, const bool write_inmodel,
								 const bool domodelfiles, const bool limit_data_range, const std::string modelname);