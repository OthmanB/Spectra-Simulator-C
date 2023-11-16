/**
 * @file artificial_spectrum.h
 *
 * Contains all kind of methods that is producing
 * a single artificial spectrum
 * 
 * @date 05 May 2016
 * @author obenomar
 */
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


/**
 * @brief Generate an artificial spectrum with asymmetric Lorentzian modes
 *
 * This function generates an artificial spectrum with asymmetric Lorentzian modes based on the given parameters.
 *
 * @param Tobs The observation time in days.
 * @param Cadence The cadence of the observations in seconds.
 * @param Nspectra Defines how much averaged spectra we should use. This changes the noise statistics to khi(2,2*Nspectra)
 * @param Nrealisation The number of realizations of noise to generate.
 * @param dir_core The directory path for the core files.
 * @param identifier The identifier for the spectrum.
 * @param doplots Flag indicating whether to generate plots of the spectrum.
 * @param write_inmodel Flag indicating whether to write the file with model information.
 */
void artificial_spectrum_act_asym(const double Tobs, const double Cadence, const double Nspectra, const long Nrealisation, const std::string dir_core, const std::string identifier, const bool doplots, const bool write_inmodel);

/**
 * @brief Generate an artificial spectrum with a1, a2, a3 coefficient and asymetrical Lorentzian modes
 *
 * This function generates an artificial spectrum with a1, a2, a3 coefficient and asymetrical Lorentzian modes based on the given parameters.
 *
 * @param Tobs The observation time in days.
 * @param Cadence The cadence of the observations in seconds.
 * @param Nspectra Defines how much averaged spectra we should use. This changes the noise statistics to khi(2,2*Nspectra)
 * @param Nrealisation The number of realizations of noise to generate.
 * @param dir_core The directory path for the core files.
 * @param identifier The identifier for the spectrum.
 * @param doplots Flag indicating whether to generate plots of the spectrum.
 * @param write_inmodel Flag indicating whether to write the file with model information.
 */
void artificial_spectrum_a1a2a3asym(const double Tobs, const double Cadence, const double Nspectra, const long Nrealisation, const std::string dir_core, const std::string identifier, const bool doplots, const bool write_inmodel);

/**
 * @brief Generate an artificial spectrum with a1, a3 coefficient and Alm activity and asymetrical Lorentzian modes
 *
 * This function generates an artificial spectrum with a1, a2, a3 coefficient and Alm activity and asymetrical Lorentzian modes based on the given parameters.
 *
 * @param Tobs The observation time in days.
 * @param Cadence The cadence of the observations in seconds.
 * @param Nspectra Defines how much averaged spectra we should use. This changes the noise statistics to khi(2,2*Nspectra)
 * @param Nrealisation The number of realizations of noise to generate.
 * @param dir_core The directory path for the core files.
 * @param identifier The identifier for the spectrum.
 * @param doplots Flag indicating whether to generate plots of the spectrum.
 * @param write_inmodel Flag indicating whether to write the file with model information.
 * @param domodelfiles Flag indicating whether to write .model files suitable for the TAMCMC code
 * @param limit_data_range Flag indicating wheter to truncate the spectrum by focusing in the region with modes. Save space if on. But not suitable for ML algorithms.
 * @param modelname If domodelfiles is true, use this to define the name of the model that has to be fitted.
 */
void artificial_spectrum_a1Alma3(const double Tobs, const double Cadence, const double Nspectra, const long Nrealisation, 
								 const std::string dir_core, const std::string identifier, const bool doplots, const bool write_inmodel,
								 const bool domodelfiles, const bool limit_data_range, const std::string modelname);

/**
 * @brief Generate an artificial spectrum with aj coefficient with j={1,2,3,4,5,6} and asymetrical Lorentzian modes
 *
 * This function generates an artificial spectrum with a1, a2, a3 coefficient and Alm activity and asymetrical Lorentzian modes based on the given parameters.
 *
 * @param Tobs The observation time in days.
 * @param Cadence The cadence of the observations in seconds.
 * @param Nspectra Defines how much averaged spectra we should use. This changes the noise statistics to khi(2,2*Nspectra)
 * @param Nrealisation The number of realizations of noise to generate.
 * @param dir_core The directory path for the core files.
 * @param identifier The identifier for the spectrum.
 * @param doplots Flag indicating whether to generate plots of the spectrum.
 * @param write_inmodel Flag indicating whether to write the file with model information.
 * @param domodelfiles Flag indicating whether to write .model files suitable for the TAMCMC code
 * @param limit_data_range Flag indicating wheter to truncate the spectrum by focusing in the region with modes. Save space if on. But not suitable for ML algorithms.
 * @param modelname If domodelfiles is true, use this to define the name of the model that has to be fitted.
 * @param noise_modelname Defines the noise model that has to be used. By default, it is "harvey_1985"
 */
void artificial_spectrum_aj(const double Tobs, const double Cadence, const double Nspectra, const long Nrealisation, 
								 const std::string dir_core, const std::string identifier, const bool doplots, const bool write_inmodel,
								 const bool domodelfiles, const bool limit_data_range, const std::string modelname, const std::string noise_modelname="harvey_1985");