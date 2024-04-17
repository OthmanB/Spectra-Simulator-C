/**
 * @file data.h
 * @brief Header file that contains all kinds of classes/structures used to process and/or encapsulate data.
 * 
 * This file contains the declarations of various classes and structures that are used to process and encapsulate data.
 * 
 * @date 22 Feb 2016
 * @author obenomar
 */

#pragma once
#include <Eigen/Dense>
#include <string>
#include <vector>
#include "gnuplot-iostream.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * @struct Data_Nd
 * @brief Structure that represents a 2D array of data.
 * 
 * This structure represents a 2D array of data, where each column contains the values of a given parameter.
 */
struct Data_Nd {
    MatrixXd data; ///< 2D array in which each column contains the values of a given parameter.
    std::vector<std::string> header; ///< Header of the data.
    std::vector<std::string> labels; ///< Labels for each column of data.
    std::vector<std::string> units; ///< Units for each column of data.
};

/**
 * @struct Model_data
 * @brief Structure that represents data for a single model.
 * 
 * This structure represents the data for a single model from a read text file. It contains the parameters, frequencies, and their labels.
 */
struct Model_data {
    VectorXd params; ///< 1D array with the data for one single model from the read text file.
    std::vector<std::string> labels_params; ///< Name of the parameters within params.
    MatrixXd freqs; ///< Table of frequencies with (l,n).
    std::vector<std::string> labels_freqs; ///< Labels for each frequency in freqs.
};

/**
 * @struct gnuplt_Data
 * @brief Structure used for data encapsulation when plotting with gnuplot-iostream.h.
 * 
 * This structure is used as an encapsulator for data when plotting with gnuplot-iostream.h. It contains x and y values for up to three series.
 */
struct gnuplt_Data {
    double x; ///< x-axis value.
    double y1; ///< y-axis value for series 1.
    double y2; ///< y-axis value for series 2.
    double y3; ///< y-axis value for series 3.
};

/**
 * @struct Star_params
 * @brief Data structure used by models that need to read parameters of a reference star from a .in file.
 * 
 * This data structure is used by models that need to read the parameters of a reference star from a .in file to operate.
 */
struct Star_params {
    VectorXd spec_params; ///< Parameters of the reference star.
    MatrixXd mode_params; ///< Mode parameters.
    MatrixXd noise_params; ///< Noise parameters.
    std::string identifier; ///< Identifier of the reference star.
    std::string noise_model="Harvey-like"; ///< Type of noise model used in the star. Could be Harvey-like or Kallinger+2014
    double cadence; ///< Cadence of the observations.
    double Tobs; ///< Total observation time.
};

/**
 * @struct Config_Data
 * @brief Encapsulates the parameters used to generate the spectra.
 * 
 * This structure encapsulates the parameters used to generate the spectra.
 */
struct Config_Data {
    bool erase_old_files; ///< Flag indicating whether to erase older files of the same name.
    bool doplots; ///< Flag indicating whether to generate plots.
    bool write_inmodel; ///< Flag indicating whether to include the input model in the output ASCII file.
    bool limit_data_range; ///< Flag indicating whether to limit the data range.
    bool do_modelfiles; ///< Flag indicating whether to generate model files.
    std::string modefile_modelname; ///< Model name written within the model file.
    std::string model_name; ///< Model name within this program.
    std::string forest_type; ///< Type of forest.
    std::string extra_params; ///< Extra parameters encoded in a string.
    std::vector<std::string> template_files; ///< List of template files used to determine height and widths.
    double Tobs; ///< Total observation time.
    double Cadence; ///< Cadence of the observations.
    double Nspectra; ///< Number of spectra.
    double Nrealisation; ///< Number of realizations.
    std::vector<std::string> labels; ///< Labels for the parameters.
    std::vector<std::string> distrib; ///< Distribution for the variables. Either Uniform, Gaussian or Fix.
    std::vector<double> val_min; ///< Minimum values for the parameters.
    std::vector<double> val_max; ///< Maximum values for the parameters.
    std::vector<double> step; ///< Step values for the parameters.
    std::vector<double> forest_params; ///< Forest parameters.
};

/**
 * @brief Noise configuration used for Kallinger+2014 noise models
 * 
 * This structure is to encapsulate the data read from a noise configuration file.
 * It is used in models that involve a Kallinger+2014 noise parametrisation. This parametrisation has specificities (a large set of parameters, with Gaussian distribution) that does not make it easy to read if put directly in the main.cfg file. Furthermore, it avoids the user to make too much errors in setting the noise, as the values are data-driven and should not be modified
 */
struct Config_Noise {
    // For the random case
    std::vector<std::string> name_random; ///< Name of the variable in the random case.
    std::vector<std::string> distrib_random; ///< Name of the distribution from which we draw random variables in the grid case.
    std::vector<double> x1_random; ///< First value associated to the distribution (eg. the mean of a Gaussian).
    std::vector<double> x2_random; ///< Second value associated to the distribution (eg. the standard deviation of a Gaussian).
    std::vector<double> kerror_random; ///< Enlargement factor in the case of a Gaussian distribution.
    // For the grid case
    std::vector<std::string> name_grid; ///< Name of the variable in the grid case.
    std::vector<std::string> distrib_grid; ///< Name of the distribution from which we draw random variables in the grid case.
    std::vector<double> x1_grid; ///< First value associated to the distribution in the grid case (eg. the min of the range of a Uniform distribution).
    std::vector<double> x2_grid; ///< Second value associated to the distribution in the grid case (eg. the max of the range of a Uniform distribution).
    std::vector<double> kerror_grid; ///< Enlargement factor or any relevant value required to compute the distribution.
};

/**
 * @typedef gnuplt_Dataset
 * @brief Typedef for a vector of gnuplt_Data structures.
 */
typedef std::vector<gnuplt_Data> gnuplt_Dataset;

namespace gnuplotio {
    /**
     * @brief Template specialization for sending gnuplt_Data to an output stream.
     */
    template<>
    struct TextSender<gnuplt_Data> {
        /**
         * @brief Sends gnuplt_Data to an output stream.
         * @param stream The output stream.
         * @param v The gnuplt_Data structure to send.
         */
        static void send(std::ostream &stream, const gnuplt_Data &v) {
            TextSender<double>::send(stream, v.x);
            stream << " ";
            TextSender<double>::send(stream, v.y1);
            stream << " ";
            TextSender<double>::send(stream, v.y2);
            stream << " ";
            TextSender<float>::send(stream, v.y3);
        }
    };
}

/**
 * @struct Deriv_out
 * @brief Structure that keeps information about derivatives.
 * 
 * This structure keeps information about derivatives, including the x values, derivative values, and error values.
 */
struct Deriv_out {
    VectorXd xderiv; ///< x values for the derivatives.
    VectorXd deriv; ///< Derivative values.
    VectorXd error; ///< Error values.
};
