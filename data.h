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
    std::vector<double> val_min; ///< Minimum values for the parameters.
    std::vector<double> val_max; ///< Maximum values for the parameters.
    std::vector<double> step; ///< Step values for the parameters.
    std::vector<double> forest_params; ///< Forest parameters.
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
