/**
 * @file data.h
 * @brief Header file that contains all kinds of classes/structures used to process and encapsulate data.
 * 
 * This file contains the declarations for various classes and structures used to process and encapsulate data. It provides functionality for storing input data, model inputs and priors, generic data with metadata, parameter headers, model files, derivatives, gnuplot data, and basic data. This is mostly for MCMC analysis.
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
using Eigen::VectorXi;
using Eigen::VectorXd;

/**
 * @struct Data
 * @brief The general structure that contains the input data that has to be analyzed.
 * 
 * This structure is used to hold the input data that needs to be analyzed. It includes vectors for x, y, sigma_y, and other metadata such as xlabel, ylabel, xunit, and yunit.
 */
struct Data{
		VectorXd x; ///< In the case of a 1D fit (e.g. fit of Mass, [Fe/H], ...), this variable will be ignored by the model function.
		VectorXd xrange; ///< The min and max of x... in practice, it is used to limit the data range if requested by the cfg file (e.g. freq_range option in the .model).
		VectorXd y; ///< The y axis values such as the power spectrum.
		VectorXd sigma_y; ///< The standard deviation for the values in the y-axis
		long Nx; ///< Length of the x vector. Must be the same as the length of the y vector.
		std::string xlabel; ///< label for x-axis. In ASCII file, the axis labels are identified by the symbol ! at the begining of the line
		std::string ylabel;  ///<  label for y-axis.
		std::string xunit; ///< Units for the x-axis.
		std::string yunit; ///< Units for the y-axis.
		std::vector<std::string> header; ///< Any kind of information that is useful. e.g. if it is a test, model info + model inputs. Any header (In ASCII, marked by #), will be put there.
};

/**
 * @struct Input_Data
 * @brief A generic structure that helps to encapsulate model inputs and priors.
 * 
 * This structure is used to encapsulate model inputs and priors. It includes vectors for inputs_names, priors_names, priors_names_switch, inputs, relax, priors, plength, and extra_priors.
 */
struct Input_Data{
	std::string model_fullname; ///< The fullname of the model that is going to be processed.
	std::vector<std::string> inputs_names; ///< Names of the parameters of the MCMC fit.
	std::vector<std::string> priors_names; ///< Names of the priors on each parameters.
	VectorXi priors_names_switch; ///< Used instead of priors_names in loops. This allows the use of the switch(x) - case statements instead of if/else.
	VectorXd inputs; ///< Vector of input parameters.
	VectorXi relax; ///< Vector indicating the relaxed (fitted) parameters.
	MatrixXd priors; ///< List of prior parameters on each parameter of the inputs.
	VectorXi plength; ///< Pointer-like vector giving indications on the structure of the inputs parameters (frequencies positions, heights postions,...).
	VectorXd extra_priors; ///< Contains extra parameters that could be used for priors, such as hyper priors.
};

/**
 * @struct Data_Nd
 * @brief A generic structure that helps to encapsulate a Matrix of information along with some metadata.
 * 
 * This structure is used to encapsulate a matrix of information along with some metadata. It includes a matrix for data and vectors for header, labels, and units.
 */
struct Data_Nd{
	MatrixXd data; ///< 2D array in which each columns are expected to contain the values of a given parameter
	std::vector<std::string> header; ///< Header of the data.
	std::vector<std::string> labels; ///< Labels for the Matrix of data.
	std::vector<std::string> units; ///< Units for the Matrix of data.
};

/**
 * @struct Params_hdr
 * @brief Header of the parameters.
 * 
 * This structure is used to hold the header of the parameters. It includes vectors for header, relax, plength, constant_names, constant_values, and variable_names.
 */
struct Params_hdr{
	std::vector<std::string> header; ///< Any comment
	int Nsamples; ///< Total number of samples (no necessarily the actual number of written samples).
	int Nchains; ///< Total number of parallel chains.
 	int Nvars; ///< Total number of variables.
	int Ncons; ///< Total number of constants.
	VectorXi relax; ///< Vector indicating the fitted parameters.
	VectorXi plength;  ///< Pointer-like vector giving indications on the structure of the inputs parameters (frequencies positions, heights postions,...).
	std::vector<std::string> constant_names; ///< The name of the constants.
	VectorXd constant_values; ///< List of constant values.
	std::vector<std::string> variable_names; ///< The name of the variables.
};

/**
 * @struct MCMC_files
 * @brief Structure that keeps information about MCMC files.
 * 
 * This structure is used to hold information about MCMC files. It includes variables for ID, Dnu, numax, C_l, els, freq_range, param_type, freqs_ref, relax_freq, relax_gamma, relax_H, hyper_priors, eigen_params, noise_params, noise_s2, common_names, common_names_priors, and modes_common.
 */
struct MCMC_files{
	std::string ID; ///< Identification number of the star.
	double Dnu; ///< Large Frequency separation. 
	double numax; ///< Frequency at maximum power.
	double C_l; ///< Ordinate at origin for the asymptotic p modes. It is equal to epsilon_p * Dnu.
	VectorXi els; ///< List of degrees
	VectorXd freq_range; ///< Range of frequencies included in the fit.
	std::vector<std::string> param_type; ///< Either p, g or m. Indicates the type of modes.
	VectorXd freqs_ref; ///< List of frequencies.
	std::vector<bool> relax_freq; ///< indicate if a specific frequency is fixed or fitted.
	std::vector<bool> relax_gamma; ///< indicate if a specific mode width is fixed or fitted.
	std::vector<bool> relax_H; ///< indicate if a specific mode height is fixed or fitted.
	VectorXd hyper_priors; ///< Hyper priors information.
	MatrixXd eigen_params; ///< Full list of degree, frequency, width, height for all fitted modes.
	VectorXd noise_params; ///< Noise parameters. Usually a sum of up to 3 Harvey Like profiles plus a white noise. Used as initial guesses.
	MatrixXd noise_s2; ///< Noise parameters obtained by the means of a Gaussian envelope fit. Used as priors.
	std::vector<std::string> common_names; ///< Names of the global parameters, such as the inclination or splitting.
	std::vector<std::string> common_names_priors; ///< Names of the priors used for each global parameter.
	MatrixXd modes_common; ///< Characteristic values for the priors in the global parameters.
};

/**
 * @struct Deriv_out
 * @brief Structure that keeps information of the derivatives.
 *
 * This structure is used to store information about derivatives. It includes vectors for x-axis derivatives, derivatives, and errors.
 */
struct Deriv_out{
	VectorXd xderiv; ///< x positions where derivatives are evaluated.
	VectorXd deriv; ///< Derivatives.
	VectorXd error; ///< Error on the derivatives.
};

/**
 * @struct gnuplt_Data
 * @brief An encapsulator for data when plotting with gnuplot-iostream.h.
 *
 * This structure is used as an encapsulator for data when plotting with gnuplot-iostream.h. It includes variables for x-axis value, y-axis series 1, y-axis series 2, and y-axis series 3.
 */
struct gnuplt_Data {
    double x;  ///< x axis value
    double y1; ///< y axis series 1
    double y2;             ///< y axis series 2
    double y3;             ///< y axis series 3
};

/**
 * @typedef gnuplt_Dataset
 * @brief A typedef for a vector of gnuplt_Data.
 */
typedef std::vector<gnuplt_Data> gnuplt_Dataset;

namespace gnuplotio {
    template<>
    struct TextSender<gnuplt_Data> {
        /**
         * @brief Sends gnuplt_Data to the output stream.
         *
         * This function sends gnuplt_Data to the output stream.
         *
         * @param stream The output stream.
         * @param v The gnuplt_Data object.
         */
        static void send(std::ostream &stream, const gnuplt_Data &v) {
            stream << " ";
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
 * @struct Data_Basic
 * @brief A basic structure that stores an array of strings and a vector of integers.
 *
 * This structure is used to store an array of strings and a vector of integers. It includes variables for the array of strings and the vector of integers.
 */
struct Data_Basic{
	std::vector<std::string> strarr; ///< Any comment
	VectorXi vecXi; ///< Case number
};
