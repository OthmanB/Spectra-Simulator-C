/**
 * @file data.h
 * @brief Declarations of the structures  used to decompose and rescale p mode frequencies.
 *
 * This file contains the declarations of the functions used to decompose into the asymptotic terms and rescale frequencies nu(n,l). 
 *
 * @date 24 Feb 2023
 * @author obenomar
 */

#pragma once 
#include <Eigen/Dense>
#include <vector>

using Eigen::VectorXd;
using Eigen::MatrixXd;

/**
 * @struct Data_asympt_p
 * @brief Structure that represents the asymptotic elements.
 *
 * This structure represents the asymptotic elements, including n and O2_term.
 */
struct Data_asympt_p{
	bool error_status; ///< Used to detect if there was a problem when writting the data
	VectorXd n; ///< Array of radial orders n.
	double Dnu; ///< The large Separation
	double epsilon; ///< The phase offset for p modes
	double d0l; ///< The small separation
	VectorXd O2_term;	 ///< Array of O2 term values.
};

/**
 * @struct Freq_modes
 * @brief Structure that represents the frequencies.
 *
 * This structure represents the frequencies, including fl0, fl1, fl2, and fl3.
 */
struct Freq_modes{
	bool error_status; ///< Bool flag used to detect if there was a problem when writting the data
    VectorXd fl0; ///< Frequencies for l=0.
    VectorXd fl1; ///< Frequencies for l=1.
    VectorXd fl2; ///< Frequencies for l=2.
    VectorXd fl3; ///< Frequencies for l=3.
	Data_asympt_p p_asymptotic; ///< Structure containing information  on the asymptotic parameters for p modes
};

/**
 * @struct Data_file
 * @brief Structure that represents a data file.
 *
 * This structure represents a data file, including the error status, header, rescaling flag, target Dnu value,
 * target epsilon value, target small separation, and raw data matrix.
 */
struct Data_file{
	bool error_status; ///< Used to detect if there was a problem when writting the data
	std::vector<std::string> header; ///< Array of string with the header (lines marked with "#" on the top of the files) of the read file
	bool do_rescale; ///< Boolean flag stating if we perform a rescaling or not
	double Dnu_target; ///< Value of the Large separation that is obtained after rescaling
	double epsilon_target; ///< Value of epsilon after the rescaling
	VectorXd d0l_target; ///< Value of the small separation obtained after rescaling
	MatrixXd data; ///< Raw matrix of data
};
