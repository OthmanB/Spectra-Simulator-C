/**
 * @file io_star_params.h
 * @brief Methods to handle IO in the program.
 * 
 * Contains all kind of methods used to process and/or read and write files for configuration and output files
 * 
 * @date  22 Feb 2016
 * @author obenomar
 */ 

# pragma once
# include <string>
# include <vector>
# include <Eigen/Dense>
# include "data.h" // contains the structure Data
# include "external/ARMM/data_solver.h" // contains the structure Data
# include <fstream>
# include <string>
# include "ioproc.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

/**
 * @brief Remove comments from a string
 *
 * This function removes anything after the specified terminator character in the input string, including the terminator itself. It also removes any leading or trailing whitespace from the string.
 *
 * @param str0 The input string
 * @param terminator The character indicating the end of the line
 * @return The modified string with comments removed and leading/trailing whitespace removed
 */
std::string rem_comments(std::string str0, std::string terminator);

/**
 * @brief Write global information to the output file
 *
 * This function writes the global information, including the identifier and spectrum parameters, to the specified output file.
 *
 * @param spec_params The spectrum parameters to be written to the file
 * @param file_out The name of the output file
 * @param identifier The identifier to be written to the file
 * @param append Boolean value indicating whether to append to an existing file or create a new file
 */
void write_global_info(VectorXd& spec_params, std::string file_out, std::string identifier, bool append=false); // Write Identifier, Cadence and observation duration


/**
 * @brief Write the mode parameters matrix to an output file in a specific format
 *
 * This function writes the mode parameters to an output file in a specific format. It allows appending the data to an existing file or creating a new file. The function includes a header indicating the purpose of the file and the format of the input mode parameters. The mode parameters matrix should contain degree, frequency, H, W, a1, eta0, epsilon_Alm, theta0_Alm, delta_Alm, a3, asymmetry, and inclination values. The function uses the Nchars and precision vectors to specify the width and precision of each value in the output file.
 *
 * @param spec_params The spectrum parameters vector
 * @param mode_params The mode parameters matrix
 * @param noise_params The noise parameters matrix
 * @param file_out The path to the output file
 * @param identifier The identifier of the star
 */
void write_star_params_act_asym(VectorXd& spec_params, MatrixXd& mode_params, MatrixXd& noise_params, std::string file_out, std::string identifier);

/**
 * @brief Write the mode parameters matrix to an output file in a specific format
 *
 * This function writes the mode parameters matrix to an output file in a specific format. It allows appending the data to an existing file or creating a new file. The function includes a header indicating the purpose of the file and the format of the input mode parameters. The mode parameters matrix should contain degree, frequency, H, W, a1, eta0, epsilon_Alm, theta0_Alm, delta_Alm, a3, asymmetry, and inclination values. The function uses the Nchars and precision vectors to specify the width and precision of each value in the output file.
 *
 * @param mode_params The mode parameters matrix
 * @param file_out The path to the output file
 * @param append Flag indicating whether to append the data to an existing file (default: false)
 */
void write_star_mode_params_a1a2a3(MatrixXd& mode_params, std::string file_out, bool append=false);

/**
 * @brief Write the star parameters (mode parameters, noise parameters, and global info) to an output file in a specific format
 *
 * This function writes the star parameters (mode parameters, noise parameters, and global info) to an output file in a specific format. It allows appending the data to an existing file or creating a new file. The function includes a header indicating the purpose of the file and the format of the input parameters. The mode parameters matrix should contain degree, frequency, H, W, a1, a2, a3, a4, a5, a6, asymmetry, and inclination values. The noise parameters matrix should contain H0, tau_0, p0, H1, tau_1, p1, and N0 values. The spec_params vector should contain the global info of the star. The function uses the Nchars and precision vectors to specify the width and precision of each value in the output file.
 *
 * @param mode_params The mode parameters matrix
 * @param file_out The path to the output file
 * @param append Flag indicating whether to append the data to an existing file (default: false)
 */
void write_star_mode_params_act_asym(MatrixXd& mode_params, std::string file_out, bool append=false);

/**
 * @brief Write the star parameters (mode parameters, noise parameters, and global info) to an output file in a specific format
 *
 * This function writes the star parameters (mode parameters, noise parameters, and global info) to an output file in a specific format. It allows appending the data to an existing file or creating a new file. The function includes a header indicating the purpose of the file and the format of the input parameters. The mode parameters matrix should contain degree, frequency, H, W, a1, eta0, epsilon_Alm, theta0_Alm, delta_Alm, a3, asymmetry, and inclination values. The noise parameters matrix should contain H0, tau_0, p0, H1, tau_1, p1, and N0 values. The spec_params vector should contain the global info of the star. The function uses the Nchars and precision vectors to specify the width and precision of each value in the output file.
 *
 * @param spec_params The global info vector
 * @param mode_params The mode parameters matrix
 * @param noise_params The noise parameters matrix
 * @param file_out The path to the output file
 * @param identifier The identifier of the star
 */
void write_star_params_Alm(VectorXd& spec_params, MatrixXd& mode_params, MatrixXd& noise_params, std::string file_out, std::string identifier);

/**
 * @brief Write the star parameters (mode parameters, noise parameters, and global info) to an output file in a specific format
 *
 * This function writes the star parameters (mode parameters, noise parameters, and global info) to an output file in a specific format. It allows appending the data to an existing file or creating a new file. The function includes a header indicating the purpose of the file and the format of the input parameters. The mode parameters matrix should contain degree, frequency, H, W, a1, a2, a3, a4, a5, a6, asymmetry, and inclination values. The noise parameters matrix should contain H0, tau_0, p0, H1, tau_1, p1, and N0 values. The spec_params vector should contain the global info of the star. The function uses the Nchars and precision vectors to specify the width and precision of each value in the output file.
 *
 * @param spec_params The global info vector
 * @param mode_params The mode parameters matrix
 * @param noise_params The noise parameters matrix
 * @param file_out The path to the output file
 * @param identifier The identifier of the star
 */
void write_star_params_aj(VectorXd& spec_params, MatrixXd& mode_params, MatrixXd& noise_params, std::string file_out, std::string identifier);

/**
 * @brief Write the mode parameters matrix to an output file in a specific format
 *
 * This function writes the mode parameters matrix to an output file in a specific format. It allows appending the data to an existing file or creating a new file. The function includes a header indicating the purpose of the file and the format of the input mode parameters. The mode parameters matrix should contain degree, frequency, H, W, a1, eta0, epsilon_Alm, theta0_Alm, delta_Alm, a3, asymmetry, and inclination values. The function uses the Nchars and precision vectors to specify the width and precision of each value in the output file.
 *
 * @param mode_params The mode parameters matrix
 * @param file_out The path to the output file
 * @param append Flag indicating whether to append the data to an existing file (default: false)
 */
void write_star_mode_params_Alm(MatrixXd& mode_params, std::string file_out, bool append=false);

/**
 * @brief Write the mode parameters matrix to an output file in a specific format
 *
 * This function writes the mode parameters matrix to an output file in a specific format. It allows appending the data to an existing file or creating a new file. The function includes a header indicating the purpose of the file and the format of the input mode parameters. The mode parameters matrix should contain degree, frequency, H, W, a1, a2, a3, a4, a5, a6, asymmetry, and inclination values. The function uses the Nchars and precision vectors to specify the width and precision of each value in the output file.
 *
 * @param mode_params The mode parameters matrix
 * @param file_out The path to the output file
 * @param append Flag indicating whether to append the data to an existing file (default: false)
 */
void write_star_mode_params_aj(MatrixXd& mode_params, std::string file_out, bool append=false);

/**
 * @brief Write the noise parameters matrix to an output file in a specific format
 *
 * This function writes the noise parameters matrix to an output file in a specific format. It allows appending the data to an existing file or creating a new file. The function includes a header indicating the purpose of the file and the format of the input noise parameters. The noise parameters matrix should contain H0, tau_0, p0, H1, tau_1, p1, and N0 values. The function uses the Nchars and precision vectors to specify the width and precision of each value in the output file.
 *
 * @param noise_params The noise parameters matrix
 * @param file_out The path to the output file
 * @param append Flag indicating whether to append the data to an existing file (default: false)
 */
void write_star_noise_params(MatrixXd& noise_params, std::string file_out, bool append=false);

/**
 * @brief Write the range of relevant modes to an output file
 *
 * This function writes the range of relevant modes to an output file. It includes the numax, minimum frequency, maximum frequency, and nmax_star values. The function uses the Nchars and precision variables to specify the width and precision of each value in the output file. The output file is created or overwritten.
 *
 * @param cfg_star The configuration of the synthetic star
 * @param params The parameters of the synthetic star
 * @param output_file The path to the output file
 */
void write_range_modes(Cfg_synthetic_star cfg_star, Params_synthetic_star params, std::string output_file);

/**
 * @brief Write the spectrum data to an output file
 *
 * This function writes the spectrum data to an output file. It includes the frequency and spectrum values. Optionally, it can also include the model spectrum values. The function allows specifying a valid range of frequencies to write to the file.
 *
 * @param x The frequency vector
 * @param y The spectrum vector
 * @param z The model spectrum vector
 * @param file_out The path to the output file
 * @param write_inmodel Flag indicating whether to include the model spectrum in the output file
 */
void write_spectrum(const VectorXd& x, const VectorXd& y, const VectorXd& z, std::string file_out, bool write_inmodel); // Alias of the write_spectrum() below, but with fmin=-1, fmax=-1 (de facto optional parameterss)

/**
 * @brief Write the spectrum data to an output file
 *
 * This function writes the spectrum data to an output file. It includes the frequency and spectrum values. Optionally, it can also include the model spectrum values. The function allows specifying a valid range of frequencies to write to the file.
 *
 * @param x The frequency vector
 * @param y The spectrum vector
 * @param z The model spectrum vector
 * @param file_out The path to the output file
 * @param write_inmodel Flag indicating whether to include the model spectrum in the output file (default: false)
 * @param fmin The minimum frequency value to write (default: -1)
 * @param fmax The maximum frequency value to write (default: -1)
 */
void write_spectrum(const VectorXd& x, const VectorXd& y, const VectorXd& z, const std::string file_out, const bool write_inmodel, const double fmin, const double fmax); 

/**
 * @brief Write the spectrum data to an output file
 *
 * This function writes the spectrum data to an output file. It adds two additional columns with values for the smoothed vector y using a boxcar with different smoothing coefficients. The output file includes the frequency, spectrum, smoothed spectrum with scoef1, smoothed spectrum with scoef2, and model spectrum.
 *
 * @param x The frequency vector
 * @param y The spectrum vector
 * @param z The model spectrum vector
 * @param scoef1 The smoothing coefficient for the first smoothed spectrum
 * @param scoef2 The smoothing coefficient for the second smoothed spectrum
 * @param file_out The path to the output file
 */
void write_spectrum(const VectorXd& x, const VectorXd& y, const VectorXd& z, const double scoef1, double scoef2, const std::string file_out);

/**
 * @brief Write the star model to an output file
 *
 * This function writes the star model to an output file in a specific format. It calculates the Dnu and C_l values from the list of l=0 frequencies. The function also includes the mode parameters, noise parameters, and hyper priors in the output file. It allows specifying a common template directory and a specific model name.
 *
 * @param mode_params The mode parameters matrix
 * @param noise_params The noise parameters matrix
 * @param file_out The path to the output file
 * @param identifier The identifier of the star
 * @param modelname The name of the model
 * @param common_template_dir The directory of the common template
 *
 * @return The mode range vector
 */
VectorXd write_star_model(const MatrixXd& mode_params, const MatrixXd& noise_params, const std::string file_out, 
                          const std::string identifier, const std::string modelname, const std::string common_template_dir);

/**
 * @brief Convert parameters to MatrixXd format
 *
 * This function converts the given parameters to a MatrixXd format. The resulting MatrixXd object contains the mode parameters for each mode specified in the parameters.
 *
 * @param params The parameters containing the mode information
 * @param inc The inclination value
 * @return The mode parameters in MatrixXd format
 */
MatrixXd bumpoutputs_2_MatrixXd(Params_synthetic_star params, double inc);

/**
 * @brief Read the main configuration file
 *
 * This function reads the main configuration file and extracts the relevant data to populate a Config_Data object.
 *
 * @param cfg_file The path to the main configuration file
 * @return The populated Config_Data object
 */
Config_Data read_main_cfg(std::string cfg_file);			


/**
 * @brief Read the noise configuration file used eg. by model with a Kallinger+2014 noise 
 *
* @param filename file name of the configuration file
 * @return a NoiseConfig structure that contains all of the noise parameters
 *
 */
Config_Noise readNoiseConfigFile(const std::string& filename);

/**
 * @brief Read the data from an ASCII file with a specified number of columns
 *
 * This function reads an input file and extracts the data, header, labels, and units. The file can contain a header, labels, and units indicated by specific characters at the beginning of the lines. The number of columns can be specified.
 *
 * @param file_in_name The path to the input file
 * @param delimiter The delimiter used to separate columns in the file
 * @param verbose_data Flag indicating whether to display additional information about the data file
 * @return The extracted data, header, labels, and units in a Data_Nd object
 */
Data_Nd read_data_ascii_Ncols(const std::string file_in_name, const std::string delimiter, const bool verbose_data);

/**
 * @brief Read the star parameters from an input file
 *
 * This function reads an input file and extracts the star parameters, including the identifier, observation duration, cadence, mode parameters, and noise parameters.
 *
 * @param file_in_name The path to the input file
 * @param verbose_data (Optional) If true (default value), then show the read data. Otherwise, only show errors  
 * @return The extracted star parameters in a Star_params object
 */
Star_params read_star_params(const std::string file_in_name, const bool verbose_data=true);

/**
 * @brief Handle file read error
 *
 * This function prints an error message when unable to open a file for reading and exits the program.
 *
 * @param file_out The name of the file that could not be opened
 */
void file_read_error(std::string);

//void convert_in_to_model(const std::string file_out, const std::string identifier);


/**
 * @brief Return a boxcar smooth of a vector
 *
 * This function performs a boxcar smoothing operation on a vector. The smoothing coefficient is in the natural unit of the vector, such as microHz instead of bins.
 *
 * @param in The input vector to be smoothed
 * @param scoef The smoothing coefficient in the natural unit of the vector
 * @return The smoothed vector
 */
VectorXd smooth(const VectorXd& in, const double scoef);

/**
 * @brief Read the entire content of a file as a single string
 *
 * This function reads the entire content of a file and returns it as a single string. It removes any leading or trailing whitespace from each line.
 *
 * @param file The path to the input file
 * @return The content of the file as a single string
 */
std::string read_allfile(const std::string file);

/**
 * @brief Read the entire content of a file as a vector of strings
 *
 * This function reads the entire content of a file and returns it as a vector of strings. Each line of the file is stored as a separate element in the vector. It removes any leading or trailing whitespace from each line.
 *
 * @param file The path to the input file
 * @return The content of the file as a vector of strings
 */
std::vector<std::string> read_allfile_vect(const std::string file);

/**
 * @brief Read the theoretical frequencies from an input file
 *
 * This function reads an input file and extracts the theoretical frequencies for a synthetic star. It returns the frequencies organized in a Cfg_synthetic_star object.
 *
 * @param file The path to the input file
 * @param critical (Optional) If set to true (default) any failure to read the file causes exit. Otherwise, it returns an empty structure
 * @return The extracted theoretical frequencies in a Cfg_synthetic_star object
 */
Cfg_synthetic_star read_theoretical_freqs(const std::string file, const bool critical=true);
