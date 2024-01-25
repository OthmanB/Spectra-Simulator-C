/**
 * @file models_database.h
 * @brief Header file that contains models
 *
 * Header file that contains all kind of methods used to generate models for the pulsation/noise
 * 
 *
 * @date 20 Apr 2016
 * @author obenomar
 */

#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/Dense>
# include <string>
#include "io_star_params.h"
#include "external/rescale/rescale_freqs.h"
#include "external/rescale/decompose_nu.h"
#include "external/rescale/data.h"



/**
 * @brief Generate a configuration file for a simulated star based on input parameters.
 *
 * This function generates a configuration file for a simulated star based on the input parameters. 
 * It initializes the random number generator, sets up constants, and deploys the parameters for the simulation.
 * It then evaluates the mixed modes and injects a Gaussian random error.
 * The function determines the isotropic inclination and generates the mode profiles and frequencies.
 * Finally, it writes the mode parameters and noise parameters to the specified output files.
 *
 * @param input_params A VectorXd object containing the input parameters for the simulation.
 * @param file_out_modes The file path to write the mode parameters of the simulated star.
 * @param file_out_noise The file path to write the noise parameters of the simulated star.
 * @param file_cfg_mm The file path for the configuration file of the simulated star.
 * @param external_path The external path for the ARMM-solver.
 * @param template_file The template file for the simulated star.
 */
void asymptotic_mm_v1(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path, std::string template_file);

/**
 * @brief Generate synthetic mode parameters and noise for an asymptotic star.
 * 
 * @param input_params The input parameters for the star model.
 * @param file_out_modes The output file path for the mode parameters.
 * @param file_out_noise The output file path for the noise parameters.
 * @param file_cfg_mm The configuration file path for the star model.
 * @param external_path The external path for the star model.
 * @param template_file The template file path for the star model.
 */
void asymptotic_mm_v2(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path, std::string template_file);

/**
 * @brief Generate synthetic mode parameters and noise for an asymptotic star.
 * 
 * @param input_params The input parameters for the star model.
 * @param file_out_modes The output file path for the mode parameters.
 * @param file_out_noise The output file path for the noise parameters.
 * @param file_cfg_mm The configuration file path for the star model.
 * @param external_path The external path for the star model.
 * @param template_file The template file path for the star model.
 */
void asymptotic_mm_v3(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path, std::string template_file);

/**
 * @brief Generate synthetic mode parameters and noise for an asymptotic star.
 *
 * @param input_params The input parameters for the star model.
 * @param file_out_modes The output file path for the mode parameters.
 * @param file_out_noise The output file path for the noise parameters.
 * @param file_cfg_mm The configuration file path for the star model.
 * @param external_path The external path for the star model.
 * @param template_file The template file path for the star model.
 */
void asymptotic_mm_freeDp_numaxspread_curvepmodes_v1(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path, std::string template_file);

/**
 * @brief Generate synthetic mode parameters and noise for an asymptotic star.
 *
 * @param input_params The input parameters for the star model.
 * @param file_out_modes The output file path for the mode parameters.
 * @param file_out_noise The output file path for the noise parameters.
 * @param file_cfg_mm The configuration file path for the star model.
 * @param external_path The external path for the star model.
 * @param template_file The template file path for the star model.
 */
void asymptotic_mm_freeDp_numaxspread_curvepmodes_v2(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path, std::string template_file);

/**
 * @brief Generate synthetic mode parameters and noise for an asymptotic star.
 *
 * @param input_params The input parameters for the star model.
 * @param file_out_modes The output file path for the mode parameters.
 * @param file_out_noise The output file path for the noise parameters.
 * @param file_cfg_mm The configuration file path for the star model.
 * @param external_path The external path for the star model.
 * @param template_file The template file path for the star model.
 */
void asymptotic_mm_freeDp_numaxspread_curvepmodes_v3(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path, std::string template_file);

/**
 * @brief Generate synthetic mode parameters and noise for an asymptotic star.
 *
 * @param input_params The input parameters for the star model.
 * @param file_out_modes The output file path for the mode parameters.
 * @param file_out_noise The output file path for the noise parameters.
 * @param file_cfg_mm The configuration file path for the star model.
 * @param external_path The external path for the star model.
 * @param template_file The template file path for the star model.
 */
void asymptotic_mm_freeDp_numaxspread_curvepmodes_v3_GRANscaled_Kallinger2014(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string file_cfg_mm, std::string external_path, std::string template_file);

/**
 * @brief Generate configuration file for an asymptotic star with active asymmetry and Gaussian height profile.
 *
 * @param input_params The input parameters for the star model.
 * @param file_out_modes The output file path for the mode parameters.
 * @param file_out_noise The output file path for the noise parameters.
 *
 * This function generates a configuration file for an asymptotic star with active asymmetry and Gaussian height profile.
 * It takes the following parameters:
 * - input_params: The input parameters for the star model.
 * - file_out_modes: The output file path for the mode parameters.
 * - file_out_noise: The output file path for the noise parameters.
 *
 * The function calculates various constants and deploys the input parameters. It then generates a list of frequencies,
 * heights, widths, splitting, centrifugal terms, latitudinal terms, and stellar inclination. Finally, it writes the mode
 * parameters and noise parameters to the output files.
 */
void generate_cfg_asymptotic_act_asym_Hgauss(VectorXd input_params, std::string file_out_modes, std::string file_out_noise);


/**
 * @brief Generate configuration file for an asymptotic star with active asymmetry and Gaussian height profile.
 *
 * @param input_params The input parameters for the star model.
 * @param file_out_modes The output file path for the mode parameters.
 * @param file_out_noise The output file path for the noise parameters.
 * @param extra The extra parameter for the star model.
 *
 * This function uses a reference star as a template to generate frequencies and width, height profiles. It can be rescaled so that you can modify the HNR but keep the same height profile. Note that the user here provides a target a1/Width so that a1 is automatically adjusted to match the requested a1/Width. The code will not change the Width so that code is not adapted to test blending between adjacent l modes, such as the l=0 and l=2 mode blending.
 *
 * The function takes the following parameters:
 * - input_params: The input parameters for the star model.
 * - file_out_modes: The output file path for the mode parameters.
 * - file_out_noise: The output file path for the noise parameters.
 * - extra: The extra parameter for the star model.
 *
 * The function calculates various constants and deploys the input parameters. It then generates a list of frequencies, heights, widths, splitting, centrifugal terms, latitudinal terms, and stellar inclination. Finally, it writes the mode parameters and noise parameters to the output files.
 */
void generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string extra);


/**
 * @brief Generate configuration file for an asymptotic star with active asymmetry and Gaussian height profile.
 *
 * @param input_params The input parameters for the star model.
 * @param file_out_modes The output file path for the mode parameters.
 * @param file_out_noise The output file path for the noise parameters.
 * @param extra The extra parameter for the star model.
 *
 * This function uses a reference star as a template to generate frequencies and width, height profiles. It can be rescaled so that you can modify the HNR but keep the same height profile. Note that the user here provides a target a1/Width so that a1 is automatically adjusted to match the requested a1/Width. The code will not change the Width so that code is not adapted to test blending between adjacent l modes, such as the l=0 and l=2 mode blending.
 *
 * The function takes the following parameters:
 * - input_params: The input parameters for the star model.
 * - file_out_modes: The output file path for the mode parameters.
 * - file_out_noise: The output file path for the noise parameters.
 * - extra: The extra parameter for the star model.
 *
 * The function calculates various constants and deploys the input parameters. It then generates a list of frequencies, heights, widths, splitting, centrifugal terms, latitudinal terms, and stellar inclination. Finally, it writes the mode parameters and noise parameters to the output files.
 */
void generate_cfg_from_synthese_file_Wscaled_a1a2a3asymovGamma(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string extra);


/**
 * @brief Generate configuration file for an asymptotic star with active asymmetry and Gaussian height profile.
 *
 * @param input_params The input parameters for the star model.
 * @param file_out_modes The output file path for the mode parameters.
 * @param file_out_noise The output file path for the noise parameters.
 * @param extra The extra parameter for the star model.
 *
 * This function uses a reference star as a template to generate frequencies and width, height profiles. It can be rescaled so that you can modify the HNR but keep the same height profile. Note that the user here provides a target a1/Width so that a1 is automatically adjusted to match the requested a1/Width. The code will not change the Width so that code is not adapted to test blending between adjacent l modes, such as the l=0 and l=2 mode blending.
 *
 * The function calculates various constants and deploys the input parameters. It then generates a list of frequencies, heights, widths, splitting, centrifugal terms, latitudinal terms, and stellar inclination. Finally, it writes the mode parameters and noise parameters to the output files.
 */
void generate_cfg_from_synthese_file_Wscaled_Alm(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string extra);


/**
 * @brief Generate configuration file for an asymptotic star with active asymmetry and Gaussian height profile.
 *
 * @param input_params The input parameters for the star model.
 * @param file_out_modes The output file path for the mode parameters.
 * @param file_out_noise The output file path for the noise parameters.
 * @param extra The extra parameter for the star model.
 *
 * This function uses a reference star as a template to generate frequencies and width, height profiles. It can be rescaled so that you can modify the HNR but keep the same height profile. Note that the user here provides a target a1/Width so that a1 is automatically adjusted to match the requested a1/Width. The code will not change the Width so that code is not adapted to test blending between adjacent l modes, such as the l=0 and l=2 mode blending.
 *
 * The function calculates various constants and deploys the input parameters. It then generates a list of frequencies, heights, widths, splitting, centrifugal terms, latitudinal terms, and stellar inclination. Finally, it writes the mode parameters and noise parameters to the output files.
 */
void generate_cfg_from_synthese_file_Wscaled_aj(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string extra);


/**
 * @brief Generate configuration file for an asymptotic star with active asymmetry and Gaussian height profile.
 *
 * @param input_params The input parameters for the star model.
 * @param file_out_modes The output file path for the mode parameters.
 * @param file_out_noise The output file path for the noise parameters.
 * @param extra The extra parameter for the star model.
 *
 * This function uses a reference star as a template to generate frequencies and width, height profiles. It can be rescaled so that you can modify the HNR but keep the same height profile. Note that the user here provides a target a1/Width so that a1 is automatically adjusted to match the requested a1/Width. The code will not change the Width so that code is not adapted to test blending between adjacent l modes, such as the l=0 and l=2 mode blending.
 *
 * The function calculates various constants and deploys the input parameters. It then generates a list of frequencies, heights, widths, splitting, centrifugal terms, latitudinal terms, and stellar inclination. Finally, it writes the mode parameters and noise parameters to the output files.
 */
void generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string extra);

/**
 * @brief Generate configuration file for an asymptotic star with modes determined from a reference template and with noise as per defined in Kallinger+2014.
 *
 * @param input_params The input parameters for the star model.
 * @param file_out_modes The output file path for the mode parameters.
 * @param file_out_noise The output file path for the noise parameters.
 * @param extra The extra parameter for the star model.
 * @date 6 Dec 2023
 *
 * This function uses a reference star as a template to generate frequencies and width, height profiles. It can be rescaled so that you can modify the HNR but keep the same height profile. Note that the user here provides a target a1/Width so that a1 is automatically adjusted to match the requested a1/Width. The code will not change the Width so that code is not adapted to test blending between adjacent l modes, such as the l=0 and l=2 mode blending.
 *
 * The function calculates various constants and deploys the input parameters. It then generates a list of frequencies, heights, widths, splitting, centrifugal terms, latitudinal terms, and stellar inclination. Finally, it writes the mode parameters and noise parameters to the output files.
 */
void generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled_Kallinger2014(VectorXd input_params, std::string file_out_modes, std::string file_out_noise, std::string extra);

/**
 * @brief Compute the effect of centrifugal distortion on mode frequencies.
 *
 * @param fl0_all The input mode frequencies.
 *
 * @return The value of eta0.
 *
 * This function calculates the effect of centrifugal distortion on mode frequencies. It uses the input mode frequencies and various constants to calculate the value of eta0. The constants used are G (gravitational constant), Dnu_sun (solar large frequency separation), R_sun (solar radius), M_sun (solar mass), and rho_sun (solar density). It then returns the calculated value of eta0.
 */
double eta0_fct(const VectorXd& fl0_all); // To compute centrifugal term
