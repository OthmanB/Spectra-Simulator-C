# Version history #

### v1.55 ###
        * Change in the definition of d0l for aj MS models so that it is consistent with RGB models that are based on the ARMM code. thus the prescription nu = (n + l/2 + epsilon)*Dnu - d0l is used now
        * Cleanup configuration directory
        * Adding new templates for MS and RGB taken from a select number of good quality MS stars
        * removal of some obselete test files
        * Consistency improvment for numax_spread between aj model and RGB model configuration file

### v1.54 ###
	* Consistency improvment: The definition of numax_spread was different between model generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled_Kallinger2014 and
          model generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled, and model generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled_Kallinger2014 and other models.
	  numax_spread had to be given in fraction of numax in one case, while it had to be in percentage in others. Now all quantities are in percentage.

### v1.53 ###
	* Edge case handling in bump.cpp + iterative_artificial_spectrum.cpp: When attempting to generate a Subgiant with no mixed modes solutions, the code was crashing.
	  As a crash without explanation is not satisfactory, I have added message specifying that the generated star has no mixed modes. A hard-coded switch is currently 
          set to a mode that lead to skipping the "failed" computation of mixed modes. This switch, called "neverfail" inside bump.cpp can be also set to a safe exit (neverfail = 0)
          or to a process continuing after replacing the solution with pure p modes (neverfail =2)
 
### v1.52 ###
	* Adding the possibility to modify the output Data directory. After compilation, type ./specsim --help for further details.
	* Adding an option to create automatically the Data directory and its subdirectories if they don't exist. 

### v1.51 ###
	* New model "asymptotic_mm_freeDp_numaxspread_curvepmodes_v3_GRANscaled_Kallinger2014": This is the equivalent to 
	 "generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled_Kallinger2014" but for Red Giants. The noise follows Kallinger2014
	 which allows for a better training. Read carefully the notes on v1.50 to understand the new way of handling the noise.
	 Eventually, these models will replace the old non-GRANscaled. 
	* The new model requiring ARMM v1.11, this one is also updated in specsim v1.51
	
### v1.50 ###
	* New model: "generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled_Kallinger2014". This model allows to generate
	noise following closely the Kallinger+2014 prescription and ensemble values derived from a sample of RGB + MS stars. 
	Using MCMC, I have verified that this model is fairly accurate. Note that this model use the new structure of parameter
	that is decomposed into two files: the main.cfg (see the example configuration main.cfg.aj_GRANscaledKallinger) and the
	noise_Kallinger2014.cfg. The program accolate the configuration present in the two, to create a full model of noise + modes. This structure is prefered here due to (1) the high number of parameters in the Kallinger+2014 model (enhanced clarity with two shorter tables that fit on screen), (2) the necessity to define many Kallinger+2014 parameters as drawn from a Gaussian distribution instead of a Uniform distribution and, (3) the necessity to avoid confusion between the variables systematically drawn from a uniform distribution (unless fixed) defined the main.cfg files and the noise variables that very often drawn from the Gaussian distributions given in the Kallinger's paper. 
	Incidentally, this has imposed an important restructuration of the way the random generation works.
	* The adoption of v1.45_alt as the main branch, instead of v1.45.
	* Bug fix: The maximum frequency of the spectrum was wrong by a factor 2. 
	
### v1.45 and v1.45_alt ###
	* Bug fix in model asymptotic_mm_freeDp_numaxspread_curvepmodes_v3: indexes for the noise parameters were pointing to the wrong parameters + indexes for a2_l3, a3_l3, a4_l3 when checking <=-9999 were wrong.

### v1.44_alt ##
        * Take the version 1.44 but change the way the Gamma_spread parameter acts.  THIS VERSION IS EXPERIMENTAL: Gamma_spread IS APPLIED ON THE GAMMA_REF INSTEAD OF GAMMA_STAR
	
### v.1.44 ##
	* For model "asymptotic_mm_freeDp_numaxspread_curvepmodes_v3", Fix the logic that controls how we can switch between using predefined frequencies in the configuration file "Configurations/MixedModes_models/star_params.theoretical" and the global main.cfg configuration files (in "Configurations/main.cfg*"). This avoids configuration mistakes that may not be obvious to the final user. In addition, I remived the "Configurations/MixedModes_models/star_params.theoretical" and kept only the "Configurations/MixedModes_models/star_params.theoretical.EXAMPLE". This again to avoid misconfigurations.

### v1.43 ##
	* Adding do_flat_noise parameter in the aj model. If do_flat_noise<=0, the noise background will the one of the template (Harvey-like). Otherwise, the noise background will be flat and set at the value given by do_flat_noise (eg. if do_flat_noise = 0.7, then N0=0.7). do_flat_noise can be a variable and does not need to be a constant (ie, if one wants to train a ML with a variable flat noise background). But only >0 do_flat_noise will behave that way. <=0 will necessarily lead to fixed background: Be careful with the min and max value of the parameters and the fix/variable flag.
 
### v1.42 ##
	* Bug fixes in the noise background of aj model

### v1.41 ##
	* Adding Gamma_spread parameter in generate_cfg_from_synthese_file_Wscaled_aj

### v1.4 ##
        * Adding H_spread and nu_spread parameters in generate_cfg_from_synthese_file_Wscaled_aj and generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled models. The idea is to perturbate the heights (H_spread) and frequencies (nu_spread) with a random quantity to avoid deterministic learning of these quantities by a ML method
	* New parallelised Alm code (version 1.3)

### v1.31 ##
	* Adding a help for functions using Doxygen library
	* Improved the algorithm for creating the combinations in combi.cpp
	* Replacing the old "integrate" by the new "Alm"
	
### v1.3 ###
	* Models asymptotic_mm_freeDp_numaxspread_curvepmodes_v2 and asymptotic_mm_freeDp_numaxspread_curvepmodes_v3 updated to work with the new parallelized ARMM solver. 
	* Models asymptotic_mm_freeDp_numaxspread_curvepmodes_v2 and asymptotic_mm_freeDp_numaxspread_curvepmodes_v3 now handle ajl coefficients. Read instructions in the example main.cfg files associated to these models

### v1.2 ###
	* New model:
		- generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled allows you to generate a simulated spectra with the noise brackground of your choice, assuming that this noise scales with numax

### v1.1.2 ###
	* Improvments:
		- Adding the possibility to enter your own main file name (instead of the default main.cfg) and your own directory for the configuration file (instead of the default 'Configurations/')
		- Improved showversion() function
		- Improved handling of options by using the boost::program_options library

### v1.1.1 ###
	* Improvments:
		- In the aj model, adding rescaling capabilities so that Dnu, epsilon, numax and d0l can be modified
### v1.1.0 ###
	* Bug Fix:
		- Removed system command that were copying .rot and .range files from a temporary directory to the data directory. This was a relic of the time when python3 was called to solve the ARMM
	* Addition:
		- In asymptotic_mm_freeDp_numaxspread_curvepmodes_v3, added the possibility to provide a frequency input file. Can be used to link the simulator with frequencies generated by a Stellar evolution code + oscillation code. The file must be placed in  Configurations/MixedModes_models. Its name must be 'star_params.theoretical'. An example of such a file is provided in that directory. Note that the syntax is strict and no change in the file structure (order of its inputs, number of comments, etc...) should be performed.
		Note also that this file will be used ONLY if the second line contains a 1
		Like this, for the first two lines:
		# use_nu_nl : 0 = false / 1 = true
		1	 
### v1.0.4-dev ###
	* Bug Fix:
		- Corrected a bug happening when reading the cfg file due to imporper syntax
	* Improvements:
		- Minimising the changes required to bump_DP.cpp amd solver_mm.cpp when migrating them from the standalone version: Only need to remove the omp section
		- Slight restructuration of the ARMM calls
	* Addition:
		- Introducing the Hfactor and the Wfactor in asymptotic models for the mixed modes: It allows to 'break' the assumption of (1) equipartition of energy (Hfactor>0) and (2) inertia scaling with the Width, ie not considering that the damping is purely due to the evanescent zone between cavity (Wfactor<1)
### v1.0.3-dev ###
	* Bug Fix:
		- Correcting eta in the *Hgauss model
		- The rescaling of HNR in the models 'generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma' , 'generate_cfg_from_synthese_file_Wscaled_a1a2a3asymovGamma', 'generate_cfg_from_synthese_file_Wscaled_Alm', 'generate_cfg_from_synthese_file_Wscaled_aj' was incorrect due to a wrong generalisation from the old IDL code. 
		The old code was considering a fix noise background of N=N0=1. The new implementation was supposed to instead keep the noise background
		profile of the reference star in the new artificial star. This was not performed as it was still considering N0=1. 
		The modification consist in replacing N0 by N(v)=local_noise 
		- The rescaling was made considering the absolute highest HNR. This implies that it was considering the HNRmax for the l=1 (due to higher visibility of that mode). Instead, it is better to use the HNRmax(l=0) for reference when defining the HNR. A change was made in that way.
		- Fix an issue with some infiles used as a template: The Cadence and Tobs was wrongly commented by #, which was leading to ignore the first l=0
		- Fix an issue with the names of the variables in the Combinations.txt file: Some parameters names could be sticked together
		
### v1.0.2-dev ###
	* Improvment: 
		- Adding the possibility to use a file for defining the common section. See explanations in the Configuration/common_modelfile directory

### v1.0.1-dev ###
	* Bug fixes and improvments:
		- generate_cfg_from_synthese_file_Wscaled_a1Alma3asymovGamma():
			* Renaming to generate_cfg_from_synthese_file_Wscaled_Alm()
			* Fixing error in the Alm computation + optimisation following the latest release of acoef_checks
			* Possibility to use of all of aj coefficient up to a6 
	* New model: generate_cfg_from_synthese_file_Wscaled_aj()
	* Improvments:
		- Grid approach now handles *Alm() model type
		- Grid approach now handles *aj() model type
		- Update of build_lorentzian.* to the latest version from the cpp and adding dependences. Large code-quality improvments

### v1.0.0-dev ###
	* Added functionalities:
		- Implement the new model generate_cfg_from_synthese_file_Wscaled_a1Alma3asymovGamma: Handling simulations with Activity effect on splittings using Alm() function
                - Adding the grid approach within this code: NOTE THAT THERE MIGHT BE SOME BUGS HERE
                - Adding new functionalities to the configuration file: 
                	* do_modelfiles : If set to 1, it will make a .model file adequate for the TAMCMC code
                	* limit_data_range: If set to 1, it will limit the data file to the range where there is modes +/- 2*Dnu. Can hugely save space for large simulations
                	NOTE: This functionality is only valid for "artificial_spectrum_a1Alma3"
                	
### v0.9.1 [100%] ###
	* Added functionalities:
		- Implement on generate_cfg_from_synthese_file_Wscaled_a1a2a3asymovGamma: Handling simulation for a2 coeficient and a3 considering polynomials O2 for those
	* Bug Fix:
		- Fixing cfg template for main.cfg.MS.scaled
		
### v0.9.0 (*Released on 20 Nov 2020*) ###
	* Improvments:
		- Compilation is now made possible using cmake
	* Added functionalities:
		- Using bump_DP.cpp and solver_mm.cpp instead of their python counter parts
	* Bug Fix:
		- Bugs on the definition of np_min and np_max corrected
		- Bugs on the solver_mm that was not handling properly the second order polynomials for the p modes

### v0.8.97 (*Released on 21 April 2020*) ###
        * Improvments:
		- Adding nmax_star used value in .range files because the random control nmax_spread is not sufficient to determine the curvature of the p-modes for each models.
### v0.8.96 (*Released on 23 Jan 2020*) ###
	* Added functionalities:
		- Adding the possibility to generate multiple noise realisation for a single model (several .ascii files for a single .in file). The new control variable should appear in the main.cfg file and is called Nrealisation (see example configuration files)
	* Improvments:
		- Adding failsafes if negative frequencies are found for p modes

### v0.8.95 (*Released on 23 Jan 2020*) ###	
	* Added functionalities:
		- Implementation of the granulation noise using a single Harvey-profile, numax driven and following the Karoff+2011 prescription

### v0.8.93 (*Released on 20 Jan 2020*) ###		
	* Added functionalities:
		- Adding Templates database for Heights
                         WARNING: A new line is required in the main.cfg to specify the template file. This line must be set to NONE if the model does not require templates
		- Allow random jumps between heights a templates population 

### v0.8.92 (*Released on 16 Jan 2020*) ###		
	* Added functionalities:
		- Adding controls for Visibilities

### v0.8.91 ###
	* Bug Fix:
		- Crash of the code due to failure in the mixed mode solver when the second order curvature is too large

### v0.8.9 (*Released on 26 Dec 2019*) ###
	* Bug fix:
		[1] numax_spread improperly doubling the value of numax instead of just adding a %age of numax [SOLVED]
		[2] nmax_spread bug fix [SOLVED]
		[3] interpolation of height/width from the sun failing due to a forbiden extrapolation [SOLVED]
		[4] Incorrect numax_star inside bump_DP.width_height_MS_sun_rescaled() due to a typo between nmax_star and numax_star [SOLVED]
		[5] The model asymptotic_mm_freeDp_numaxspread_curvepmodes_v3 was not properly linked to the function of the same name into models_database.cpp [SOLVED]
		[6] Fix of Nmax issue: At high HNR, the fix number of generated modes can make the spectrum look 'cut' at high frequency. To solve this temporarily
		  I generate +2 radial order at high frequency (6 below numax, 8 above number) instead of having the same number of modes. The apparent issue is 
		  mostly due to the fact that high degree modes have extremely large FWHM, making them very visible.      [SOLVED]

	* Added functionalities:
		[1] on models with mixed modes, add second order effects on frequencies (curvature for p modes):
		 	D0 was a fix parameter while it should be a variable. Instead of sticking to D0, I implemented the second order effect on the
		  	asymptotic relation for p modes. Thus now, we have Eq. 22 from Mosser+2018 (https://www.aanda.org/articles/aa/pdf/2018/10/aa32777-18.pdf)
		  	encoded. This means that we have 2 new variable (+ 1 variable change for D0):
		  		- alpha_p : The curvature coefficient which should be around 0.076/nmax according to Mosser+2018, but it seems that 0.0076/nmax or 0.076/nmax^2 is more reasonable (typo in their paper?)
	  			- nmax_spread : The location of the inflexion point for the parabola is described by alpha_p * (np - nmax)^2 / 2. nmax is in principle given by 
	  			  				nmax_th= nu_max/Dnu - epsilon. But we don't necessarily want to strictly enforce that. Instead, nmax_spread allows you to depart from numax_th. The value is given in %
	  			- delta0l : The small separation. It relates to D0 and replaces it. The relationship is: delta0l=-l(l+1) D0 /Dnu
	  					    Note that the actual used parameter is delta0l_percent, defined in unit of l(l+1) because by definition:
	  					    delta0l = - l(l+1) delta0l_percent / 100.
	  		Modifications required for that implementation:
	  			[1] Changes in solver_mm.py: compatibility changes in the test functions, update of asympt_nu_p() and introduction of solve_mm_asymptotic_O2p() to 	handle second order effects. It introduced a new test function as well specifically for O2p tests test_asymptotic()  ===> [100%] [TESTED]
		  		[2] Changes in bump_DP.py: 
		  			[a] Write a new test function test_asymptotic_star_O2p() to evaluate the impacts for the changes in solver_mm.py  				 [100%]
					[b] Update make_synthetic_asymptotic_star to handle the new syntax in the file_cfg_mm created by the mixed modes models in C++   [100%]
					[c] Update main_star_generator																									 [100%]
	  			[3] In iterative_articial_spectrum.cpp: Add for all previous mixed modes models the new parameters    							 
	  					- asymptotic_mm_v1 																											[100%] [TESTED] 
	  					- asymptotic_mm_v2  																										[100%] [TESTED]
	  					- asymptotic_mm_v3  																										[100%] [TESTED]
	  					- asymptotic_mm_freeDp_numaxspread_curvepmodes_v1																			[100%] [TESTED] 
	  					- asymptotic_mm_freeDp_numaxspread_curvepmodes_v2																	  		[100%] [TESTED]
	  					- asymptotic_mm_freeDp_numaxspread_curvepmodes_v3																			[100%] 
	  			[4] In models_database.cpp: In all mixed modes models, replace D0 by the new parameters for the exchange file file_cfg_mm			[100%]
	  					- asymptotic_mm_v1 																											[100%] [TESTED]
	  					- asymptotic_mm_v2  																										[100%] [TESTED]
	  					- asymptotic_mm_v3  																										[100%] [TESTED]
	  					- asymptotic_mm_freeDp_numaxspread_curvepmodes_v1																			[100%] [TESTED] 
	  					- asymptotic_mm_freeDp_numaxspread_curvepmodes_v2																			[100%] [TESTED] 
	  					- asymptotic_mm_freeDp_numaxspread_curvepmodes_v3																			[100%] 
	  			[5] Change the main.cfg.XX for all templates according to the new format 															[100%]

		[2] At each iteration, remove the modes_tmp.cfg and noise_tmp.cfg in order to have a program termination if those files are not created properly [100%] [TESTED]



### v0.8.8 (*Released 19 Dec 2019*)###
	* Added functionalities:
		* Merging with several functions present in the IDL version 1.4. The model that is imported the followings:
			- generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma
		  This model has the particularity to allow you to use a reference star as a template. On IDL, the output was on a sav binary file.
		  Here the inputs must have the same format as the .in files
		Development status [FINISHED]: 
			[1] Changes in main.cfg  [100%]
			[2] Templates from main.cfg.synthese_file [0%]
			[3] Rename write_star_params.cpp and .h into io_star_params.cpp and .h 							 [100%]
			[4] Update write_star_params::read_main_cfg() to accept extra string argument (e.g. filenames)   [100%]
			[5] Changes in models_database.cpp: 
					- generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma							 [100%]
			[6] Prerequisite to [3]: 
				[a] Change to allow reading in file format:
					- Implement a reading function for in files: read_star_params() 						 [100%]
					- Implement a structure star_params in data.h   										[100%]
	                  The format of the outputs will be in a structure star_params with those parameters: 
			          VectorXd spec_params, MatrixXd mode_params, MatrixXd noise_params, std::string identifier
			    [b] Import noise_models.cpp and .h from TAMCMC-CPP 											[100%] 
			    [c] Update noise_modes.cpp and .h to handle MatrixXd shape inputs 							[100%] 
			    [d] Import string_handler.cpp and .h from TACMMC-CPP										[100%]
			    		- WARNING: During the merging process, I found duplicate but slightly different strsplit() function
			    		  Will need to investigate the effect of the changes on the code                         =======>   [TESTED]
			    		  At the moment, the original strtrim was put in strplit2()
			[7] Changes in iterative_artificial_spectrum 													[100%]
			[8] Write a function that convert sav files into .in format. 									[100%]
			[9] Update the README.md regarding the compilation line: 
				- noise_models.cpp must be added in the list of input files 								[100%]
				- string_handler.cpp and .h also 															[100%]
		Test status:
			- Compilation errors: No
			- Quick test running the program: Tested
			- Stress test using multiple configuration: Not performed

### v0.8.3 (*Released on 2 Dec 2019*)####
	* Added functionalities:
		* new model:asymptotic_mm_freeDp_numaxspread_curvepmodes_v1, v2 and v3. These are similar to asymptotic_mm_vx, 
		but I put the DP as a free parameter and the numax(Dnu) relation is not strictly enforced. The user has the possibility 
		To add a uniform spread around the numax(Dnu) relation. Recommended value is 0.2 (20%).
		Development status: 
			[1] Changes in main.cfg  [100%]
			[2] Templates from main.cfg.freeDP_curvepmodes.v1, v2 and v3 [100%]
			[3] Changes in models_database.cpp:
					- asymptotic_mm_freeDp_numaxspread_curvepmodes_v1 [100%]
					- asymptotic_mm_freeDp_numaxspread_curvepmodes_v2 [100%]
					- asymptotic_mm_freeDp_numaxspread_curvepmodes_v3 [100%]
			[4] Changes in iterative_artificial_spectrum [100%]
			[5] Changes in bump_DP [100%]
			
### v0.8.2 (*Released on 15 Nov 2019*)###
	* Added functionalities:
		* new model: asymptotic_mm_v2. This generates a uniform population of rotation rate in the envelope and 
		  compute the rotation rate in the core considering a uniform population of rotation in the core
		  The envelope rotation and the envelope-to-core ratio are therefore two new variables in the main.cfg, that replaces the variable Teff (effective temperature) of asymptotic_mm_v1
		* new model: asymptotic_mm_v3. This generates a uniform population of rotation rate in the envelope and in the core considering a uniform population in both case
		  The envelope rotation and the core rotation are therefore two new variables in the main.cfg, that replaces the variable Teff (effective temperature) of asymptotic_mm_v1

### v0.8.1: Starting public revision history (*Released on 21 Oct 2019*) ### 

	* Added functionalities:
        * Added this changelog.md file
		* Show the version and author information at start, or using the argument 'version' on the command line
		* adding a xxxx.python.rot file that gives the average core and average envelope rotation rate in microHz used for generating splittings

### v0.8.0: First release with rotation (*Released on 16 Oct 2019*) ### 

	* Radial differential rotation Rotation implemented for RGB and SG following Deheuvels+2015 and Mosser+2015
	
### v0.1.0: Basic initial release (* Released on August 2019*) ###
	
	* Added functionalities:
		* Use of asymptotic relations for p and mixed modes to generate spectra of Red Giants
	
	* Known issues: 
		* No support of SG or MS stars at the moment
		* Python routines are blended into C++... the python code is for test purpose and will have to be converted to C++ when reaching version v1.0.0
