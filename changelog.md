# Version history #

### v0.9.0 (*Released on 2 Dec 2019*)####
	* Added functionalities:
		* Merging with several functions present in the IDL version 1.4. The models that have been imported are the followings:
			- generate_cfg_from_synthese_file_Wscaled_act_asym_cosi
			- generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma
		  These models have the particularity to allow you to use a reference star as a template. On IDL, the output was on a sav binary file.
		  Here the inputs must have the same format as the .in files
		Development status: 
			[1] Changes in main.cfg  [30%]
			[2] Templates from main.cfg.synthese_file [0%]
			[3] Rename write_star_params.cpp and .h into io_star_params.cpp and .h 							 [0%]
			[4] Changes in models_database.cpp: [START WITH THIS AND CREATE EMPTY TEMPLATES FOR PREREQUISITES]
					- generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma							 [70%]
			[5] Prerequisite to [3]: 
				[a] Change to allow reading in file format:
					- Implement a reading function for in files: read_star_params() 						 [0%]
					- Implement a structure star_params in data.h   										[100%] [Need testing]
	                  The format of the outputs will be in a structure star_params with those parameters: 
			          VectorXd spec_params, MatrixXd mode_params, MatrixXd noise_params, std::string identifier
			    [b] Import noise_models.cpp and .h from TAMCMC-CPP 											[100%]
			    [c] Update noise_modes.cpp and .h to handle MatrixXd shape inputs 							[100%] [Need Testing]
			    [d] Import string_handler.cpp and .h from TACMMC-CPP										[100%]
			[6] Changes in iterative_artificial_spectrum 													[0%]
			[7] Write a function that convert sav files into .in format. 									[0%]
					- Corrolary: Will need to update the IDLpostMCMC code to create a native .in output file [0%]
			[8] Code the grid capability 																	[0%]
			[9] Update the README.md regarding the compilation line: 
				- noise_models.cpp must be added in the list of input files 								[0%]
				- string_handler.cpp and .h also 															[0%]
			
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
