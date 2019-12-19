# Version history #

### ROADMAP v0.9.0 (*Expected Release: 30 Dec 2019*) ####
	* Added functionalities:
		- Code the grid capability
		- Implement numax variable on generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma

### ROADMAP v0.8.9 (*Expected Release: 23 Dec 2019*) ####
	* Added functionalities:
		- on models with mixed modes, add second order effects on frequencies: 
			[a] curvature for p modes
			[b] random error around that curvature

### IN DEV v0.8.8 (*Released 19 Dec 2019*)####
	* Added functionalities:
		* Merging with several functions present in the IDL version 1.4. The model that is imported the followings:
			- generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma
		  This model has the particularity to allow you to use a reference star as a template. On IDL, the output was on a sav binary file.
		  Here the inputs must have the same format as the .in files
		Development status: 
			[1] Changes in main.cfg  [100%]
			[2] Templates from main.cfg.synthese_file [0%]
			[3] Rename write_star_params.cpp and .h into io_star_params.cpp and .h 							 [0%]
			[4] update write_star_params::read_main_cfg() to accept extra string argument (e.g. filenames)   [100%]
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
					- Corrolary: Will need to update the IDLpostMCMC code to create a native .in output file [0%]
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
