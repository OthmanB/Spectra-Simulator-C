# Version history #

### ROADMAP v0.9.0 (*Expected Release: Jan 2020*) ####
	* Added functionalities:
		- Code the grid capability
		- Implement numax variable on generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma

### ROADMAP v0.8.95 (*Expected Release: Jan 2020*) ####
	* Added functionalities:
		- Adding controls for Visibilities
		- Addning controls for Nmax (?): Probably no need since the fix [6] in v0.8.9
		- Adding random errors around the second order asymptotic relations. Try to make them spline-consistent

### v0.8.9 (*Release on 26 Dec 2019*) ####
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



### v0.8.8 (*Released 19 Dec 2019*)####
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
