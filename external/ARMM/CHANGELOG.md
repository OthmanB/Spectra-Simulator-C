### 0.75 ##
	- Adding new elements in cfg_star in order to allow the user to provide Stellar model-based input frequencies instead of asymptotic. Note that the asymptotic is still used to compute Widths, Heights and rotation for the two zone case. The new parameters are:
			- cfg_star.nu_nl : Input frequencies to be used for l=0, l=1, l=2 and l=3. Must be a matrix type
			- cfg_star.use_nu_nl : If true, uses the frequencies provided in cfg_star.nu_nl
	Note that the following (old) parameters MUST also be provided in order for the function work:
			- cfg_star.Dnu_star : Can be calculated eg from a stellar model sound speed
			- cfg_star.DPl_star : Can also be calculated theoretically
			- cfg_star.q_star
			- cfg_star.alpha_g_star
	epsilon_p_star delta_0l will be calculated within make_synthetic_asymptotic_star, if cfg_star.use_nu_nl = true

### 0.7 ###
        - Adding parallelisation in solve_mm_asymptotic_O2from_l0() and solve_mm_asymptotic_O2from_nupl() with fixed 4 threads. Not added to solve_mm_asymptotic_O2p()
        - Adding a filter to avoid any computation when nu_g is beyond the numin / numax range provided in solver_mm(). Might be faster for large number of nu_g
	- Added a solve_mm_asymptotic_O2from_l0_DEV() for tests... exclusively linked to the test_asymptotic_sg_O2from_l0_DEV() function which is used when generating the test.prg through the test.sh script

### 0.6  ###
	- Code rewriten using proper headers that include only constants so that there is no clashes. Vectors from Eigen are as well passed now by reference
	- Added a main.cpp for tests
	Bug Fix: 
		- Correction of ksi_fct1() to be sure that in precise mode, you get something smaller than 1
		- Correction of an error that was imposing too often a constant derivative of nu_l (hence constant Dnu_p). 

### 0.55 [100%]
	Including a new function 'solve_mm_asymptotic_O2from_l0()' inside solver_mm.cpp that can be used to generate l>0 p modes
	using a series of provided l=0 modes. This differs from solve_mm_asymptotic_O2p() from the fact that Dnu_p, epsilon, fmin, fmax
	alpha and nmax are not required because the new introduced approach treat the Oj(n), with j>1 terms as descrived by the l=0 p modes. 
	This new function is mostly useful for fitting the asymptotic relation of the mixed modes over a spectra.

### 0.5 ###
	Fully operational version with exhaustive testing.	
### 0.35alpha [100%]
	Take back 0.3alpha and test the new implementations regarding random number generators for randomness of nu_p and nu_m 

### 0.3alpha [100%]
	Implementation of all relevant function from bump_DP.py into bump_DP.cpp [100%]
	Testing all the functions internally in the c++ and check the behavior of of the master code creating values [100%]
	Testing comparatively with python code the master functions... BEWARE THAT SMALL BUGS WERE FOUND THAT MAY CHANGE THE RESULTS [100%]
	Note:  I found out that the 'fast' method used to compute ksi_fct2() is not accurate enough. While this inaccuracy could be acceptable
		   for RGB (as the spectrum resolution does not allow us to resolve the modes that are mostly affected by the approximation), it is not
		   the case for SG. After some thought, I decided to just use the 'precise' method, assuming a resolutino of 4 years in the data for defining
		   the normalisation constant. It slower than the fast method, but remains managable in the C++ implementation.
	
### 0.2alpha [100%]
	Implementation of solve_mm_asymptotic_O2p() [100%]
	Implementation of a testing function for solve_mm_asymptotic_O2p() [100%]
	Testing comparatively with the python code: [100%] 

	Code Improvements: The python version of solver_mm had an error: The curvature of the p modes is not properly handled in the python code. The global large separation is used to determine the solutions of tan(\thetap) = tan(\thetag) while it should be the local one, see Eq. 3.31 in Charlotte Gehan thesis (https://tel.archives-ouvertes.fr/tel-02128409/document). Frequencies are barely unchanged during my test but this might be important if large curvatures exists in p modes. 

	Performance Improvements: Limiting the use of conservativeResize by using static arrays.

### 0.1alpha [DONE]
	Implementation of solver_mm()
	Implementation of a testing function for sg
	Implementation of a testing function for RGB
	Testing comparatively with the python code: Works well (actually may be even better in terms of precision)
