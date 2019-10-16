# Spectra-Simulator-C
Spectrum Simulator in C++. It is intended to create synthetic spectra that follows the chi(2,2) noise statistics.

### Requirements ###
* python3 with numpy and scipy
* gnuplot-iostream (provided)
* Eigen Library (provided)
* Boost Library
* openMP
* g++ compiler

### How to compile? ###

```
g++ -O3 -I eigen -fopenmp -lutil -lboost_iostreams -lboost_system -lboost_filesystem -lgsl -lgslcblas artificial_spectrum.cpp write_star_params.cpp build_lorentzian.cpp function_rot.cpp plots_diags.cpp iterative_artificial_spectrum.cpp models_database.cpp random_JB.cpp -o ./sim.out
```

or without openmp (MacOS compiler does not support openmp):
```
g++ -O3 -I eigen -lutil -lboost_iostreams -lboost_system -lboost_filesystem -lgsl -lgslcblas artificial_spectrum.cpp write_star_params.cpp build_lorentzian.cpp function_rot.cpp plots_diags.cpp iterative_artificial_spectrum.cpp models_database.cpp  random_JB.cpp -o ./sim.out
```

### The program ###

The Program 'Spectrum Simulator' creates synthetic spectra for asteroseismology.
It is composed of different functions and procedures, enumerated and explained here.

* Main procedures 
   
  1. **iterative_artificial_spectrum.cpp**
	    * Create an ensemble of artificial spectra by reading a configuration file (main.cfg). These are for the moment created using a random generator (uniform).
	  You can create up to 9 999 999 spectra with a single execution and formated properly (with 0s in front of the name). 
	    * Two possible way of using the code:
		     - parameter erase_old_files=1: Any older summary file (Combinations.txt) will be OVERWRITTEN (data loss then). If the Data subdirectories are not
		  files of same name as those processed (e.g. 0000010.ascii) will also be overwritten. If this mode is used, it is STRONGLY recommended to empty the 
		  Data directory and its subdirectories before running the program.
          - parameter erase_old_files=0: Any older summary file (Combinations.txt) will be APPENDED (data kept). Any data within the Data subdirectories will
 		    not be erased, ASSUMING that their ID number is matching the ID number within the Combinations.txt file. 

  2. **artificial_spectrum.cpp**
	  * Create a single artificial spectra by reading a configuration files (modes_tmp.cfg and noise_tmp.cfg)
	  provided into the Configuration directory. See these files to learn about their syntax
	  * Save the final spectrum parameters (noise parameters and mode parameters) into files within the Spectra_info directory
	  * Save the model parameters into a [identifier].ASCII file. Allows a easy access to the configuration
	  * Data are saved into the Data directory

  3. **models_database.cpp**
	  * Different kind of models that can be used in order to generate the modes_tmp.cfg and noise_tmp.cfg files

### Main Dependences ###
	
  * **build_lorentzian.cpp**: is in charge of creating the lorentzian modes
  * **plot_diags.cpp**: is in charge of graphical plots. Possible only if gnuplot is installed (MacOS may have issues to generate plots)
  * **write_star_params.cpp** contains function intended to read/write 
  			- summary configuration files in the format [identifier].ASCII 
  			- configuration files for the modes (modes.cfg file)
  			- configuration files for the noise (noise.cfg file)
			- other formating subroutines for processing/converting strings and numbers (strsplit(), strtrim(),...)
   * **BETA** python external functions **bump_DP.py** and **solver_mm.py** that handle the generation of the model for evolved stars
             THIS WILL NEED TO BE IMPLEMENTED IN PURE C++ AT LATER STAGE FOR EFFICIENCY REASONS

### The model ###

As a general description,
* Each mode is described by a lorentzian function, with parameters:

	- H(n,l,m) : The mode maximum Height (ppm^2/uHz)    
	- nu(n,l,m): The mode central frequency  (uHz)
	- W(n,l,m) : The mode width (uHz)

* The noise is described using two Harvey-like profile + White noise, with parameters:
The Harvey-like profiles have for parameters:
	 - H : Maximum noise level (ppm^2/uHz)
	 - tau: Characteristic timescale (of convection) (kilo-sec)
	 - p  : Characteristic power law (of convection) (no unit)
	 - N0 : White noise level (ppm^2/uHz).

### Assumptions for the model 'generate_cfg_asymptotic_act_asym_Hgauss' ###

**Noise:** 
  Follows the general description given above.

**Modes:**

* Heights: Variations of heights are modeled by gaussian (a more accurate description would be a Voigt Profile)

    - ```H(n,l) = V(l).maxH exp( -0.5 * (nu - numax)^2/sigma^2), with sigma= 2.Dnu```
		- ```H(m | n,l) = Constant```
		- ```V(0) = 1, V(1) = 1.5, V(2) =0.5 , V(3) =0.07.``` This is the solar values (and verified empirically for Main sequence stars)

* Widths: Modes widths are assumed to be constant (a more accurate description would use the Eq.1 of Appourchaux et al. 2014 (http://adsabs.harvard.edu/abs/2014A%26A...566A..20A)
		- ```W(n,l,m) = Constant = G```
        
* Frequency: Strictly follow the asymptotic relation for the p modes:
	- ```nu(n,l) = ( n + epsilon + l/2) Dnu - l(l+1)D0```
	
* Rotation and inclination: Different effect on rotation are incoroporated,
  - a1: rotational splitting as ```nu(n,l, m) = nu(n,l) + m a1```
  - eta: Centrifugal force (fixed parameter) as,
		    '''
        eta=(4./3)*PI * a1^2 / (G * rho_sun) * (Dnu_sun/Dnu)^2
        '''
			This because the star density if determined at <1% (from theory) by rho=rho_sun * (Dnu_sun/Dnu)^2

  - a3: Effect of latitudinal rotation on the frequencies. See Gizon & Solanki 2004 
		      'Measuring Stellar Differential rotation with asteroseismology' for further details
  - alfa and b: Describes (empirically) the frequency shift due to stellar activity
				 The functional form is: ```b.nu^alfa```
				 alfa~3
				 b<<1
  - beta: Lorentzian asymmetry expressed as defined by Eq.8 of Gizon.L, CEAB 2006 (http://adsabs.harvard.edu/abs/2006CEAB...30....1G)
 
  - i : Stellar inclination, assumed constant for all modes (this is a very weak assumption and thus, very reasonable)
		  See Gizon & Solanki 2003 (http://adsabs.harvard.edu/abs/2003ApJ...589.1009G) for further details


### Assumptions for the model 'asymptotic_mm_v1' ###

**Noise:** 

Follows the general description given above

**Modes:**

* Heights: 
	- Variations of heights are rescaled using Solar Heights for l=0 modes
	
	- *[Code update on 27 August 2019]* Due to the fact that we impose Inertia(l=1).Width(l=1) = Inertia(l=0).Width(l=0), the Heights of l=1 mixed modes are function of sqrt(1-ksi) and scale with the height of l=0 modes, modulo the bolometric visbility. 
	This is the result from discussion with Kevin Belkacem. Note however that this assumption might only be valid for not-too-evolved RGB stars and for subgiants (see Belkacem+2018, 'Angular momentum redistribution by mixed modes in evolved low-mass stars')
	
	- Bolometric visibilities in height are assumed to be V(l=0) = 1, V(l=1) = 1.5, V(l=2) =0.5 , V(l=3) =0.07

* Widths:
	- l=0, 2, 3 modes are rescaled using the synthetic relation from Appourchaux+2014 applied to the solar profile
	- l=1 mixed modes are defined using the ksi function, scaled using l=0 modes. Thus l=0 modes fixes the upper limit of the mode width
        
* Frequency: 
  - l=0,2,3 Strictly follow the asymptotic relation for the p modes:
		```nu(n,l) = ( n + epsilon + l/2) Dnu - l(l+1)D0```
  - The frequencies of the l=1 modes follow exactly the asymptotitc relation for the mixed modes	
	
* Rotation and inclination: *Splitting and inclination implementation* 
  - *[Code update on 16 Oct 2019]* The splitting assumptions:
	* We generate a population of stars with a surface rotation randomly distributed (truncated Gaussian) between 30 and 90 days. There is the possibility to enter your own rotation at the surface (rot_env_input variable in bump_DP.make_synthetic_asymptotic_star()), but this requires manual intervention or extra coding.

	* For non-mixed modes, the splitting is assumed to be equal to the surface rotation (solid body rotation)

	* For l=1 mixed modes, the splitting account for the degree of mixture of the modes using the ksi function and assuming a core-to-envelope rotation ratio function of the log(g) dependence that follows the Fig 13 of Deheuvels et al. 2014 for 3.55 <=log(g)< 3.8825 (early RGB and SG regime).
        
        * It is unity (solid body) if log(g)>=3.8825 (MS star).

        * It is constant for log(g) < 3.55 (evolved RGB), at a value witch is in continuity with Fig 13 of Deheuvels et al. 2014. The flat ratio in agreement with Charlote Gehand's results. BUT THIS IS OBVIOUSLY NOT VALID FOR CLUMP STARS. 
	
        * The core-to-envelope rotation contrasts from Deheuvels et al. 2014 are parametrised with a 2nd order polynomial fit, function of log(g). log(g) is determined are a pure function of the temperature: ```g = (numax/numax_sun) ( Teff / Teff_sun)^0.5```. Uncertainties on the log(g) are fixed to 0.1. Randomisation of the results are also performed using the uncertainties of the polynomial fit, but in a simplified manner. As shown in the code (bump_DP.py), weuse the quadratic mean of the relative error for stars B-E, not accounting for the uncertainties given by stars A and F. The relative uncertainty on the ratio is then fixed to 44.8%.
        
	*  The stellar inclination follows an isotropic distribution (uniform in cos(i)).
    
**Noise:***
	White noise only. No frequency-dependent noise so far implemented.
	

### QUICK START ###

  1. Compile and verfity that the python3 program is properly installed with its dependencies
  2. check the configuration of the main.cfg file. The default setup should be good to go, but check:
	    * The binary variable handling the outputs at the end of the main.cfg file
	    * forest type parameters. The default is ```random 4``` (for test purpose mostly). If you need ```100``` models, created in a row, you need to set ```random 100```.
  3. run: ```./sim.out```
  4. Check the outputs in the Data directory. Check also the Combinations.txt file generated in the Configuration directory.
 
 
### Contribution guidelines ###

No external contribution is expected. This project is constantly improved, so please contact me if you need to see some worthy functionnality implemented. 

### Who do I talk to? ###

* Owner: Othman Benomar (NAOJ Research Fellow, Visiting Scientist at NYU Abu Dhabi)
* Contact: othman.benomar@nao.ac.jp   ob19@nyu.edu

