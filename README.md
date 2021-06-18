# Spectra-Simulator-C
Spectrum Simulator in C++. It is intended to create synthetic spectra that follows the chi(2,2) noise statistics.

### Requirements ###
* gnuplot-iostream (provided)
* Eigen Library
* Boost Library
* openMP
* g++ compiler
* cmake (optional)

### How to compile? ###

The best way to compile is to use cmake as it will handle automatically user-specific configuration and platforms. To do so, you need to:

1. Make a **build** directory

2. Enter in this new directory and run **cmake ..** 

3. Transfer the created binary executable **specsim** file into the base directory of the program. 
 
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

  3. **models_database.cpp and model_database_grid.cpp**
	  * Different kind of models that can be used in order to generate the modes_tmp.cfg and noise_tmp.cfg files
    * The grid model database is exclusively for grid. The other one contains what is available for both grid and random case

### Main Dependences ###
	
  * **build_lorentzian.cpp**: is in charge of creating the lorentzian modes
  * **plot_diags.cpp**: is in charge of graphical plots. Possible only if gnuplot is installed (MacOS may have issues to generate plots)
  * **write_star_params.cpp** contains function intended to read/write 
  			- summary configuration files in the format [identifier].ASCII 
  			- configuration files for the modes (modes.cfg file)
  			- configuration files for the noise (noise.cfg file)
			- other formating subroutines for processing/converting strings and numbers (strsplit(), strtrim(),...)
   * **bump_DP.cpp** and **solver_mm.cpp** that handle the generation of the model for evolved stars
         

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


###  Assumptions for the model 'generate_cfg_from_refstar_HWscaled' ###

This is a generic model that use knowledge of a reference star and of stellar models in order to determine a relastic 
simulated power spectrum

**Noise:** 
	- 2 Harvey profiles + White noise

**Modes:**
	* Heights and Widths profiles: Are given by a reference star. All of the parameters of that reference star should be in the path: [program dir]/Configurations/ref_spectra.params
		The default reference star is the Sun. BEWARE THAT THE REFERENCE STAR CONTAINS ENOUGH MODES

	* Widths: Modes widths are rescaled using a reparametrisation of the frequency axis ```\nu``` into  ```(\nu - \nu_{max})/\Delta\nu```. Then we simply rescale the y-axis such as
		```W(n,l) = W_ref(n,l) * W_{maxHNR}/W_ref{maxHNR}```. This is similar to what was proposed in Figures of Appourchaux et al. 2014 (http://adsabs.harvard.edu/abs/2014A%26A...566A..20A)
        
	* Frequency: 
		Use a set of model parameters contained into ```[program dir]/external/MESA_grid/models.params``` For test-purpose, you can use instead models_samples.params (save the original and rename 'models_samples.params' in 'models.params')
	
	* Rotation: The first order effect of rotation is incoroporated,
	
	* a1: rotational splitting as nu(n,l, m) = nu(n,l) + m a1
  
  * Lorentzian asymmetry: No asymmetry

	* i : Stellar inclination, assumed constant for all modes (this is a very weak assumption and thus, very reasonable)
		  See Gizon & Solanki 2003 (http://adsabs.harvard.edu/abs/2003ApJ...589.1009G) for further details



###  Assumptions for the model 'generate_cfg_from_refstar_HWscaled_GRANscaled' ### 

This is a generic model that use knowledge of a reference star and of stellar models in order to determine a relastic 
simulated power spectrum.
By allowing to use a granulation noise that scales with nu_{max}, this model is drastically reducing the parameter space.
Typically, 1 to 2 free parameters for the noise are required. Instead of 6 (assuming the N0 is normalized to 1 in both cases).

**Noise:**
	- 1 Harvey profile that has its parameters scaled with nu_{max}. Indeed, several studies, both theoretical and observational (e.g. Mathur et al. 2011, 741:119) suggest a strong correlation of the granulation noise properties with nu_{max}. P, the maximum power of the Harvey profile scales approximately with (nu_{max})^(-2). tau, the timescale of convection scales with (nu_{max})^(-1). From this, the user of this code can set the functionnal form for P and tau as follow,
		```P = Ap * (nu_{max})^(Bp) + Cp``` and ```tau=At * (nu_{max})^(Bt) + Ct```
	  With the following recommendation: 
		```Ap=At=1```
		```Bp between 1.8 and 2.2```
		```Bt between 0.85 and 1.15```
		```Cp=Ct=1```
              
	- White noise

Modes (same as for 'generate_cfg_from_refstar_HWscaled') :
	- Heights and Widths profiles: Are given by a reference star. All of the parameters of that reference star should be in the path: [program dir]/Configurations/ref_spectra.params
		The default reference star is the Sun. BEWARE THAT THE REFERENCE STAR CONTAINS ENOUGH MODES

	- Widths: Modes widths are rescaled using a reparametrisation of the frequency axis \nu into  (\nu - \nu_{max})/\Delta\nu. Then we simply rescale the y-axis such as
		W(n,l) = W_ref(n,l) * W_{maxHNR}/W_ref{maxHNR}. This is similar to what was proposed in Figures of Appourchaux et al. 2014 (http://adsabs.harvard.edu/abs/2014A%26A...566A..20A)
        
	- Frequency: 
		Use a set of model parameters contained into [program dir]/external/MESA_grid/models.params
		For test-purpose, you can use instead models_samples.params (save the original and rename 'models_samples.params' in 'models.params')
	
	- Rotation: The first order effect of rotation is incoroporated,
		- a1: rotational splitting as nu(n,l, m) = nu(n,l) + m a1
        - Lorentzian asymmetry: No asymmetry
	- i : Stellar inclination, assumed constant for all modes (this is a very weak assumption and thus, very reasonable)
		  See Gizon & Solanki 2003 (http://adsabs.harvard.edu/abs/2003ApJ...589.1009G) for further details



### Assumptions for the model 'asymptotic_mm_v1' ###

**Noise:** 

Follows the general description given above

**Modes:**

* Heights: 
	- Variations of heights are rescaled using Solar Heights for l=0 modes
	
	- Due to the fact that we impose Inertia(l=1).Width(l=1) = Inertia(l=0).Width(l=0), the Heights of l=1 mixed modes are function of sqrt(1-ksi) and scale with the height of l=0 modes, modulo the bolometric visbility. 
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
  -  The splitting assumptions:
	* We generate a population of stars with a surface rotation randomly distributed (truncated Gaussian) between 30 and 90 days. 

	* For non-mixed modes, the splitting is assumed to be equal to the surface rotation (solid body rotation)

	* For l=1 mixed modes, the splitting account for the degree of mixture of the modes using the ksi function and assuming a core-to-envelope rotation ratio function of the log(g) dependence that follows the Fig 13 of Deheuvels et al. 2014 for 3.55 <=log(g)< 3.8825 (early RGB and SG regime).
        
        * It is unity (solid body) if log(g)>=3.8825 (MS star).

        * It is constant for log(g) < 3.55 (evolved RGB), at a value witch is in continuity with Fig 13 of Deheuvels et al. 2014. The flat ratio in agreement with Charlote Gehand's results. BUT THIS IS OBVIOUSLY NOT VALID FOR CLUMP STARS. 
	
        * The core-to-envelope rotation contrasts from Deheuvels et al. 2014 are parametrised with a 2nd order polynomial fit, function of log(g). log(g) is determined are a pure function of the temperature: ```g = (numax/numax_sun) ( Teff / Teff_sun)^0.5```. Uncertainties on the log(g) are fixed to 0.1. Randomisation of the results are also performed using the uncertainties of the polynomial fit, but in a simplified manner. As shown in the code (bump_DP.py), weuse the quadratic mean of the relative error for stars B-E, not accounting for the uncertainties given by stars A and F. The relative uncertainty on the ratio is then fixed to 44.8%.
        
	*  The stellar inclination follows an isotropic distribution (uniform in cos(i)).
    
**Noise:***
	White noise only. No frequency-dependent noise so far implemented.
	

### Assumptions for the model 'asymptotic_mm_v2' ###

Same as for asymptotic_mm_v1 but with the following changes:

**Modes:**

* Rotation and inclination: *Splitting and inclination implementation* 
  -  The splitting assumptions:
	* The user can generate a population of stars with a surface rotation uniformly and randomly distributed two values (min and max). 
	* For l=1 mixed modes, instead of imposing some relationship with the evolution of the star regarding the core-to-envelope ratio, we let the user define two bundaries for an uniform sampling of the core-to-envelope ratio (see main.cfg.v2 for a configuration example)
    
### Assumptions for the model 'asymptotic_mm_v3' ###

Same as for asymptotic_mm_v1 but with the following changes:

**Modes:**

* Rotation and inclination: *Splitting and inclination implementation* 
  -  The splitting assumptions:
	* The user can generate a population of stars with a surface rotation uniformly and randomly distributed two values (min and max). 
	* For l=1 mixed modes, instead of imposing some relationship with the evolution of the star regarding the core rotation rate, we let the user define two bundaries for an uniform sampling of the core rotation rate (see main.cfg.v3 for a configuration example)
	
### Assumptions for the models 'asymptotic_mm_freeDp_numaxspread_curvepmodes_v1' / 'asymptotic_mm_freeDp_numaxspread_curvepmodes_v2' / 'asymptotic_mm_freeDp_numaxspread_curvepmodes_v3' ###

Same as for asymptotic_mm_v1, asymptotic_mm_v2 and asymptotic_mm_v3 but with the following changes:

**Modes:**

* Period Spacing: Uniformly generated over a range defined by the user into the cfg file
* Possibility to introduce some extra random spreads on numax

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

