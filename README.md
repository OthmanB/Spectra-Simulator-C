# Spectra-Simulator-C

> New here? Start with **`docs/User_Guide.md`**, then browse the model manuals in **`docs/models/README.md`**. Doxygen remains the API reference (`docs/Readme.md`).

Spectrum Simulator in C++. It is intended to create synthetic spectra that follows the chi(2,2) noise statistics.

Tests for MacOS and Linux:
 [![CI](https://github.com/OthmanB/Spectra-Simulator-C/actions/workflows/cmake_UbuntuMac.yml/badge.svg?branch=master)](https://github.com/OthmanB/Spectra-Simulator-C/actions/workflows/cmake_UbuntuMac.yml)

### Documentation ###

* Online docs (MkDocs): https://othmanb.github.io/Spectra-Simulator-C/
* User Guide: `docs/User_Guide.md`
* Model manuals (non-obsolete): `docs/models/README.md`
* Legacy model notes (obsolete): `docs/legacy_models.md`
* API reference (Doxygen): `docs/Readme.md`
  

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



### 1. iterative_artificial_spectrum.cpp

  

- Create an ensemble of artificial spectra by reading a configuration file (main.cfg). These are currently created using a random generator (uniform).

- You can create up to 9,999,999 spectra with a single execution and format them properly (with leading zeros in the name).

- There are two possible ways of using the code:

	- If `erase_old_files=1` parameter is used: Any older summary file (Combinations.txt) will be overwritten (data loss). If the Data subdirectories contain files with the same names as those processed (e.g., 0000010.ascii), they will also be overwritten. If this mode is used, it is strongly recommended to empty the Data directory and its subdirectories before running the program.

  

	- If `erase_old_files=0` parameter is used: Any older summary file (Combinations.txt) will be appended (data kept). Any data within the Data subdirectories will not be erased, assuming that their ID number matches the ID number within the Combinations.txt file.

  
  
2.  **artificial_spectrum.cpp**

* Create a single artificial spectra by reading a configuration files (modes_tmp.cfg and noise_tmp.cfg)

provided into the Configuration directory. See these files to learn about their syntax

* Save the final spectrum parameters (noise parameters and mode parameters) into files within the Spectra_info directory

* Save the model parameters into a [identifier].ASCII file. Allows a easy access to the configuration

* Data are saved into the Data directory

  

3.  **models_database.cpp and model_database_grid.cpp**

* Different kind of models that can be used in order to generate the modes_tmp.cfg and noise_tmp.cfg files

* The grid model database is exclusively for grid. The other one contains what is available for both grid and random case

  

### Main Dependences ###

*  **build_lorentzian.cpp**: is in charge of creating the lorentzian modes

*  **plot_diags.cpp**: is in charge of graphical plots. Possible only if gnuplot is installed (MacOS may have issues to generate plots)

*  **write_star_params.cpp** contains function intended to read/write

- summary configuration files in the format [identifier].ASCII

- configuration files for the modes (modes.cfg file)

- configuration files for the noise (noise.cfg file)

- other formating subroutines for processing/converting strings and numbers (strsplit(), strtrim(),...)

*  **bump_DP.cpp** and **solver_mm.cpp** that handle the generation of the model for evolved stars

  

### The model ###

  

As a general description,

* Each mode is described by a lorentzian function, with parameters:

  

- H(n,l,m) : The mode maximum Height (ppm^2/uHz)

- nu(n,l,m): The mode central frequency (uHz)

- W(n,l,m) : The mode width (uHz)

  

* The noise is described using two Harvey-like profile + White noise, with parameters:

The Harvey-like profiles have for parameters:

- H : Maximum noise level (ppm^2/uHz)

- tau: Characteristic timescale (of convection) (kilo-sec)

- p : Characteristic power law (of convection) (no unit)

- N0 : White noise level (ppm^2/uHz).

Model-specific manuals are in `docs/models/README.md`. Legacy model notes were moved to `docs/legacy_models.md`.

  

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

  

* Owner: Othman Benomar (Visiting researcher at NAOJ; Craftsman Software)

* Contact: othman.benomar.pro@gmail.com othman.benomar@craftsman-software.com
