# Configuration in order to generate an ensemble of spectra automatically using [val_min] and [val_max] as the limits of the parameter space
# Seven inputs are required: 
# 		- model_name : Name of the model, as defined in models_database.pro
#		- The name of the parameters as defined in models_database.pro, for the subfunction [model_name]
#		- val_min : vector listing the initial parameters of the model [model_name]. See models_database.pro for more information about those parameters
# 		- val_max:  vector of same size as [val_min] with the final parameters of the model.
#		- Tobs and Cadence: Observation duration (in days) and the Cadence of observation (in seconds)
#		- forest_type: either grid or random. Currently, only random (uniform) is implemented.
#               - erase_old_file: If 1, then (1) the combination file is overwritten and (2) the model number (identifier) is reset to 0. 
#                                 If 0, then (1) append the combination file and (2) model number = last model number + 1
grid 10 # forest_type, followed by the forest_params. If forest_type=random, then forest_params is a single value that corresponds to the Number of samples
generate_cfg_from_synthese_file_Wscaled_aj   Configurations/infiles/8379927.in      
NONE		        # Used template(s) name(s). If several, randomly select one/iteration. If set to 'all', will use all *.template files in Configuration/templates 
Dnu  epsilon    delta0l_percent     HNR     a1ovGamma  Gamma_at_numax     a2      a3         a4      a5    a6      beta_asym     i       	H_spread   nu_spread   Gamma_spread    do_flat_noise     #Variable names 
70    0.5          1.               10        0.6         1.              0.1     -0.1       0.15    0.2   0.05      10          0.0    	5          0.2         5                    1            #val_min
70    0.5          1.               10        0.6         1.     	      0.1     -0.1       0.3     0.2   0.10     10           20         5          0.2         5                    1            #val_max
0      0           0                0.        0.0         0.     	      0       0.         0.15    0.    0.05     0            10         0          0           0                    0            #If forest_type="random" ==> 1=Variable OR 0=Constant. If forest_type="grid" then must be the stepsize of the grid
Tobs   Cadence  Naverage    Nrealisation
100    120       1            1
1     # It is erase_old_files. If set to 1, will remove old Combination.txt and restart counting from 1. Otherwise append Combination.txt
0     # Do you want plots ? 0 = No, 1 = Yes
1     # Do you list values of the input model in the output ascii file?
0     # Limit Data to mode range?
0     model_MS_Global_aj_HarveyLike  # Do .model files? If yes, use the name provided here to define the used model name within the .model file
