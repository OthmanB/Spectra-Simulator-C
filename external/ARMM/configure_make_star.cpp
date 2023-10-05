#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "data_solver.h"
#include "configure_make_star.h"
#include "readparams_job.h"
#include "string_handler.h"
#include "solver_mm.h"
#include "bump_DP.h"


Cfg_synthetic_star configure_make_star(std::unordered_map<std::string, std::string> input_params){
    double xmin, xmax, nmax_spread, inc_star;

	std::random_device rd;
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<double> distrib(0 , 1);
    Cfg_synthetic_star cfg_star;
    // ---- Deploy parameters and check them ----
     //------------------    Un-used parameters  ------------------   
    cfg_star.Teff_star=-1;
    cfg_star.sigma_m=0;
	cfg_star.sigma_p=0;
    // ------------------   Used parameters  ------------------   
    // 	-----------------  p/g modes   ------------------
    cfg_star.Dnu_star=str_to_dbl(input_params["Dnu_star"]);
	cfg_star.DPl_star=str_to_dbl(input_params["DPl"]);                
	cfg_star.q_star=str_to_dbl(input_params["q_star"]);
	cfg_star.alpha_g_star=str_to_dbl(input_params["epsilon_g"]);
	cfg_star.epsilon_star=str_to_dbl(input_params["epsilon_p"]);    
	cfg_star.delta0l_percent_star=str_to_dbl(input_params["delta0l_percent"]);
    // ------------------   Rotation  ------------------   
	cfg_star.rot_env_input=str_to_dbl(input_params["rot_env"]);
    if (input_params["rot_core"] == "None" && input_params["rot_ratio"] == "None"){
        std::cerr << "Error: rot_core and rot_ratio cannot be both set to None" << std::endl;
        std::cerr << "Set one of them to >=0 and the other to None" << std::endl;
        exit(EXIT_FAILURE);
    }    
    if (input_params["rot_core"] != "None" && input_params["rot_ratio"] != "None"){
        std::cerr << "Error: rot_core and rot_ratio cannot be both have a value" << std::endl;
        std::cerr << "Set only one of them to None." << std::endl;
        exit(EXIT_FAILURE);
    }    
    if (input_params["rot_core"] != "None" && input_params["rot_ratio"] == "None"){
        cfg_star.rot_core_input=str_to_dbl(input_params["rot_core"]);
        cfg_star.rot_ratio_input=-1;
    }
    if (input_params["rot_core"] == "None" && input_params["rot_ratio"] != "None"){
        cfg_star.rot_ratio_input=str_to_dbl(input_params["rot_ratio"]);
        cfg_star.rot_core_input=-1;
    }
	cfg_star.rot_ratio_input=str_to_dbl(input_params["rot_ratio_input"]);
	cfg_star.rot_core_input=str_to_dbl(input_params["rot_core_input"]);
    // latitudinal effects
    cfg_star.env_lat_dif_rot.a3_l2=str_to_dbl(input_params["a3_l2"]);
    cfg_star.env_lat_dif_rot.a3_l3=str_to_dbl(input_params["a3_l3"]);
    cfg_star.env_lat_dif_rot.a5_l3=str_to_dbl(input_params["a5_l3"]);
    // apshericity effects
    cfg_star.env_aspher.a2_l1=str_to_dbl(input_params["a2_l1"]);
    cfg_star.env_aspher.a2_l2=str_to_dbl(input_params["a2_l2"]); 
    cfg_star.env_aspher.a2_l3=str_to_dbl(input_params["a2_l3"]);
    cfg_star.env_aspher.a4_l2=str_to_dbl(input_params["a4_l2"]);
    cfg_star.env_aspher.a4_l3=str_to_dbl(input_params["a4_l3"]);
    cfg_star.env_aspher.a6_l3=str_to_dbl(input_params["a6_l3"]);

    // ------------------   Global ------------------   
    cfg_star.maxHNR_l0=str_to_dbl(input_params["max_HNR"]);
	cfg_star.H0_spread=str_to_dbl(input_params["H0_spread"]);
	cfg_star.Gamma_max_l0=str_to_dbl(input_params["Gamma_max_l0"]);
	cfg_star.Hfactor=str_to_dbl(input_params["Hfactor"]);
	cfg_star.Wfactor=str_to_dbl(input_params["Wfactor"]);
    if (input_params["numax_star"] == "Auto"){
        if (input_params["numax_spread"] != "None"){ // numax_from_stello2009() will ingore cases with numax_spread <= 0
            cfg_star.numax_star=numax_from_stello2009(cfg_star.Dnu_star, str_to_dbl(input_params["numax_spread"])); // Second argument is the random spread on numax
        } 
    } else{
        if (str_to_dbl(input_params["numax_star"]) > 0){
            cfg_star.numax_star=str_to_dbl(input_params["numax_star"]);
        } else{
            std::cerr << "Error : numax_star must be > 0" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    if (str_to_dbl(input_params["fmin_in_Dnu"]) > 0 && str_to_dbl(input_params["fmax_in_Dnu"]) > 0){
	    cfg_star.fmin=cfg_star.numax_star - str_to_dbl(input_params["fmin_in_Dnu"])*cfg_star.Dnu_star;
	    cfg_star.fmax=cfg_star.numax_star + str_to_dbl(input_params["fmax_in_Dnu"])*cfg_star.Dnu_star;
    } else{
        std::cerr << "Error : fmin_in_Dnu and fmax_in_Dnu must be > 0" << std::endl;
        std::cerr << "        Recommended values are fmin_in_Dnu=6 and fmax_in_Dnu=8" << std::endl;
        exit(EXIT_FAILURE);
    }
	cfg_star.output_file_rot=input_params["output_file_rot"];
	cfg_star.filetemplate = input_params["file_template"];
	//cfg_star.Vl.resize(3);
    cfg_star.Vl = str_to_Xdarr(input_params["Vl"], " \t");
    if (cfg_star.Vl.size() !=3){
        std::cerr << "Error while reading the Visibility parameters" << std::endl;
        std::cerr << "You must provide 3 parameters" << std::endl;
        exit(EXIT_FAILURE);
    }
	if (input_params["nmax_star"] == "Auto"){
        cfg_star.nmax_star=cfg_star.numax_star/cfg_star.Dnu_star - cfg_star.epsilon_star;
    } else{
        if (str_to_dbl(input_params["nmax_star"]) < 0){
            std::cerr << "Error for nmax_star : Set it either to a positive value or to 'Auto'" << std::endl;
            exit(EXIT_FAILURE);
        } else{
            cfg_star.nmax_star=str_to_dbl(input_params["nmax_star"]);
        }
    }
    if (input_params["nmax_spread"] == "None"){ 
        nmax_spread=0;
    } else{ // If the parameter is set to "None", it will go to the catch part, which is no spread.        
        if (nmax_spread < 0){
            std::cerr << "Error: Please enter nmax_spread >=0" << std::endl;
            exit(EXIT_FAILURE);
        } else{
            str_to_dbl(input_params["nmax_spread"]);
        }
    }
	if (std::abs(nmax_spread) > 0)
	{
		try
		{
			xmin=cfg_star.nmax_star*(1. - std::abs(nmax_spread)/100.);
			xmax=cfg_star.nmax_star*(1. + std::abs(nmax_spread)/100.);
			cfg_star.nmax_star=xmin + (xmax-xmin)*distrib(gen);
			
		}
		catch(...)
		{
			std::cout << "Error debug info:" << std::endl;
			std::cout << "cfg_star.nmax: " << cfg_star.nmax_star << std::endl;
			std::cout << "nmax_spread: " << nmax_spread << std::endl;
			exit(EXIT_FAILURE);
		}
	}
    if (input_params["beta_p"] == "None" && input_params["alpha_p"] == "None"){
        std::cerr << "Error: alpha_p and beta_p cannot be both set to None" << std::endl;
        std::cerr << "Set one of them to 0 and the other to None if you do not want a p mode curvature" << std::endl;
        exit(EXIT_FAILURE);
    }    
    if (input_params["beta_p"] != "None" && input_params["alpha_p"] != "None"){
        std::cerr << "Error: alpha_p and beta_p cannot be both have a value" << std::endl;
        std::cerr << "Set one of them to None. Please read explanations on the parameters to choose which one." << std::endl;
        exit(EXIT_FAILURE);
    }    
    if (input_params["beta_p"] != "None" && input_params["alpha_p"] == "None"){
        cfg_star.beta_p_star=str_to_dbl(input_params["beta_p"]);
        cfg_star.alpha_p_star=cfg_star.beta_p_star/cfg_star.nmax_star;
    }
    if (input_params["beta_p"] == "None" && input_params["alpha_p"] != "None"){
        cfg_star.alpha_p_star=str_to_dbl(input_params["alpha_p"]);
        cfg_star.beta_p_star=cfg_star.alpha_p_star*cfg_star.nmax_star;
        std::cout << " INSIDE " << std::endl;
    }
    // ------------------  Noise ------------------   
	cfg_star.noise_params_harvey_like =  str_to_Xdarr(input_params["params_harvey_like"], " \t");    //[A_Pgran ,  B_Pgran , C_Pgran   ,  A_taugran ,  B_taugran  , C_taugran    , p      N0]
    if(cfg_star.noise_params_harvey_like.size() != 8) {
        std::cerr << "Error while reading the harvey like parameters" << std::endl;
        std::cerr << "You must provide 8 parameters" << std::endl;
        exit(EXIT_FAILURE);
    }    

    if (input_params["inclination"] == "Auto"){
        // Determination of an isotropic inclination
        double inc_rad=distrib(gen)*M_PI_2;
        double inc_y=distrib(gen);
        double inc_c=std::cos(inc_rad);
        while (inc_y <=inc_c){
            inc_rad=distrib(gen)*M_PI_2;
            inc_y=distrib(gen);
            inc_c=std::cos(inc_rad);
        }
        inc_star=inc_rad*180./M_PI;
    } else{
        if (str_to_dbl(input_params["inclination"]) > 0 && str_to_dbl(input_params["inclination"]) < 90){
            inc_star = str_to_dbl(input_params["inclination"]);
        } else{
            std::cerr << "Error : The stellar inclination must be either set to Auto or to a value between 0 and 90" << std::endl;
            std::cerr << "        The Auto mode will generate a random inclination distributed isotropically over a sphere" << std::endl;
        }
    }
    cfg_star.inclination=inc_star;
    return cfg_star;
 }