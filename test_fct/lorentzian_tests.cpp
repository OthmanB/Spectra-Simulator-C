//#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include "../build_lorentzian.h"
#include "../function_rot.h"
//#include "../string_handler.h" // REPLACED BY ioproc.h on 17/06/2021
#include "../ioproc.h"
#include "../version.h"
#include <unistd.h>

using Eigen::VectorXd;

/* To generate tests

*/

void showversion();
int options(int argc, char* argv[]);
void usage(int argc, char* argv[]);
void write_model(const VectorXd& x, const VectorXd& y, std::string file_out);

int main(int argc, char* argv[]){

		bool verbose=0, pass=false;
		int i, testval;
		std::string filename_params;
		std::string char0, line0;
		std::ifstream cfg_session;

 		int el, Ndata;
		double fc_l, H_l, gamma_l, a1, a2, a3, a4, a5, a6, eta0, epsilon, asym, inc;
		VectorXd range, V, x_l, model, thetas(2);
	
		//Model_def model_list;
		// Set the current path variables
		std::string cpath=getcwd(NULL, 0);
		std::string file_out = cpath + "/model.out";
		std::string model_name;

		testval=options(argc, argv);
	
		model_name=argv[1];
		filename_params=argv[2];
		
		
		std::cout << "  0. Configuration: " << std::endl;
		std::cout << "      - Parameters file: " << filename_params << std::endl;
		std::cout << "      - Output file: " << file_out << std::endl;
		std::cout << "      - Model Name: " << model_name << std::endl;
		std::cout << "  2. Reading the file with the parameters for the build_lorentzian::" << model_name << " function..." << std::endl;
		if(model_name == "build_l_mode_a1a2a3") {
			std::cout << "     expecting the following paramters: [fmin, fmax, resol], l, fc_l, H_l, gamma_l, a1, a2, a3, asym, inclination" << std::endl;
			pass=true;
		}
		if(model_name == "build_l_mode_aj") {
			std::cout << "     expecting the following paramters: [fmin, fmax, resol], l, fc_l, H_l, gamma_l, a1, a2, a3, a4, a5, a6, eta0, asym, inclination" << std::endl;
			pass=true;
		}
		if(model_name == "build_l_mode_a1etaAlma3") {
			std::cout << "     expecting the following paramters: [fmin, fmax, resol], l, fc_l, H_l, gamma_l, a1, eta0, epsilon_nl, theta, delta, a3, asym, inclination" << std::endl;
			pass=true;
		}
		if (pass == false){
			std::cout << " Error: Provided model_name is not recognized: " << model_name << std::endl;
			std::cout << "        Usage instructions:" << std::endl;
			usage(argc, argv);
			exit(EXIT_FAILURE);
		}
		cfg_session.open(filename_params.c_str());
		if (cfg_session.is_open()) {
			char0="#";
			std::getline(cfg_session, line0);	
			while(!cfg_session.eof() && char0 == "#"){ // Jump comments lines in the header
				line0=strtrim(line0);
				char0=strtrim(line0.substr(0, 1));
				if (char0 == "#"){
					std::getline(cfg_session, line0);
				}
			}
			// After all the comments, the first line must contain the frequency range fmin,fmax
			range=str_to_Xdarr(line0, " \t");
			std::getline(cfg_session, line0);
			el=str_to_int(line0);
			std::getline(cfg_session, line0);
			fc_l=str_to_dbl(line0);
			std::getline(cfg_session, line0);
			H_l=str_to_dbl(line0);
			std::getline(cfg_session, line0);
			gamma_l=str_to_dbl(line0);
			std::getline(cfg_session, line0);
			a1=str_to_dbl(line0);
			if(model_name == "build_l_mode_a1a2a3"){
				std::getline(cfg_session, line0);
				a2=str_to_dbl(line0);
				std::getline(cfg_session, line0);
				a3=str_to_dbl(line0);
			}
			if(model_name == "build_l_mode_aj"){
				std::getline(cfg_session, line0);
				a2=str_to_dbl(line0);
				std::getline(cfg_session, line0);
				a3=str_to_dbl(line0);
				std::getline(cfg_session, line0);
				a4=str_to_dbl(line0);
				std::getline(cfg_session, line0);
				a5=str_to_dbl(line0);
				std::getline(cfg_session, line0);
				a6=str_to_dbl(line0);
				std::getline(cfg_session, line0);
				eta0=str_to_dbl(line0);
			}
			if(model_name == "build_l_mode_a1etaAlma3"){
				std::getline(cfg_session, line0);
				eta0=str_to_dbl(line0);
				std::getline(cfg_session, line0);
				epsilon=str_to_dbl(line0);
				std::getline(cfg_session, line0);
				thetas[0]=str_to_dbl(line0);
				std::getline(cfg_session, line0);
				thetas[1]=str_to_dbl(line0);
				std::getline(cfg_session, line0);
				a3=str_to_dbl(line0);
			}
			std::getline(cfg_session, line0);
			asym=str_to_dbl(line0);
			std::getline(cfg_session, line0);
			inc=str_to_dbl(line0);

			cfg_session.close();
 	 	} else {
   			std::cout << "Unable to open the file: " << filename_params << std::endl;
   			std::cout << "Check that the file exist and that the path is correct" << std::endl;
   			std::cout << "Cannot proceed" << std::endl;
   			std::cout << "The program will exit now" << std::endl;
   			exit(EXIT_FAILURE);
		}

		std::cout << "  3. Computing..." << std::endl;
		//std::cout << "el = " << el << std::endl;
		//std::cout << "inc: " << inc << std::endl;
		V=amplitude_ratio(el, inc);	
		//Delta=1e6/Cadence/2;
		//std::cout << "V = " << V.transpose() << std::endl;
		//std::cout << "range = " << range.transpose() << std::endl;
		Ndata=(range[1]-range[0])/range[2]; // (fmax - fmin)/resol
		//std::cout << "Ndata =" << Ndata << std::endl;
		x_l.setLinSpaced(Ndata, range[0], range[1]);
		if(model_name == "build_l_mode_a1a2a3"){
			model=build_l_mode_a1a2a3(x_l, H_l, fc_l, a1, a2, a3, asym, gamma_l, el, V);
		}
		if(model_name == "build_l_mode_aj"){
			model=build_l_mode_aj(x_l, H_l, fc_l, a1, a2, a3, a4, a5, a6, eta0, asym, gamma_l, el, V);
		}
		if(model_name == "build_l_mode_a1etaAlma3"){
			model=build_l_mode_a1etaAlma3(x_l, H_l, fc_l, a1, eta0, epsilon, thetas, a3, asym, gamma_l, el, V);
		}
		//std::cout << " Model done" << std::endl;
		write_model(x_l, model, file_out);
}

void write_model(const VectorXd& x, const VectorXd& y, std::string file_out){

	VectorXi Nchars(3), precision(3);

	std::ofstream outfile;

	Nchars << 20, 20, 20;
	precision << 10, 10, 7;

	outfile.open(file_out.c_str());
	if(outfile.is_open()){
        outfile << "#       freq (microHz)    model spectrum (ppm2/microHz)" << std::endl;
        
        for(int i=0; i<x.size(); i++){
			outfile << std::setw(Nchars[0]) << std::setprecision(precision[0]) << x[i];
			outfile << std::setw(Nchars[2]) << std::setprecision(precision[2]) << y[i];
            outfile << std::endl;
		}
		outfile.close();
	}  
	else {
		std::cout << " Unable to open file " << file_out << std::endl;	
		std::cout << " Check that the full path exists" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

}


int options(int argc, char* argv[]){

	std::string arg1;//, arg2;
	int val;
	
	val=-1;
	arg1="";
	
	if(argc == 2){
		arg1=argv[1];
		if(arg1 == "version"){
			 val=0;
		} 
	}
	if (argc == 3){
		//arg1=argv[1];
		//arg2=argv[2];
		val=1; // 
	}
	if (val == -1){ // Error code
		usage(argc, argv);
	} 
	if (val == 0){ // Version code
		showversion(); 
		exit(EXIT_SUCCESS);
	}
	if (val > 0 ){
		return val; // Execution code val
	} else{
		return -1; // Default value is to return an error code
	}
}

void showversion()
{
    std::cout << APP_NAME " - Lorentzian with a1a2a3-asym tool - " APP_VERSION "\n built on " __DATE__ << std::endl;

#   if defined(__clang__)
    	printf(" with clang " __clang_version__);
#   elif defined(__GNUC__)
    	printf(" with GCC");
    	printf(" %d.%d.%d", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#   elif defined(_MSC_VER)
    	printf(" with MSVC");
    	printf(" %d", MSVC_VERSION);
#   else
    printf(" unknown compiler");
#   endif

    std::cout << "\n features:";
#   if defined(__i386__) || defined(_M_IX86)
    std::cout << " i386" << std::endl;
#   elif defined(__x86_64__) || defined(_M_AMD64)
    std::cout << " x86_64" << std::endl;
#   endif
    std::cout << " Author: " << APP_COPYRIGHT << std::endl;

}


void usage(int argc, char* argv[]){

		if( argc < 2){
			std::cout << " This program is mostly a debug program for some of the functions inside the build_lorentzian.cpp files (the core function making Lorentzian models)" << std::endl;
			std::cout << " You need to provide one argument to that function. The available argument is: " << std::endl;
			std::cout << "     [1] The model_name. Availabe model names are:" << std::endl;
			std::cout << "           - build_l_mode_a1a2a3" << std::endl;
			std::cout << "           - build_l_mode_aj" << std::endl;
			std::cout << "           - build_l_mode_a1etaAlma3" << std::endl;
			std::cout << "     [2] The filename for the parameters to read. After comments ('#'), this files contains one parameter per line. List and order of the parameters: DP, alpha_g, delta0l, q_star, fl0_all, step" <<  std::endl;
			std::cout << "         THIS FILE IS NOT A DIRECT OUTPUT OF THE TAMCMC program. You must make it yourself" <<std::endl;
			std::cout << "         Parameters should appear in the following order (one type per line): [fmin, fmax, resol], l, fc_l, H_l, gamma_l, a1, a2, a3, asym, inclination" << std::endl;
			std::cout << " Call sequence: " << std::endl;
			std::cout << "     " << argv[0] << " [model_name]   [parameter filename]" << std::endl;
			std::cout << " Alternatively you can check the code version by typing " << argv[0] << " version" << std::endl;
			exit(EXIT_FAILURE);
		}
	
}

