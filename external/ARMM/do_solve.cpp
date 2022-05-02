/*
   A small function that compute the results of the solver as it is implemented into models.cpp for
   mixed modes models. Mostly for debug purpose
*/

#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <iomanip>
#include "string_handler.h"
#include "noise_models.h"
#include "version_solver.h"
#include "solver_mm.h"
#include "bump_DP.h"


void showversion();
int options(int argc, char* argv[]);
void usage(int argc, char* argv[]);

//Data Data_Nd2Data(Data_Nd dat);

int main(int argc, char* argv[]){

		bool verbose=0;
		int i, testval;
		Data_Nd data;
		//Config cfg;
		std::string filename_data, filename_params;
		std::string char0, line0;
		std::ifstream cfg_session;
		
		const long double sigma_p_l1=0;
		const long double sigma_m_l1=0;
		const long double sigma_limit=100; // Dummy value for test purpose only here
    	
    	double r; // Dummy value for the random qty for test purpose only here
		int el;
		double DPl, alpha_g, delta0l, q_star, step, fmin ,fmax;
		VectorXd fl0_all, h1_h0_ratio, ksi_pg;
	
		Data_eigensols freqs_l1;
    	
    	int i_dbg;
    	VectorXi posOK;  
    

		//Model_def model_list;
		// Set the current path variables
		std::string cpath=getcwd(NULL, 0);

		testval=options(argc, argv);

		//filename_data=argv[1];
		//filename_params=argv[2];
		//modelname=argv[3]; 
	
		filename_params=argv[1];
		
		
		std::cout << "  0. Configuration: " << std::endl;
//		std::cout << "      - Data file: " << filename_data << std::endl;
		std::cout << "      - Parameters file: " << filename_params << std::endl;
	
		/*std::cout << "  1. Reading the data file..." << std::endl;
		verbose=0;
		data=cfg.read_data_ascii_Ncols(filename_data, " \t", verbose); // the data in a matricial form				
		data_model=Data_Nd2Data(data);
		*/

		std::cout << "  2. Reading the file with the parameters for the solver_mm.cpp::solve_mm_asymptotic_O2from_l0() function..." << std::endl;
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
			// After all the comments, the first line must contain el
			el=str_to_int(line0);
			std::getline(cfg_session, line0);
			DPl=str_to_dbl(line0);
			std::getline(cfg_session, line0);
			alpha_g=str_to_dbl(line0);
			std::getline(cfg_session, line0);
			delta0l=str_to_dbl(line0);
			std::getline(cfg_session, line0);
			q_star=str_to_dbl(line0);
			std::getline(cfg_session, line0);
			fl0_all=str_to_Xdarr(line0, " \t"); // The separator is either a space of a tabulation
			std::getline(cfg_session, line0);
			step=str_to_dbl(line0);
			std::getline(cfg_session, line0);
			fmin=str_to_dbl(line0);
			std::getline(cfg_session, line0);
			fmax=str_to_dbl(line0);

			cfg_session.close();
 	 	} else {
   			std::cout << "Unable to open the file: " << filename_params << std::endl;
   			std::cout << "Check that the file exist and that the path is correct" << std::endl;
   			std::cout << "Cannot proceed" << std::endl;
   			std::cout << "The program will exit now" << std::endl;
   			exit(EXIT_FAILURE);
		}
		std::cout << "     Summary of the requested configuration: " << std::endl;
		std::cout << " step = " << step << std::endl; 
		std::cout << " fmin = " << fmin << std::endl;
		std::cout << " fmax = " << fmax << std::endl;
		std::cout << " el = " << el << std::endl;
		std::cout << " DPl = " << DPl << std::endl;
		std::cout << " alpha_g = " << alpha_g << std::endl;
		std::cout << " delta0l = " << delta0l << std::endl;
		std::cout << " q_star = " << q_star << std::endl;
		std::cout << " fl0_all = " << fl0_all.transpose() << std::endl;
		std::cout << "---" << std::endl;
		std::cout << "  3. Computing l=" << el << " frequencies..." << std::endl;
		freqs_l1=solve_mm_asymptotic_O2from_l0(fl0_all, el, delta0l, DPl, alpha_g, q_star, sigma_p_l1, step, true, false, fmin, fmax);

		// Generating widths profiles for l modes using the ksi function
    	std::cout << "  4. Computing l=" << el << " heights through the zeta function (see code for assumptions)..." << std::endl;
    	ksi_pg=ksi_fct2(freqs_l1.nu_m, freqs_l1.nu_p, freqs_l1.nu_g, freqs_l1.dnup, freqs_l1.dPg, q_star, "precise"); //"precise" // assume Dnu_p, DPl and q constant   	
    	h1_h0_ratio=h_l_rgb(ksi_pg); // WARNING: Valid assummption only not too evolved RGB stars (below the bump, see Kevin mail 10 August 2019)
    	
		std::cout << "  5. Writing Results on-screen... " << std::endl;
		std::cout << "     - Pure p modes:" << std::endl;
		std::cout << freqs_l1.nu_p.transpose() << std::endl;
		std::cout << "     - Derivative of pure p modes dnup/dn:" << std::endl;
		std::cout << freqs_l1.dnup.transpose() << std::endl;
		std::cout << "     - Pure g modes:" << std::endl;
		std::cout << freqs_l1.nu_g.transpose() << std::endl;
		std::cout << "     - Derivative of pure g modes dPg/dn:" << std::endl;
		std::cout << freqs_l1.dPg.transpose() << std::endl;
		std::cout << "     - Mixed modes:" << std::endl;
		std::cout << freqs_l1.nu_m.transpose() << std::endl;
		std::cout << "     - zeta_pg:" << std::endl;
		std::cout << ksi_pg.transpose() << std::endl;
		std::cout << "     - hl_h0 (does not account of mode visibility) :" << std::endl;
		std::cout << h1_h0_ratio.transpose() << std::endl;
		
		std::cout << " ---- " << std::endl;

}

/*
Data Data_Nd2Data(Data_Nd dat){

		Data dat_out;
		dat_out.x=dat.data.col(0);
		dat_out.y=dat.data.col(1);
		if (dat.data.cols() == 3){
			dat_out.sigma_y=dat.data.col(2);
		}
		dat_out.Nx=dat.data.rows();
	return dat_out;
}
*/

void showversion()
{
    std::cout << APP_NAME " - getmodel tool - " APP_VERSION "\n built on " __DATE__ << std::endl;

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

int options(int argc, char* argv[]){

	std::string arg1, arg2;
	int val;
	
	val=-1;
	arg1="";
	
	if(argc == 2){
		arg1=argv[1];
		if(arg1 == "version"){
			 val=0;
		} 
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

void usage(int argc, char* argv[]){

		if( argc < 2){
			std::cout << " This program allows you to use the solver as stand-alone program (in command line) " << std::endl;
			std::cout << " You need to provide an argument to the program, which is : " << std::endl;
			std::cout << "     - filename for the parameters to read (Text file). After comments ('#'), this files contains one parameter per line. List and order of the parameters: " << std::endl;
			std::cout << " (1) l" << std::endl;
			std::cout << " (2) DP" <<  std::endl;
			std::cout << " (3) alpha_g" << std::endl;
			std::cout << " (4) delta0l" << std::endl;
			std::cout << " (5) q_star" << std::endl;
			std::cout << " (6) fl0, a list of l=0 p mode frequencies that will be shift by l*Dnu/2 + delta0l to produce l=1 p modes" << std::endl;
			std::cout << " (7) step / resolution of the spectrum to define the precision on the retrieve modes" <<  std::endl;
			std::cout << " (8) fmin: Remove solutions found below fmin. Note that fl0 defines the total number of calculated frequencies. Those may/may not be below fmin" <<  std::endl;
			std::cout << " (9) fmax: Remove solutions found above fmax. Note that fl0 defines the total number of calculated frequencies. Those may/may not be above fmax" <<  std::endl;
			std::cout << " Call sequence: " << std::endl;
			//std::cout << "     " << argv[0] << " [data filename] [parameter filename]" << std::endl;
			std::cout << "     " << argv[0] << " [filename]" << std::endl;
			exit(EXIT_FAILURE);
		}
	
}
