
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "version_solver.h"
#include "data_solver.h"
#include "string_handler.h"
#include "interpol.h"
#include "noise_models.h" // get the harvey_1985 function
#include "solver_mm.h"
#include "bump_DP.h"

int main(void)
{

	// --------------------------------------------------
	// ---- Content to be modified for testing ----
	// --------------------------------------------------
	const int el=1;
	//const long double Dnu_p=60.;
	const long double beta_p=0.0076;
	//const long double delta0l_percent=2;
	const long double epsilon=0.4;
	//const long double DPl=400;
	const long double alpha_g=0.;
	const long double q=0.15;

	const int npmin=12;
	const int npmax=20;
	// --------------------------------------------

	//VectorXd nu_m_from_l0, nu_m_from_O2p;
	//Data_eigensols nu_from_l0, nu_from_O2p;
	Data_eigensols nu_m_from_l0;

	//std::cout << " Testing solver_mm() for the case of an SubGiant..." << std::endl;
	//test_sg_solver_mm();
	//std::cout << " Testing solver_mm() the case of a RedGiant..." << std::endl;
	//test_rgb_solver_mm();
	
	//std::cout << " Testing solve_mm_asymptotic_O2p() the case of a SugGiant..." << std::endl;
	//nu_from_O2p=test_asymptotic_sg_O2();
	std::cout << " Testing solve_mm_asymptotic_O2from_l0() the case of a SugGiant..." << std::endl;
	for(double DPl=70; DPl<200; DPl=DPl+5){
		for(double Dnu_p=15; Dnu_p<40;Dnu_p=Dnu_p+2){
			nu_m_from_l0=test_asymptotic_sg_O2from_l0_DEV(el, Dnu_p, beta_p, epsilon, DPl, alpha_g, q, npmin, npmax);
			std::cout << "DP = " << DPl << "    ,     Dnu_p = " << Dnu_p << ".... PASSED" << std::endl;
		}
	}
	//std::cout << "Comparison O2p and O2from_l0:" << std::endl;
	//std::cout << "nu_m(O2p) : " << nu_from_O2p.nu_m.transpose() << std::endl;
	//std::cout << "nu_m(l0) : " << nu_m_from_l0.nu_m.transpose() << std::endl;

	//std::cout << "nu_m_from_O2p.size()" << nu_from_O2p.nu_m.size() << std::endl;
	//std::cout << "nu_m_from_l0.size()" << nu_from_O2p.nu_m.size() << std::endl;
	std::cout << "-----" << std::endl;
	std::cout << "-----" << std::endl;
	std::cout << "-----" << std::endl;
	std::cout << "-----" << std::endl;
	std::cout << "-----" << std::endl;
	std::cout << "-----" << std::endl;
	std::cout << "-----" << std::endl;
/*
	std::cout << " Testing Full build of a synthetic Sugbiant star ..." << std::endl;
	Cfg_synthetic_star cfg_star;
	cfg_star=test_make_synthetic_asymptotic_star_sg();
	std::cout << "-----" << std::endl;
	std::cout << "-----" << std::endl;
	std::cout << " Testing Full build of a synthetic RGB star ..." << std::endl;
	cfg_star=test_make_synthetic_asymptotic_star_rgb();
*/
	std::cout << "All tests completed " << std::endl;
}


