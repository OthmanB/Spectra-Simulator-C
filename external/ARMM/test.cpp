
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "version_solver.h"
#include "data.h"
#include "string_handler.h"
#include "interpol.h"
#include "noise_models.h" // get the harvey_1985 function
#include "solver_mm.h"
#include "bump_DP.h"

int main(void)
{
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
	nu_m_from_l0=test_asymptotic_sg_O2from_l0();

	//std::cout << "Comparison O2p and O2from_l0:" << std::endl;
	//std::cout << "nu_m(O2p) : " << nu_from_O2p.nu_m.transpose() << std::endl;
	std::cout << "nu_m(l0) : " << nu_m_from_l0.nu_m.transpose() << std::endl;

	//std::cout << "nu_m_from_O2p.size()" << nu_from_O2p.nu_m.size() << std::endl;
	//std::cout << "nu_m_from_l0.size()" << nu_from_O2p.nu_m.size() << std::endl;
	std::cout << "-----" << std::endl;
	std::cout << "-----" << std::endl;
	std::cout << "-----" << std::endl;
	std::cout << "-----" << std::endl;
	std::cout << "-----" << std::endl;
	std::cout << "-----" << std::endl;
	std::cout << "-----" << std::endl;

	std::cout << " Testing Full build of a synthetic Sugbiant star ..." << std::endl;
	Cfg_synthetic_star cfg_star;
	cfg_star=test_make_synthetic_asymptotic_star_sg();
	std::cout << "-----" << std::endl;
	std::cout << "-----" << std::endl;
	std::cout << " Testing Full build of a synthetic RGB star ..." << std::endl;
	cfg_star=test_make_synthetic_asymptotic_star_rgb();
	
	std::cout << "All tests completed " << std::endl;
}


