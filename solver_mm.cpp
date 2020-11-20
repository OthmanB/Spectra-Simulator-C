// -------------------
// ---- Functions adapted from the solver_mm.py function ----
/* This contains all the functions for solving the asymptotic relation
# of the mixed modes, as they have been tested during their development
# All this arise from reading few papers from Benoit Mosser and 
# The PhD thesis from Charlotte Gehand:
# https://arxiv.org/pdf/1203.0689.pdf (Mosser paper on mixed modes)
# https://arxiv.org/pdf/1004.0449.pdf (older Mosser paper on scaling relations for gaussian_width, Amp etc.. - 2010 -)
# https://arxiv.org/pdf/1011.1928.pdf (The universal pattern introduced with the curvature - Fig. 3 - )
# https://arxiv.org/pdf/1411.1082.pdf
# https://tel.archives-ouvertes.fr/tel-02128409/document

# Examples and tests function have been built using asymptotic relations in the python original code.
# But note that they should be applicable to ANY set of value of:
#	 nu_p(nu), nu_g(nu), Dnu_p(nu) and DPl(nu) (meaning handling glitches)
*/
// ------------------
//#pragma once
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>

#include "version_solver.h"
#include "data.h"
#include "string_handler.h"
#include "interpol.h"
#include "derivatives_handler.h"

using Eigen::MatrixXd;
using Eigen::VectorXi;
using Eigen::VectorXd;

/*
// Taken from https://stackoverflow.com/questions/27028226/python-linspace-in-c
template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{

  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}
// -----
*/

// My linspace function
VectorXd linspace(const long double start_in, const long double end_in, const long num_in)
{
	if (num_in == 0) {
		std::cout << " num_in in linspace is 0. Cannot create a linspace vector. Returning -1." << std::endl;
		VectorXd linspaced(1);
		linspaced[0]=-1;
		return linspaced;
	}
	VectorXd linspaced(num_in);

	const long double delta = (end_in - start_in) / (num_in - 1);
	for(long i=0 ; i< num_in ; i++){
		linspaced[i]=start_in + delta*i;
	}

	/*std::cout << "start_in =" << start_in << std::endl;
	std::cout << "end_in =" << end_in << std::endl;
	std::cout << "num_in =" << num_in << std::endl;
	std::cout << "delta =" << delta << std::endl;
	std::cout << "linspaced.size() =" << linspaced.size() << std::endl;
	
	for(long i=0; i<linspaced.size(); i++)
	{
		std::cout << linspaced[i] << std::endl;
	}
	std::cout << "linspace needs a thorough test" << std::endl;
	std::cout << "exiting now"<< std::endl;
	exit(EXIT_SUCCESS);
	*/
	return linspaced;
}

// Function that detects sign changes
// If the sign went from + to - tag it with a -1
// If the sign went from - to + tag it with a +1
// If there is no change of sign tag it with a 0
// 0 is dealt as a zone of change of sign as well. eg. if we pass from 0 to 2, then the result is +1
// Inputs:
//    - x: input vector for which we want to know sign changes
//    - return_indices: if true (default), the function returns positions at which the sign changed
//						if false, it returns a vector of size(x)-1 with the tags for the sign changes (or not)
VectorXi sign_change(const VectorXd x, bool return_indices=true)
{
	bool bool_tmp;
	VectorXi s(x.size()-1), pos_s(x.size()-1);
	s.setConstant(0); // Vector of tags for sign changes
	pos_s.setConstant(-1); // Vector of indices
	long i, j=0; 
	for (i = 0; i<x.size()-1; i++)
	{
		if ((   x[i+1]>=0 && x[i] >=0   )|| (  x[i+1]<=0 && x[i]   )) // No sign change case or values are constant at 0
		{
			s[i]=0;
		}
		if (  (x[i+1]>=0 && x[i] <0) || (x[i+1]>0 && x[i] <=0)  ) // Sign change case from - to +
		{
			s[i]=1;
			pos_s[j]=i;
			j=j+1;
		}
		if (  (x[i+1]<=0 && x[i] >0) || (x[i+1]<=0 && x[i] >=0)   ) // Sign change case from + to -
		{
			s[i]=-1;
			pos_s[j]=i;
			j=j+1;
		}
	}
	pos_s.conservativeResize(j); // Discard slots that where not used in the vector
	/*
	std::cout << "sign_change needs a thorough test" << std::endl;
	std::cout << " x.size() = " << x.size() << std::endl; 
	for (i=0; i<x.size()-1; i++)
	{
		std::cout << "[" << i << "] " << x[i] << "  " << x[i+1] << "  " << s[i] << std::endl;
	}
	std::cout << "-----" << std::endl;
	for (i=0; i<pos_s.size(); i++)
	{
		std::cout << pos_s[i] << std::endl;
	}
	std::cout << "exiting now"<< std::endl;
	exit(EXIT_SUCCESS);
	*/
	if (return_indices == true)
	{
		return pos_s;
	} else
	{
		return s;
	}
}


VectorXd pnu_fct(const VectorXd nu, const long double nu_p)
{
	if (nu.size() == 0){
		std::cout << "Vector nu in pnu_fct is of size 0. Cannot pursue." << std::endl;
		std::cout << "The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	} else{
		VectorXd tmp(nu.size()), pnu(nu.size());
		tmp.setConstant(nu_p);
		pnu=nu - tmp;
		return pnu;
	}
}

long double pnu_fct(const long double nu, const long double nu_p)
{
	long double pnu=nu-nu_p;
	return pnu;
}

VectorXd gnu_fct(const VectorXd nu, const long double nu_g, const long double Dnu_p, const long double DPl, const long double q)
{
	const long double pi = 3.141592653589793238L;
	VectorXd X(nu.size()), gnu(nu.size()), tmp(nu.size());
	tmp.setConstant(nu_g);

	X=pi * (nu.cwiseInverse() - tmp.cwiseInverse())*1e6 / DPl;
	tmp=q*X.array().tan(); //.tan();
	gnu=Dnu_p * tmp.array().atan()/pi; //.atan() / pi;
	return gnu;
}

long double gnu_fct(const long double nu, const long double nu_g, const long double Dnu_p, const long double DPl, const long double q)
{
	const long double pi = 3.141592653589793238L;
	long double X, gnu;
	
	X=pi * (1./nu - 1./nu_g)*1e6 / DPl;
	gnu=Dnu_p * atan(q*tan(X))/pi; //.atan() / pi;
	return gnu;
}

/* A small function that generate a serie of p modes using the asymptotic relation
# at the second order as per defined in Mosser et al. 2018, equation 22 (https://www.aanda.org/articles/aa/pdf/2018/10/aa32777-18.pdf)
# delta0l and alpha and nmax must be set, a
# Note that we have the following relationship between D0 and delta0l:
#			delta0l=-l(l+1) D0 / Dnu_p
# Such that delta0l=-l(l+1) gamma / 100, if gamma is in % of Dnu_p
# r: Allows you to add an extra term to Dnu_o
*/
long double asympt_nu_p(const long double Dnu_p, const int np, const long double epsilon, const int l, 
	const long double delta0l, const long double alpha, const long double nmax, long double r=0)
{

/*
	std::cout << " np=" << np << std::endl;
	std::cout << " epsilon=" << epsilon << std::endl;
	std::cout << " l=" << l << std::endl;
	std::cout << " delta0l=" << delta0l << std::endl;
	std::cout << " alpha=" << alpha << std::endl;
	std::cout << " nmax=" << nmax << std::endl;
	std::cout << " Dnu_p=" << Dnu_p << std::endl;
*/	
	long double nu_p=(np + epsilon + l/2. + delta0l + alpha*std::pow(np - nmax, 2) / 2)*Dnu_p;
	if (nu_p < 0.0)
	{
		std::cout << " WARNING: NEGATIVE FREQUENCIES DETECTED: IMPOSING POSITIVITY" << std::endl;
		std::cout << " nu_p: " << nu_p << std::endl;
		std::cout << " Cannot pursue " << std::endl;
		exit(EXIT_FAILURE);
	}
	return nu_p+r;
}

/* 
	Compute the asymptotic relation for the g modes.
	r: an optional parameter that can be added to the Period (e.g. a random quantity)	
*/
long double asympt_nu_g(const long double DPl, const int ng, const long double alpha, long double r=0)
{
	const long double Pl=(ng + alpha)*DPl;
	return 1e6/(Pl+r);
}

/*
This the main function that solves the mixed mode asymptotic relation
which is of the type p(nu) = g(nu)
This solver specifically solve the case:
      nu - nu_p = Dnu*arctan(q tan(1/(nu DPl) - 1/(nu_g*DPl)))
      It tries to find the intersect between p(nu) and g(nu)
      using an interpolation and over a range of frequency such that nu is in [numin, numax]
Parameters:
	- Mandatory: 
	     nu_p (double) : frequency of a given p-mode (in microHz)
	     nu_g (double): frequency of a given g-mode (in microHz)
	     Dnu_p (double): Large separation for p modes (in microHz)
	     DP1 (double): Period spacing for g modes (in seconds)
	     q (double): Coupling term (no unit, should be between 0 and 1)
	- Optional:
		numin (double): Minimum frequency considered for the solution (in microHz)
		numax (double): Maximum frequency considered for the solution (in microHz)
		resol (double): Base resolution for the interpolated base function. The interpolation may miss solutions 
		       if this is set too low. Typically, the resolution parameter should be higher than the
		       spectral resolution of your spectrum so that all resolved modes should be found.
		       This is also used for creating the nu axis for visualisation (in microHz).
		factor (double): Define how fine will be the new tiny grid used for performing the interpolation. This is important
				to avoid extrapolation (which is forbiden and will result in crash of the code). Typically, the default
				value factor=0.05 can compute mixed modes for frequency down to 80microHz. Going below requires a smaller factor

		returns_axis: If True, returns nu, pnu and gnu (see optional reutrns below). Mainly for debug
Returns a structure with:
	nu_m: An array with all solutions that match p(nu) = g(nu)
	nu (optional): The frequency axis used as reference for finding the intersection
	pnu (optional): The curve for p(nu)
	gnu (optional): The curve g(nu)
*/
Data_coresolver solver_mm(const long double nu_p, const long double nu_g, const long double Dnu_p, const long double DPl, const long double q, 
	const long double numin, const long double numax, const long double resol, const bool returns_axis=false, const bool verbose=false, const long double factor=0.05)
{
	//const int Nmmax=500; // Number of maximum mixed modes solutions that can be found 

	int i, s_ok, Nsize=0;
	long double range_min, range_max, nu_m_proposed, ratio,  ysol_pnu, ysol_gnu;
	Data_coresolver results;
	VectorXi idx;
	VectorXd nu, pnu, gnu, nu_local, pnu_local, gnu_local, nu_m, ysol_all;

	// Generate a frequency axis that has a fixed resolution and that span from numin to numax
	nu=linspace(numin, numax, long((numax-numin)/resol));
	/*
	std::cout << " numin =" << numin << std::endl;
	std::cout << " numax =" << numax << std::endl;
	std::cout << " resol =" << resol << std::endl;
	std::cout << "long((numax-numin)/resol) = " << long((numax-numin)/resol) << std::endl;
	*/
	// Function p(nu) describing the p modes
	pnu=pnu_fct(nu, nu_p);
	// Function g(nu) describing the g modes 
	gnu=gnu_fct(nu, nu_g, Dnu_p, DPl, q);
	
	/* Find when p(nu) = g(nu) by looking for solution of p(nu) - g(nu) = 0
	#     Method 1: Direct Interpolation... Works only for single solutions ==> Not used here
	#int_fct = interpolate.interp1d(pnu - gnu, nu)
	#nu_m=int_fct(0)
	#     Method 2: (a) Find indices close to sign changes for p(nu) - g(nu)
	#               (b) Then perform an iterative interpolation in narrow ranges
	#                   near the approximate solutions. How narrow is the range is defined
	#					by the resolution parameter resol, which in this case can be view
	#					as the minimum precision.
	*/
	idx=sign_change(pnu-gnu);
	//std::cout << "inside solver_mm: passed idx,   idx=" << idx << std::endl;	
	/*
	std::cout << " === C++ === " << std::endl;
	std::cout << " nu_p =" << nu_p << std::endl;
	std::cout << " nu_g =" << nu_g << std::endl;
	std::cout << " len(nu)= " << nu.size() << std::endl;
	std::cout << "' len(pnu)= " << pnu.size() << std::endl;
	std::cout << "' len(gnu)= " << gnu.size() << std::endl;
	std::cout << "idx =" << idx << std::endl;
	std::cout << "idx.size() = " << idx.size() << std::endl;
	*/
	s_ok=0;
	nu_m.resize(idx.size());
	for (long ind=0; ind<idx.size();ind++)
	{
		//std::cout << "idx[" << ind << "] =" << idx[ind] << std::endl;
		// Define a small local range around each of the best solutions
		range_min=nu[idx[ind]] - 2*resol;
		range_max=nu[idx[ind]] + 2*resol;
		// Redefine nu, pnu and gnu for that local range
		nu_local=linspace(range_min, range_max, long((range_max-range_min)/(resol*factor)));
		pnu_local=pnu_fct(nu_local, nu_p);
		gnu_local=gnu_fct(nu_local, nu_g, Dnu_p, DPl, q);	

		// Perform the interpolation on the local range and append the solution to the nu_m list
		nu_m_proposed=lin_interpol(pnu_local - gnu_local, nu_local, 0);
		try
		{	
			ysol_gnu=gnu_fct(nu_m_proposed, nu_g, Dnu_p, DPl, q);
			ysol_pnu=pnu_fct(nu_m_proposed, nu_p);
		}
		catch(...)
		{
			std::cout << "Interpolation issue detected. Debuging information:" << std::endl;
			std::cout << "    nu_p: " <<  nu_p << std::endl;
			std::cout << "    nu_g: " <<  nu_g << std::endl;
			std::cout << "    Dnu_p: "<< Dnu_p << std::endl;
			std::cout << "    DPl: "<< DPl << std::endl;
			std::cout << "    q: " << q << std::endl;
			std::cout << "    numin: "<< numin << std::endl;
			std::cout << "    numax: "<< numax << std::endl;
			std::cout << "    resol:"<< resol << std::endl;
			std::cout << "    factor:"<< factor << std::endl;
			std::cout << " ------------" << std::endl;
			std::cout << "range_min/max: "<< range_min << range_max << std::endl;
			std::cout << "  nu_local: "<< nu_local << std::endl;
			std::cout << "  pnu_local: "<< pnu_local << std::endl;
			std::cout << "  gnu_local: "<< gnu_local << std::endl;
			std::cout << " ------------" << std::endl;
			std::cout << " int_fct  ==>  nu_local      /   pnu_local - gnu_local : " << std::endl;
			for (i=0; i<nu_local.size(); i++)
			{
				std::cout << "    " <<  nu_local[i]<< pnu_local[i]-gnu_local[i] << std::endl;
			}
			exit(EXIT_FAILURE);
		}
		ratio=ysol_gnu/ysol_pnu;
		if (verbose == true)
		{	std::cout << "-------"<< std::endl;
			std::cout << "nu_m:"<<  nu_m_proposed<< std::endl;
			std::cout << "Ratio:"<<  ratio << std::endl;
		}
		// Sometimes, the interpolator mess up due to the limits of validity for the atan function
		// The way to keep real intersection is to verify after interpolation that we really
		// have p(nu_m_proposed) = g(nu_m_proposed). We then only keeps solutions that satisfy
		// a precision criteria of 0.1%.
		/*if ((ratio >= 0.999) && (ratio <= 1.001))
		{
			std::cout << "   solver_mm : Adding solution only if it is really a solution" << std::endl;
			Nsize=nu_m.size() +1;
			nu_m.conservativeResize(Nsize);
			nu_m[Nsize-1]=nu_m_proposed;
			std::cout << "   solver_mm : Solution added" << std::endl;
		}
		*/
		if ((ratio >= 0.999) && (ratio <= 1.001))
		{
			nu_m[s_ok]=nu_m_proposed;
			s_ok=s_ok+1;
		}
	}
	nu_m.conservativeResize(s_ok);

	//ysol_gnu=gnu_fct(nu_m, nu_g, Dnu_p, DPl, q);
	ysol_all=gnu_fct(nu_m, nu_g, Dnu_p, DPl, q);
	
	if (returns_axis == true){
		results.nu_m=nu_m;
		results.ysol=ysol_all;
		results.nu=nu;
		results.pnu=pnu;
		results.gnu=gnu;
		return results;
	}
	else{
		results.nu_m=nu_m;
		results.ysol=ysol_all;
		return results;
	}
}

// Function to test solver_mm()
// This is a typical SG case, with  density of g modes << density of p modes
void test_sg_solver_mm()
{

	Data_coresolver sols;

	const long double Dnu_p=60; //microHz
	const long double DP1= 400; //microHz, typical for a RGB with Dnu_p=10 (see Fig. 1 of Mosser+2014, https://arxiv.org/pdf/1411.1082.pdf)

	// Generate a p-mode that follow exactly the asymptotic relation of p modes
	const long double D0=Dnu_p/100. ;
	const long double epsilon=0.4;
	const long double np=10.;
	const long double nu_p=(np + epsilon + 1./2.)*Dnu_p - 2*D0;

	// Generate a g-mode that follow exactly the asymptotic relation of g modes for a star with radiative core
	const long double ng=10;
	const long double alpha=0;
	const long double nu_g=1e6/(ng*DP1);

	// Use the solver
	const long double q=0.2; // Fix the coupling term
	const long double resol=0.01;
	sols=solver_mm(nu_p, nu_g, Dnu_p, DP1, q, nu_p - Dnu_p/2, nu_p + Dnu_p/2, resol, true, true, 0.05);
	//std::cout << "nu_m: "  <<  sols.nu_m << std::endl;
}


// This function uses solver_mm to find solutions from a spectrum
// of pure p modes and pure g modes following the asymptotic relations at the second order for p modes and the first order for g modes
//
//	Dnu_p: Average large separation for the p modes
//	epsilon: phase offset for the p modes
//	el: Degree of the mode
//	delta0l: first order shift related to core structure (and to D0)
//	alpha_p: Second order shift relate to the mode curvature
//	nmax: radial order at numax
//	DPl: average Period spacing of the g modes
//	alpha: phase offset for the g modes
//	q: coupling strength
//	sigma_p: standard deviation controling the randomisation of individual p modes. Set it to 0 for no spread
//	sigma_g: standard deviation controling the randomisation of individial g modes. Set it to 0 for no spread
//  fmin: minimum frequency to consider for p modes. Note that mixed mode solutions may be of lower frequency than this
//  fmax: maximum frequency to consider for p modes. Note that mixed mode solutions may be of lower frequency than this
//  resol: Control the grid resolution. Might be set to the resolution of the spectrum
//  returns_pg_freqs: If true, returns the values for calculated p and g modes
//  verbose: If true, print the solution on screen
Data_eigensols solve_mm_asymptotic_O2p(const long double Dnu_p, const long double epsilon, const int el, const long double delta0l, const long double alpha_p, 
	const long double nmax, const long double DPl, const long double alpha, const long double q, const long double sigma_p, 
	const long double fmin, const long double fmax, const long double resol, bool returns_pg_freqs=true, bool verbose=false)
{

	const bool returns_axis=true;
	//const int Npmax=100;
	//const int Ngmax=2000;
	const int Nmmax=5000; //Ngmax+Npmax;
	const int Nmax_attempts=4;
	const double tol=2*resol; // Tolerance while searching for double solutions of mixed modes
	const double sigma_g=0;  // FOR SOME REASON, WE FIND MULTIPLE SOLUTIONS WITHIN A NARROW RANGE WHEN USING SIGMA_G... SO IT TURNED IT OFF

	unsigned seed_p = std::chrono::system_clock::now().time_since_epoch().count();
	//unsigned seed_g = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine gen_p(seed_p); //, gen_g(seed_g);
	std::normal_distribution<double> distrib_p(0.,sigma_p);
	//std::normal_distribution<double> distrib_g(0.,sigma_g);

	bool success;
	int s0m, i, attempts, np_min, np_max, ng_min, ng_max;
	double nu_p, nu_g, Dnu_p_local, DPl_local; // Dnu_p_local and DPl_local are important if modes does not follow exactly the asymptotic relation.
	double fact=0.04;  // Default factor
	double r;

	VectorXi test;
	VectorXd nu_p_all, nu_g_all, nu_m_all(Nmmax), results(Nmmax);	
	//VectorXd nu_p_all(Nmmax), nu_g_all(Nmmax), nu_m_all(Nmmax), results(Nmmax);

	Data_coresolver sols_iter;
	Data_eigensols nu_sols;
	Deriv_out deriv_p, deriv_g;

	// Use fmin and fmax to define the number of pure p modes and pure g modes to be considered
	np_min=int(floor(fmin/Dnu_p - epsilon - el/2 - delta0l));
	np_max=int(ceil(fmax/Dnu_p - epsilon - el/2 - delta0l));

	np_min=int(floor(np_min - alpha*std::pow(np_min - nmax, 2) /2.));
	np_max=int(ceil(np_max + alpha*std::pow(np_max - nmax, 2) /2.)); // CHECK THIS DUE TO - -

	ng_min=int(floor(1e6/(fmax*DPl) - alpha));
	ng_max=int(ceil(1e6/(fmin*DPl) - alpha));

	if (np_min <= 0)
	{
		np_min=1;
	}
	if (fmin <= 150) // overrides of the default factor in case fmin is low
	{
		fact=0.01;
	}
	if (fmin <= 50)
	{
		fact=0.005;
	}

	// Handling the p and g modes, randomized or not
	nu_p_all.resize(np_max-np_min);
	nu_g_all.resize(ng_max-ng_min);
	for (int np=np_min; np<np_max; np++)
	{
		if (sigma_p == 0)
		{
			nu_p=asympt_nu_p(Dnu_p, np, epsilon, el, delta0l, alpha_p, nmax);
		} else{
			r = distrib_p(gen_p);
			nu_p=asympt_nu_p(Dnu_p, np, epsilon, el, delta0l, alpha_p, nmax, r);
		}
		nu_p_all[np-np_min]=nu_p;
	}

	for (int ng=ng_min; ng<ng_max;ng++)
	{
		//if (sigma_g == 0)
		//{
			nu_g=asympt_nu_g(DPl, ng, alpha);
		//} else{
		//	r = distrib_g(gen_g);
		//	nu_g=asympt_nu_g(DPl, ng, alpha, r);
		//}
		nu_g_all[ng-ng_min]=nu_g;
	}

	if (sigma_p != 0)
	{
		deriv_p=Frstder_adaptive_reggrid(nu_p_all);
	} else{
		deriv_p.deriv.resize(nu_p_all.size());
		deriv_p.deriv.setConstant(Dnu_p);
	}
	deriv_g.deriv.resize(nu_g_all.size());
	deriv_g.deriv.setConstant(DPl);

	//if (sigma_g != 0)
	//{
	//	deriv_g=Frstder_adaptive_reggrid(1e6 * nu_g_all.cwiseInverse());
	//}

	//std::cout << " np_min = " << np_min << std::endl;
	//std::cout << " np_max = " << np_max << std::endl;
	//std::cout << " ng_min = " << ng_min << std::endl;
	//std::cout << " ng_max = " << ng_max << std::endl;
	//std::cout << "fmin=" << fmin << std::endl;
	//std::cout << "fmax=" << fmax << std::endl;
	s0m=0;
	for (int np=np_min; np<np_max; np++)
	{
		for (int ng=ng_min; ng<ng_max;ng++)
		{
			//nu_p=asympt_nu_p(Dnu_p, np, epsilon, el, delta0l, alpha_p, nmax);
			//nu_g=asympt_nu_g(DPl, ng, alpha);
			nu_p=nu_p_all[np-np_min];
			nu_g=nu_g_all[ng-ng_min];
			
			// This is the local Dnu_p which differs from the average Dnu_p because of the curvature. The solver needs basically d(nu_p)/dnp , which is Dnu if O2 terms are 0.
			if (sigma_p == 0){
				Dnu_p_local=Dnu_p*(1. + alpha_p*(np - nmax));
			} else{ // Due to the randomisation from sigma_p, we need to evaluate the local Dnu by the means of the first derivative
				Dnu_p_local=deriv_p.deriv[np-np_min];
			} 
			//if (sigma_g == 0){
				DPl_local=DPl; // The solver needs here d(nu_g)/dng. Here we assume no core glitches so that it is the same as DPl. 	
			//} else{ // Due to the randomisation from sigma_g, we need to evaluate the local DPl by the means of the first derivative
			//	DPl_local=deriv_g.deriv[ng-ng_min];
			//}
			try
			{
				sols_iter=solver_mm(nu_p, nu_g, Dnu_p_local, DPl_local, q, nu_p - 3.*Dnu_p/4, nu_p + 3.*Dnu_p/4, resol, returns_axis, verbose, fact);
			}
			catch (...){
				success=false;
				attempts=0;
				try{
					while (success ==false && attempts < Nmax_attempts){
						try{
							fact=fact/2;
							sols_iter=solver_mm(nu_p, nu_g, Dnu_p_local, DPl_local, q, nu_p - 3.*Dnu_p/4, nu_p + 3.*Dnu_p/4, resol, returns_axis, verbose, fact);
							success=true;
						}
						catch(...){
							std::cout << " Problem with the fine grid when searching for a solution... attempting to reduce factor to " << fact << "..." << std::endl;
						}
					}
				}
				catch(...){
						std::cout << "ValueError in solver_mm... Debug information:"<< std::endl;
						std::cout << " We excedeed the number of attempts to refine the grid by reducing factor" << std::endl;
						std::cout << " np_min = " << np_min << std::endl;
						std::cout << " np_max = " << np_max << std::endl;
						std::cout << " ng_min = " << ng_min << std::endl;
						std::cout << " ng_max = " << ng_max << std::endl;
						std::cout << " ---------- " << std::endl;			
						std::cout << " Dnu_p = " << Dnu_p << std::endl;
						std::cout << " np = " << np << std::endl;
						std::cout << " epsilon= " << epsilon << std::endl;
						std::cout << " delta0l= " << delta0l << std::endl;
						std::cout << " alpha_p= " << alpha_p << std::endl;
						std::cout << " nmax= " << nmax << std::endl;
						std::cout << " ---------- " << std::endl;
						std::cout << "   nu_p: " << nu_p << std::endl;
						std::cout << "   nu_g: " << nu_g << std::endl;
						std::cout << "   Dnu_p: " << Dnu_p << std::endl;
						std::cout << "   DPl: " << DPl << std::endl;
						std::cout << "   q: " << q << std::endl;
						std::cout << "   numin=nu_p - Dnu_p: " << nu_p - Dnu_p << std::endl;
						std::cout << "   numax=nu_p + Dnu_p: " << nu_p + Dnu_p << std::endl;
						std::cout << "   resol: " << resol << std::endl;
						std::cout << "   factor: " << fact << std::endl;
						exit(EXIT_FAILURE);
				}
			}
			if (verbose == true)
			{
				std::cout << "=========================================="  << std::endl;
				std::cout << "nu_p: " << nu_p << std::endl;
				std::cout << "nu_g: " << nu_g << std::endl;
				std::cout << "solutions nu_m: " << sols_iter.nu_m << std::endl;
			}
			for (int s=0;s<sols_iter.nu_m.size();s++)
			{
				// Cleaning doubles: Assuming exact matches or within a tolerance range
				test=where_dbl(nu_m_all, sols_iter.nu_m[s], tol);
				if (test[0] == -1)
				{
					/*std::cout << " adding a solution" << std::endl;
					std::cout << "    nu_m_all.size()= " << nu_m_all.size() << std::endl;
					std::cout << "    nu_p_all.size()= " << nu_p_all.size() << std::endl;
					std::cout << "    nu_g_all.size()= " << nu_g_all.size() << std::endl;
					std::cout << "    s0m               =" << s0m << std::endl;
					std::cout << "    s                =" << s << std::endl;
					std::cout << "    sols_iter.nu_m[s]=" << sols_iter.nu_m[s] << std::endl;
					std::cout << "    nu_p             =" << nu_p << std::endl;
					std::cout << "    nu_g             =" << nu_g << std::endl;
					*/
					nu_m_all[s0m]=sols_iter.nu_m[s];
					s0m=s0m+1;
				}
			}
		}
	}
	nu_m_all.conservativeResize(s0m);	
	//nu_p_all.conservativeResize(s0m);	
	//nu_g_all.conservativeResize(s0m);	

	if (returns_pg_freqs == true)
	{
		nu_sols.nu_m=nu_m_all;
		nu_sols.nu_p=nu_p_all;
		nu_sols.nu_g=nu_g_all;
		nu_sols.dnup=deriv_p.deriv;
		nu_sols.dPg=deriv_g.deriv;
		return nu_sols;
	} else
	{
		nu_sols.nu_m=nu_m_all;
		return nu_sols;
	}
}

// Function to test solver_mm()
// This is a typical RGB case, with  density of g modes >> density of p modes
void test_rgb_solver_mm()
{
	Data_coresolver sols;

	const long double Dnu_p=15; // microHz
	const long double DP1= 80; // microHz, typical for a RGB with Dnu_p=10 (see Fig. 1 of Mosser+2014, https://arxiv.org/pdf/1411.1082.pdf)

	// Generate a p-mode that follow exactly the asymptotic relation of p modes
	const long double D0=Dnu_p/100.; 
	const long double epsilon=0.4;
	const long double np=10.;
	const long double nu_p=(np + epsilon + 1./2.)*Dnu_p - 2*D0;
	// Generate a g-mode that follow exactly the asymptotic relation of g modes for a star with radiative core
	const long double ng=50;
	const long double alpha=0.01;
	const long double nu_g=1e6/(ng*DP1);

	// Use the solver
	const long double q=0.1; // Fix the coupling term
	const long double resol=0.005;
	sols=solver_mm(nu_p, nu_g, Dnu_p, DP1, q, nu_p - Dnu_p/2, nu_p + Dnu_p/2, resol, true, true, 0.05);
	//std::cout << "nu_m: "  <<  sols.nu_m << std::endl;

}


/* Function to test solve_mm_asymptotic
# The parameters are typical for a RGB in the g mode asymptotic regime
# Default parameters are for an early SG... The asymptotic is not accurate then
# consider: test_asymptotic(el=1, Dnu_p=30, beta_p=0.01, gamma0l=2., epsilon=0.4, DPl=110, alpha_g=0., q=0.15)
# for a RGB
*/
void test_asymptotic_rgb()
{

	// --------------------------------------------------
	// ---- Content to be modified for testing a RGB ----
	// --------------------------------------------------
	
	const int el=1;
	const long double Dnu_p=10.;
	const long double beta_p=0.076;
	const long double delta0l_percent=2;
	const long double epsilon=0.4;
	const long double DPl=70;
	const long double alpha_g=0.;
	const long double q=0.15;
	
	// --------------------------------------------

	// --------------------------------------------------
	// ---- Content to be modified for testing a SG ----
	// --------------------------------------------------
	/*
	const int el=1;
	const long double Dnu_p=60.;
	const long double beta_p=0.0076;
	const long double delta0l_percent=2;
	const long double epsilon=0.4;
	const long double DPl=400;
	const long double alpha_g=0.;
	const long double q=0.15;
	*/
	// --------------------------------------------

	// Define global Pulsation parameters
	// Parameters for p modes that follow exactly the asymptotic relation of p modes
	const long double D0=Dnu_p/100.;
	const long double delta0l=-el*(el + 1) * delta0l_percent / 100.;

	// Parameters for g modes that follow exactly the asymptotic relation of g modes for a star with radiative core
	const long double alpha=0.;

	// Define the frequency range for the calculation by (1) getting numax from Dnu and (2) fixing a range around numax
	const long double beta0=0.263; // according to Stello+2009, we have Dnu_p ~ 0.263*numax^0.77 (https://arxiv.org/pdf/0909.5193.pdf)
	const long double beta1=0.77; // according to Stello+2009, we have Dnu_p ~ 0.263*numax^0.77 (https://arxiv.org/pdf/0909.5193.pdf)
	const long double nu_max=std::pow(10, log10(Dnu_p/beta0)/beta1);

	const long double fmin=nu_max - 6*Dnu_p;
	const long double fmax=nu_max + 4*Dnu_p;

	const long double nmax=nu_max/Dnu_p - epsilon;
	const long double alpha_p=beta_p/nmax;

	// Fix the resolution to 4 years (converted into microHz)
	const long double data_resol=1e6/(4.*365.*86400.);

	const long double sigma_p=0.10*Dnu_p;
	//const long double sigma_g=0.0*DPl;
	Data_eigensols freqs;

	//std::cout << "nmax= " << nmax << std::endl;

	// Use the solver
	freqs=solve_mm_asymptotic_O2p(Dnu_p, epsilon, el, delta0l, alpha_p, nmax, DPl, alpha_g, q, sigma_p, fmin, fmax, data_resol, true, false);

	std::cout << " --- Lenghts ----"  << std::endl;
	std::cout << " L(nu_m): " << freqs.nu_m.size()  << std::endl;

	std::cout << " nu_m =" << freqs.nu_m << std::endl;
}

/* Function to test solve_mm_asymptotic
# The parameters are typical for a RGB in the g mode asymptotic regime
# Default parameters are for an early SG... The asymptotic is not accurate then
# consider: test_asymptotic(el=1, Dnu_p=30, beta_p=0.01, gamma0l=2., epsilon=0.4, DPl=110, alpha_g=0., q=0.15)
# for a RGB
*/
void test_asymptotic_sg()
{

	// --------------------------------------------------
	// ---- Content to be modified for testing a SG ----
	// --------------------------------------------------
	const int el=1;
	const long double Dnu_p=60.;
	const long double beta_p=0.0076;
	const long double delta0l_percent=2;
	const long double epsilon=0.4;
	const long double DPl=400;
	const long double alpha_g=0.;
	const long double q=0.15;
	// --------------------------------------------

	// Define global Pulsation parameters
	// Parameters for p modes that follow exactly the asymptotic relation of p modes
	const long double D0=Dnu_p/100.;
	const long double delta0l=-el*(el + 1) * delta0l_percent / 100.;

	// Parameters for g modes that follow exactly the asymptotic relation of g modes for a star with radiative core
	const long double alpha=0.;

	// Define the frequency range for the calculation by (1) getting numax from Dnu and (2) fixing a range around numax
	const long double beta0=0.263; // according to Stello+2009, we have Dnu_p ~ 0.263*numax^0.77 (https://arxiv.org/pdf/0909.5193.pdf)
	const long double beta1=0.77; // according to Stello+2009, we have Dnu_p ~ 0.263*numax^0.77 (https://arxiv.org/pdf/0909.5193.pdf)
	const long double nu_max=std::pow(10, log10(Dnu_p/beta0)/beta1);

	const long double fmin=nu_max - 6*Dnu_p;
	const long double fmax=nu_max + 4*Dnu_p;

	const long double nmax=nu_max/Dnu_p - epsilon;
	const long double alpha_p=beta_p/nmax;

	// Fix the resolution to 4 years (converted into microHz)
	const long double data_resol=1e6/(4.*365.*86400.);

	const long double sigma_p=0.1*Dnu_p;
	//const long double sigma_g=0.01*DPl;
	Data_eigensols freqs;

	//std::cout << "nmax= " << nmax << std::endl;

	// Use the solver
	freqs=solve_mm_asymptotic_O2p(Dnu_p, epsilon, el, delta0l, alpha_p, nmax, DPl, alpha_g, q, sigma_p, fmin, fmax, data_resol, true, true);

	std::cout << " --- Lenghts ----"  << std::endl;
	std::cout << " L(nu_m): " << freqs.nu_m.size()  << std::endl;
	std::cout << " L(nu_p): " << freqs.nu_p.size()  << std::endl;
	std::cout << " L(nu_g): " << freqs.nu_g.size()  << std::endl;

	std::cout << " nu_m =" << freqs.nu_m << std::endl;
}

/*
int main(void)
{
	//std::cout << " Testing solver_mm() for the case of an SubGiant..." << std::endl;
	//test_sg_solver_mm();
	//std::cout << " Testing solver_mm() the case of a RedGiant..." << std::endl;
	//test_rgb_solver_mm();
	//std::cout << " Testing solve_mm_asymptotic_O2p() the case of a RedGiant..." << std::endl;
	//test_asymptotic_rgb();
	std::cout << " Testing solve_mm_asymptotic_O2p() the case of a SugGiant..." << std::endl;
	test_asymptotic_sg();
}
*/